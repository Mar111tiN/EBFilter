#! /usr/bin/env python
import sys
import re
import os
import pandas as pd
import numpy as np
import math
import subprocess
import multiprocessing
from multiprocessing import Pool
from functools import partial
from . import utils
from . import anno
from . import vcf
from .cache import generate_cache, get_EB_from_cache
from datetime import datetime as dt


def main(args, config):
    '''
    validates files and refers to respective functions
    '''

    # ############## ARGUMENTS and CONFIG ##################################
    config['cache_mode'] = False
    config['generate_cache'] = args['generate_cache'] if 'generate_cache' in args.keys() else False

    threads = config['threads']
    debug_mode = config['debug_mode']
    # store list and pon_df in pon_dict, store pon chroms in config['pon_chr']
    pon_dict = utils.validate_pon(args['pon_list'], config)

    # ##### set chromosome ###########################
    # set config['chr'] to chromosome if provided in makeEBcache
    config['chr'] = [args['chrom']] if 'chrom' in args.keys() else config['pon_chr']
    # set config['chr'] to region if provided in EBscore
    region = args['region'] if 'region' in args.keys() else ''
    if region:
        config['chr'] == region.split(':')[0]
    if not set(config['chr']).issubset(set(config['pon_chr'])):
        sys.stderr.write('''
        Provided chromosome name is not found in bam files! Exiting..
        ''')
        sys.exit(1)
    # else config['chr'] is the one selected chromosome

    # ###########################################################
    # #################### --> GENERATE CACHE ###################
    if config['generate_cache']:

        # ##### set cache folder #########################
        if 'cache_folder' in args.keys():
            config['cache_folder'] = cache_folder = args['cache_folder']
        else:
            # if no cache folder is provided, pon_list folder is used
            config['cache_folder'] = os.path.join(os.path.dirname(pon_dict['list']), f"EBcache_{os.path.splitext(pon_dict['list'][0])}")
        if not os.path.isdir(cache_folder):
            os.mkdir(cache_folder)

        # ##### set bed file #############################
        # validate bed file and get the chromosomes needed for caching
        bed_file = args['bed_file'] if 'bed_file' in args.keys() else None
        if bed_file:
            config['bed_file'], config['bed_chr'] = utils.validate_bed(bed_file, config)
            valid_chrs = list(set(config['bed_chr']) & set(config['chr']))
            if len(valid_chrs):
                # load the valid chroms in the bed file into the active chroms
                config['chr'] = valid_chrs
            # if restricted chrom is not in bed file
            else:
                print(f'Chromosome {config["chr"]} is not in bed file. Nothing to do here.')
                return
        else:   # no bed_file
            # write empty
            if args['force_caching']:
                config['bed_file'], config['bed_chr'] = None, []
            else:
                sys.stderr.write('''
                    Please provide a bed_file or use -force caching option!\n
                    Generating a cache file for the entire genomic stretch covered by the PoN potentially takes forever!!
                    ''')
                sys.exit(1)
        success = generate_cache(pon_dict, config)
        print(success)
        return
    # ###########################################################
    # #################### EBscore ###################
    # implicit else:
    cache_folder = args['use_cache'] if 'use_cache' in args.keys() else None
    if cache_folder:
        # validate cache returns the cache_folder (same as input) and the chrom list of the..
        # ..existing cache files --> stored in config
        config['cache_folder'], config['cache_chr'] = utils.validate_cache(cache_folder, config)
        config['cache_mode'] = True
        print('Running EBscore in EBcache mode...')

    # get arguments for EBscore
    mut_file = utils.validate(args['mut_file'], "No target mutation file")
    if not config['sep'] in [',', '\t']:
        print(f'Separator \" {config["sep"]} \" cannot be used. Trying to open mutation file with separator \" \\t \"..')

    # check if tumor_bam and bai exists and whether it has the same chrom set as pon_file
    tumor_bam = utils.validate_bam(args['tumor_bam'])
    output_path = args['output_path']
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')

    # #################### EBscore ANNO ###################
    if is_anno:
        print(f'Loading annotation file {mut_file}')
        # create anno_df, store original other info and get anno_chr list
        anno_df, original_columns, config['anno_chr'] = utils.read_anno_csv(mut_file, config)
        # create small copy for working with
        mut_df = anno_df[anno_df.columns[:5]].copy()

        # ############## CACHE MODE ###############
        EB_dfs = []
        if config['cache_mode']:
            # ################ CACHE RETRIEVAL ####################################
            # do multi-threaded merging of the different AB_files per chromosome

            # get the SNPs for the annotation
            snp_df = mut_df.query('not (Ref == "-" or Alt == "-")')
            # get the valid chromosomes
            cache_chromes = list(set(config['chr']) & set(config['cache_chr']))
            # init multicore
            getEB_pool = Pool(threads)
            # threads are mapped to the pool of chromosomes
            print('Drawing BB parameters from cache and piling up the target bam..')
            mut_EBs = getEB_pool.map(partial(get_EB_from_cache, snp_df, tumor_bam, config), cache_chromes)
            # EBcached_outputs will be concatenated with the annoworker dfs
            getEB_pool.close()
            getEB_pool.join()
            # store EBscored dfs in EB_dfs
            EB_dfs += mut_EBs

        # ################ NOT CACHE MODE // INDELS IN CACHE MODE ###############
            # reduce mut_df to the remaining indel occurrences and uncached chroms for anno worker
            mut_df = mut_df.query(f'(Ref == "-" or Alt == "-") and Chr not in {config["cache_chr"]}')
        # if all are snps, mut_df is now empty after caching and would raise error in if clause
        mut_len = len(mut_df.index)
        if mut_len > 0:
            #  setup the processor pool
            anno_pool = Pool(threads)
            # set the minimum number of lines for one thread to 2000
            split_factor = min(math.ceil(mut_len / 2000), threads)

            # mut_split is the argument pool for anno_pool
            mut_split = np.array_split(mut_df, split_factor)
            print('Piling up target and PoN for EBscore computation')
            # run partial function anno_partial with mut_df as remaining argument to iterate over for multiprocessing
            out_dfs = anno_pool.map(partial(anno.worker, tumor_bam, pon_dict, output_path, region, config), mut_split)  # mut_split is the iterable df_pool
            anno_pool.close()
            anno_pool.join()
            # add the anno.worker output to EB_dfs
            EB_dfs += out_dfs
        # concat EB_dfs to one
        EB_df = pd.concat(EB_dfs)
        EB_df = EB_df.sort_values([EB_df.columns[0], EB_df.columns[1]])

        # add the original columns back to the annotation file
        final_df = pd.merge(left=anno_df, right=EB_df, how='outer', on=['Chr', 'Start'])
        if original_columns is not None:
            final_df.columns = list(original_columns) + ['EB_score']
        print(f'Writing annotation file {output_path} with EBscores to disc..')
        final_df.to_csv(output_path, sep='\t', index=False)

        # ############### DEBUG #####################################
        if config['debug_mode']:
            out_file = output_path.replace('eb', "eb_only")
            EB_df.to_csv(out_file, sep=config['sep'], index=False)
        # ############# ##############################################
        print('EBscore is finished!')

    else:   # is_anno is False
        print('EBFilter is not YET compatible with .vcf files. ')
