#! /usr/bin/env python
import sys
import re
import os
import pandas as pd
import numpy as np
import subprocess
import multiprocessing
from multiprocessing import Pool
from functools import partial

from .utils import validate, read_anno_csv, validate_bam, validate_pon, validate_cache
from . import anno
from . import vcf
from .cache import generate_cache, get_EB_from_cache



def main(args, config):
    '''
    validates files and refers to respective functions
    '''
    

    ############### ARGUMENTS and CONFIG ##################################
    config['cache_mode'] = False
    generate_cache_mode = args['generate_cache'] if 'generate_cache' in args.keys() else False

    threads = config['threads']
    debug_mode = config['debug_mode']
    pon_dict = validate_pon(args['pon_list'])  


    ##################### --> GENERATE CACHE ###################
    if generate_cache_mode:

        ###### set cache folder #########################    
        if 'cache_folder' in args.keys():
            config['cache_folder'] = cache_folder = args['cache_folder']
        else:
        # if no cache folder is provided, pon_list folder is used
            config['cache_folder'] = os.path.join(os.path.dirname(pon_dict['list']), '.cache')
        if not os.path.isdir(cache_folder):
            os.mkdir(cache_folder)

        ###### set chromosome ########################### 
        config['chr'] = args['chromosome'] if 'chromosome' in args.keys() else 'all'
        if not config['chr'] in (pon_dict['chr_list'] + ['all']):
            sys.stderr.write('''
                    Provided chromosome name is not found in bam files! Exiting..
                    ''')
            sys.exit(1)

        ###### set bed file #############################
        config['bed_file'] = args['bed_file'] if 'bed_file' in args.keys() else None
        if not config['bed_file']:
            if not args['force_caching']:
                sys.stderr.write('''
                    Please provide a bed_file or use -force caching option!\n
                    Generating a cache file for the entire genomic stretch covered by the PoN potentially takes forever!!
                    ''')
                sys.exit(1)
        success = generate_cache(pon_dict, config)
        print(success)
        return

    ##################### EBscore ###################
    else: # EBscore mode
        cache_folder = args['use_cache'] if 'use_cache' in args.keys() else None
        if cache_folder:
            config['cache_folder'] = validate_cache(cache_folder, pon_dict)
            config['cache_mode'] = True
            print('Running EBscore in Cache mode...')

    # returns a pon dict: {'df': pon_df, 'list': pon_list}
          
        
    # get arguments for EBscore 
    sep = config['sep']
    mut_file = validate(args['mut_file'], "No target mutation file")
    tumor_bam = validate_bam(args['tumor_bam'])
    output_path = args['output_path']   
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')
    region = args['region'] if 'region' in args.keys() else ''
    
    ##################### EBscore ANNO ###################   
    if is_anno:
        # create dataframe (maybe more options for other input file formats)
        anno_df, original_columns = read_anno_csv(mut_file, config)
        # create small copy for working with
        mut_df = anno_df[anno_df.columns[:5]].copy()

        ############### CACHE MODE ###############
        EBcached_outputs = []
        if config['cache_mode']:
            # do multi-threaded merging of the different AB_files per chromosome
            # get the SNPs for the annotation
            snp_df = mut_df.query('not (Ref == "-" or Alt == "-")')

            # reduce mut_df to the remaining indel occurrences
            mut_df = mut_df.query('Ref == "-" or Alt == "-"')
            # get the valid chromosomes
            chromosomes = pon_dict['chr_list']
            # init multicore
            getEB_pool = Pool(threads)
            # threads are mapped to the pool of chromosomes
            mut_EBs = getEB_pool.map(partial(get_EB_from_cache, snp_df, tumor_bam, config), chromosomes)
            # EBcached_outputs will be concatenated with the annoworker dfs
            EBcached_outputs = mut_EBs
        # if all are snps, mut_df is now empty after caching and will raise error
        if len(mut_df.index) > 0:
            #  setup the processor pool
            anno_pool = Pool(threads)
            # mut_split is the argument pool for anno_pool
            mut_split = np.array_split(mut_df, threads)
            # run partial function anno_partial with mut_df as remaining argument to iterate over for multiprocessing
            out_dfs = anno_pool.map(partial(anno.worker, tumor_bam, pon_dict, output_path, region, config), mut_split)  # mut_split is the iterable df_pool
            anno_pool.close()
            anno_pool.join()
            # add the output from EBcached to the other output
            out_dfs += EBcached_outputs
        else:
            out_dfs = EBcached_outputs
        out_df = pd.concat(out_dfs)
        out_df = out_df.sort_values([out_df.columns[0], out_df.columns[1]])

        # add the original columns back to the annotation file
        final_df = pd.merge(left=anno_df, right=out_df, how='outer', on=['Chr', 'Start'])
        if original_columns is not None:
            final_df.columns = list(original_columns) + ['EB_score']
        final_df.to_csv(output_path, sep='\t', index=False)

        ################ DEBUG #####################################
        if config['debug_mode']:
            out_file = output_path.replace('eb', "eb_only")
            out_df.to_csv(out_file, sep=sep, index=False)
        ############################################################     
        
    else: 
        print('EBFilter is not YET compatible with .vcf files. ')
