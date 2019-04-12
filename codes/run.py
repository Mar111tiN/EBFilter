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
from .cache import generate_cache



def main(args, config):
    '''
    validates files and refers to respective functions
    '''
    

    ############### ARGUMENTS and CONFIG ##################################
    config['cache_mode'] = False
    generate_cache = args['generate_cache'] if 'generate_cache' in args.keys() else False
    config['cache_folder'] = cache_folder = args['use_cache'] if 'use_cache' in args.keys() else None
    config['bed_file'] = args['bed_file'] if 'bed_files' in args.keys() else None
    threads = config['threads']
    debug_mode = config['debug_mode']


    ##################### --> GENERATE CACHE ###################
    if args['generate_cache']:
        ###### set cache folder #########################
        
        if not os.path.isdir(cache_folder):
            os.mkdir(cache_folder)

        ###### set bed file #########################
        if not args['bed_file']:
            if not args['force_caching']:
                sys.stderr.write('''
                    Please provide a bed_file or use -force caching option!\n
                    Generating a cache file for the entire genomic stretch covered by the PoN potentially takes forever!!
                    ''')
                sys.exit(1)
        pon_list = validate_pon(args['pon_list'])
        success = generate_cache(pon_list, config)
        print(success)
        return

    else: # EBscore mode
        if config['cache_folder']:
            cache_folder = validate_cache(config['cache_folder'])
            config['cache_mode'] = True

    pon_list = validate_pon(args['pon_list'])            
        
    # get arguments for EBscore 
    sep = config['sep']
    mut_file = validate(args['mut_file'], "No target mutation file")
    print('mutfile: ', mut_file )
    tumor_bam = validate_bam(args['tumor_bam'])
    output_path = args['output_path']   
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')
    region = args['region'] if 'region' in args.keys() else ''

    # create log directory - remove in snakemake
    log_folder = os.path.split(config['log'])[0]
    if not os.path.exists(log_folder) or os.path.isfile(log_folder):
        os.makedirs(log_folder)
    

    if is_anno:
        # create dataframe (maybe more options for other input file formats)
        anno_df, original_columns = read_anno_csv(mut_file, config)
        # create small copy for working with
        mut_df = anno_df[anno_df.columns[:5]].copy()

        if threads == 1:
        # non multi-threading mode
            if config['cache_mode']:
                # do multi-threaded merging of the different AB_files per chromosome
                #AB_columns = pd.MultiIndex.from_product([['A','C','T','G'],['+', '-'],['a','b']], names=['var', 'strand', 'param'])
                # read in the AB_df for the different chromosomes or total
                #AB_df = pd.read_csv('')

                out_df = anno.worker(tumor_bam, pon_list, output_path, region, config, mut_df) # -1 means single-threaded


        else: # multi-threading mode
            if config['cache_mode']:
                # do multi-threaded merging of the different AB_files per chromosome
                AB_columns = pd.MultiIndex.from_product([['A','C','T','G'],['+', '-'],['a','b']], names=['var', 'strand', 'param'])
                # read in the AB_df for the different chromosomes or total
                AB_df = pd.read_csv('')

            #  setup the processor pool
            anno_pool = Pool(threads)
            # mut_split is the argument pool for anno_pool
            mut_split = np.array_split(mut_df, threads)
            # run partial function anno_partial with mut_df as remaining argument to iterate over for multiprocessing
            out_dfs = anno_pool.map(partial(anno.worker, tumor_bam, pon_list, output_path, region, config), mut_split)  # mut_split is the iterable df_pool
            anno_pool.close()
            anno_pool.join()
            out_df = pd.concat(out_dfs)
            out_df = out_df.sort_values([out_df.columns[0], out_df.columns[1]])

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
        if threads == 1:
            vcf.worker(mut_file, tumor_bam, pon_list, output_path, region, config)
            # delete intermediate files
        else:
            # partition vcf files
            vcf.partition(mut_file, f"{output_path}.sub.vcf", threads)

            jobs = []
            for i in range(threads):
                process = multiprocessing.Process(target = EBFilter_worker_vcf, args = \
                    (output_path + ".sub.vcf." + str(i), tumor_bam, pon_list, output_path + "." + str(i), region))
                jobs.append(process)
                process.start()

            # wait all the jobs to be done
            for i in range(threads):
                jobs[i].join()

            # merge the individual results
            vcf.merge(output_path + ".", output_path, threads)

            # delete intermediate files
            if debug_mode == False:
                for i in range(threads):
                    subprocess.check_call(["rm", output_path + ".sub.vcf." + str(i)])
                    subprocess.check_call(["rm", output_path + "." + str(i)])
