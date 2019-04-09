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

from .utils import validate, read_anno_csv, validate_bam, validate_pon
from . import anno
from . import vcf
from .cache import generate_cache



def main(args, state):
    '''
    validates files and refers to respective functions
    '''

    state['cache_mode'] = False
    ############### ARGUMENTS #######################
    threads = state['threads']
    debug_mode = state['debug_mode']

    # generate cache
    if args['generate_cache']:
        if 'cache_path' in args.keys():
            cache_file = args['cache_path']
        else:
            # if no path to cache file is given, it will be generated at pon_list destination
            cache_file = os.path.join(os.path.splitext(args['pon_list'])[0], '.ABcache')
        state['cache_dir'] = os.path.dirname(cache_file)
        pon_list = validate_pon(args['pon_list'])
        return generate_cache(pon_list, state)
    else: # EBscore mode
        if 'cache_path' in args.keys():
            cache_file = validate(args['cache_path'], "No ABcache file found")
            state['cache_mode'] = True
        else:
            pon_list = validate_pon(args['pon_list'])            
        
    # get arguments for EBscore 
    sep = state['sep']
    mut_file = validate(args['mut_file'], "No target mutation file")
    tumor_bam = validate_bam(args['tumor_bam'])
    output_path = args['output_path']   
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')
    region = args['region']

    # create log directory - remove in snakemake
    log_folder = os.path.split(state['log'])[0]
    if not os.path.exists(log_folder) or os.path.isfile(log_folder):
        os.makedirs(log_folder)
    

    if is_anno:
        # create dataframe (maybe more options for other input file formats)
        anno_df, original_columns = read_anno_csv(mut_file, state)
        # create small copy for working with
        mut_df = anno_df[anno_df.columns[:5]].copy()

        if threads == 1:
        # non multi-threading mode
            out_df = anno.worker(tumor_bam, pon_list, output_path, region, state, mut_df) # -1 means single-threaded

        else: # multi-threading mode
            mut_split = np.array_split(mut_df, threads)
            # create partial function anno_partial with mut_df as remaining argument to iterate over for multiprocessing
            anno_partial = partial(anno.worker, tumor_bam, pon_list, output_path, region, state)
            anno_pool = Pool(threads)
            out_dfs = anno_pool.map(anno_partial, mut_split)  # mut_split is the iterable df_pool
            anno_pool.close()
            anno_pool.join()
            out_df = pd.concat(out_dfs)
            out_df = out_df.sort_values([out_df.columns[0], out_df.columns[1]])

        final_df = pd.merge(left=anno_df, right=out_df, how='outer', on=['Chr', 'Start'])
        if original_columns is not None:
            final_df.columns = list(original_columns) + ['EB_score']
        final_df.to_csv(output_path, sep='\t', index=False)

        ################ DEBUG #####################################
        if state['debug_mode']:
            out_file = output_path.replace('eb', "eb_only")
            out_df.to_csv(out_file, sep=sep, index=False)
        ############################################################

        
        
    else: 
        if threads == 1:
            vcf.worker(mut_file, tumor_bam, pon_list, output_path, region,state)
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
