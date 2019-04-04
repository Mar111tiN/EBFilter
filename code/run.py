#! /usr/bin/env python
import pysam
import sys
import re
import os
import pandas as pd
import numpy as np
import subprocess
import multiprocessing
from multiprocessing import Pool
from functools import partial

from .utils import validate
from . import anno
from . import vcf



def main(args, state):
    '''
    validates files and refers to respective functions
    '''

    ############### ARGUMENTS #######################
    mut_file = args['mut_file']
    tumor_bam = args['tumor_bam']
    pon_list = args['pon_list']
    output_path = args['output_path']
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')
    region = args['region']
    threads = state['threads']
    debug_mode = state['debug_mode']
    sep = state['sep']
    # create log directory - remove in snakemake
    log_folder = os.path.split(state['log'])[0]
    if not os.path.exists(log_folder) or os.path.isfile(log_folder):
        os.makedirs(log_folder)
    
    ############### ARGUMENTS #######################
    # file existence check for files and bams in pon_list
    validate(mut_file, tumor_bam, pon_list) 

    if is_anno:
        # create dataframe (maybe more options for other input file formats)
        anno_df = pd.read_csv(args['mut_file'], sep=sep).sort_values(['Chr', 'Start'])
        # create small copy for working with
        mut_df = anno_df[anno_df.columns[:5]].copy()

        if threads == 1:
        # non multi-threading mode
            out_df = anno.worker(tumor_bam, pon_list, output_path, region, state, mut_df) # -1 means single-threaded
            if not state['debug_mode']:
                subprocess.check_call(["rm", output_path + '.target.pileup'])
                subprocess.check_call(["rm", output_path + '.control.pileup'])
                subprocess.check_call(["rm", "-f", f"{mut_file}.region_list.bed"])
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

            # jobs = []
            # out_dfs = [0] * threads # out_dfs as placeholder for the sub_dataframes
            # for i in range(threads):
            #     worker_args = (mut_split[i], tumor_bam, pon_list, out_dfs[i], region, state, i)
            #     process = multiprocessing.Process(target=anno.worker, args=worker_args)                    
            #     jobs.append(process)
            #     process.start()       
            # # wait for jobs to finish
            # for i in range(threads):
            #     jobs[i].join()      
            # # merge the individual results

            # anno.merge(output_path, threads)      
            # # delete intermediate files
            # if not debug_mode:
            #     for i in range(threads):
            #         rm_list = [f"{output_path}.{i}"]
            #         rm_list += [f"{output_path}.sub.{i}"]
            #         rm_list += [f"{output_path}.sub.{i}.control.pileup", f"{output_path}.sub.{i}.target.pileup"]
            #         rm_list += [f"{output_path}.{i}.region_list.bed"]
            #         subprocess.check_call(["rm", *rm_list])


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
