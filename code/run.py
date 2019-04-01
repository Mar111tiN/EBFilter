#! /usr/bin/env python
import pysam
import sys
import re
import os
import pandas as pd
import subprocess
import multiprocessing

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
    # create log directory - remove in snakemake
    log_folder = os.path.split(state['log'])[0]
    if not os.path.exists(log_folder) or os.path.isfile(log_folder):
        os.makedirs(log_folder)
    

    # file existence check for files and bams in pon_list
    validate(mut_file, tumor_bam, pon_list) 
    if threads == 1:
        # non multi-threading mode
        if is_anno:
            anno.worker(mut_file, tumor_bam, pon_list, output_path, region,state)
        else: 
            vcf.worker(mut_file, tumor_bam, pon_list, output_path, region,state)
            # delete intermediate files
        if not state['debug_mode']:
            subprocess.check_call(["rm", output_path + '.target.pileup'])
            subprocess.check_call(["rm", output_path + '.control.pileup'])
            subprocess.check_call(["rm", "-f", f"{mut_file}.region_list.bed"])
    else:
        # multi-threading mode
        ##########
        if is_anno:
            # partition anno files
            anno.partition(mut_file, output_path, threads)
            jobs = []
            for i in range(threads):
                worker_args = (f"{output_path}.{i}", tumor_bam, pon_list, f"{output_path}.sub.{i}", region, state)
                process = multiprocessing.Process(target=anno.worker, args=worker_args)                    
                jobs.append(process)
                process.start()       
            # wait all the jobs to be done
            for i in range(threads):
                jobs[i].join()      
            # merge the individual results
            anno.merge(output_path, threads)      
            # delete intermediate files
            if not debug_mode:
                for i in range(threads):
                    rm_list = [f"{output_path}.{i}"]
                    rm_list += [f"{output_path}.sub.{i}"]
                    rm_list += [f"{output_path}.sub.{i}.control.pileup", f"{output_path}.sub.{i}.target.pileup"]
                    rm_list += [f"{output_path}.{i}.region_list.bed"]
                    subprocess.check_call(["rm", *rm_list])

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
