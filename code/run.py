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

region_exp = re.compile(r'^([^ \t\n\r\f\v,]+):(\d+)\-(\d+)')
indel_re = re.compile(r'([\+\-])([0-9]+)([ACGTNacgtn]+)') # +23ATTTNNGC or -34TTCCAAG
sign_re = re.compile(r'\^\]|\$')


def main(args,_):
    '''
    validates files and refers to respective functions
    '''

    # should add validity check for os.path.exists(os.path.splitext(tumor_bam)[0] + '.bai')
    mut_file = args['mut_file']
    tumor_bam = args['tumor_bam']
    pon_list = args['pon_list']
    output_path = args['output_path']
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')
    region = args['region']
    # create log directory - remove in snakemake
    log_folder = os.path.split(_['log'])[0]
    if not os.path.exists(log_folder) or os.path.isfile(log_folder):
        os.makedirs(log_folder)
    

    # file existence check for files and bams in pon_list
    validate(mut_file, tumor_bam, pon_list) 
    if _['threads'] == 1:
        # non multi-threading mode
        if is_anno:
            anno.worker(mut_file, tumor_bam, pon_list, output_path, region,_)
        else: 
            vcf.worker(mut_file, tumor_bam, pon_list, output_path, region,_)
    else:
        # multi-threading mode
        ##########
        if is_anno:
            # partition anno files
            anno.partition(mut_file, output_path, threads)
            jobs = []
            for i in range(threads):
                worker_args = (f"{output_path}.{i}", tumor_bam, pon_list, f"{output_path}.sub.{i}", region, _)
                process = multiprocessing.Process(target=EBFilter_worker_anno, args=worker_args)                    
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
                    print('delete')
                    subprocess.check_call(["rm", f"{output_path}.{i}", f"{output_path}.{i}.control.pileup", f"{output_path}.{i}.target.pileup"])

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
