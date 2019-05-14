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
from .cache import generate_cache
from .ebscore import get_EB_from_cache
from datetime import datetime as dt


def main(args, config):
    '''
    validates files and refers to respective functions
    '''

    # ############## ARGUMENTS and CONFIG ##################################
    config['cache_mode'] = False
    config['generate_cache'] = args['generate_cache'] if 'generate_cache' in args.keys() else False
    # store list and pon_df in pon_dict, store pon chroms in config['pon_chr']
    pon_dict = utils.validate_pon(args['pon_list'], config)
    threads = config['threads']
    debug_mode = config['debug_mode']

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

    # ######################################################################
    # #################### --> GENERATE CACHE ##############################
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

    # #########################################################################
    # #################### EBscore ############################################
    # implicit else:
    now = dt.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"{now} Starting EBscore for file {os.path.basename(args['mut_file'])}!")
    if cache_folder:
        cache_folder = args['use_cache'] if 'use_cache' in args.keys() else None
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
    config['output_path'] = args['output_path']
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')

    # #################### EBscore ANNO ###################
    if is_anno:
        success = generate_EBscore(mut_file, config)
        print(success)
        return

    else:   # is_anno is False
        print('EBFilter is not YET compatible with .vcf files. ')
        return
