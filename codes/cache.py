from multiprocessing import Pool
from functools import partial
import subprocess
from subprocess import Popen, PIPE
from io import StringIO
import os
import pandas as pd
import numpy as np
from .anno import anno2pileup
from . import utils
from .eb import get_count_df_snp
from .beta_binomial import fit_bb, bb_pvalues, fisher_combination
import re
import math
from datetime import datetime as dt

sign_re = re.compile(r'\^.|\$')
acgt = ['A', 'C', 'T', 'G']


# ############ PILEUP2BAM ############################################
def pon2pileup(pon_dict, config, chromosome):
    '''
    create the dataframe of AB parameters per region per mismatch base
    '''

    # get the number of files in the pon_list
    pon_count = len(pon_dict['df'].index)
    pon_folder = config['pon_folder']
    # create a chromosome-bound bam and bai for each bam in the pon_list
    # write the created bams to a dataframe and output as pon_list
    pon_sub_df = pd.DataFrame()
    print(f"{dt.now().strftime('%H:%M:%S')} Splitting bam files of PoN for chromosome {chromosome}..")
    pon_sub_df['bam'] = pon_dict['df'].apply(partial(utils.split_bam, chromosome, pon_folder), axis=1)
    # use pon_list_chr?.txt instead of the global pon_list.txt
    pon_sub_list = os.path.join(pon_folder, f"pon_list_{chromosome}.txt")
    print(f"{dt.now().strftime('%H:%M:%S')} Writing pon list for chromosome {chromosome}..")
    # write the pon_list_{chr#} to file in output/pon for access by pon2pileup
    pon_sub_df.to_csv(pon_sub_list, header=None, index=False)

    # ############ PON2PILEUP ############################################
    # from the pon_list return a pileup_df
    # mpileup configs
    mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(config['q']), "-Q", str(config['Q']), "--ff", config['ff']]
    # mpileup PoN list
    mpileup_cmd += ["-b", pon_sub_list]
    # optional (but highly recommended:-) bed_file for -l option
    if config['bed_file']:
        mpileup_cmd += ["-l", config['bed_file']]

    # pileup_file = os.path.join(config['pileup_folder'], f"cache_{chromosome}.pileup")
    # mpileup_cmd += ['-o', pileup_file]
    print(f"{dt.now().strftime('%H:%M:%S')} Generating pileup for chromosome {chromosome}..")

    utils.show_command(mpileup_cmd, config, multi=False)
    pileup_stream = Popen(mpileup_cmd, stdout=PIPE)
    # !!!!!!!!!! MEMORY ERROR FOR LARGE FILES !!!!!!!!!!!!!!!!!!!!! --> provide more memory or write to file
    pileup_file = StringIO(pileup_stream.communicate()[0].decode('utf-8'))
    ###############################################################################

    names = ['Chr', 'Start', 'Ref']
    for i in range(pon_count):
        names += [f"depth{i}", f"read{i}", f"Q{i}"]
    # read StringIO object from stream into pileup_df in memory
    pileup_df = pd.read_csv(pileup_file, sep='\t', header=None, names=names)
    pileup_len = len(pileup_df.index)
    if pileup_len == 0:
        print(f"Pileup for chromosome {chromosome} is empty and will be dropped..")
        return {'file': 'empty', 'chr': chromosome, 'pileup_len': 0}

    # ############# CLEANUP PILEUP #######################################
    # get the count from the number of columns of pileup_df
    pon_count = int(pileup_df.shape[1] / 3 - 1)
    # get the read columns for the clean_columns apply
    read_columns = [f'read{i}' for i in range(pon_count)]
    # apply the cleaning function on the read rows of the pileup
    pileup_df[read_columns] = pileup_df[read_columns].apply(utils.clean_read_column)

    # ############# WRITE TO FILE #############################################
    # write the pileup df to csv as <cache_folder>/ab_cache[_chr1].pileup

    # fragment the pileup at this stage using numpy array to only get chunks into memory
    pileup_file = os.path.join(config['pileup_folder'], f"cache_{chromosome}.pileup")
    print(f"{dt.now().strftime('%H:%M:%S')}: Writing pileup of Chr {chromosome} to file {pileup_file}")
    pileup_df.to_csv(pileup_file, sep='\t', index=False)

    # ############################ DEBUG ######################################
    if not config['debug_mode']:
        # delete the split_bam pon list
        subprocess.check_call(['rm', pon_sub_list])
        pon_sub_df['bam'].apply(utils.delete_pom_bams)
    ##########################################################################

    pileup_file_dict = {'file': pileup_file, 'chr': chromosome, 'pileup_len': pileup_len}
    return pileup_file_dict


def pileup2AB(config, chromosome, chr_len, pileup_df):
    '''
    creates the AB_df for a pileup
    '''

    # ###################### AB FITTING ###############################################
    # create a copy of just the Chr and Start coords to store the AB values in
    AB_df = pileup_df.iloc[:, [0, 1]].copy()

    def get_AB(penalty, row):
        '''
        returns the AB parameters (A+a A+b A-a A-b G+a G+b....) for each pileup row
        main computational load
        '''

        bb_s = pd.Series()
        # ########## get count matrix ###########################################
        for var in acgt:
            # get the count matrix
            count_df = get_count_df_snp(row, var, 4)
            # get the AB parameters for 
            bb_params = fit_bb(count_df, penalty)
        # dump the different parameters into bb_s
        # keys have to fit with the var_columns for the receiving AB_df
            bb_s[f'{var}+a'] = bb_params['p'][0]
            bb_s[f'{var}+b'] = bb_params['p'][1]
            bb_s[f'{var}-a'] = bb_params['n'][0]
            bb_s[f'{var}-b'] = bb_params['n'][1]
        return bb_s

    # ################# Store AB data into df ####################################
    # create the columns (A+a A+b A-a A-b G+a G+b....) for the recipient df of the pileup_df apply function
    var_columns = [f'{var}{strand}{param}' for var in acgt for strand in ['+', '-'] for param in ['a', 'b']]

    pileup_start = pileup_df.iloc[0].name
    frac = pileup_start / chr_len
    bar_size = 25
    progress_bar = '|' + '.' * math.ceil(frac * bar_size) + ' ' * math.floor(bar_size * (1 - frac)) + '|'

    print(f'Process {os.getpid()}: {pileup_start } lines ({round(frac * 100, 1)}%) of Chr {chromosome}\t{progress_bar}')

    AB_df[var_columns] = pileup_df.apply(partial(get_AB, config['fitting_penalty']), axis=1)

    # ##################### DEBUG #######################################################
    if config['debug_mode']:
        # for multithreading also output the sub_files in debug_mode
        chr_cache = os.path.join(config['cache_folder'], f"{chromosome}_{os.getpid()}.cache")
        i = 1
        while os.path.isfile(chr_cache):
            chr_cache = os.path.join(config['cache_folder'], f"{chromosome}_{os.getpid()}-{i}.cache")
            i += 1
    ####################################################################################

    # ####################### OUTPUT ####################################################
    else:
        chr_cache = os.path.join(config['cache_folder'], f"{chromosome}.cache")
    AB_df.to_csv(chr_cache, sep=',', index=False)
    return AB_df


def generate_cache(pon_dict, config):
    '''
    create a cache file for ab-parameters
    '''

    # #####################INIT + VALIDATION #########################
    # get the list of desired chroms without existing cache file from config['chr']
    config['chr'] = utils.check_cache_files(config)
    if len(config['chr']):
        now = dt.now().strftime('%Y-%m-%d %H:%M:%S')
        suf = 's'
        if len(config['chr']) == 1:
            suf = ''
        print(f"{now} Generating Cache for chromosome{suf} {' '.join(config['chr'])}.. ")
    else:
        return "Everything is there. No need for computation."

    threads = config['threads']

    # ######################## PON2PILEUP ###################################
    config['pileup_folder'] = pileup_folder = os.path.join(config['cache_folder'], 'cache_pileups')
    config['pon_folder'] = pon_folder = os.path.join(config['cache_folder'], 'pon')

    # get the list of chroms without existing pileup_files from config['chr'] 
    # ...and store the existing pileups in pileup_file_dicts

    pileup_chrs, pileup_file_dicts = utils.check_pileup_files(config)

    # if still some pileups remain
    if len(pileup_chrs):
        # make directory for temporary bam and pileup files(only if pileups are missing)
        if not os.path.isdir(pon_folder):
            os.mkdir(pon_folder)
        if not os.path.isdir(pileup_folder):
            os.mkdir(pileup_folder)

        # in order to use the threads, no multiprocessing is used here
        # otherwise, memory peaks would kill the process
        for pileup_chr in pileup_chrs:
            pileup_file_dicts.append(pon2pileup(pon_dict, config, pileup_chr))

    # reading the pileup files into df generators
    pileup_dicts = []
    for pileup_file_dict in pileup_file_dicts:
        # pileup df is loaded into memory as generator of chunksize 10000 for reducing memory peak
        pileup_dfgen = pd.read_csv(pileup_file_dict['file'], sep='\t', chunksize=10000)

        print(f"{dt.now().strftime('%H:%M:%S')}: Reading {pileup_file_dict['pileup_len']} lines of Chr {pileup_file_dict['chr']} pileup for AB computation")

        pileup_dict = {'df': pileup_dfgen, 'chr': pileup_file_dict['chr'], 'pileup_len': pileup_file_dict['pileup_len']}
        pileup_dicts.append(pileup_dict)

        # ######### DEBUG ####################################
        if not config['debug_mode']:
            subprocess.check_call(['rm', '-f', pileup_file_dict['file']])
        #####################################################

    # account for empty pileups with filter (None..
    pileup_dicts = sorted(filter(None, pileup_dicts), key=utils.sort_chr)

    # ######################## PILEUP2AB ###################################
    # split the pileups for each Chr into threaded chunks to even out different chromosome sizes
    # go through each chromosome with required number of threads
    # create the job pool for the threads
    for pileup_dict in pileup_dicts:
        chromosome = pileup_dict['chr']
        # why not working every round
        # processed = 0
        chr_cache = os.path.join(config['cache_folder'], f"{chromosome}.cache")
        # accomodate for empty pileups (should not be neccessary because of validation)
        if not pileup_dict['pileup_len']:
            print(f"Writing empty cache for Chr {chromosome} to file {chr_cache}.")
            open(chr_cache, 'a').close()
            continue
        chr_len = pileup_dict['pileup_len']  # get length for progress info
        pileup2AB_pool = Pool(threads)
        # pass in the chunked-up df generator to the pool
        AB_chr_dfs = pileup2AB_pool.imap(partial(pileup2AB, config, chromosome, chr_len), pileup_dict['df'])
        pileup2AB_pool.close()
        pileup2AB_pool.join()
        # concatenate the AB_dfs for each chromosome to AB_chr_df
        AB_chr_df = pd.concat(AB_chr_dfs)
        AB_chr_df = AB_chr_df.sort_values([AB_chr_df.columns[0], AB_chr_df.columns[1]])

        # ############## OUTPUT ###########################################################
        chr_cache = os.path.join(config['cache_folder'], f"{chromosome}.cache")
        print(f"{dt.now().strftime('%H:%M:%S')}: Writing ABcache for Chr {chromosome} to file {chr_cache}.")
        AB_chr_df.to_csv(chr_cache, sep=',', index=False, compression='gzip')

    #########################################################################################

    return f"{dt.now().strftime('%H:%M:%S')}: Generation of AB file finished."
