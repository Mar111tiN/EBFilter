from multiprocessing import Pool
from functools import partial
import subprocess
from subprocess import Popen, PIPE
from io import StringIO
import os
import pandas as pd
import numpy as np
from .utils import clean_up_df, cleanup_badQ, bam_to_chr_list, clean_read_column, sort_chr
from .eb import get_count_df_snp
from .beta_binomial import fit_beta_binomial
import re
import math

sign_re = re.compile(r'\^.|\$')
acgt = ['A','C','T','G']

############# PON2SPLITBAM ############################################
def split_bam(chromosome, pon_folder, pon_row):
    '''
    creates a sub bam file (+bai) of the specified chromosomal region per pon_list row
    returns the name of the bam file for creating a sub pon_list
    '''

    bam_file = pon_row[0]
    # create the name for the sub-bam file as org_bam_name_chr?.bam (and .bai)
    bam_out = os.path.join(pon_folder, f"{os.path.splitext(os.path.basename(bam_file))[0]}_{str(chromosome)}.bam")
    split_bam_cmd = ["samtools", "view", "-b", "-o", bam_out, bam_file, str(chromosome)]
    bam_index_cmd = ["samtools", "index", bam_out]
    subprocess.check_call(split_bam_cmd)
    subprocess.check_call(bam_index_cmd)
    return bam_out

############# PILEUP2BAM ############################################
def pon2pileup(pon_list, config, pon_folder, chromosome):
    '''
    create the dataframe of AB parameters per region per mismatch base
    '''

    # get the number of files in the pon_list
    pon_count = sum(1 for line in open(pon_list, 'r'))

    if chromosome != 'all_chromosomes': # multithreading
        # create a chromosome-bound bam and bai for each bam in the pon_list
        # write the created bams to a dataframe and output as pon_list
        print(f'Generating  pileup for chromosome {chromosome}..')
        pon_df = pd.read_csv(pon_list, header=None)
        pon_sub_df = pd.DataFrame()
        print(f'Splitting bam files for chromosome {chromosome}..')
        pon_sub_df['bam'] = pon_df.apply(partial(split_bam, chromosome, pon_folder), axis=1)
        # use pon_list_chr?.txt instead of the global pon_list.txt
        pon_sub_list = os.path.join(pon_folder, f"pon_list_{chromosome}.txt")
        # write the pon_list_{chr#} to file in output/pon for access by pon2pileup
        pon_sub_df.to_csv(pon_sub_list, header=None, index=False)
    else:
        print(f'Generating  pileup..')

    ############# PON2PILEUP ############################################
    # from the pon_list return a pileup_df
    with open(config['log'], 'w+') as log:
        # mpileup configs
        mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(config['q']), "-Q", str(config['Q']), "--ff", config['ff']]
        # mpileup PoN list
        mpileup_cmd += ["-b", pon_sub_list]
        # optional (but highly recommended:-) bed_file for -l option
        if config['bed_file']:
            mpileup_cmd += ["-l", config['bed_file']]

        pileup_stream = Popen(mpileup_cmd, stdout=PIPE, stderr=log)
        pileup_string = StringIO(pileup_stream.communicate()[0].decode('utf-8'))
        names = ['Chr', 'Start', 'Ref']
        for i in range(pon_count):
            names += [f"depth{i}", f"read{i}", f"Q{i}"]
    pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names)
    if len(pileup_df.index) == 0:
        print(f'Pileup for chromosome {chromosome} is empty and will be dropped..')
        return


    ############## CLEANUP PILEUP #######################################
    # get the count from the number of columns of pileup_df
    pon_count = int(pileup_df.shape[1] / 3 - 1)
    # get the read columns for the clean_columns apply
    read_columns = [f'read{i}' for i in range(pon_count)]
    # apply the cleaning function on the read rows of the pileup
    pileup_df[read_columns] = pileup_df[read_columns].apply(clean_read_column)

    ########## DEBUG ####################################
    if config['debug_mode']:
        # 
        pileup_file = 'cache'
        if config['threads'] > 1:
            pileup_file += f"_{chromosome}"
        pileup_file += '.pileup'
        pileup_folder = os.path.join(config['cache_folder'], 'cache_pileups')

        if not os.path.isdir(pileup_folder):
            os.mkdir(pileup_folder)
        pileup_path = os.path.join(pileup_folder, pileup_file)
        # write the pileup df to csv as <cache_folder>/ab_cache[_chr1].pileup
        pileup_df.to_csv(pileup_path, sep='\t', index=False)
    #####################################################
    return {'df': pileup_df, 'chr': chromosome}


def pileup2AB(config, chromosome, chr_len, pileup_df):
    '''
    creates the AB_df for a pileup
    '''

    ####################### AB FITTING ###############################################
    # create a copy of just the Chr and Start coords to store the AB values in
    AB_df = pileup_df.iloc[:,[0,1]].copy()

    def get_AB(penalty, length, start, chromosome, row):
        '''
        returns the AB parameters (A+a A+b A-a A-b G+a G+b....) for each pileup row
        '''
        bb_s = pd.Series()
        ########### Progress Bar ###############################################
        if (row.name - start) % 1000 == 0 and (row.name - start) > 0:
            bar_size = 25
            perc = math.floor((int(row.name) - start)/ length * bar_size)
            progress = '|' + '.' * math.floor(perc / 100 * bar_size) + ' ' * (bar_size - perc) + '|'
            print(f"P{os.getpid()}: {row.name - start} lines (total {perc}% of Chr {chromosome.replace('chr', '').replace('Chr', '')}) processed\t {progress}")

        ########### get count matrix ###########################################    
        for var in acgt:
            # get the count matrix
            count_df = get_count_df_snp(row, var, 4)
            # get the AB parameters for 
            bb_params = fit_beta_binomial(count_df, penalty)
        # dump the different parameters into bb_s
        # keys have to fit with the var_columns for the receiving AB_df
            bb_s[f'{var}+a'] = bb_params['p'][0]
            bb_s[f'{var}+b'] = bb_params['p'][1]
            bb_s[f'{var}-a'] = bb_params['n'][0]
            bb_s[f'{var}-b'] = bb_params['n'][1]
        return bb_s

    ################## Store AB data into df ####################################   
    # create the columns (A+a A+b A-a A-b G+a G+b....) for the recipient df of the pileup_df apply function 
    var_columns = [f'{var}{strand}{param}' for var in acgt for strand in ['+', '-'] for param in ['a','b']]
    pileup_length = len(pileup_df.index)
    pileup_start = pileup_df.iloc[0].name
    print(f'Process {os.getpid()}: Computing ABs for Chr {chromosome} \t(lines {pileup_start }\t to \t{pileup_length + pileup_start})')
    AB_df[var_columns] = pileup_df.apply(partial(get_AB, config['fitting_penalty'], pileup_length, pileup_start, chromosome), axis=1)
    print(f'Process {os.getpid()}: Computing ABs for chromosome {chromosome} \t(lines {pileup_start }\t to \t{pileup_length + pileup_start}) finished.')

    ###################### DEBUG #######################################################
    if chromosome != 'all_chromosomes' and config['debug_mode']: 
        # for multithreading also output the sub_files in debug_mode
        chr_cache = os.path.join(config['cache_folder'], f"{chromosome}_{os.getpid()}.cache")
        i = 1
        while os.path.isfile(chr_cache):
            chr_cache = os.path.join(config['cache_folder'], f"{chromosome}_{os.getpid()}-{i}.cache")
            i += 1
    ####################################################################################

    ######################## OUTPUT ####################################################
    else:
        chr_cache = os.path.join(config['cache_folder'], f"{chromosome}.cache")
    AB_df.to_csv(chr_cache, sep=',', index=False)
    return AB_df

def generate_cache(pon_list, config):
    '''
    create a cache file for ab-parameters
    '''

    print('Generating Cache...')
    threads = config['threads']
    pon_df = pd.read_csv(pon_list)  
    if threads == 1:
        pileup_dict = pon2pileup(pon_list, config, 'all_chromosomes')
        AB_df = pileup2AB(config, 'all_chromosomes', pileup_dict['df'])
    else: # multithreading
        ####################### SET TO FINAL OUTPUT ########################
        # make directory for temporary bam files
        pon_folder = os.path.join(os.path.dirname(config['cache_folder']), 'pon')
        if not os.path.isdir(pon_folder):
            os.mkdir(pon_folder)

        ######################### PON2PILEUP ###################################
        # get the list of valid chromosomes from the first bam in the pon_list
        chromosomes = bam_to_chr_list(pon_df.iloc[0,0])
        cache_pool = Pool(threads)
        # threads are mapped to the pool of chromosomes
        # returns a dict {'df':pileup_df, 'chr': 'chr1'}
        pileup_dicts = cache_pool.map(partial(pon2pileup, pon_list, config, pon_folder), chromosomes)

        # remove Nones and sort for chromosome name
        pileup_dicts = sorted(filter(None, pileup_dicts), key=sort_chr)
        cache_pool.close()
        cache_pool.join()

        ######################### PILEUP2AB ###################################
        
        AB_dfs = []
        # split the pileups for each Chr into threaded chunks to even out different chromosome sizes
        # go through each chromosome with required number of threads
        # create the job pool for the threads
        for pileup_dict in pileup_dicts:  # account for empty pileups with filter(None..
            chromosome = pileup_dict['chr']
            chr_len = len(pileup_dict['df'].index)      # get length for progress info
            # set the minimum number of lines for one thread to 10000
            split_factor = min(math.ceil(chr_len / 2000), threads)
            # split the arrays into litte fractions for computation
            pileup_split = np.array_split(pileup_dict['df'], split_factor)
            pileup2AB_pool = Pool(threads)
            AB_chr_dfs = pileup2AB_pool.map(partial(pileup2AB, config, chromosome, chr_len), pileup_split)
            pileup2AB_pool.close()
            pileup2AB_pool.join()
            # concatenate the AB_dfs for each chromosome to AB_chr_df
            AB_chr_df = pd.concat(AB_chr_dfs)
            AB_dfs.append(AB_chr_df.sort_values([AB_chr_df.columns[0], AB_chr_df.columns[1]]))

            ############### OUTPUT ###########################################################
            chr_cache = os.path.join(config['cache_folder'], f"{chromosome}.cache")
            print(f"Writing ABcache for Chr {chromosome} to file {chr_cache}.")
            AB_chr_df.to_csv(chr_cache, sep=',', index=False)

        AB_df = pd.concat(AB_dfs)
        AB_df = AB_df.sort_values([AB_df.columns[0], AB_df.columns[1]])
        # AB_df = out_df.sort_values([out_df.columns[0], out_df.columns[1]])

    ############# DEBUG #####################################################################
    # remove the pon_list and all bam files (by removing the whole pon folder)
    if not config['debug_mode'] and threads > 1:    
        subprocess.check_call(['rm', '-r', pon_folder])
    #########################################################################################   

    all_cache = os.path.join(config['cache_folder'], "all.cache")
    print(f"Writing final ABcache for covered regions to file {all_cache}.")
    AB_df.to_csv(all_cache, sep=',', index=False)

    return 'Generation of AB file finished.'
