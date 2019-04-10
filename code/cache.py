from multiprocessing import Pool
from functools import partial
import subprocess
from subprocess import Popen, PIPE
from io import StringIO
import os
import pandas as pd
import numpy as np
from .utils import clean_up_df, cleanup_badQ, bam_to_chr_list, clean_read_column
from .eb import get_count_df_snp
from .beta_binomial import fit_beta_binomial
import re

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
def pon2pileup(pon_list, state, pon_folder, chromosome):
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
    with open(state['log'], 'w+') as log:
        mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(state['q']), "-Q", str(state['Q']), "--ff", state['ff']]
        mpileup_cmd += ["-b", pon_sub_list]
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
    if state['debug_mode']:
        # 
        pileup_file = 'cache'
        if state['threads'] > 1:
            pileup_file += f"_{chromosome}"
        pileup_file += '.pileup'
        pileup_folder = os.path.join(os.path.dirname(state['cache_name']), 'cache_pileups')

        if not os.path.isdir(pileup_folder):
            os.mkdir(pileup_folder)
        pileup_path = os.path.join(pileup_folder, pileup_file)
        # write the pileup df to csv as <cache_dir>/ab_cache[_chr1].pileup
        pileup_df.to_csv(pileup_path, sep='\t', index=False)
    #####################################################
    return {'df': pileup_df, 'chr': chromosome}


def pileup2AB(state, chromosome, pileup_df):
    '''
    creates the AB_df for a pileup
    '''

    ####################### AB FITTING ###############################################
    # create a copy of just the Chr and Start coords to store the AB values in
    AB_df = pileup_df.iloc[:,[0,1]].copy()

    def get_AB(penalty, pileup_length, row):
        bb_s = pd.Series()
        if row.name % 5000 == 0 and row.name > 0:
            print(f"{math.round(int(row.name) / pileup_length, 1)}% ({row.name} lines) processed..")
        for var in acgt:
            # get the count matrix
            count_df = get_count_df_snp(row, var, 3)
            # get the AB parameters for 
            bb_params = fit_beta_binomial(count_df, penalty)
        # dump the different parameters into bb_s
        # keys have to fit with the var_columns for the receiving AB_df
            bb_s[f'{var}+a'] = bb_params['p'][0]
            bb_s[f'{var}+b'] = bb_params['p'][1]
            bb_s[f'{var}-a'] = bb_params['n'][0]
            bb_s[f'{var}-b'] = bb_params['n'][1]
        return bb_s

    var_columns = [f'{var}{strand}{param}' for var in acgt for strand in ['+', '-'] for param in ['a','b']]
    pileup_length = len(pileup_df.index)
    print(f'Computing ABs for chromosome: {chromosome}\n{pileup_length} rows to go..')
    AB_df[var_columns] = pileup_df.apply(partial(get_AB, state['fitting_penalty'], pileup_length), axis=1)

    ######################## OUTPUT ####################################################
    if chromosome != 'all_chromosomes' and state['debug_mode']: 
        # for multithreading also output the sub_files
        chr_cache = f"{os.path.splitext(state['cache_name'])[0]}_{chromosome}_{os.getpid()}.{os.path.splitext(state['cache_name'])[1]}"
        AB_df.to_csv(chr_cache, sep=',', index=False)

    return AB_df

def generate_cache(pon_list, state):
    '''
    create a cache file for ab-parameters
    '''

    print('Generating Cache...')
    threads = state['threads']
    pon_df = pd.read_csv(pon_list)  
    if threads == 1:
        pileup_dict = pon2pileup(pon_list, state, 'all_chromosomes')
        AB_df = pileup2AB(state, 'all_chromosomes', pileup_dict['df'])
    else: # multithreading
        ####################### SET TO FINAL OUTPUT ########################
        # make directory for temporary bam files
        pon_folder = os.path.join(os.path.dirname(state['cache_name']), 'pon')
        if not os.path.isdir(pon_folder):
            os.mkdir(pon_folder)

        # get the list of valid chromosomes from the first bam in the pon_list
        chromosomes = bam_to_chr_list(pon_df.iloc[0,0])

        cache_pool = Pool(threads)
        # threads are mapped to the pool of chromosomes
        pileup_dicts = cache_pool.map(partial(pon2pileup, pon_list, state, pon_folder), chromosomes)
        cache_pool.close()
        cache_pool.join()
        # split the pileups in multiprocessing chunks to even out different chromosome sizes
        AB_dfs = []
        for pileup_dict in filter(None, pileup_dicts):  # account for empty pileups with filter(None..
            chromosome = pileup_dict['chr']
            pileup_split = np.array_split(pileup_dict['df'], threads * 3)
            pileup2AB_pool = Pool(threads)
            AB_chr_dfs = pileup2AB_pool.map(partial(pileup2AB, state, chromosome), pileup_split)
            pileup2AB_pool.close()
            pileup2AB_pool.join()
            # concatenate the AB_dfs for each chromosome
            AB_chr_df = pd.concat(AB_chr_dfs)
            AB_dfs.append(AB_chr_df.sort_values([out_df.columns[0], out_df.columns[1]]))

            ############### OUTPUT ###########################################################
            chr_cache = f"{os.path.splitext(state['cache_name'])[0]}_{chromosome}.{os.path.splitext(state['cache_name'])[1]}"
            AB_df.to_csv(chr_cache, sep=',', index=False)

        AB_df = pd.concat(AB_dfs).sort_values([out_df.columns[0], out_df.columns[1]])
        AB_df.to_csv(state['cache_name'], sep=',', index=False)

        # AB_df = out_df.sort_values([out_df.columns[0], out_df.columns[1]])

    ############# DEBUG ################################
    # remove the pon_list and all bam files (by removing the whole pon folder)
    if not state['debug_mode']:    
        subprocess.check_call('rm', '-r', pon_folder)

    AB_df.to_csv(state['cache_name'], sep=',', index=False)
    #   AB_columns = pd.MultiIndex.from_product([acgt,['+', '-'],['a','b']], names=['var', 'strand', 'param'])


