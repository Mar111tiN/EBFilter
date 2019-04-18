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

sign_re = re.compile(r'\^.|\$')
acgt = ['A','C','T','G']


############# PILEUP2BAM ############################################
def pon2pileup(pon_dict, config, chromosome):
    '''
    create the dataframe of AB parameters per region per mismatch base
    '''

    # get the number of files in the pon_list
    pon_count = len(pon_dict['df'].index)
    pon_folder = config['pon_folder']
    # create a chromosome-bound bam and bai for each bam in the pon_list
    # write the created bams to a dataframe and output as pon_list
    print(f'Generating  pileup for chromosome {chromosome}..')
    pon_sub_df = pd.DataFrame()
    print(f'Splitting bam files for chromosome {chromosome}..')
    pon_sub_df['bam'] = pon_dict['df'].apply(partial(utils.split_bam, chromosome, pon_folder), axis=1)
    # use pon_list_chr?.txt instead of the global pon_list.txt
    pon_sub_list = os.path.join(pon_folder, f"pon_list_{chromosome}.txt")
    # write the pon_list_{chr#} to file in output/pon for access by pon2pileup
    pon_sub_df.to_csv(pon_sub_list, header=None, index=False)


    ############# PON2PILEUP ############################################
    # from the pon_list return a pileup_df
    # mpileup configs
    mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(config['q']), "-Q", str(config['Q']), "--ff", config['ff']]
    # mpileup PoN list
    mpileup_cmd += ["-b", pon_sub_list]
    # optional (but highly recommended:-) bed_file for -l option
    if config['bed_file']:
        mpileup_cmd += ["-l", config['bed_file']]

    # create the pileup folder
    config['pileup_folder'] = pileup_folder = os.path.join(config['cache_folder'], 'cache_pileups')
    if not os.path.isdir(pileup_folder):
        os.mkdir(pileup_folder)

    pileup_file = os.path.join(pileup_folder, f"cache_{chromosome}.pileup")
    mpileup_cmd += ['-o', pileup_file]
    subprocess.check_call(mpileup_cmd)

###################### THE BETTER (but not working) WAY #######################       
    # # use StringIO for small files (not implemented)
    # else:
    #     pileup_stream = Popen(mpileup_cmd, stdout=PIPE)
    #     pileup_file = StringIO(pileup_stream.communicate()[0].decode('utf-8'))
###############################################################################

    names = ['Chr', 'Start', 'Ref']
    for i in range(pon_count):
        names += [f"depth{i}", f"read{i}", f"Q{i}"]
    pileup_df = pd.read_csv(pileup_file, sep='\t', header=None, names=names)
    if len(pileup_df.index) == 0:
        print(f"Pileup for chromosome {chromosome} is empty and will be dropped..")
        return {'file': 'empty', 'chr': chromosome}

    ############## CLEANUP PILEUP #######################################
    # get the count from the number of columns of pileup_df
    pon_count = int(pileup_df.shape[1] / 3 - 1)
    # get the read columns for the clean_columns apply
    read_columns = [f'read{i}' for i in range(pon_count)]
    # apply the cleaning function on the read rows of the pileup
    pileup_df[read_columns] = pileup_df[read_columns].apply(utils.clean_read_column)

    ############## WRITE TO FILE #############################################
    # write the pileup df to csv as <cache_folder>/ab_cache[_chr1].pileup
    print(f"Writing pileup of Chr {chromosome} to file {pileup_file}")
    pileup_df.to_csv(pileup_file, sep='\t', index=False)

    return {'file': pileup_file, 'chr': chromosome}


def pileup2AB(config, chromosome, pileup_df):
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
        if (row.name - start) % 2500 == 0 and (row.name - start) > 0:
            bar_size = 25
            frac = (int(row.name) - start) / length
            progress = '|' + '.' * math.floor(frac * bar_size) + ' ' * math.ceil(bar_size * (1 - frac)) + '|'
            print(f"P{os.getpid()}: {row.name - start} lines (total {round(frac * 100, 1)}% of Chr {chromosome.replace('chr', '').replace('Chr', '')}) processed\t {progress}")

        ########### get count matrix ###########################################    
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

    ################## Store AB data into df ####################################   
    # create the columns (A+a A+b A-a A-b G+a G+b....) for the recipient df of the pileup_df apply function 
    var_columns = [f'{var}{strand}{param}' for var in acgt for strand in ['+', '-'] for param in ['a','b']]
    pileup_length = len(pileup_df.index)
    pileup_start = pileup_df.iloc[0].name
    print(f'Process {os.getpid()}: Computing ABs for Chr {chromosome} \t(lines {pileup_start }\t to \t{pileup_length + pileup_start})')
    AB_df[var_columns] = pileup_df.apply(partial(get_AB, config['fitting_penalty'], pileup_length, pileup_start, chromosome), axis=1)
    print(f'Process {os.getpid()}: Computing ABs for chromosome {chromosome} (lines {pileup_start } to {pileup_length + pileup_start}) finished.')

    ###################### DEBUG #######################################################
    if config['debug_mode']: 
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

def generate_cache(pon_dict, config):
    '''
    create a cache file for ab-parameters
    '''

    print('Generating Cache...')
    threads = config['threads']

    # make directory for temporary bam files
    config['pon_folder'] = pon_folder = os.path.join(config['cache_folder'], 'pon')
    if not os.path.isdir(pon_folder):
        os.mkdir(pon_folder)

    ######################### PON2PILEUP ###################################
    # get the list of chromswithout existing cache file from config['chr']
    config['chr'] = utils.check_cache_files(config)
    # get the list of chroms without existing pileup_files from config['chr']
    config['chr'] = utils.check_pileup_files(config)

    # init the processor pool
    cache_pool = Pool(threads)
    # threads are mapped to the pool of chroms

    ###################### !!!! multiprocessing.Pool can only transfer data up to a certain size
    # pileup_file_dicts stores the list of pileup file dictionaries [{'file': 'chr11.pileup', 'chr': 'chr11'}, {'file': 'chr2.pileup', 'chr': 'chr2'} ]
    pileup_file_dicts = []
    success_messages = []
    # OR send in the path to the dict and transfer df within pon2pileup to global storage list (does not work so far)

    pileup_file_dicts = cache_pool.map(partial(pon2pileup, pon_dict, config), config['chr'])
    # remove Nones and sort for chromosome name
    cache_pool.close()
    cache_pool.join()

    # reading the pileup files into dfs
    pileup_dicts = []
    for pileup_file_dict in pileup_file_dicts:
        if pileup_file_dict['file'] == 'empty':
           pileup_dicts.append({'chr': pileup_file_dict['chr'], 'empty': True})
           continue
        pileup_df = pd.read_csv(pileup_file_dict['file'], sep='\t')

        print(f"Reading pileup {pileup_file_dict['file']} for AB computation")
        pileup_dict = {'df': pileup_df, 'chr': pileup_file_dict['chr'], 'empty': False}
        pileup_dicts.append(pileup_dict)

        ########## DEBUG ####################################
        if not config['debug_mode']:
            subprocess.check_call(['rm', '-f', pileup_file_dict['file']])
        #####################################################

    # account for empty pileups with filter(None..
    pileup_dicts = sorted(filter(None, pileup_dicts), key=utils.sort_chr)

    ######################### PILEUP2AB ###################################
    AB_dfs = []
    # split the pileups for each Chr into threaded chunks to even out different chromosome sizes
    # go through each chromosome with required number of threads
    # create the job pool for the threads
    for pileup_dict in pileup_dicts:  
        chromosome = pileup_dict['chr']
        chr_cache = os.path.join(config['cache_folder'], f"{chromosome}.cache")
        if pileup_dict['empty']:
            print(f"Writing empty cache for Chr {chromosome} to file {chr_cache}.")
            open(chr_cache, 'a').close()
            continue
        chr_len = len(pileup_dict['df'].index)      # get length for progress info     
        # set the minimum number of lines for one thread to 2000
        split_factor = min(math.ceil(chr_len / 2000), threads)
        if split_factor > 1:
            print(f"Splitting the {chr_len} lines of {chromosome}.pileup into {threads} chunks for multithreaded computation..")
        # split the arrays into litte fractions for computation
        pileup_split = np.array_split(pileup_dict['df'], split_factor)
        pileup2AB_pool = Pool(threads)
        AB_chr_dfs = pileup2AB_pool.map(partial(pileup2AB, config, chromosome), pileup_split)
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

    ############# DEBUG #####################################################################
    # remove the pon_list and all bam files (by removing the whole pon folder)
    if not config['debug_mode'] and threads > 1:    
        subprocess.check_call(['rm', '-r', pon_folder])
        subprocess.check_call(['rm', '-r', config['pileup_folder']])
    #########################################################################################   

    return 'Generation of AB file finished.'


def get_EBscore_from_AB(pen, row):
    '''
    get the EBscore from cached AB parameters
    no fitting is needed as parameters are precomputed and stored in row[5:9]
    '''

    # we only get the snp count_df, using the mut_df 'Alt' as var and adjust for AB_df with column 9
    count_df = get_count_df_snp(row, row['Alt'].upper(), 10)
    bb_params = {}

    # feed-in the AB params coming with the row
    bb_params['p'] = [row['a+'], row['b+']]
    bb_params['n'] = [row['a-'], row['b-']]

    # get the dataSeries for read0
    target_df = count_df.loc['read0']
    p_values = bb_pvalues(bb_params, target_df)

    ############ FISHER COMBINATION #########################
    # perform Fisher's combination methods for integrating two p-values of positive and negative strands
    EB_pvalue = fisher_combination(p_values)
    EB_score = 0
    if EB_pvalue < 1e-60:
        EB_score = 60
    elif EB_pvalue > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = -round(math.log10(EB_pvalue), 3)
    return EB_score


def get_EB_from_cache(snp_df, tumor_bam, config, chrom):
    '''
    takes a shortened annotation file (snps only) and loads the valid AB parameters and gets the EB score for these lines
    '''

    ############################## GET THE AB parameters
    cache_file = os.path.join(config['cache_folder'], f"{chrom}.cache")
    # the merger columns
    cols = ['Chr', 'Start', 'Alt']
    # load the file and set Chr and Start to index for proper setting of multi-index
    print(f"Loading cache {cache_file}..")
    AB_df = pd.read_csv(cache_file, sep=',').set_index(cols[:2]) 
    AB_columns = pd.MultiIndex.from_product([['A','C','T','G'],['+', '-'],['a','b']], names=['var', 'strand', 'param'])
    # set multi-indexed columns
    AB_df.columns = AB_columns
    # stack the var column level for merge with the snp_df
    AB_df = AB_df.stack('var')
    # reduce the column index level to 1
    AB_df.columns = AB_df.columns.droplevel(0)
    # unset the indices and transfer to columns
    AB_df = AB_df.reset_index()
    # rename columns for merge
    AB_df.columns = cols + ['a+', 'b+', 'a-', 'b-']
    snpAB_df = pd.merge(snp_df, AB_df, on=cols)
    
    # set region list file_path
    region_list = os.path.join(config['cache_folder'], chrom)
    snpAB_df = anno2pileup(snpAB_df, region_list, tumor_bam, None, None, config)
    # remove start/end signs
    utils.cleanup_df(snpAB_df, 0, config)
    snpAB_df['EB_score'] = snpAB_df.apply(partial(get_EBscore_from_AB, config['fitting_penalty']), axis=1)

    return snpAB_df.loc[:,['Chr', 'Start', 'EB_score']]
