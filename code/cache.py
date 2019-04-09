from multiprocessing import Pool
from functools import partial
import subprocess
import os
import pandas as pd
from .utils import clean_up_df, cleanup_badQ, bam_to_chr_list
from .eb import get_count_df_snp
from .beta_binomial import fit_beta_binomial

def generate_cache(pon_list, state):
    '''
    create a cache file for ab-parameters
    '''

    print('Generating Cache...')
    threads = state['threads']
    pon_df = pd.read_csv(pon_list)  
    if threads == 1:
        out_df = pile2AB(pon_list, state, 'all_chromosomes')

    else: # multithreading

        ####################### SET TO FINAL OUTPUT ########################
        # make directory for temporary bam files
        pon_sub_folder = os.path.join(state['cache_dir'], 'pon')
        if not os.path.isdir(pon_sub_folder):
            os.mkdir(pon_sub_folder)

        ##################### TODO ##########################################

        chromosomes = bam_to_chr_list(pon_df.iloc[0,0])
        print(chromosomes)
        #####################################################################

        cache_pool = Pool(threads)
        # threads are mapped to the pool of chromosomes
        out_dfs = cache_pool.map(partial(pile2AB, pon_list, state), chromosomes)  # mut_split is the iterable df_pool
        anno_pool.close()
        anno_pool.join()
        out_df = pd.concat(out_dfs)
        out_df = out_df.sort_values([out_df.columns[0], out_df.columns[1]])


    def pile2AB(pon_list, state, chromosome):
        '''
        create the dataframe of AB parameters per region per mismatch base
        '''
        def get_pileup_df(pon_file):
            '''
            from the pon_list return a pileup_df
            '''

            with open(state['log'], 'w+') as log:
                mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q",str(state['q']), "-Q",str(state['Q']), "--ff",state['ff']]
                mpileup_cmd += ["-b", pon_file]
                pileup_stream = Popen(mpileup_cmd, stdout=PIPE, stderr=log)
                pileup_string = StringIO(pileup_stream.communicate()[0].decode('utf-8'))
                pileup_stream = Popen(mpileup_cmd, stdout=PIPE, stderr=log)
                pileup_string = StringIO(pileup_stream.communicate()[0].decode('utf-8'))
                names = ['Chr', 'Start', 'Ref']
                for i in range(pon_count):
                    names += [f"depth{i}", f"read{i}", f"Q{i}"]
            pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names)

            ########## DEBUG ####################################
            if state['debug_mode']:
                pileup_path = os.path.join(state['cache_dir'], 'ab_cache')
                if state['threads'] > 1:
                    pileup_path += f"_{os.getpid()}"
                pileup_path += '.pileup'
                pileup_df.to_csv(out_file, sep='\t', index=False)
            else:
                subprocess.check_call('rm', '-r', pon_sub_folder)
            #####################################################

            return pileup_df


        def split_bam(chromosome, pon_row):
            '''
            creates a sub bam file (+bai) of the specified chromosomal region per pon_list row
            returns the name of the bam file for creating a sub pon_list
            '''

            bam_file = pon_row[0]
            bam_out = os.path.join(pon_sub_folder, os.path.splitext(os.path.basename(bam_file))[0], chromosome, '.bam')
            split_bam_cmd = ["samtools", "view", "-b", "-o", "bam_out", bam_file, str(chromosome)]
            bam_index_cmd = ["samtools", "index", bam_out]
            subprocess.check_call(split_bam_cmd)
            subprocess.check_call(bam_index_cmd)
            return bam_out


        pon_count = sum(1 for line in open(pon_list, 'r'))
        mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(state['q']), "-Q",str(state['Q']), "--ff",state['ff']]
        if chromosome == 'all_chromosomes': # single threaded
            pileup_df = get_pileup_df(pon_list)
        else:  # multithreading
            # create a chromosome-bound bam and bai for each bam in the pon_list
            pon_sub_df = pd.DataFrame()
            pon_sub_df['bam'] = pon_df.apply(partial(split_bam, chromosome), axis=1)
            # write the pon_list_{chr#} to file in output/pon

            ############# DEBUG ################################
            if state['debug_mode']:
                pon_sub_file = os.path.join(pon_sub_folder, f"pon_list_{chromosome}.txt")
                pon_sub_df.to_csv(pon_sub_file, header=None, index=False)
            # get pileup_df from the pon file
        pileup_df = get_pileup_df(pon_sub_file)
        
        ################### ALL THREADS #######################################
        # clean up the dataframe pon_count -1 because there is no target


        ############# FOR DEBUGGING #######################
        if state['debug_mode']:
            pileup_df.to_csv(out_file, sep='\t', index=False)
        AB_df = pileup_df.iloc[:,[0,1]].copy()

        def get_AB(var, penalty, pon_count, row):
            for i in range(pon_count):
                read = f"read{i}"
                mut_df[read] = mut_df[read].str.replace(sign_re, '')
                count_df = get_count_df_snp(row, 3)
                bb_params = fit_beta_binomial(control_df, pen)
                

        for var in ['A', 'C', 'T', 'G']:
            AB_df[var] = pileup_df.apply(partial(get_AB, var, state['fitting_penalty']), axis=1)

    
        
        


