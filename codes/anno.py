import os
import pandas as pd
import subprocess
from subprocess import Popen, PIPE
from io import StringIO
from functools import partial
from .utils import make_region_list, clean_up_df, cleanup_badQ
import re
from .eb import get_EB_score

def worker(tumor_bam, pon_dict, output_path, region, config, mut_df):

    pon_count = len(pon_dict['df'].index)

    ########### PANDAS IMPORT ################
    # mut_pd = pd.read_csv(mutfile, sep=',')
    # generate pileup files and store data in mut_df as 
    mut_df = anno2pileup(mut_df, output_path, tumor_bam, pon_dict, region, config)

    # in_place removal of indel traces and start/end signs in pileup data
    clean_up_df(mut_df, pon_count)

    # cleanup_badQ should not be necessary because these bases have been removed using mpileup -Q option (?)  
    #cleanup_badQ(mut_df, pon_count, config['filter_quals'])

    ############# FOR DEBUGGING #######################
    if config['debug_mode']:
        out_file = output_path.replace('eb', "clean")
        mut_df.to_csv(out_file, sep='\t', index=False)

    ########### EB score ###############################
    mut_df['EB_score'] = mut_df.apply(partial(get_EB_score, config['fitting_penalty']), axis=1)

    # return minimal consensus df for merge with annotated df
    return mut_df.loc[:,['Chr', 'Start', 'EB_score']]

def anno2pileup(mut_df, out_path, bam, pon_dict, region, config):
    '''
    creates a pileup from all the entries in mut_df (the mutation dataframe) and stores the pileup data in mut_df
    '''

    # make region list for use in l_option of mpileup
    bed_file = make_region_list(mut_df, out_path, config['threads']) # in utils --> out_1.region_list.bed
    # get the numbers of control bams from the pon_list
    pon_count = len(pon_dict['df'].index)


    # determine wether it is bam or pon      
    mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q",str(config['q']), "-Q",str(config['Q']), "--ff",config['ff'], "-l", bed_file]
    if region:
        mpileup_cmd = mpileup_cmd + ["-r", region]

    # adjust for pileup positions for correct merge
    mut_df['Start'] -= (mut_df['Alt'] == '-')

    # execute for both pileup and target:
    #   - create the respective mpileup_command
    #   - execute subprocess with Popen and store stream into pileup_df
    #   - merge to mut_df (everything in one big dataframe)
    #   - 
    for out in ['target', 'control']:
        if 'target' in out:                     # pileup from bam file
            mpileup_cmd += [bam]
        else:
            mpileup_cmd += ["-b", pon_dict['list']]   # pileup from pon_list
        pileup_stream = Popen(mpileup_cmd, stdout=PIPE)
        pileup_string = StringIO(pileup_stream.communicate()[0].decode('utf-8'))
        # the columns needed in the dataframe
        names = ['Chr', 'Start', 'Ref']
        # target pileup
        if 'target' in out:   
            names += ['depth0', 'read0', 'Q0']
            pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names).drop(columns='Ref')
        # control pileup
        else:
            # create the columns for the control pileup data: depth0 read0 Q0 depth1 read1 Q1 depth2 ....
            for i in range(pon_count):
                names += [f"depth{i+1}", f"read{i+1}", f"Q{i+1}"]
            pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names).drop(columns='Ref')

        ############# FOR DEBUGGING #######################
        if config['debug_mode']:
            out_file = bed_file.replace('region_list.bed', f"{out}.pileup")
            pileup_df.to_csv(out_file, sep='\t', index=False)

        ############## MERGE DFs ###########################
        mut_df = pd.merge(left=mut_df, right=pileup_df, how='outer', on=['Chr', 'Start'], left_index=True)

    ############## FOR DEBUGGING #######################
    if not config['debug_mode']:
        subprocess.check_call(["rm", bed_file])
    else:
        out_file = bed_file.replace('region_list.bed', f"{out}.merged.csv")
        mut_df.to_csv(out_file, sep='\t', index=False)
    ####################################################

    # revert to the original Start positions
    mut_df['Start'] += (mut_df['Alt'] == '-')

    return mut_df
