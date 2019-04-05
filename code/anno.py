import os
import pandas as pd
from subprocess import Popen, PIPE
from io import StringIO
from .utils import make_region_list, clean_up_df, cleanup_badQ
import re
from .eb import get_eb_score

def worker(tumor_bam, pon_list, output_path, region, state, mut_df):

    pon_count = sum(1 for line in open(pon_list, 'r'))

    ########### PANDAS IMPORT ################
    # mut_pd = pd.read_csv(mutfile, sep=',')
    # generate pileup files and store data in mut_df as 
    mut_df = anno2pileup(mut_df, output_path, tumor_bam, pon_list, region, state)

    # in_place removal of indel traces and start/end signs in pileup data
    clean_up_df(mut_df, pon_count)

    # cleanup_badQ should not be necessary because these bases have been removed using mpileup -Q option (?)  
    # cleanup_badQ(mut_df, pon_count, state['filter_quals'])

    ############# FOR DEBUGGING #######################
    if state['debug_mode']:
        out_file = output_path.replace('eb', "clean")
        mut_df.to_csv(out_file, sep='\t', index=False)

    ########### EB score ############
    # mut_pd['EB_score'] = mut_pd.apply(EB_score, axis=1)

    #                 # get_eb_score('+A', [depth, reads, rQ], [depth1, reads1, rQ1, depth2, reads2, rQ2, depth3, reads3, rQ3], 3, state)
    #                 EB_score = get_eb_score(var, field_target, field_control, pon_count, state['filter_quals'])
                
                
    #             # add the score and write the record
    #             print(sep.join(field + [str(EB_score)]), file=file_out)
    return mut_df

def anno2pileup(mut_df, out_path, bam, pon_list, region, state):
    '''
    creates a pileup from all the entries in mut_df (the mutation dataframe) and stores the pileup data in mut_df
    '''

    # make region list for use in l_option of mpileup
    bed_file = make_region_list(mut_df, out_path, state['threads']) # in utils --> out_1.region_list.bed
    # get the numbers of control bams from the pon_list
    pon_count = sum(1 for line in open(pon_list, 'r'))

    with open(state['log'], 'w+') as log:
        # determine wether it is bam or pon      
        mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q",str(state['q']), "-Q",str(state['Q']), "--ff",state['ff'], "-l", bed_file]
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
                mpileup_cmd += ["-b", pon_list]   # pileup from pon_list
            pileup_stream = Popen(mpileup_cmd, stdout=PIPE, stderr=log)
            pileup_string = StringIO(pileup_stream.communicate()[0].decode('utf-8'))

            names = ['Chr', 'Start', 'Ref']

            # target pileup
            if 'target' in out:   
                names += ['depth0', 'read0', 'Q0']
                pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names, dtype = {'Chr':int, 'Start':int, 'Ref':str, 'depth':int, 'read':str, 'Q':str}).drop(columns='Ref')
            # control pileup
            else:
                # create the columns for the control pileup data: depth0 read0 Q0 depth1 read1 Q1 depth2 ....
                for i in range(10):
                    names += [f"depth{i+1}", f"read{i+1}", f"Q{i+1}"]
                pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names, dtype = {'Chr':int, 'Start':int, 'Ref':str}).drop(columns='Ref')

            ############# FOR DEBUGGING #######################
            if state['debug_mode']:
                out_file = bed_file.replace('region_list.bed', f"{out}.pileup")
                pileup_df.to_csv(out_file, sep='\t', index=False)

            ############## MERGE DFs ###########################
            mut_df = pd.merge(left=mut_df, right=pileup_df, how='outer', on=['Chr', 'Start'], left_index=True)

    ############## FOR DEBUGGING #######################
    if not state['debug_mode']:
        subprocess.check_call(["rm", bedfile])
    else:
        out_file = bed_file.replace('region_list.bed', f"{out}.merged.csv")
        mut_df.to_csv(out_file, sep='\t', index=False)
    ####################################################

    # revert to the original Start positions
    mut_df['Start'] += (mut_df['Alt'] == '-')

    return mut_df
