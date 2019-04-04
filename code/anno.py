import os
import pandas as pd
from subprocess import Popen, PIPE
from io import StringIO
from .utils import make_region_list
import re
from .eb import get_eb_score

def worker(tumor_bam, pon_list, output_path, region, state, mut_df):

    pon_count = sum(1 for line in open(pon_list, 'r'))

    ########### PANDAS IMPORT ################
    # mut_pd = pd.read_csv(mutfile, sep=',')
    # generate pileup files and store data in mut_df as 
    mut_df = anno2pileup(mut_df, output_path, tumor_bam, pon_list, region, state)
     
    ##########
    # load pileup files into dictionaries pos2pileup_target['chr1:123453'] = "depth \t reads \t rQ"
    # pos2pileup_target = {}
    # pos2pileup_control = {}
 
    # with open(f"{output_path}.target.pileup", 'r') as file_in:
    #     for line in file_in:
    #         field = line.rstrip('\n').split('\t')
    #         # store all but the first 3 colummns
    #         pos2pileup_target[f"{field[0]}:{field[1]}"] = '\t'.join(field[3:])

    # with open(f"{output_path}.control.pileup", 'r') as file_in:
    #     for line in file_in:
    #         field = line.rstrip('\n').split('\t')
    #         pos2pileup_control[f"{field[0]}:{field[1]}"] = '\t'.join(field[3:])
    # ##########

    #  ##########
    # # get restricted region if not None
    # if region:
    #     region_match = region_exp.match(region)
    #     reg_chr = region_match.group(1)
    #     reg_start = int(region_match.group(2))
    #     reg_end = int(region_match.group(3))

    # ##########

    # ############ EB score with pandas ############
    # # mut_pd['EB_score'] = mut_pd.apply(EB_score, axis=1)
 
    # with open(mut_df, 'r') as file_in:
    #     with open(output_path, 'w') as file_out:
    #         for line in file_in:
    #             field = line.rstrip('\n').split(sep)
    #             # header alarm
    #             if field[0] == 'Chr':
    #                 continue
    #             chr, pos, pos2, ref, alt = field[0], int(field[1]), field[2], field[3], field[4]
    #             # adjust pos for deletion
    #             if alt == "-":
    #                 pos -= 1
    #             if region:
    #                 if reg_chr != chr:
    #                     continue
    #                 if (int(pos) < reg_start) or (int(pos) > reg_end):
    #                     continue

    #             # pileup dicts are read into field_target as arrays
    #             field_target = pos2pileup_target[f"{chr}:{pos}"].split('\t') if f"{chr}:{pos}" in pos2pileup_target else []
    #             field_control = pos2pileup_control[f"{chr}:{pos}"].split('\t') if f"{chr}:{pos}" in pos2pileup_control else [] 

    #             # set the variance
    #             # ref   alt    var
    #             #  A     T      T
    #             #  -     T     +T
    #             #  A     -     -A
    #             var = ""
    #             if ref != "-" and alt != "-":
    #                 var = alt
    #             else:
    #                 if ref == "-":
    #                     var = "+" + alt
    #                 elif alt == "-":
    #                     var = "-" + ref
    #             EB_score = "." # if the variant is complex, we ignore that
    #             if var:
    #                 # get_eb_score('+A', [depth, reads, rQ], [depth1, reads1, rQ1, depth2, reads2, rQ2, depth3, reads3, rQ3], 3, state)
    #                 EB_score = get_eb_score(var, field_target, field_control, pon_count, state['filter_quals'])
                
                
    #             # add the score and write the record
    #             print(sep.join(field + [str(EB_score)]), file=file_out)

def anno2pileup(mut_df, out_path, bam, pon_list, region, state):
    '''
    creates a pileup from all the entries in the anno file
    --> out_path.target.pileup
    --> out_path.control.pileup
    '''
    # make region list for use in l_option of mpileup
    bed_file = make_region_list(mut_df, out_path, state['threads']) # in utils --> out_1.region_list.bed

    with open(state['log'], 'w+') as log:
        # determine wether it is bam or pon      
        mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q",str(state['q']), "-Q",str(state['Q']), "--ff",state['ff'], "-l", bed_file]
        if region:
            mpileup_cmd = mpileup_cmd + ["-r", region]

        # adjust for pileup positions for correct merge
        mut_df['Start'] -= (mut_df['Alt'] == '-')
        for out in ['target', 'control']:
            if 'target' in out:                     # pileup from bam file
                mpileup_cmd += [bam]
            else:
                mpileup_cmd += ["-b", pon_list]   # pileup from pon_list
            pileup_stream = Popen(mpileup_cmd, stdout=PIPE, stderr=log)
            pileup_string = StringIO(pileup_stream.communicate()[0].decode('utf-8'))

            names = ['Chr', 'Start', 'Ref']
            if 'target' in out:   
                names += ['depth', 'read', 'Q']
                pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names, dtype = {'Chr':int, 'Start':int, 'Ref':str, 'depth':int, 'read':str, 'Q':str}).drop(columns='Ref')
            else:
                for i in range(10):
                    names += [f"depth{i}", f"read{i}", f"Q{i}"]
                pileup_df = pd.read_csv(pileup_string, sep='\t', header=None, names=names, dtype = {'Chr':int, 'Start':int, 'Ref':str}).drop(columns='Ref')
            if state['debug_mode']:
                out_file = bed_file.replace('region_list.bed', f"{out}.pileup")
                pileup_df.to_csv(out_file, sep='\t')
            ############## MERGE DFs ###########################
            mut_df = pd.merge(left=mut_df, right=pileup_df, how='inner', on=['Chr', 'Start'], how='outer', left_index=True)
    ############## FOR DEBUGGING #######################
    if not state['debug_mode']:
        subprocess.check_call(["rm", bedfile])
    mut_df['Start'] += (mut_df['Alt'] == '-')
    return mut_df
