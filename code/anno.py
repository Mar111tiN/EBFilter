import os
import subprocess
from .utils import make_region_list
import re
from .eb import get_eb_score

def worker(mut_file, tumor_bam, pon_list, output_path, region,state):

    pon_count = sum(1 for line in open(pon_list, 'r'))

    ########### PANDAS IMPORT ################
    # mut_pd = pd.read_csv(mutfile, sep=',')
    # generate pileup files
    anno2pileup(mut_file, output_path, tumor_bam, region,state)
    anno2pileup(mut_file, output_path, pon_list, region,state)
    ##########
    # load pileup files into dictionaries pos2pileup_target['chr1:123453'] = "depth \t reads \t rQ"
    pos2pileup_target = {}
    pos2pileup_control = {}
 
    with open(f"{output_path}.target.pileup", 'r') as file_in:
        for line in file_in:
            field = line.rstrip('\n').split('\t')
            # store all but the first 3 colummns
            pos2pileup_target[f"{field[0]}:{field[1]}"] = '\t'.join(field[3:])

    with open(f"{output_path}.control.pileup", 'r') as file_in:
        for line in file_in:
            field = line.rstrip('\n').split('\t')
            pos2pileup_control[f"{field[0]}:{field[1]}"] = '\t'.join(field[3:])
    ##########

     ##########
    # get restricted region if not None
    if region:
        region_match = region_exp.match(region)
        reg_chr = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))

    ##########

    ############ EB score with pandas ############
    # mut_pd['EB_score'] = mut_pd.apply(EB_score, axis=1)
 
    with open(mut_file, 'r') as file_in:
        with open(output_path, 'w') as file_out:
            for line in file_in:
                field = line.rstrip('\n').split(',')
                # header alarm
                if field[0] == 'Chr':
                    continue
                chr, pos, pos2, ref, alt = field[0], int(field[1]), field[2], field[3], field[4]
                # adjust pos for deletion
                if alt == "-":
                    pos -= 1
                if region:
                    if reg_chr != chr:
                        continue
                    if (int(pos) < reg_start) or (int(pos) > reg_end):
                        continue

                # pileup dicts are read into field_target as arrays
                field_target = pos2pileup_target[f"{chr}:{pos}"].split('\t') if f"{chr}:{pos}" in pos2pileup_target else []
                field_control = pos2pileup_control[f"{chr}:{pos}"].split('\t') if f"{chr}:{pos}" in pos2pileup_control else [] 

                # set the variance
                # ref   alt    var
                #  A     T      T
                #  -     T     +T
                #  A     -     -A
                var = ""
                if ref != "-" and alt != "-":
                    var = alt
                else:
                    if ref == "-":
                        var = "+" + alt
                    elif alt == "-":
                        var = "-" + ref
                EB_score = "." # if the variant is complex, we ignore that
                if var:
                    # get_eb_score('+A', [depth, reads, rQ], [depth1, reads1, rQ1, depth2, reads2, rQ2, depth3, reads3, rQ3], 3, state)
                    EB_score = get_eb_score(var, field_target, field_control, pon_count, state['filter_quals'])
                
                
                # add the score and write the vcf record
                print('\t'.join(field + [str(EB_score)]), file=file_out)

def anno2pileup(mut_file, out_path, bam_or_pon, region, state):
    '''
    creates a pileup from all the entries in the anno file
    --> out_path.target.pileup
    --> out_path.control.pileup
    '''
    # make region list for use in l_option of mpileup
    make_region_list(mut_file) # in utils --> mut_file.region_list.bed
    with open(state['log'], 'w') as log:
        with open(mut_file, 'r') as file_in:
            # determine wether it is bam or pon
            is_bam = (os.path.splitext(bam_or_pon)[-1] == '.bam')
            if is_bam:
                out_file = f"{out_path}.target.pileup"
            else:
                out_file = f"{out_path}.control.pileup"
            with open(out_file, 'w') as file_out:
                mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q",state['q'], "-Q",state['Q'], "--ff",state['ff'], "-l", f"{mut_file}.region_list.bed"]

                # add tumor_bam or pon_list of bam files depending on file extension of bam_or_pon
                if is_bam:
                    mpileup_cmd += [bam_or_pon]
                else:
                    mpileup_cmd += ["-b", bam_or_pon]
                if region:
                    mpileup_cmd = mpileup_cmd + ["-r", region]
                subprocess.check_call([str(command) for command in mpileup_cmd], stdout=file_out, stderr=log)

            # direct sed to directly remove the $ and ^] signs seems to break certain things in the quality
                # samtools = subprocess.Popen([str(command) for command in mpileup_cmd], stdout=subprocess.PIPE, stderr=log)
                # sed_cmd = ['sed', '-E', 's/\^\]|\$//g'] # maybe this breaks the whole thing
                # sed = subprocess.Popen(sed_cmd, stdin=samtools.stdout, stdout=file_out)
                # samtools.stdout.close()
                # success = sed.communicate()


def partition(anno_path, out_path, threads):

    
    with open(anno_path, 'r') as file_in:
        # get line number
        record_num = sum(1 for line in file_in)
        file_in.seek(0,0)
        threads = min(record_num, threads)
        # get lines per subprocess
        frac_lines = record_num / threads

        current_sub = current_line = 0
        file_out = open(f"{out_path}.{current_sub}", 'w')
        for line in file_in:
            print(line.rstrip("\n"), file=file_out) 
            current_line += 1
            if (current_line >= frac_lines) and (current_sub < threads - 1):
                current_sub += 1
                current_line = 0
                file_out.close()
                file_out = open(f"{out_path}.{current_sub}", 'w')
        file_out.close()

    return threads


def merge(out_path, threads):

    file_out = open(out_path, 'w')
    for i in range(threads):
        file_in = open(f"{out_path}.sub.{i}", 'r')
        for line in file_in:
            print(line.rstrip('\n'), file=file_out)