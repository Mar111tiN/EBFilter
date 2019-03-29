#! /usr/bin/env python

import process_vcf
import process_anno
import get_eb_score
import sys, os, subprocess, math, re, multiprocessing 
import vcf, pysam, numpy

# ___:123124-124124


def EBFilter_worker_vcf(mut_file, tumor_bam, pon_list, output_path, region):

    pon_count = sum(1 for line in open(pon_list, 'r'))

    ##########
    # generate pileup files
    process_vcf.vcf2pileup(mut_file, output_path + '.target.pileup', tumor_bam, False, region)
    process_vcf.vcf2pileup(mut_file, output_path + '.control.pileup', pon_list, True, region)
    ##########

    ##########
    # load pileup files
    pos2pileup_target = {}
    pos2pileup_control = {}

    hIN = open(output_path + '.target.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_target[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()

    hIN = open(output_path + '.control.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_control[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()
    ##########

    ##########
    # get restricted region if not None
    if is_loption == True and region != "":
        region_match = region_exp.match(region)
        reg_chr = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))
    ##########

    vcf_reader = vcf.Reader(open(mut_file, 'r'))
    vcf_reader.infos['EB'] = vcf.parser._Info('EB', 1, 'Float', "EBCall Score", "EBCall", "ver0.2.0")
    vcf_writer =vcf.Writer(open(output_path, 'w'), vcf_reader)


    for vcf_record in vcf_reader:
        current_pos = str(vcf_record.CHROM) + '\t' + str(vcf_record.POS) 

        if is_loption == True and region != "":
            if reg_chr != vcf_record.CHROM: continue
            if int(vcf_record.POS) < reg_start or int(vcf_record.POS) > reg_end: continue

        F_target = pos2pileup_target[current_pos].split('\t') if current_pos in pos2pileup_target else []
        F_control = pos2pileup_control[current_pos].split('\t') if current_pos in pos2pileup_control else []

        current_ref = str(vcf_record.REF)
        current_alt = str(vcf_record.ALT[0])
        var = ""
        if len(current_ref) == 1 and len(current_alt) == 1:
            var = current_alt
        else:
            if len(current_ref) == 1:
                var = "+" + current_alt[1:]
            elif len(current_alt) == 1:
                var = "-" + current_ref[1:]

        EB_score = "." # if the variant is complex, we ignore that
        if not var == "":
            EB_score = get_eb_score.get_eb_score(var, F_target, F_control, pon_count)

        # add the score and write the vcf record
        vcf_record.INFO['EB'] = EB_score
        vcf_writer.write_record(vcf_record)

    vcf_writer.close()


    # delete intermediate files
    if debug_mode == False:
        subprocess.check_call(["rm", output_path + '.target.pileup'])
        subprocess.check_call(["rm", output_path + '.control.pileup'])



def EBFilter_worker_anno(mut_file, tumor_bam, pon_list, output_path, region):

    pon_count = sum(1 for line in open(pon_list, 'r'))

    # --> process_anno
    if is_loption:
        make_region_list(mut_file) # in utils

    # generate pileup files
    anno2pileup(mut_file, output_path, tumor_bam, region)
    anno2pileup(mut_file, output_path, pon_list, region)
    ##########

    # delete region_list.bed
    if is_loption and not debug_mode:
        subprocess.check_call(["rm", "-f", f"{mut_file}.region_list.bed"])

    ##########
    # load pileup files into dictionaries pos2pileup_target['chr1:123453'] = "depth \t reads \t rQ"
    pos2pileup_target = {}
    pos2pileup_control = {}
 
    with open(f"{output_path}.target.pileup", 'r') as file_in:
        for line in file_in:
            field = line.rstrip('\n').split('\t')
            pos2pileup_target[f"{field[0]}:{field[1]}"] = '\t'.join(field[3:])

    with open(f"{output_path}.control.pileup", 'r') as file_in:
        for line in file_in:
            field = line.rstrip('\n').split('\t')
            pos2pileup_control[f"{field[0]}:{field[1]}"] = '\t'.join(field[3:])
    ##########

     ##########
    # get restricted region if not None
    if is_loption and region:
        region_match = region_exp.match(region)
        reg_chr = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))

    ##########

    with open(mut_file, 'r') as file_in:
        with open(output_path, 'w') as file_out:

            for line in file_in:

                field = line.rstrip('\n').split('\t')
                chr, pos, pos2, ref, alt = field[0], field[1], field[2], field[3], field[4]
                # adjust pos for deletion
                if alt == "-":
                    pos -= 1

                if is_loption and region:
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
                #  -     T      +T
                #  A     -      -A
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
                    EB_score = get_eb_score(var, F_target, F_control, pon_count)

                # add the score and write the vcf record
                print('\t'.join(F + [str(EB_score)]), file=file_out)


    # delete intermediate files
    if debug_mode == False:
        subprocess.check_call(["rm", output_path + '.target.pileup'])
        subprocess.check_call(["rm", output_path + '.control.pileup'])


############ GLOBAL STATE ########################
debug_mode = True
params = config['EB']['params']
threads = config['EB']['threads']
# mapping quality
_q = str(params['map_quality'])  # 20
# base quality
_Q = str(params['base_quality'])
filter_quals = ''
for qual in range( 33, 33 + _Q ):
    filter_quals += chr(qual)

_ff = params['filter_flags'] # 'UNMAP,SECONDARY,QCFAIL,DUP'
is_loption = params['loption'] # False

def main(args):
    '''
    validates files and refers to respective functions
    '''

    # should add validity check for arguments
    mut_file = args['mut_file']
    tumor_bam = args['tumor_bam']
    pon_list = args['pon_list']
    output_path = args['output_path']
    is_anno = not(os.path.splitext(mut_file)[-1] == '.vcf')
    region = args['region']

    # file existence check
    validate(mut_file, tumor_bam, pon_list) 
    if threads == 1:
        # non multi-threading mode
        if is_anno:
            EBFilter_worker_anno(mut_file, tumor_bam, pon_list, output_path, region)
        else: 
            EBFilter_worker_vcf(mut_file, tumor_bam, pon_list, output_path, region)
    else:
        # multi-threading mode
        ##########

        if is_anno:
            # partition anno files
            partition_anno(mut_file, output_path, threads)

            jobs = []

            for i in range(threads):
                worker_args = (f"{output_path}.{i}", tumor_bam, pon_list, f"{output_path}.{i}", region)
                process = multiprocessing.Process(target=EBFilter_worker_anno, args=worker_args)
                    
                jobs.append(process)
                process.start()
        
            # wait all the jobs to be done
            for i in range(threads):
                jobs[i].join()
        
            # merge the individual results
            merge_anno(output_path, threads)
        
            # delete intermediate files
            if debug_mode == False:
                for i in range(threads):
                    subprocess.check_call(["rm", f"{output_path}.{i}", f"{output_path}.{i}.control.pileup", f"{output_path}.{i}.target.pileup"])

        else:
            # partition vcf files
            process_vcf.partition_vcf(mut_file, f"{output_path}.sub.vcf", threads)

            jobs = []
            for i in range(threads):
                process = multiprocessing.Process(target = EBFilter_worker_vcf, args = \
                    (output_path + ".sub.vcf." + str(i), tumor_bam, pon_list, output_path + "." + str(i), region))
                jobs.append(process)
                process.start()

            # wait all the jobs to be done
            for i in range(threads):
                jobs[i].join()

            # merge the individual results
            process_vcf.merge_vcf(output_path + ".", output_path, threads)

            # delete intermediate files
            if debug_mode == False:
                for i in range(threads):
                    subprocess.check_call(["rm", output_path + ".sub.vcf." + str(i)])
                    subprocess.check_call(["rm", output_path + "." + str(i)])



