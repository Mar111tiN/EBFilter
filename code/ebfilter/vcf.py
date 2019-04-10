#! /usr/bin/env python
import vcf
import os
import subprocess
import re

def worker(mut_file, tumor_bam, pon_list, output_path, region):

    pon_count = sum(1 for line in open(pon_list, 'r'))

    ##########
    # generate pileup files
    vcf2pileup(mut_file, output_path + '.target.pileup', tumor_bam, False, region)
    vcf2pileup(mut_file, output_path + '.control.pileup', pon_list, True, region)
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


def vcf2pileup(inputFilePath, outputFilePath, bamPath, mapping_qual_thres, base_qual_thres, filter_flags, is_multi, is_loption, region):

    vcf_reader = vcf.Reader(open(inputFilePath, 'r'))
    hOUT = open(outputFilePath, 'w')
    FNULL = open(os.devnull, 'w')

    if is_loption == True:

        # make bed file for mpileup
        hOUT2 = open(outputFilePath + ".region_list.bed", 'w')
        for record in vcf_reader:
            print >> hOUT2, record.CHROM + '\t' + str(record.POS - 1) + '\t' + str(record.POS)

        hOUT2.close()

        samtools_mpileup_commands = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", \
         str(mapping_qual_thres), "-Q", str(base_qual_thres), "--ff", filter_flags, "-l", outputFilePath + ".region_list.bed"]

        if region != "":
            samtools_mpileup_commands = samtools_mpileup_commands + ["-r", region]

        if is_multi == True:
            samtools_mpileup_commands = samtools_mpileup_commands + ["-b", bamPath]
        else:
            samtools_mpileup_commands = samtools_mpileup_commands + [bamPath]

        subprocess.check_call(samtools_mpileup_commands, stdout = hOUT, stderr = FNULL)
        subprocess.check_call(["rm", "-f", outputFilePath + ".region_list.bed"])

    else:

        for record in vcf_reader:

            mutReg = record.CHROM + ":" + str(record.POS) + "-" + str(record.POS)
        
            samtools_mpileup_commands = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "--ff", filter_flags, "-r", mutReg]

            if is_multi == True:
                samtools_mpileup_commands = samtools_mpileup_commands + ["-b", bamPath]
            else:
                samtools_mpileup_commands = samtools_mpileup_commands + [bamPath]

            # print ' '.join(samtools_mpileup_commands)

            subprocess.check_call(samtools_mpileup_commands, stdout = hOUT, stderr = FNULL)


    FNULL.close()
    hOUT.close()

def partition(inputFilePath, outputFilePrefix, partitionNum):

    vcf_reader1 = vcf.Reader(filename = inputFilePath)
    recordNum = 0
    for record in vcf_reader1:
        recordNum += 1

    partitionNum_mod = min(recordNum, partitionNum)
    eachPartitionNum = recordNum / partitionNum_mod

    currentPartition = 0
    currentRecordNum = 0
    
    vcf_reader2 = vcf.Reader(filename = inputFilePath)
    vcf_writer = vcf.Writer(open(outputFilePrefix + "0", 'w'), vcf_reader2)
    for record in vcf_reader2:
        vcf_writer.write_record(record)
        currentRecordNum += 1
        if currentRecordNum >= eachPartitionNum and currentPartition < partitionNum_mod - 1:
            currentPartition += 1
            currentRecordNum = 0
            vcf_writer.close()
            vcf_writer = vcf.Writer(open(outputFilePrefix + str(currentPartition), 'w'), vcf_reader2) 

    vcf_writer.close()

    return partitionNum_mod


def merge(inputFilePrefix, outputFilePath, partitionNum):

    vcf_reader = vcf.Reader(filename = inputFilePrefix + '0')
    vcf_writer = vcf.Writer(open(outputFilePath, 'w'), vcf_reader)    

    for i in range(partitionNum):
        vcf_reader = vcf.Reader(filename = inputFilePrefix + str(i))
        for record in vcf_reader:
            vcf_writer.write_record(record)
    
    vcf_writer.close()
    


if __name__ == "__main__":
    import sys
    vcf2bed(sys.argv[1], sys.argv[2])

 
