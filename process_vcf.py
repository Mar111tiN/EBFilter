#! /usr/bin/env python

import vcf, os, subprocess


def partition_vcf(inputFilePath, outputFilePrefix, partitionNum):

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


def merge_vcf(inputFilePrefix, outputFilePath, partitionNum):

    vcf_reader = vcf.Reader(filename = inputFilePrefix + '0')
    vcf_writer = vcf.Writer(open(outputFilePath, 'w'), vcf_reader)    

    for i in range(partitionNum):
        vcf_reader = vcf.Reader(filename = inputFilePrefix + str(i))
        for record in vcf_reader:
            vcf_writer.write_record(record)
    
    vcf_writer.close()



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


if __name__ == "__main__":
    import sys
    vcf2bed(sys.argv[1], sys.argv[2])
 
