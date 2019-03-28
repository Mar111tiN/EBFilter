#! /usr/bin/env python

import sys, vcf, os, subprocess

def partition_anno(inputFilePath, outputFilePrefix, partitionNum):

    hIN = open(inputFilePath, 'r')
    recordNum = 0 
    for line in hIN:
        recordNum += 1
    hIN.seek(0, 0)

    partitionNum_mod = min(recordNum, partitionNum)
    eachPartitionNum = recordNum / partitionNum_mod

    currentPartition = 0
    currentRecordNum = 0


    hOUT = open(outputFilePrefix + "0", 'w')
    for line in hIN:
        print >> hOUT, line.rstrip("\n")
        currentRecordNum += 1
        if currentRecordNum >= eachPartitionNum and currentPartition < partitionNum_mod - 1:
            currentPartition += 1
            currentRecordNum = 0
            hOUT.close()
            hOUT = open(outputFilePrefix + str(currentPartition), 'w')

    hIN.close()
    hOUT.close()

    return partitionNum_mod


def merge_anno(inputFilePrefix, outputFilePath, partitionNum):

    hIN = open(inputFilePrefix + "0", 'r')
    hOUT = open(outputFilePath, 'w')

    for i in range(partitionNum):
        hIN = open(inputFilePrefix + str(i), 'r')
        for line in hIN:
            print >> hOUT, line.rstrip('\n')
        hIN.close()

    hOUT.close()



def anno2pileup(inputFilePath, outputFilePath, bamPath, mapping_qual_thres, base_qual_thres, filter_flags, is_multi, is_loption, region):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
    FNULL = open(os.devnull, 'w')

    if is_loption == True:

        # make bed file for mpileup
        hOUT2 = open(outputFilePath + ".region_list.bed", 'w')
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            if F[4] == "-": # for deletion in anno format
                print >> hOUT2, F[0] + '\t' + str(int(F[1]) - 2) + '\t' + str(int(F[1]) - 1) 
                # print >> hOUT2, F[0] + '\t' + str(int(F[1]) - 1)
            else:
                print >> hOUT2, F[0] + '\t' + str(int(F[1]) - 1) + '\t' + F[1]
                # print >> hOUT2, F[0] + '\t' + F[1]

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
        for line in hIN:

            F = line.rstrip('\n').split('\t')

            if F[4] == "-": # for deletion in anno format
                mutReg = F[0] + ":" + str(int(F[1]) - 1)  + "-" + str(int(F[1]) - 1)  
            else:
                mutReg = F[0] + ":" + F[1] + "-" + F[1]
    
            samtools_mpileup_commands = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "--ff", filter_flags, "-r", mutReg]

            if is_multi == True:
                samtools_mpileup_commands = samtools_mpileup_commands + ["-b", bamPath]
            else:
                samtools_mpileup_commands = samtools_mpileup_commands + [bamPath]

            subprocess.check_call(samtools_mpileup_commands, stdout = hOUT, stderr = FNULL)

    FNULL.close()
    hIN.close()
    hOUT.close()


