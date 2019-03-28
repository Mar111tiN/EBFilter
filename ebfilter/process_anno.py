def partition_anno(anno_path, outputFilePrefix, partitionNum):

    file_in = open(anno_path, 'r')
    recordNum = 0 
    for line in file_in:
        recordNum += 1
    file_in.seek(0, 0)

    partitionNum_mod = min(recordNum, partitionNum)
    eachPartitionNum = recordNum / partitionNum_mod

    currentPartition = 0
    currentRecordNum = 0


    file_out = open(outputFilePrefix + "0", 'w')
    for line in file_in:
        print >> file_out, line.rstrip("\n")
        currentRecordNum += 1
        if currentRecordNum >= eachPartitionNum and currentPartition < partitionNum_mod - 1:
            currentPartition += 1
            currentRecordNum = 0
            file_out.close()
            file_out = open(outputFilePrefix + str(currentPartition), 'w')

    file_in.close()
    file_out.close()

    return partitionNum_mod


def merge_anno(inputFilePrefix, out_path, partitionNum):

    file_in = open(inputFilePrefix + "0", 'r')
    file_out = open(out_path, 'w')

    for i in range(partitionNum):
        file_in = open(inputFilePrefix + str(i), 'r')
        for line in file_in:
            print >> file_out, line.rstrip('\n')
        file_in.close()

    file_out.close()



def anno2pileup(anno_path, out_path, bam_or_pon, region):
    '''
    creates a pileup from all the entries in the anno file
    '''
    file_in = open(anno_path, 'r')
    file_out = open(out_path, 'w')
    FNULL = open(os.devnull, 'w')

    mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", _q, "-Q", _Q, "--ff", _ff]

    if is_loption:
        # make bed file for mpileup
        with open(out_path + ".region_list.bed", 'w') as file_out2:
            for line in file_in:
                field = line.rstrip('\n').split('\t')
                loc = int(fields[1]) - (fields[4] == "-")  # -1 if fields 4 == '-' eg. deletion 
                print(field[0], (loc - 1), loc, file=file_out2, sep='\t') 

        mpileup_cmd += ["-l", f"{out_path}.region_list.bed"]

        if region:
            mpileup_cmd = mpileup_cmd + ["-r", region]

    else: # no loption 
        # get lines of anno file
        for line in file_in:
            fields = line.rstrip('\n').split('\t')
            loc = int(fields[1]) - (fields[4] == "-") # -1 if fields 4 == '-' eg. deletion
            mutReg = f"{fields[0]}:{loc}-{loc}"
    
            # set region for mpileup
            mpileup_cmd += ["-r", mutReg]

    # add tumor_bam or pon_list of bam files
    if (os.path.splitext(bam_or_pon)[-1] == '.bam'):
        mpileup_cmd += [bam_or_pon]
    else:
        mpileup_cmd += ["-b", bam_or_pon]
    mpileup_cmd = [str(command) for command in mpileup_cmd]
    subprocess.check_call(mpileup_cmd, stdout=file_out, stderr=FNULL) # maybe logging
    # delete region_list
    if is_loption:
            subprocess.check_call(["rm", "-f", out_path + ".region_list.bed"])

    FNULL.close()
    file_in.close()
    file_out.close()


