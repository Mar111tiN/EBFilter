#! /usr/bin/env python

import pysam, re, sys

indel_re = re.compile(r'([\+\-])([0-9]+)([ACGTNacgtn]+)') # +23ATTTNNGC or -34TTCCAAG
sign_re = re.compile(r'\^\]|\$')


def varCountCheck(var, depth, reads, rQ, is_verbose):
    '''
    per anno entry: outputs 
    '''
    # var   = '+A'
    # depth = 20
    # reads = 'AAATTCCGGG^]ACGTA$CCT'
    # rQ = 'IIIIIIIFFCDDD'

    if var[0].upper() not in "+-ATGCN":
            print(var + ": input var has wrong format!", file=sys.stderr)
            sys.exit(1)
    if len(reads) != len(rQ):
        print("f{reads}\n{rQ}", file=sys.stderr)
        print("lengths of bases and qualities are different!", file=sys.stderr)
        sys.exit(1)

    if depth == 0:
        return [0,0,0,0]

    # delete the start and end signs
    reads = sign_re.sub('', reads)

    # init
    ins_p, ins_n, del_p, del_n = 0
    ins_vb_p, ins_vb_n, del_vb_p, del_vb_n = 0

    #######################       INDELS   #########################################
    #######################       ??????   #########################################
    deleted = 0

    for m in indel_re.finditer(reads):  # match object generator
        site = m.start()
        type = m.group(1)
        indel_size = int(m.group(2))
        varChar = m.group(3)[0:indel_size]

        """
        # just leave these codes for the case where I want to evaluate indels in more detail....
        if not (type in indel and varChar.upper() in indel[type]):
            indel[type][varChar.upper()]['+'] = 0
            indel[type][varChar.upper()]['-'] = 0
       
        strand = '+' if varChar.islower() else '-' 
        indel[type][varChar.upper()][strand] += 1
        """
        # checking if size of del is OK
        var_match = False
        if var[0] == "-":
            if len(var[1:]) == len(varChar):
                var_match = True
        elif var[1:].upper() == varChar.upper():
            var_match = True


        if type == "+": # ins
            if varChar.isupper():
                ins_vb_p += 1
                if var_match:
                    ins_p += 1
            else:
                ins_vb_n += 1
                if var_match:
                    ins_n += 1
        else:   # del
            if varChar.isupper():
                del_vb_p += 1
                if var_match:
                    del_p += 1
            else:
                del_vb_n += 1
                if var_match:
                    del_n += 1

        print("Indels: {type}\t}{var}\t'{varChar.upper()}\t{strand}")

        reads = reads[0:(site - deleted)] + reads[(site + indel_size + len(indel_size) + 1 - deleted):]
        deleted += 1 + len(indel_size) + indel_size

    #############################################################################

    base_count = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "a": 0, "c": 0, "g": 0, "t": 0, "n": 0}

    # count all the bases in the read and allocate to base_count dict
    if var.upper() in "ACGT":
        for base, qual in zip(reads, rQ):
            if not (qual in filter_quals):
                if base in "ATGCNatgcn":
                    base_count[base] += 1
    # for indel check we ignore base qualities
    else:
        for base in reads:
            if base in "ATGCNatgcn":
                base_count[base] += 1


    # sum up the forward and reverse bases
    depth_p = base_count["A"] + base_count["C"] + base_count["G"] + base_count["T"] + base_count["N"]
    depth_n = base_count["a"] + base_count["c"] + base_count["g"] + base_count["t"] + base_count["n"]

    mismatch_p = 0
    mismatch_n = 0
   
    if var.upper() in "ACGT":
        mismatch_p = base_count[var.upper()]
        mismatch_n = base_count[var.lower()]    
    else:
        if var[0] == "+":
            if is_verbose:
                mismatch_p, mismatch_n = ins_vb_p, ins_vb_n
            else:
                mismatch_p, mismatch_n = ins_p, ins_n
        elif var[0] == "-":
            if is_verbose:
                mismatch_p, mismatch_n = del_vb_p, del_vb_n
            else:
                mismatch_p, mismatch_n = del_p, del_n

    return [mismatch_p, depth_p, mismatch_n, depth_n]




