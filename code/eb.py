import pysam
import re
import sys
import math
from .beta_binomial import fit_beta_binomial, beta_binom_pvalue, fisher_combination
import numpy as np

sign_re = re.compile(r'\^\]|\$')
region_exp = re.compile(r'^([^ \t\n\r\f\v,]+):(\d+)\-(\d+)')
indel_re = re.compile(r'([\+\-])([0-9]+)([ACGTNacgtn]+)') # +23ATTTNNGC or -34TTCCAAG


def var_count_check(var, depth, reads, rQ, is_verbose, filter_quals):
    '''
    per anno entry: outputs 
    '''
    # var   = '+A'
    # depth = 20
    # reads = 'AAATTCCGGG^]ACGTA$CCT'
    # rQ = 'IIIIIIIFFCDDD'
    # delete the start and end signs
    reads = sign_re.sub('', reads)
    
    if var[0].upper() not in "+-ATGCN":
            print(var + ": input var has wrong format!", file=sys.stderr)
            sys.exit(1)
    if len(reads) != len(rQ):
        print(f"{reads}\n{rQ}", file=sys.stderr)
        print("lengths of bases and qualities are different!", file=sys.stderr)
        sys.exit(1)

    if depth == 0:
        return [0,0,0,0]

    # init
    ins_p = ins_n = del_p = del_n = 0
    ins_vb_p = ins_vb_n = del_vb_p = del_vb_n = 0

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



def get_eb_score(var, F_target, F_control, pon_count, filter_quals):
    """
    calculate the EBCall score from pileup bases of tumor and control samples
    """

    # var = '+A'
    # F_target = [depth, reads, rQ]
    # F_control = [depth1, reads1, rQ1, depth2, reads2, rQ2, depth3, reads3, rQ3]
    # pon_count = 3)

    # obtain the mismatch numbers and depths of target sequence data for positive and negative strands
    if len(F_target) > 0:
        vars_target_p, depth_target_p, vars_target_n, depth_target_n = var_count_check(var, *F_target, False, filter_quals)
    else:
        vars_target_p, depth_target_p, vars_target_n, depth_target_n = 0

    # create [0,0,0,0,0,...,0] arrays for the 4 parameters
    vars_control_p = [0] * pon_count
    vars_control_n = [0] * pon_count
    depth_control_p = [0] * pon_count
    depth_control_n = [0] * pon_count

    # obtain the mismatch numbers and depths (for positive and negative strands) of control sequence data
    # for i in range(len(F_control) / 3):
    for i in range(pon_count):
        vars_control_p[i], depth_control_p[i], vars_control_n[i], depth_control_n[i] = var_count_check(var, *F_control[3*i:3*i+3], True, filter_quals)

    # estimate the beta-binomial parameters for positive and negative strands
    alpha_p, beta_p = fit_beta_binomial(np.array(depth_control_p), np.array(vars_control_p))
    alpha_n, beta_n = fit_beta_binomial(np.array(depth_control_n), np.array(vars_control_n))

    # evaluate the p-values of target mismatch numbers for positive and negative strands
    pvalue_p = beta_binom_pvalue([alpha_p, beta_p], depth_target_p, vars_target_p)
    pvalue_n = beta_binom_pvalue([alpha_n, beta_n], depth_target_n, vars_target_n)

    # perform Fisher's combination methods for integrating two p-values of positive and negative strands
    EB_pvalue = fisher_combination([pvalue_p, pvalue_n])
    EB_score = 0
    if EB_pvalue < 1e-60:
        EB_score = 60
    elif EB_pvalue > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = -round(math.log10(EB_pvalue), 3)
    return EB_score


