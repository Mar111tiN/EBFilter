import pysam
import re
import sys
import math
from .beta_binomial import fit_beta_binomial, beta_binom_pvalue, fisher_combination
import numpy as np

sign_re = re.compile(r'\^.|\$')
indel_simple = re.compile(r'[\+\-]([0-9]+)')
region_exp = re.compile(r'^([^ \t\n\r\f\v,]+):([0-9]+)\-([0-9]+)')


def var_count_check(var, depth, read, rQ, is_verbose, filter_quals):
    '''
    per anno entry: outputs 
    '''
    # var   = '+A'
    # depth = 20
    # read = 'AAATTCCGGG^]ACGTA$CCT'
    # rQ = 'IIIIIIIFFCDDD'
    ####### VALIDATE ################################################
    if var[0].upper() not in "+-ATGCN":
            print(var + ": input var has wrong format!", file=sys.stderr)
            sys.exit(1)
    if depth == 0:
        return [0,0,0,0]
    # remove $ and ^] and ^I and so on
    read = sign_re.sub('', read)

    #######################       INDELS   #########################################
    if indel_simple.search(read):
        if var[0] in "+-" :  # is indel
            l = len(var) - 1
        else:
        # indels in reads that where not called
            l = indel_simple.search(read).group(1)
        re_string = r"([ACGTNacgtn])([\+\-])([0-9]+)([ACGTNacgtn]{" + str(l) + "})"
        indel_re = re.compile(re_string)
        read = indel_re.sub(r'\1', read)

    
        re_string = r"([ACGTNacgtn])([\+\-])([0-9]+)([ACGTNacgtn]{" + str(l) + "})"  # for exact length of insert
        indel_re = re.compile(re_string)
        read = indel_re.sub(r'\1', read)



    #############################################################################
    if len(read) != len(rQ):
        print(f"{read}\n{rQ}", var, file=sys.stderr)
        print("lengths of bases and qualities are different!", file=sys.stderr)
        sys.exit(1)


    base_count = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "a": 0, "c": 0, "g": 0, "t": 0, "n": 0}

    # count all the bases in the read and allocate to base_count dict
    if var.upper() in "ACGT":
        for base, qual in zip(read, rQ):
            if not (qual in filter_quals):
                if base in "ATGCNatgcn":
                    base_count[base] += 1
    # for indel check we ignore base qualities
    else:
        for base in read:
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


    ################## DEBUG ###########################################            
    # print(f'var_count_check:\nvar {var}:{read}\nmismatch_p: {mismatch_p}\ndepth_p: {depth_p}\nmismatch_n: {mismatch_n}\ndepth_n: {depth_n}')
    return [mismatch_p, depth_p, mismatch_n, depth_n]



def get_eb_score(var, F_target, F_control, pon_count, filter_quals):
    """
    calculate the EBCall score from pileup bases of tumor and control samples
    """
    # var = '+A'
    # F_target = [depth, read, rQ]
    # F_control = [depth1, read1, rQ1, depth2, read2, rQ2, depth3, read3, rQ3]
    # pon_count = 3)
    # obtain the mismatch numbers and depths of target sequence data for positive and negative strands

    if len(F_target) > 0:
        vars_target_p, depth_target_p, vars_target_n, depth_target_n = var_count_check(var, *F_target, False, filter_quals)

    else:
        vars_target_p = depth_target_p = vars_target_n = depth_target_n = 0

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
    ################# Debugging #######################

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


