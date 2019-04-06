import pysam
import re
import sys
import math
from .beta_binomial import fit_beta_binomial, beta_binom_pvalue, fisher_combination
import numpy as np

sign_re = re.compile(r'\^.|\$')
indel_simple = re.compile(r'[\+\-]([0-9]+)')
region_exp = re.compile(r'^([^ \t\n\r\f\v,]+):([0-9]+)\-([0-9]+)')


def get_read_matrix(var, depth, read, rQ, is_verbose, filter_quals):
    '''
    per anno entry: outputs 
    '''

    base_count = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "a": 0, "c": 0, "g": 0, "t": 0, "n": 0}


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



def get_eb_score(pon_count, row):
    """
    calculate the EBCall score from pileup bases of tumor and control samples
    """

    if len(F_target) > 0:
        mut_matrix = get_read_matrix(var, *F_target, False, filter_quals)

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
    return 5


