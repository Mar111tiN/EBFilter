import pandas as pd
import numpy as np
import re
import sys
import math
from .beta_binomial import fit_beta_binomial, beta_binom_pvalues, fisher_combination


sign_re = re.compile(r'\^.|\$')
indel_simple = re.compile(r'[\+\-]([0-9]+)')
region_exp = re.compile(r'^([^ \t\n\r\f\v,]+):([0-9]+)\-([0-9]+)')


def get_read_df(row):
    '''
    extracts the read_count matrix from each compound data frame row
    ___[0] is target
    ___[1:] are control counts
    '''
    var = row['Alt'].upper()
    matrix = pd.DataFrame()
    matrix['depth_p'] = row.iloc[6::3].str.count(r'[ACTG]')
    matrix['mm_p'] = row.iloc[6::3].str.count(var)
    matrix['depth_n'] = row.iloc[6::3].str.count(r'[actg]')
    matrix['mm_n'] = row.iloc[6::3].str.count(var.lower())
    return matrix


def get_EB_score(pen, row):
    """
    calculate the EBCall score from pileup bases of tumor and control samples
    p is penality to use as soft constraints during parameter fitting
    """
    # convert the row into a count matrix for depth and mismatch over reads
    count_df = get_read_df(row)
    # print(count_df)
    ############ FITTING ####################################
    # get the respective control matrices (as dataframe) for positive and negative strands
    control_df = count_df.loc['read1':]
    # estimate the beta-binomial parameters for positive and negative strands  
    bb_params = fit_beta_binomial(control_df, pen)

    ##############################
    # print(f"ab_pn: {bb_params}")
    ##############################

    ########## USING FIT ON TARGET ##########################
    # get the respective target matrix (as dataframe) for positive and negative strands
    target_df = count_df.loc['read0']
    # evaluate the p-values of target mismatch numbers for positive and negative strands
    p_values = beta_binom_pvalues(bb_params, target_df)

    ############################################
    # print(p_values)
    ############################################

    ############ FISHER COMBINATION #########################
    # perform Fisher's combination methods for integrating two p-values of positive and negative strands
    EB_pvalue = fisher_combination(p_values)
    EB_score = 0
    if EB_pvalue < 1e-60:
        EB_score = 60
    elif EB_pvalue > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = -round(math.log10(EB_pvalue), 3)
    return EB_score


