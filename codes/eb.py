import pandas as pd
import numpy as np
import re
import sys
import math
from .beta_binomial import fit_bb, bb_pvalues, fisher_combination


sign_re = re.compile(r'\^.|\$')
indel_simple = re.compile(r'[\+\-]([0-9]+)')
region_exp = re.compile(r'^([^ \t\n\r\f\v,]+):([0-9]+)\-([0-9]+)')


def get_count_df_snp(row, var, start_col):
    '''
    extracts the read_count matrix from each compound data frame row
    read/Q/[0] is target
    ___[1:] are control counts
    '''
    matrix = pd.DataFrame()
    matrix['depth_p'] = row.iloc[start_col::3].str.count(r'[ACTG]')
    matrix['mm_p'] = row.iloc[start_col::3].str.count(var)
    matrix['depth_n'] = row.iloc[start_col::3].str.count(r'[actg]')
    matrix['mm_n'] = row.iloc[start_col::3].str.count(var.lower())
    return matrix


def get_count_df_indels(row, start_col):
    '''
    extracts the read_count matrix from each compound data frame row
    read/Q/[0] is target
    ___[1:] are control counts
    '''
    matrix = pd.DataFrae()

    matrix['depth_p'] = row.iloc[start_col::3].str.count(r'[ACTG\-]')
    matrix['mm_p'] = row.iloc[start_col::3].str.count('-')
    matrix['depth_n'] = row.iloc[start_col::3].str.count(r'[actg_]')
    matrix['mm_n'] = row.iloc[start_col::3].str.count('_')
    return matrix



def get_EB_score(pen, row):
    """
    calculate the EBCall score from pileup bases of tumor and control samples
    p is penality to use as soft constraints during parameter fitting
    """
    # convert the row into a count matrix for depth and mismatch over reads
    
    if (row['Ref'] == '-') or row['Alt'] == '-':
        count_df = get_count_df_indels(row, 6)
    else:
        var = row['Alt'].upper()
        count_df = get_count_df_snp(row, var, 6)

    ############ FITTING ####################################
    # get the respective control matrices (as dataframe) for positive and negative strands
    control_df = count_df.loc['read1':]
    # estimate the beta-binomial parameters for positive and negative strands  
    bb_params = fit_bb(control_df, pen)

    ########## USING FIT ON TARGET ##########################
    # get the respective target matrix (as dataframe) for positive and negative strands
    target_df = count_df.loc['read0']
    # evaluate the p-values of target mismatch numbers for positive and negative strands
    p_values = bb_pvalues(bb_params, target_df)

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
