#! /usr/bin/env python
import control_count
import beta_binomial
import utils
import math, numpy

def get_eb_score(var, F_target, F_control, base_qual_thres, controlFileNum):

    """calculate the EBCall score from pileup bases of tumor and control samples"""

    # obtain the mismatch numbers and depths of target sequence data for positive and negative strands
    if len(F_target) > 0:
        varCounts_target_p, depthCounts_target_p, varCounts_target_n, depthCounts_target_n = control_count.varCountCheck(var, F_target[0], F_target[1], F_target[2], base_qual_thres, False)
    else:
        varCounts_target_p, depthCounts_target_p, varCounts_target_n, depthCounts_target_n = 0, 0, 0, 0

    varCounts_control_p = [0] * controlFileNum
    varCounts_control_n = [0] * controlFileNum
    depthCounts_control_p = [0] * controlFileNum
    depthCounts_control_n = [0] * controlFileNum

    # obtain the mismatch numbers and depths (for positive and negative strands) of control sequence data
    # for i in range(controlFileNum):
    for i in range(len(F_control) / 3):
        varCounts_control_p[i], depthCounts_control_p[i], varCounts_control_n[i], depthCounts_control_n[i] = control_count.varCountCheck(var, F_control[3 * i], F_control[1 + 3 * i], F_control[2 + 3 * i], base_qual_thres, True)

    # estimate the beta-binomial parameters for positive and negative strands
    alpha_p, beta_p = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_control_p), numpy.array(varCounts_control_p))
    alpha_n, beta_n = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_control_n), numpy.array(varCounts_control_n))

    # evaluate the p-values of target mismatch numbers for positive and negative strands
    pvalue_p = beta_binomial.beta_binom_pvalue([alpha_p, beta_p], depthCounts_target_p, varCounts_target_p)
    pvalue_n = beta_binomial.beta_binom_pvalue([alpha_n, beta_n], depthCounts_target_n, varCounts_target_n)

    # perform Fisher's combination methods for integrating two p-values of positive and negative strands
    EB_pvalue = utils.fisher_combination([pvalue_p, pvalue_n])
    EB_score = 0
    if EB_pvalue < 1e-60:
        EB_score = 60
    elif EB_pvalue > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = - round(math.log10(EB_pvalue), 3)

    return EB_score

