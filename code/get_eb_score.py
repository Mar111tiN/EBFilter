def get_eb_score(var, F_target, F_control, pon_count):
    """
    calculate the EBCall score from pileup bases of tumor and control samples
    """

    # var = '+A'
    # F_target = [depth, reads, rQ]
    # F_control = [depth1, reads1, rQ1, depth2, reads2, rQ2, depth3, reads3, rQ3]
    # pon_count = 3)

    # obtain the mismatch numbers and depths of target sequence data for positive and negative strands
    if len(F_target) > 0:
        vars_target_p, depth_target_p, vars_target_n, depth_target_n = var_count_check(var, *F_target, False)
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
        vars_control_p[i], depth_control_p[i], vars_control_n[i], depth_control_n[i] = var_count_check(var, *F_control[3*i:3*i+3], True)

    # estimate the beta-binomial parameters for positive and negative strands
    alpha_p, beta_p = fit_beta_binomial(np.array(depth_control_p), np.array(vars_control_p))
    alpha_n, beta_n = fit_beta_binomial(np.array(depth_control_n), np.array(vars_control_n))

    # evaluate the p-values of target mismatch numbers for positive and negative strands
    pvalue_p = beta_binom_pvalue([alpha_p, beta_p], depth_target_p, vars_target_p)
    pvalue_n = beta_binom_pvalue([alpha_n, beta_n], depth_target_n, vars_target_n)

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