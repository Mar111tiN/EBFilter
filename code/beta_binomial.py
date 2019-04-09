import numpy as np
import math
from scipy.optimize import fmin_l_bfgs_b as minimize_func
from scipy.stats import chi2
from scipy.special import gammaln


def fisher_combination(p_values):

    if 0 in p_values.values():
        return 0
    else:
        return 1 - chi2.cdf(sum([-2 * math.log(x) for x in p_values.values()]), 2 * len(p_values.values()))


def beta_binom_pvalues(params, target_df):
    '''
    accumulate p_value of target observation falling in fitted bb_distribution (not a variant)
    p_values are computed per strand (pvalue_p and pvalue_n)
    p_value: exponential sum of loglikelihooks of successes greater or equal than observed
    [n, k] --> sum of density (exp of loglikelihood) [n, k] to [n, n]
    '''

    def beta_binom_pvalue(params, target_df):
        n_minus_k = target_df[0] - target_df[1]
        # get the list of observations [n, k] to [n, n]
        obs_list = [target_df + np.array([0,i]) for i in range(0, n_minus_k + 1)]
        # get the list of loglikelihoods per observation
        ll_list = [bb_loglikelihood(params, obs, True) for obs in obs_list]

        #######################################################
        # print(f'ab: {params}\n observations: {obs_list} ll {ll_list}\n')
        #######################################################

        # get the sum of exponentials of loglikelihoods (densities) per observation

        p_value = sum([math.exp(ll) for ll in ll_list])

        return p_value

    target_p = target_df.loc[['depth_p', 'mm_p']]
    target_n = target_df.loc[['depth_n', 'mm_n']]
    p_values = {}
    p_values['p'] = beta_binom_pvalue(params['p'], target_p)
    p_values['n'] = beta_binom_pvalue(params['n'], target_n)
    return p_values

# the matrices for beta-binomial calculation
KS_matrix = np.array([[1,0,1,1,0,1,0,0,0],[0,1,-1,0,1,-1,0,0,0]])
gamma_reduce = np.array([1,-1,-1,-1,1,1,1,-1,-1])

def bb_loglikelihood(params, count_df, is_1d):
    [a, b] = params
    ab_matrix = np.array([1,1,1,a+b,a,b,a+b,a,b])
    # convert df into matrix for np.array operations that change dims
    count_matrix = count_df.values
    # perform matrix multiplication to get inputs to log-gamma
    input_matrix = np.matmul(count_matrix,KS_matrix) + ab_matrix
    # get corresponding log-gamma values and reduce over pon-values
    if is_1d: # check whether gammatrix is 2-dim - otherwise sum aggregation over axis 0 is faulty
        gamma_matrix = gammaln(input_matrix)
    else:  
        gamma_matrix = np.sum(gammaln(input_matrix), axis=0)
    # add or subtract using gamma_reduce matrix and sum to loglikelihood (scalar)
    log_likelihood = np.sum(gamma_matrix * gamma_reduce)
    return log_likelihood
 

def fit_beta_binomial(count_df, pen):
    '''
    Obtaining maximum likelihood estimator of beta-binomial distribution
    count_df is the array of depth-mismatch (trials, success) pairs over the PoN list for either strand
    during minimization of fitting function (max for loglikelihood) penalty term is applied to constrain alpha and beta
        Ref for L-BFGS-B algorithm:
        A Limited Memory Algorithm for Bound Constrained Optimization
        R. H. Byrd, P. Lu and J. Nocedal. , (1995), 
        SIAM Journal on Scientific and Statistical Computing, 16, 5, pp. 1190-1208.
    '''

    def bb_loglikelihood_fitting(params, count_df, penalty):
        '''
        Fitting params [alpha, beta] to maximize loglikelihood
        '''

        # Here, we apply the penalty term of alpha and beta (default 0.5 is slightly arbitray...)
        result = 0.5 * math.log(sum(params)) - bb_loglikelihood(params, count_df, False) # matrix is dim2
        return result
    # get the respective control matrices (as dataframe) for positive and negative strands
    count_p = count_df.loc[:, ['depth_p', 'mm_p']]
    count_n = count_df.loc[:, ['depth_n', 'mm_n']]
    # minimize loglikelihood using L-BFGS-B algorithm
    ab_p = minimize_func(
                           bb_loglikelihood_fitting, [20, 20],
                           args = (count_p, pen), approx_grad = True,
                           bounds = [(0.1, 10000000), (1, 10000000)]
                          )[0]
    ab_n = minimize_func(
                           bb_loglikelihood_fitting, [20, 20],
                           args = (count_n, pen), approx_grad = True,
                           bounds = [(0.1, 10000000), (1, 10000000)]
                          )[0]
    # print(f'abP: {ab_p} - abN: {ab_n}')
    return {'p':ab_p, 'n':ab_n}
