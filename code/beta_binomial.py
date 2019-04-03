import math
import scipy.optimize
from scipy.stats import chi2
from scipy import optimize
import numpy as np
import re
from functools import reduce



def fisher_combination(pvalues):

    if 0 in pvalues:
        return 0
    else:
        return 1 - chi2.cdf(sum([-2 * math.log(x) for x in pvalues]), 2 * len(pvalues))


def beta_binomial_density(params, n, k):

    alpha = params[0]
    beta = params[1]

    tempD = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    tempD = tempD - math.lgamma(n + alpha + beta) + math.lgamma(k + alpha) + math.lgamma(n - k + beta)
    tempD = tempD + math.lgamma(alpha + beta) - math.lgamma(alpha) - math.lgamma(beta) 

    return math.exp(tempD)

def beta_binom_pvalue(params, n, k):

    tempPV = 0
    for kk in range(k, n + 1):

        currentValue = beta_binomial_density(params, n, kk)
        tempPV = tempPV + currentValue


    return tempPV


def beta_binomial_loglikelihood(params, Ns, Ks, var):

    """
    Calculating log-likelihood of beta-binomial distribution
    Args:
        params (List[float]): the parameter of beta distribution ([alpha, beta])  
        As (numpy.array([int])): the counts for success      
        Bs (numpy.array([int])): the counts of trials
    """

    alpha = params[0]    
    beta = params[1]

    ML = 0
    ML += reduce(lambda a, b: a + math.lgamma(b), np.r_[0, Ns + 1])
    ML -= reduce(lambda a, b: a + math.lgamma(b), np.r_[0, Ks + 1])
    ML -= reduce(lambda a, b: a + math.lgamma(b), np.r_[0, Ns - Ks + 1])
    
    ML -= reduce(lambda a, b: a + math.lgamma(b), np.r_[0, Ns + alpha + beta])
    ML += reduce(lambda a, b: a + math.lgamma(b), np.r_[0, Ks + alpha])
    ML += reduce(lambda a, b: a + math.lgamma(b), np.r_[0, Ns - Ks + beta])

    ML += len(Ns) * (math.lgamma(alpha + beta) - math.lgamma(alpha) - math.lgamma(beta))


    # Here, we set the penalty term of alpha and beta (0.5 is slightly arbitray...)
    ML -= 0.5 * math.log(alpha + beta)
    return -ML

 

def fit_beta_binomial(As, Bs, var):

    """
    Obtaining maximum likelihood estimator of beta-binomial distribution
    Args:
        As (numpy.array([int])): the counts for success      
        Bs (numpy.array([int])): the counts of trials
    """
    result = optimize.fmin_l_bfgs_b(beta_binomial_loglikelihood,
                                          [20, 20],
                                          args = (As, Bs, var),
                                          approx_grad = True,
                                          bounds = [(0.1, 10000000), (1, 10000000)])

    return result[0]
