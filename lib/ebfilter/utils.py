#! /usr/bin/env python

import scipy.stats
import math

def fisher_combination(pvalues):

    if 0 in pvalues:
        return 0
    else:
        return 1 - scipy.stats.chi2.cdf(sum([-2 * math.log(x) for x in pvalues]), 2 * len(pvalues))

