#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 23:48:28 2020

@author: swallan
"""
import numpy as np
from itertools import permutations



from scipy.stats import norm
from scipy.integrate import quad
from numpy.testing import assert_allclose
from scipy.special import gamma

phi = norm.pdf
Phi = norm.cdf

def F_R(q, k, nu=np.inf):
    """CDF of studentized range distribution"""
    # Direct implementation of formulas from Wikipedia
    # https://en.wikipedia.org/wiki/Studentized_range_distribution
    if nu > 40:
        # asymptotic version
        # results may not be very accurate for 40 < nu < 120; 
        def integrand(z):
            return phi(z)*(Phi(z+q)-Phi(z))**(k-1)
        F, err = quad(integrand, -np.inf, np.inf, epsabs=1e-12, limit=1000)
        return k*F

    def outer(s):
        def inner(z):
            return phi(z)*(Phi(z+q*s)-Phi(z))**(k-1)
        inner_int = quad(inner, -np.inf, np.inf, epsabs=1e-12, limit=1000)[0]
        return s**(nu-1)*phi(np.sqrt(nu)*s)*inner_int
    
    outer_int = quad(outer, 0, np.inf, epsabs=1e-12, limit=1000)[0]
    return np.sqrt(2*np.pi)*k*nu**(nu/2) / (gamma(nu/2)*2**(nu/2-1))*outer_int

def wrapper_cdf(q_crit, treatments, ddof, alpha):
    return alpha - (1 - F_R(q_crit, treatments, ddof))

def stud_range_crit(alpha, n, ddof):
    soln = root(wrapper_cdf, 2, args=(n, ddof, alpha))
    if soln.success:
        return soln.x
    else: 
        raise ValueError(str(soln.message))
fromwrapper = np.zeros([31,8])

for i in range(31):
    for j in range(8):
        x = root(wrapper_cdf, 3, args=(i, j, .05))
        print(x)
        fromwrapper[31:8] = x
        

from scipy.optimize import root

studentized_range = [
     [26.976, 32.819, 37.081, 40.407, 43.118, 45.397, 47.356, 49.07],
     [8.331, 9.798, 10.881, 11.734, 12.434, 13.027, 13.538, 13.987],
     [5.91, 6.825, 7.502, 8.037, 8.478, 8.852, 9.177, 9.462],
     [5.04, 5.757, 6.287, 6.706, 7.053, 7.347, 7.602, 7.826],
     [4.602, 5.218, 5.673, 6.033, 6.33, 6.582, 6.801, 6.995],
     [4.339, 4.896, 5.305, 5.629, 5.895, 6.122, 6.319, 6.493],
     [4.165, 4.681, 5.06, 5.359, 5.606, 5.815, 5.997, 6.158],
     [4.041, 4.529, 4.886, 5.167, 5.399, 5.596, 5.767, 5.918],
     [3.948, 4.415, 4.755, 5.024, 5.244, 5.432, 5.595, 5.738],
     [3.877, 4.327, 4.654, 4.912, 5.124, 5.304, 5.46, 5.598],
     [3.82, 4.256, 4.574, 4.823, 5.028, 5.202, 5.353, 5.486],
     [3.773, 4.199, 4.508, 4.748, 4.947, 5.116, 5.262, 5.395],
     [3.734, 4.151, 4.453, 4.69, 4.884, 5.049, 5.192, 5.318],
     [3.701, 4.111, 4.407, 4.639, 4.829, 4.99, 5.13, 5.253],
     [3.673, 4.076, 4.367, 4.595, 4.782, 4.94, 5.077, 5.198],
     [3.649, 4.046, 4.333, 4.557, 4.741, 4.896, 5.031, 5.15],
     [3.628, 4.02, 4.303, 4.524, 4.705, 4.858, 4.991, 5.108],
     [3.609, 3.997, 4.276, 4.494, 4.673, 4.824, 4.955, 5.071],
     [3.593, 3.977, 4.253, 4.468, 4.645, 4.794, 4.924, 5.037],
     [3.578, 3.958, 4.232, 4.445, 4.62, 4.768, 4.895, 5.008],
     [3.565, 3.942, 4.213, 4.424, 4.597, 4.743, 4.87, 4.981],
     [3.553, 3.927, 4.196, 4.405, 4.577, 4.722, 4.847, 4.957],
     [3.542, 3.914, 4.18, 4.388, 4.558, 4.702, 4.826, 4.935],
     [3.532, 3.901, 4.166, 4.373, 4.541, 4.684, 4.807, 4.915],
     [3.523, 3.89, 4.153, 4.358, 4.526, 4.667, 4.789, 4.897],
     [3.514, 3.88, 4.141, 4.345, 4.511, 4.652, 4.773, 4.88],
     [3.506, 3.87, 4.13, 4.333, 4.498, 4.638, 4.758, 4.864],
     [3.499, 3.861, 4.12, 4.322, 4.486, 4.625, 4.745, 4.85],
     [3.493, 3.853, 4.111, 4.311, 4.475, 4.613, 4.732, 4.837],
     [3.487, 3.845, 4.102, 4.301, 4.464, 4.601, 4.72, 4.824],
     [3.481, 3.838, 4.094, 4.292, 4.454, 4.591, 4.709, 4.813]
     ]

# computation of the Tukey Criterion
def stud_range(sig_lev, num_treatments, ddof):
    print(num_treatments, ddof)
    """
    for studentized range distribution, we could try to do the integral,
    or use table temporarily.

    # -1 for row since it is zero indexed to get ddof = 1
    # -3 for column since it starts at 3 groups
    """
    return studentized_range[ddof - 1][num_treatments - 3]

class TukeyKramerResult():
    def __init__(self):
        self.comp = dict()
        
    def __repr__(self):
        s =f"{'group' : <7}{'min' : ^10}{'mean' : ^10}{'max' : >5}{'sig' : >9}\n"
        for key in self.comp.keys():
            s = s + f"{key[0]} - {key[1]}  " \
                f"{self.comp[key][0] :>7.3f}{self.comp[key][1] : >10.3f}"\
                f"{self.comp[key][2] :>9.3f}{self.comp[key][3] : >7}\n"
        return s
    

def tukeykramer(*args, sig_level=.05):
    """
    Perform Tukey-Kramer. The order of the input arrays is

    Under the assumption of a rejected null hypothesis that two or more samples
    have the same population mean, the Tukey Kramer test compares the absolute
    mean difference between each input group to the Tukey criterion for
    significance. For inputs of differing sample sizes, the Tukey Kramer
    method is used, albeit with a higher confidence coefficient.

    Assumptions:
    1. There are equal variances between the samples.
    2. The samples are equally distributed


    Parameters
    ----------
    sample1, sample2, ... : array_like
        The sample measurements for each group.  There must be at least
        two arguments.
    sig_level : float, optional
        Significance level for which to determine the appropriate studentized
        range distribution value.
        Default is .05.


    SAS example, adapted from
    https://www.stattutorials.com/SAS/TUTORIAL-PROC-GLM.htm
    ```
    DATA ACHE;
    INPUT BRAND RELIEF;
    CARDS;
    1 24.5
    1 23.5
    1 26.4
    1 27.1
    1 29.9
    2 28.4
    2 34.2
    2 29.5
    2 32.2
    2 30.1
    3 26.1
    3 28.3
    3 24.3
    3 26.2
    3 27.8
    ;
    ODS RTF;ODS LISTING CLOSE;
    PROC ANOVA DATA=ACHE;
        CLASS BRAND;
        MODEL RELIEF=BRAND;
        MEANS BRAND/TUKEY CLDIFF;
    TITLE 'COMPARE RELIEF ACROSS MEDICINES  - ANOVA EXAMPLE';
    ```

    Additional Info: the critical value of the Studentized Range is 3.77289.

    Yeilds the following results :
    (*** denotes significance at the .05 significance level)
    Brand Comparison    Mean Difference     Simultanious 95% Confidence Limits
    2 - 3               4.340               0.691           7.989   ***
    2 - 1               4.600               0.951           8.249   ***
    3 - 2              -4.340              -7.989          -0.691   ***
    3 - 1               0.260              -3.389           3.909
    1 - 2              -4.600              -8.249          -0.951   ***
    1 - 3              -0.260              -3.909           3.389

    Let's test our function against the SAS result:
     group1 = [24.5, 23.5, 26.4, 27.1, 29.9]
     group2 = [28.4, 34.2, 29.5, 32.2, 30.1]
     group3 = [26.1, 28.3, 24.3, 26.2, 27.8]
    
    
    # second just different sizes
    group1 = [24.5, 23.5, 26.28, 26.4, 27.1, 29.9, 30.1, 30.1]
    group2 = [28.4, 34.2, 29.5, 32.2, 30.1]
    group3 = [26.1, 28.3, 24.3, 26.2, 27.8]
    args=[group1, group2, group3]

    from https://www.youtube.com/watch?v=zQr190cacC0&t=124s
    a = [77, 79, 87, 85, 78]
    b = [83, 91, 94, 88, 85]
    c = [80, 82,  86, 85, 80]

    """
    #%%
    args = [np.asarray(arg) for arg in args]
    means = [np.mean(arg) for arg in args]
    lengths = [len(arg) for arg in args]

    """
    Critical Values of Studentized Range Distribution(q) for
    Familywise ALPHA = .05.
    Rows: DF (index 0 = 1 DF)
    Cols: Number of Groups (a.k.a. Treatments) (index 0 = 3 groups)
    source:
    https://www.stat.purdue.edu/~lingsong/teaching/2018fall/q-table.pdf
    """


    # determine the studentized range distribution value
    # sig_level: always .05 right now.
    # Number of treatments is number of args.
    # DF: Number of datapoints minus number of treatments
    srd = stud_range(sig_level, len(args),
                     len(args[0]) * len(args) - len(args))

    # determine mean square error
    mse = np.mean([np.var(arg, ddof=1) for arg in args])

    # also called maxmimum critical value, tukey criterion is the studentized
    # range value * the square root of mean square error over the sample size.
    # This only applies for treatments of equal sizes. The criterion must be
    # calculated for each comparison when treatments differ in size, with th
    permutations_lengths = list(permutations(lengths, 2))
    size_relation = np.asarray([(1/n_i + 1/n_j) for (n_i, n_j) in permutations_lengths])
    tukey_criterion = srd * ((mse * size_relation / 2) ** .5)

    # create permutations of input group means along with keys to keep track
    # of which is which
    permutations_key = list(permutations(np.arange(1, len(means) + 1), 2))
    permutations_means = list(permutations(means, 2))

    # determine the mean difference
    mean_differences = np.asarray([x[0] - x[1] for x in permutations_means])

    # the confidence levels are determined by the mean diff +- tukey_criterion
    conf_levels = (mean_differences - tukey_criterion,
                   mean_differences + tukey_criterion)

    # The simultaneous pairwise comparisons are not significantly different
    # from 0 if their confidence intervals include 0. 
    # ("Conclusions")[https://www.itl.nist.gov/div898/handbook/prc/section4/prc471.htm]
    is_significant = np.abs(mean_differences) > tukey_criterion

    res = TukeyKramerResult()
    for i in range(len(mean_differences)):
        res.comp[permutations_key[i]] = [conf_levels[0][i],
                                         mean_differences[i],
                                         conf_levels[1][i], is_significant[i]]
    #return res

    #%%
# group1 = [24.5, 23.5, 26.4, 27.1, 29.9]
# group2 = [28.4, 34.2, 29.5, 32.2, 30.1]
# group3 = [26.1, 28.3, 24.3, 26.2, 27.8]
# print("""
# expected: 
# Brand Comparison    Mean Difference     Simultanious 95% Confidence Limits
# 2 - 3               4.340               0.691           7.989   ***
# 2 - 1               4.600               0.951           8.249   ***
# 3 - 2              -4.340              -7.989          -0.691   ***
# 3 - 1               0.260              -3.389           3.909
# 1 - 2              -4.600              -8.249          -0.951   ***
# 1 - 3              -0.260              -3.909           3.389
      
#       """)
# print(tukeykramer(group1, group2, group3))

# # second just different sizes
# group1 = [24.5, 23.5, 26.28, 26.4, 27.1, 29.9, 30.1, 30.1]
# group2 = [28.4, 34.2, 29.5, 32.2, 30.1]
# group3 = [26.1, 28.3, 24.3, 26.2, 27.8]

# print("""
# expected: 
# Brand Comparison    Mean Difference     Simultanious 95% Confidence Limits
# 2 - 3               3.645               0.268        	7.022   ***
# 2 - 1               4.34                0.593        	8.087   ***
# 3 - 2              -3.645              -7.022      	   -0.268   ***
# 3 - 1               0.695              -2.682         	4.072
# 1 - 2              -4.34               -8.087	       -0.593   ***
# 1 - 3              -0.695              -4.072	        2.682
      
#       """)

# print(tukeykramer(group1, group2, group3))