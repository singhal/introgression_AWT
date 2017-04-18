import numpy as np
import random
import matplotlib.pyplot as plt
import re
import math
import pandas as pd
import collections
import scipy.optimize

avg_cov = 100
num_sim = 1000
sample_sizes = [10] + [32] * 9 + [10]
out = "sweep_simulations.csv"

def false_negative(af_diff, in_btn, sample_sizes):
    af1 = 0.5 - (af_diff / 2.0)
    af2 = 0.5 + (af_diff / 2.0)

    exp_af = [af1] + [in_btn] * 9 + [af2]

    real_af = []

    for ix, (af, ss) in enumerate(zip(exp_af, sample_sizes)):
        # sampling error 1
        est_af = np.random.binomial(ss, af, 1)[0] / float(ss)

        # only get sampling error 2 on ends
        if ix > 0 and ix < (len(exp_af) - 1):
            est_af = np.random.binomial(avg_cov, est_af, 1)[0] / float(avg_cov)

        real_af.append(est_af)

    end_diff = False
    if (real_af[-1] - real_af[0]) >= 0.5:
        end_diff = True

    in_between = True
    for i in range(1,10):
        if real_af[i] >= 0.2:
            in_between = False

    if in_between and end_diff:
        return(True)
    else:
        return(False)


def false_positive_constant(af, af_diff, sample_sizes):
    af1 = 0.5 - (af_diff / 2.0)
    af2 = 0.5 + (af_diff / 2.0)

    exp_af = [af1] + [af] * 9 + [af2]

    real_af = []

    for ix, (af, ss) in enumerate(zip(exp_af, sample_sizes)):
        # sampling error 1
        est_af = np.random.binomial(ss, af, 1)[0] / float(ss)

        # only get sampling error 2 on ends
        if ix > 0 and ix < (len(exp_af) - 1):
            est_af = np.random.binomial(avg_cov, est_af, 1)[0] / float(avg_cov)

        real_af.append(est_af)

    end_diff = False
    if (real_af[-1] - real_af[0]) >= 0.5:
        end_diff = True

    in_between = True
    for i in range(1,10):
        if real_af[i] >= 0.2:
            in_between = False

    if in_between and end_diff:
        return(True)
    else:
        return(False)


def af_cline(c, w, x):
    af = (1 + math.tanh(2 * (x - c) / float(w))) / 2.0
    return af


def unscale(af, pmin, pmax):
    '''
    unscale a cline
    based on pmin and pmax
    '''

    new_af = [ freq * (pmax - pmin) + pmin for freq in af]
    return new_af


def rescale(af, pmin, pmax):
    '''
    rescale a cline
    based on pmin and pmax
    '''

    new_af = [(freq - pmin) / float(pmax - pmin) for freq in af]
    return new_af


def false_positive_slope(af_diff, pmin, slope, sample_sizes):
    af1 = pmin
    af2 = pmin + af_diff

    dists = [-2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2]
    cline_af = [af_cline(0, slope, x) for x in dists]
    cline_af = unscale(cline_af, af1, af2)

    exp_af = [af1] + cline_af + [af2]

    real_af = []

    for ix, (af, ss) in enumerate(zip(exp_af, sample_sizes)):
        # sampling error 1
        est_af = np.random.binomial(ss, af, 1)[0] / float(ss)

        # only get sampling error 2 on ends
        if ix > 0 and ix < (len(exp_af) - 1):
            est_af = np.random.binomial(avg_cov, est_af, 1)[0] / float(avg_cov)

        real_af.append(est_af)

    end_diff = False
    if (real_af[-1] - real_af[0]) >= 0.5:
        end_diff = True

    in_between = True
    for i in range(1,10):
        if real_af[i] >= 0.2:
            in_between = False

    if in_between and end_diff:
        return(True)
    else:
        return(False)


def run_true_cline(num_sim):
    for diff in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        for pmin in [0, 0.1, 0.2, 0.3, 0.4, 0.5]:
            if (pmin + diff) <= 1.0:
                for slope in [0.1, 0.5, 1, 2, 5, 10]:
                    sims = [false_positive_slope(diff, pmin, slope, sample_sizes) for i in range(0, num_sim)]
                    rate = sims.count(True) / float(len(sims))
                    o.write('false_positive_cline,%s,%s,%s,%.3f\n' % (diff, pmin, slope, rate))

def run_false(num_sim):
    for diff1 in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
     for diff2 in [0.01, 0.05, 0.1]:
         sims = [false_negative(diff1, diff2, sample_sizes) for i in range(0, num_sim)]
         rate = sims.count(False) / float(len(sims))
         o.write('false_negative,%s,NA,%s,%.3f\n' % (diff1, diff2, rate))

def run_true_constant(num_sim):
    for diff in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:
        for af_diff in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            sims = [false_positive_constant(diff, af_diff, sample_sizes) for i in range(0, num_sim)]
            rate = sims.count(True) / float(len(sims))
            o.write('false_positive_constant,%s,NA,%s,%.3f\n' % (diff, af_diff, rate))


o = open(out, 'w')
o.write('sim_type,diff1,pmin,diff2,rate\n')
run_true_constant(num_sim)
run_true_cline(num_sim)
run_false(num_sim)
o.close()
