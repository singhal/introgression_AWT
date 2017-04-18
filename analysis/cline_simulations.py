import numpy as np
import random
import matplotlib.pyplot as plt
import re
import math
import pandas as pd
import collections
import scipy.optimize

avg_cov = 100
num_sim = 500
out = open("cline_simulations.csv", 'w')

# ## SIMULATION APPROACH
# 
# 1. Generate a sigmoidal cline under the model - vary width & center
# 2. Sample populations along the cline - vary number of demes & number of inds per deme
# 3. Get allele frequency for a given population 
# 4. Get pooled allele frequency for populations with 50x coverage
# 5. Fit cline
# 
# -- repeat 1000 times for each width & center?
# 
# ## ASSUMPTIONS
# 
# 1. Clines fit a sigmoidal model
# 2. DNA is pooled in equimolar amounts
# 3. Pmin and Pmax are inferred accurately

# In[101]:

def af_cline(c, w, x):
    af = (1 + math.tanh(2 * (x - c) / float(w))) / 2.0
    return af

def get_rid_of_extremes(af):
    '''
    lnl equations can't handle 0 or 1 because of
    log factor, so get rid of those get_rid_of_extreme
    '''

    new = []
    for freq in af:
        if freq < 0.001:
            freq = 0.001
        if freq > 0.999:
            freq = 0.999
        new.append(freq)
    return new


def cline(eval, *params):
    '''
    cline fitting equation to use 
    '''

    # parameters to evaluate
    c = eval[0]
    w = eval[1]

    # fixed parameters
    _af = params[0]
    _dists = params[1]
    _sample_sizes = params[2]
    
    # based on standard sigmoidal cline
    exp_af = [((1 + math.tanh(2* (dist - c) / float(w))) / 2.0) for dist in _dists]
    # get rid of values very close to 0 or 1
    # leads to weird edge behavior
    exp_af = get_rid_of_extremes(exp_af)
    
    # calculate likelihood
    # this is based on Porter from ClineFit
    # Brumfield 2001 appears to be wrong (or I'm misinterpreting)
    pop_lnls = []
    for n, pobs, pexp in zip(_sample_sizes, _af, exp_af):
        pop_lnl =  n * pobs * math.log(pexp) +                     (n * (1 - pobs)) * math.log(1 - pexp)
        pop_lnls.append(pop_lnl)
        
    # needs to return positive so minimization func works
    return -1 *np.sum(pop_lnls)

def rescale(af, pmin, pmax):
    '''
    rescale a cline
    based on pmin and pmax
    '''

    new_af = [(freq - pmin) / float(pmax - pmin) for freq in af]
    return new_af


d = pd.read_csv("/Users/Sonal/thesisWork/introgression/distances/distances", sep="\t")

contacts = ['gillies', 'sjo', 'carlia', 'nBMB']
cinfo = {}

for contact in contacts:
    cfile = '/Users/sonal/thesisWork/introgression/clines/%s.fitting2.csv' % contact
    c = pd.read_csv(cfile)
    c = c[c.transect == 'cline']

    meanw = c.w.mean()
    percents = [meanw] + list(np.percentile(c.w, [25, 50, 75]))
    center = c.c.median()
    width = percents[2]
   
    ancN = -40e3
    center = 0
    kS1 = 1e3
    kN1 = -1e3
    kS25 = 2.5e3
    kN25 = -2.5e3
    kS10 = 10e3
    kN10 = -10e3
    nTail = (-0.1 * width) - center
    sTail = (0.1 * width) - center
    ancS = 40e3
    
    locs = [ancN, kN10, kN25, kN1, nTail, center, sTail,
            kS1, kS25, kS10, ancS]
    
    cinfo[contact] = {'dist': locs, 's': [10] + [32] * 9 + [10], 
                      'widths': percents}


wtypes = ['mean', 'per25', 'median', 'per75']
out.write("cline,width_type,width,num_sim,center,width\n")
for c in cinfo:
    for ix, width in enumerate(cinfo[c]['widths']):

        for i in range(num_sim):
            dists = cinfo[c]['dist']
            samps = cinfo[c]['s']
            af = []
            for dist, samp in zip(dists, samps):
                exp_af = af_cline(0, width, dist)
                # sampling error 1
                est_af1 = np.random.binomial(samp, exp_af, 1)[0] / float(samp)
                
                # sampling error 2
                est_af2 = np.random.binomial(avg_cov, est_af1, 1)[0] / float(avg_cov)
                
                af.append(est_af2)
            
            af1 = af[2:9]
            pmin = np.min(af[0:2])
            pmax = np.min(af[10:12])
            dist1 = dists[2:9]
            
            d = dict(zip(dist1, af1))
            od = collections.OrderedDict(sorted(d.items()))
            cline_af = [y for x, y in od.items()]
            cline_dists = [x for x, y in od.items()]
            
            cline_af = get_rid_of_extremes(cline_af)
            cline_af = rescale(cline_af, pmin, pmax)
            sample_sizes = samps[2:9]
            
            w_start = 100
            w_end = 20000
            w_slice = int((w_end - w_start) / 100.)
            
            ranges = (slice(100, 5000, 25), slice(w_start, w_end, w_slice))
            
            fixed_params = (cline_af, cline_dists, sample_sizes)
            res1 = scipy.optimize.brute(cline, ranges, args=fixed_params, full_output=True)
            
            out.write('%s,%s,%.2f,%s,%.2f,%.2f\n' % (c, wtypes[ix], width, i, res1[0][0], res1[0][1]))
out.close()


