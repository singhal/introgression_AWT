import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.optimize
import math
from scipy.stats import pearsonr
import argparse


parser = argparse.ArgumentParser(description='cline-fitting')
parser.add_argument('--contact', help='the contact for which to run the analysis')
args = parser.parse_args()
contact = args.contact

# suggested starting parameters
start = pd.read_csv('/Users/sonal/thesisWork/introgression/clineAF/starting', sep='\t')
# allele frequency data
d = pd.read_csv('/Users/sonal/thesisWork/introgression/clineAF/%s.cline.out' % contact, 
                sep='\t', names=['locus', 'position', 'population', 'af'])
# distances on transect
dist = pd.read_csv('/Users/sonal/thesisWork/introgression/distances/distances', sep='\t')
out = open('/Users/sonal/thesisWork/introgression/clines/%s.fitting_new.csv' % contact, 'w')


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
        pop_lnl =  n * pobs * math.log(pexp) + \
                    (n * (1 - pobs)) * math.log(1 - pexp)
        pop_lnls.append(pop_lnl)
        
    # needs to return positive so minimization func works
    return -1 *np.sum(pop_lnls)

  
def get_cline(c, w, dists, pmin, pmax):
    '''
    get cline values to use for plotting
    '''

    # sigmoidal cline
    exp_af = [((1 + math.tanh(2* (dist - c) / float(w))) / 2.0) for dist in dists]
    # scale it to pmin and pmax
    exp_af = [pmin + x * (pmax - pmin) for x in exp_af]  
    
    return exp_af


def rescale(af, pmin, pmax):
    '''
    rescale a cline
    based on pmin and pmax
    '''

    new_af = [(freq - pmin) / float(pmax - pmin) for freq in af]
    return new_af


# first get rid of control loci
d = d.drop(d[d.locus.str.contains('16s') == True].index)
d = d.drop(d[d.locus.str.contains('utr') == True].index)
d = d.drop(d[d.locus.str.contains('ND4') == True].index)

out.write('contact,transect,locus,pos,num_pops,num_cline_pops,pmin,pmax,anc_diff,tenk_diff,trans_diff,c,w,lnl\n')
# go through each snp one by one
d = d.groupby(['locus', 'position'])
for (locus, pos), cur_d in d:
    # note keeping
    printed = False
    anc_diff = np.nan
    trans_diff = np.nan
    tenk_diff = np.nan
    pmin = np.nan
    pmax = np.nan
    
    # merge to get distance values for pops
    cur_d = cur_d.merge(dist[dist.contact == contact], left_on='population', right_on='library')
    cur_d = cur_d.dropna()
    # some weird behavior if not sorted, so sort
    cur_d = cur_d.sort(['distance'])
    # create a data frame with transect pops only
    trans = cur_d[(cur_d.population != 'ancN') & (cur_d.population != 'ancS')]
   
    # identify gradient in allele frequencies
    # across full transect, across ancestors, across 10k pops
    pops = dict([(x, y) for x, y in zip(cur_d.population, cur_d.af)])
    if ('ancN' in pops) and ('ancS' in pops):
        anc_diff = pops['ancS'] - pops['ancN']
    if ('10kS' in pops) and ('10kN' in pops):
        tenk_diff =  pops['10kS'] - pops['10kN']

    # calculate max transition between two ends
    endN = [ pops['10kN'], pops['2kN'], pops['1kN'], pops['nTail'] ] 
    endS = [ pops['10kS'], pops['2kS'], pops['1kS'], pops['sTail'] ]
    endN = [x for x in endN if np.isfinite(x)]
    endS = [x for x in endS if np.isfinite(x)]
    if len(endN) > 1 and len(endS) > 1:
        # don't know direction yet so don't know up from down
        trans_diff1 = abs(min(endN) - max(endS))
        trans_diff2 = abs(max(endN) - min(endS))
        trans_diff = max([trans_diff1, trans_diff2])
    else:
        trans_diff = np.nan        
    
    diffs = pd.Series([anc_diff, trans_diff, tenk_diff])
    diffs = diffs[~np.isnan(diffs)]
    clines = trans[(trans.population != '10kN') & (trans.population != '10kS')]
             
    # can't calculate differences
    # and / or there are only 4 sampled
    # af in the center
    # well what's the point then?
    if len(diffs) > 0 and trans.shape[0] >= 5:
        max_abs_diff = max([abs(x) for x in diffs])
        
        # not going to try to characterize something
        # if af difference across all pops is less than 0.4
        if max_abs_diff >= 0.5:
            # need to flip?
            flip = False
            if math.isnan(tenk_diff):
                slope = np.polyfit(trans.distance, trans.af, 1)[0]
                if slope < 0:
                    flip = True
            else:
                # these need to be flipped
                if tenk_diff < 0:
                      flip = True
            
            if flip:
                trans.af = 1 - trans.af
                cur_d.af = 1 - cur_d.af
            
            # in a true sweep
            # no one should have big pops    
            sweep = True
            for ix, x in trans.iterrows():
                if x['af'] > 0.2 and x['af'] < 0.8:
                    sweep = False          
    
            # sweep
            if anc_diff >= 0.5 and trans_diff < 0.2 and sweep:
                printed = True
                out.write('%s,sweep,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % \
                    (contact, locus, pos, trans.shape[0], clines.shape[0], pmin, pmax, abs(anc_diff), 
                    abs(tenk_diff), abs(trans_diff), np.nan, np.nan, np.nan))
            # try to fit cline!
            # only fit it if at least 5 pops in cline
            elif trans_diff >= 0.5 and clines.shape[0] >= 5:
                pmin = trans.af.min()
                pmax = trans.af.max()
                
                # don't want to fit 10k pops
                # they are too far off the transect
                # so from here on out use clines
                
                # get rid of extreme values
                af = get_rid_of_extremes(clines.af.tolist())
                af = rescale(af, pmin, pmax)
                
                dists = clines.distance.tolist()
                # get sample sizes
                sample_sizes = [16] * len(dists)

                # fitting procedure prone to getting stuck in local minima
                # so will use brute force approach
                # but to minimize time will set parameters
                # for search based on a priori sense 
                # of starting values
                w_start = int(start[(start.contact == contact) & (start.measure == 'w')].value)
                w_end = w_start * 40
                if w_end > 20000:
                    w_end = 20000
                w_start = int(w_start / 2.0)
                w_slice = int((w_end - w_start) / 500.)

                # define search space for brute
                # note that final values can be outside of this
                # because of "finishing"
                ranges = (slice(100, 5000, 25), slice(w_start, w_end, w_slice))
                fixed_params = (af, dists, sample_sizes)

                # fit ML cline
                res1 = scipy.optimize.brute(cline, ranges, args=fixed_params, full_output=True)
               
                printed = True
                out.write('%s,cline,%s,%s,%s,%s,%s,%s,%s,%s,%s,%.1f,%.1f,%.1f\n' % \
                                (contact, locus, pos, trans.shape[0], clines.shape[0], pmin, pmax, abs(anc_diff), 
                                abs(tenk_diff), abs(trans_diff), res1[0][0], res1[0][1], res1[1]))
                    
    if not printed:
        out.write('%s,no_fit,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % \
                    (contact, locus, pos, trans.shape[0], clines.shape[0], pmin, pmax, 
                    abs(anc_diff), abs(tenk_diff), abs(trans_diff), np.nan, np.nan, np.nan))
out.close()
