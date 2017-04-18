import argparse
import re
import subprocess
import os
import glob
import numpy as np

dir = '/Users/sonal/thesisWork/introgression/'

parser = argparse.ArgumentParser(description='cline-fitting')
parser.add_argument('--contact', help='the contact for which to run the analysis')
args = parser.parse_args()
contact = args.contact

def get_seq(contact):
	seqfile = os.path.join(dir, 'fullTranscripts', '%s.fa' % contact)
	f = open(seqfile, 'r')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			d = re.split('\t', l.rstrip())
			id = re.search('>(\S+)', d[0]).group(1)
			start = re.search('gs(\d+)', l).group(1)
			end = re.search('ge(\d+)', l).group(1)
			seq[id] = {'loc': [int(start), int(end)], 'seq': '',
			           'ann': d[2]}
		else:
			seq[id]['seq'] += l.rstrip()

	return seq


def fst_estimator(counts, sample_sizes):
	'''
	modified from G. Bradford's R code in bedassle
	'calculate.pairwise.Fst'

	both inputs are arrays where each row is an individual
	and each column is a SNP
	'''

	counts = np.array(counts)
	sample_sizes = np.array(sample_sizes).astype('float')

	pop_af = counts / sample_sizes
	mean_af = np.sum(counts, axis=0) / np.sum(sample_sizes, axis = 0)

	MSP = np.sum((pop_af - mean_af) ** 2 * sample_sizes, axis=0)
	MSG = np.sum((1 - pop_af) * pop_af * sample_sizes, axis=0) \
                * (1 / np.sum(sample_sizes - 1, axis=0))
	n_c = np.sum(sample_sizes, axis = 0) - np.sum(sample_sizes ** 2, axis=0) \
	                / np.sum(sample_sizes, axis=0)

	fst = np.sum(MSP - MSG) / np.sum(MSP + (n_c - 1) * MSG)

	return fst


def get_div_data(contact):
	vcffile = os.path.join(dir, 'vcf', '%s_full.vcf' % contact)
	f = open(vcffile, 'r')

	bases = ['A', 'T', 'G', 'C']

	div = {}
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l)
			# only consider non-indels, biallelics
			if d[4] in bases and d[3] in bases:
				genos = []
				for x in d[9:]:
					geno = re.search('^(\S/\\S)', x).group(1)
					geno = re.split('/', geno)
					genos += geno
				pop1 = [x for x in genos[0:10] if x != '.']
				pop2 = [x for x in genos[10:] if x != '.']

				if len(pop1) > 0 and len(pop2) > 0:
					if d[0] not in div:
						div[d[0]] = {}
					div[d[0]][int(d[1])] = {'geno': [pop1.count('1'), pop2.count('1')],
					                        'ss': [len(pop1), len(pop2)] 
					                        }

	return div


def summarize(seq, div, contact):
	outfile = os.path.join(dir, 'summaryStatistics2', '%s.full.fst.csv' % contact)
	o = open(outfile, 'w')
	o.write('loci,contig,seq_len,cds_len,seq_snps,'
		    'cds_snps,seq_fst,cds_fst\n')

	for c in seq:
		seq_len = len(seq[c]['seq'])
		cds = seq[c]['loc'][1] - seq[c]['loc'][0]
		if c in div:
			all_snps = len(div[c])
			counts = [[],[]]
			ss = [[],[]]
			for pos in div[c]:
				counts[0].append(div[c][pos]['geno'][0])
				counts[1].append(div[c][pos]['geno'][1])
				ss[0].append(div[c][pos]['ss'][0])
				ss[1].append(div[c][pos]['ss'][1])
			if len(counts[0]) > 0:
				seq_fst = fst_estimator(counts, ss)
			else:
				seq_fst = 0

			cds_pos = [x for x in div[c] if (x >= seq[c]['loc'][0] and x <= seq[c]['loc'][1])]
			cds_snps = len(cds_pos)
			counts = [[],[]]
			ss = [[],[]]
			for pos in cds_pos:
				counts[0].append(div[c][pos]['geno'][0])
				counts[1].append(div[c][pos]['geno'][1])
				ss[0].append(div[c][pos]['ss'][0])
				ss[1].append(div[c][pos]['ss'][1])
			if len(counts[0]) > 0:
				cds_fst = fst_estimator(counts, ss)
			else:
				cds_fst = 0
		else:
			# ack, this is bad
			# there could be no variants
			# because of low coverage
			# or because it is not polymorphic
			# treat as unknown
			cds_fst = seq_fst = np.nan
			all_snps = cds_snps = np.nan
		vals = [seq_len, cds, all_snps, cds_snps, 
		        seq_fst, cds_fst]
		vals = [str(round(x, 4)) for x in vals]
		o.write('%s,%s,%s\n' % (seq[c]['ann'], c, ','.join(vals)))
	o.close()


seq = get_seq(contact)
div = get_div_data(contact)
summarize(seq, div, contact)
