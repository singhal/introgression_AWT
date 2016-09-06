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
	seqfile = os.path.join(dir, 'targetSequences/final', '%s_targets.fa.annotated' % contact)
	f = open(seqfile, 'r')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			d = re.split('\t', l.rstrip())
			id = re.search('>(\S+)', d[0]).group(1)
			start = re.search('gs(\d+)', l).group(1)
			end = re.search('ge(\d+)', l).group(1)
			seq[id] = {'loc': [int(start), int(end)], 'seq': ''}
		else:
			seq[id]['seq'] += l.rstrip()

	return seq

def get_pi(alleles):
	alleles = dict([(x, alleles.count(x)) for x in set(alleles)])
	if len(alleles) > 1:
		# total alleles
		n = float(np.sum(list(alleles.values())))
		# minor count
		j = float(np.min(list(alleles.values())))
		pi_prop = (2 * j * (n - j)) / (n * (n - 1))
	elif len(alleles) == 0:
		pi_prop = np.nan
	else:
		pi_prop = 0

	return pi_prop


def get_div(pop1, pop2):
	diff = 0
	num_c = 0

	for i in pop1:
		for j in pop2:
			num_c += 1
			if i != j:
				diff += 1

	if num_c > 0:
		div = diff / float(num_c)
	else:
		div = np.nan

	return div

def get_div_data(contact):
	vcffile = os.path.join(dir, 'vcf', '%s.vcf' % contact)
	f = open(vcffile, 'r')

	bases = ['A', 'T', 'G', 'C']

	div = {}
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l)
			# only consider indels
			if d[4] in bases and d[3] in bases:
				genos = []
				for x in d[9:]:
					geno = re.search('^(\S/\\S)', x).group(1)
					geno = re.split('/', geno)
					genos += geno
				pop1 = [x for x in genos[0:10] if x != '.']
				pop2 = [x for x in genos[10:] if x != '.']

				pi1 = get_pi(pop1)
				pi2 = get_pi(pop2)

				div_seq = get_div(pop1, pop2)

				if d[0] not in div:
					div[d[0]] = {}
				div[d[0]][int(d[1])] = [pi1, pi2, div_seq]

	return div

def summarize(seq, div, contact):
	outfile = os.path.join(dir, 'summaryStatistics2', '%s.contig.divpoly.csv' % contact)
	o = open(outfile, 'w')
	o.write('loci,seq_len,cds_len,seq_snps,'
		    'cds_snps,seq_pi1,seq_pi2,seq_div,cds_pi1,'
		    'cds_pi2,cds_div\n')

	for c in seq:
		seq_len = len(seq[c]['seq'])
		cds = seq[c]['loc'][1] - seq[c]['loc'][0]
		if c in div:
			pi_full1 = np.sum([div[c][x][0] for x in div[c]]) / float(seq_len)
			pi_full2 = np.sum([div[c][x][1] for x in div[c]]) / float(seq_len)
			div_full = np.sum([div[c][x][2] for x in div[c]]) / float(seq_len)
			all_snps = len(div[c])

			pos = [x for x in div[c] if (x >= seq[c]['loc'][0] and x <= seq[c]['loc'][1])]
			cds_snps = len(pos)
			pi1 = np.sum([div[c][x][0] for x in pos]) / float(cds)
			pi2 = np.sum([div[c][x][1] for x in pos]) / float(cds)
			divseq = np.sum([div[c][x][2] for x in pos]) / float(cds)
		else:
			# ack, this is bad
			# there could be no variants
			# because of low coverage
			# or because it is not polymorphic
			# treat as unknown
			pi_full1 = pi1 = np.nan
			pi_full2 = pi2 = np.nan
			div_full = divseq = np.nan
			all_snps = cds_snps = np.nan
		vals = [seq_len, cds, all_snps, cds_snps, 
		        pi_full1, pi_full2, div_full,
		        pi1, pi2, divseq]
		vals = [str(round(x, 4)) for x in vals]
		o.write('%s,%s\n' % (c, ','.join(vals)))
	o.close()

seq = get_seq(contact)
div = get_div_data(contact)
summarize(seq, div, contact)
