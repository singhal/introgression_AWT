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


def get_fst(pop1, pop2):
	a1_1 = pop1.count('0') / float(len(pop1))
	a2_1 = pop1.count('1') / float(len(pop1))
	a1_2 = pop2.count('0') / float(len(pop2))
	a2_2 = pop2.count('1') / float(len(pop2))

	het1 = 2 * a1_1 * a2_1
	het2 = 2 * a1_2 * a2_2
	avg_het = (het1 + het2) / 2.0

	a1 = (pop1.count('0') + pop2.count('0')) / float(len(pop1) + len(pop2))
	a2 = (pop1.count('1') + pop2.count('1')) / float(len(pop1) + len(pop2))
	het = 2 * a1 * a2

	if het > 0:
		fst = (het - avg_het) / het
	else:
		fst = 0

	return fst


def get_div_data(contact):
	vcffile = os.path.join(dir, 'vcf', '%s_full.vcf' % contact)
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

				fst = get_fst(pop1, pop2)

				if fst > 0:
					if d[0] not in div:
						div[d[0]] = {}
					div[d[0]][int(d[1])] = fst

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
			seq_fst = np.mean([div[c][x] for x in div[c]])
			all_snps = len(div[c])

			pos = [x for x in div[c] if (x >= seq[c]['loc'][0] and x <= seq[c]['loc'][1])]
			cds_snps = len(pos)
			fst = [div[c][x] for x in pos]
			if len(fst) > 0:
				cds_fst = np.mean(fst)
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
