import argparse
import re
import subprocess
import os
import glob
import numpy as np
import sys

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


def translate(seq):
	gencode = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

	seq = [gencode[''.join(seq[i:i+3])] if ''.join(seq[i:i+3]) in gencode else 'X' for i in range(0, int(len(seq)/3)*3, 3)]
	seq = ''.join(seq)

	return seq


def get_frame(seq, c):
	trans = {}

	winner = 0
	start = None
	end = None

	for frame in [0, 1, 2]:
		cds = seq[c]['seq'][(seq[c]['loc'][0] + frame - 1):(seq[c]['loc'][1])]
		aa = translate(cds)

		matches = re.findall('([^\*]+)', aa)
		longest_aa = ''
		for m in matches:
			if len(m) > len(longest_aa):
				longest_aa = m

		if len(longest_aa) > winner:
			winner = len(longest_aa)
			cds_span = re.search(longest_aa, aa).span()
			start = cds_span[0] * 3
			end = cds_span[1] * 3
			
	seq[c]['loc'] = [(seq[c]['loc'][0] + start), (seq[c]['loc'][0] + end)]

	return seq


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
				pop1 = [int(x) for x in genos[0:10] if x != '.']
				pop2 = [int(x) for x in genos[10:] if x != '.']
				alleles = [d[3], d[4]]

				pop1 = alleles[max(set(pop1), key=pop1.count)]
				pop2 = alleles[max(set(pop2), key=pop2.count)]

				if d[0] not in div:
					div[d[0]] = {}

				div[d[0]][int(d[1])] = [pop1, pop2]

	return div


def make_paml_ctl(contact):
	ctlfile = '%s.ctl' % contact
	o = open(ctlfile, 'w')
	o.write("seqfile = %s.nuc\n" % contact)
	o.write("outfile = %s_out\n" % contact)
	o.write("verbose = 0\n")
	o.write("icode = 0\n")
	o.write("weighting = 0\n")
	o.write("commonf3x4 = 0\n")
	o.close()
	return ctlfile

def run_paml(seq, div):
	final = os.path.join(dir, 'summaryStatistics2', '%s.full.dnds.csv' % contact)
	final = open(final, 'w')
	final.write('contig,loci,dn,ds,dnds\n')
	for c in seq:
		if c in div:
			seq = get_frame(seq, c)
			seq1 = list(seq[c]['seq'])
			seq2 = list(seq[c]['seq'])
			for pos in div[c]:
				seq1[pos - 1] = div[c][pos][0]
				seq2[pos - 1] = div[c][pos][1]
			seq1 = seq1[(seq[c]['loc'][0] - 1):(seq[c]['loc'][1] - 1)]
			seq2 = seq2[(seq[c]['loc'][0] - 1):(seq[c]['loc'][1] - 1)]

			o = open('%s.nuc' % contact, 'w')
			o.write(' 2 %s\n' % (len(seq1)))
			o.write('seq1  %s\n' % (''.join(seq1)))
			o.write('seq2  %s\n' % (''.join(seq2)))
			o.close()

			subprocess.call("/usr/local/bin/paml4.8/bin/yn00 %s.ctl" % contact, shell=True)
			i = open('%s_out' % contact, 'r')
			i = i.readlines()
			wrote = False
			for ix, l in enumerate(i):
				if re.search('omega', l):
					d = i[ix + 2]
					d = re.split('\s+', d)
					dN = abs(float(d[8]))
					dS = float(d[11])
					if dS > 0:
						ratio = dN / dS
					else:
						ratio = 0
					final.write('%s,%s,%s,%s,%s\n' % (c, seq[c]['ann'], dN, dS, ratio))
					wrote = True
			if not wrote:
				final.write('%s,%s,%s,%s,%s\n' % (c, seq[c]['ann'], np.nan, np.nan, np.nan))
		else:
			final.write('%s,%s,%s,%s,%s\n' % (c, seq[c]['ann'], np.nan, np.nan, np.nan))
	final.close()

seq = get_seq(contact)
div = get_div_data(contact)
ctlfile = make_paml_ctl(contact)
run_paml(seq, div)