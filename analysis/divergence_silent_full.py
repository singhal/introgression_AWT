import argparse
import re
import subprocess
import os
import glob
import numpy as np

#######################
# doesn't account for unequal base frequences
# doesn't account for unequal codon use (i.e., codon bias)
# doesn't account for multiple mutations at the same codon

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


def site_count(gencode):

	sites = {}
	for codon in gencode:
		codon_sites = []
		for ix, base in enumerate(codon):
			tmp = list(codon)
			aa = []
			for mut in ['A', 'T', 'C', 'G']:
				tmp[ix] = mut
				aa.append(gencode[''.join(tmp)])
			num = len(set(aa))
			codon_sites.append( ( 4 - num) / 3.0)
		sites[codon] = codon_sites

	return sites


def translate(seq, gencode):
	seq = [gencode[''.join(seq[i:i+3])] 
			if ''.join(seq[i:i+3]) in gencode else 'X' 
			for i in range(0, int(len(seq)/3)*3, 3)]
	seq = ''.join(seq)
	return(seq)


def get_frame(seq, c):
	trans = {}

	winner = 0
	start = None
	end = None

	for frame in [0, 1, 2]:
		cds = seq[c]['seq'][(seq[c]['loc'][0] + frame - 1):(seq[c]['loc'][1])]
		aa = translate(cds, gencode)

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
	# vcffile = os.path.join(dir, 'vcf', '%s_test_full.vcf' % contact)
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

				if d[0] not in div:
					div[d[0]] = {}
				div[d[0]][int(d[1])] = {'alleles': [d[3], d[4]], 'geno': [pop1, pop2]}

	return div


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
	'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
	}
sites =  site_count(gencode)
seq = get_seq(contact)
div = get_div_data(contact)
out = os.path.join('/Users/Sonal/thesisWork/introgression/summaryStatistics', '%s.divpoly_silent.csv' % contact)
o = open(out, 'w')
keys = ['pi1', 'pi2', 'btn', 'denom']
o.write('contig,loci,%s\n' % (','.join(keys)))
res = {}

for c in div:
	res[c] = {'denom': 0, 'pi1': 0, 'pi2': 0, 'btn': 0}
	seq = get_frame(seq, c)
	start = seq[c]['loc'][0]
	end = seq[c]['loc'][1]

	cds = seq[c]['seq'][start - 1:end - 1]
	for i in range(start, end, 3):
		rel_pos = i - start
		codon = cds[rel_pos: rel_pos + 3]
		pos = [i, i+1, i+2]

		if codon in gencode:
			res[c]['denom'] += np.sum(sites[codon])
			for ix, p in enumerate(pos):
				if p in div[c]:
					syn = True
					mut = list(codon)
					mut[ix] = div[c][p]['alleles'][1]
					aa1 = gencode[''.join(codon)]
					aa2 = gencode[''.join(mut)]

					if aa1 != aa2:
						syn = False
					
					if syn:
						res[c]['pi1'] += get_pi(div[c][p]['geno'][0])
						res[c]['pi2'] += get_pi(div[c][p]['geno'][1])
						res[c]['btn'] += get_div(div[c][p]['geno'][0], div[c][p]['geno'][1])
	vals = [res[c][key] for key in keys]
	vals = [str(round(val, 4)) for val in vals]				
	o.write('%s,%s,%s\n' % (c, seq[c]['ann'], ','.join(vals)))
o.close()