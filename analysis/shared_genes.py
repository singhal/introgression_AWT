import re
import glob
import os

contacts = ['carlia', 'gillies', 'nBMB', 'sjo']

exons = {}
genes = {}

for c in contacts:
	file = os.path.join('/Users/Sonal/thesisWork/introgression/targetSequences/final',
		'%s_targets.fa.annotated' % c)
	f = open(file, 'r')

	for l in f:
		if re.search('>', l.rstrip()):
			id = re.search('>(\S+)', l).group(1)
			exon = re.sub('_\d+$', '', id)
			gene = re.sub('_exon\d+$', '', exon)

			if exon not in exons:
				exons[exon] = []
			if gene not in genes:
				genes[gene] = []

			if c not in exons[exon]:
				exons[exon].append(c)
			if c not in genes[gene]:
				genes[gene].append(c)

	f.close()

for exon in exons:
	if len(exons[exon]) > 3:
		print(exon)

#for gene in genes:
#	if len(genes[gene]) > 3:
#		print(gene)