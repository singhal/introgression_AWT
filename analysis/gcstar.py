import re
import glob
import subprocess
import multiprocessing as mp
from collections import Counter
import numpy as np

sps = {	'Carlia_rubrigularis_N': 'Carlia_N', 
		'Carlia_rubrigularis_S': 'Carlia_S', 
		'Lampropholis_coggeri_N': 'Lampro_N',
		'Lampropholis_coggeri_C': 'Lampro_C', 
		'Lampropholis_coggeri_S': 'Lampro_S', 
		'Saproscincus_basiliscus_N': 'Sapro_C', 
		'Saproscincus_basiliscus_C': 'Sapro_S', 
		'Pseudemoia_entrecasteauxii': 'Pseudemoia_entrecasteauxii'
}

seqdir = '/Users/Sonal/thesisWork/sutureGenomics/seqfiles/'
outdir = '/Users/Sonal/thesisWork/introgression/gcstar/alns/'

'''
def get_seq(seqfile):
	id = ''
	seq = {}

	f = open(seqfile, 'r')
	for l in f:
		if re.search('>.*', l):
			d = re.split('\t', l.rstrip())
			if len(d) > 1:
				id = d[2]
				info = d[1]
			else:
				id = d[0]
				info = 'NA'
			seq[id] = {'seq': '', 'info': info}
		else:
			seq[id]['seq'] += l.strip()
	f.close()

	return seq

for sp in sps:
	seqfile = '%s%s_trinity.fa.final.annotated' % (seqdir, sps[sp])
	seq = get_seq(seqfile)

	sps[sp] = seq
	print(sp)

ids = []
for sp in sps:
	seq = sps[sp]

	for id in seq:
		if re.search('ENS', id):
			if id not in ids:
				ids.append(id)

for id in ids:
	out = '%s%s.fasta' % (outdir, id)
	o = open(out, 'w')

	for sp in sps:
		if id in sps[sp]:
			seq = list(sps[sp][id]['seq'])

			start = re.search('gs(\d+)', sps[sp][id]['info']).group(1)
			start = int(start) - 1

			end = re.search('ge(\d+)', sps[sp][id]['info']).group(1)
			end = int(end)

			o.write('>%s\n%s\n' % (sp, ''.join(seq[start:end])))
	o.close()
'''


def align(params):
        file, outdir = params

        aln_out = file.replace('.fasta', '.fasta.aln')
        proc = subprocess.call("muscle -in %s -out %s -quiet" %
                               (file, aln_out), shell=True)

        # os.remove(file)
        return aln_out


def run_alignments(outdir):
        files = glob.glob(outdir + '/*fasta')
        
        params = zip(files, [outdir] * len(files))
        
        if cpus > 1:
                pool = mp.Pool(cpus)
                alns = pool.map(align, params)
        
        return alns

def read_aln(seqfile):
	id = ''
	seq = {}

	f = open(seqfile, 'r')
	for l in f:
		if re.search('>.*', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seq[id] = ''
		else:
			seq[id] += l.strip()
	f.close()

	for id, s in seq.items():
		seq[id] = list(s)

	return seq


def get_ancestral(seq):
	loc_len = len(list(seq.values())[0])
	anc_seq = ['N'] * loc_len
	
	legal = ['A', 'G', 'T', 'C', 'a', 't', 'g', 'c' ]

	for i in range(0, loc_len):
		bases = [seq[id][i] for id in seq]
		bases = [base for base in bases if base in legal]
		if len(set(bases)) == 1:
			anc_seq[i] = bases[0]
		else:
			bases = Counter(bases)

			# don't pursue if triallelic
			if len(bases) == 2:
				cbases = sorted(bases, key=bases.get)
	
				mina = bases[cbases[0]]
				maxa = bases[cbases[1]]

				if mina < 3:
					anc_seq[i] = cbases[1]
				elif (maxa - mina) >= 2:
					anc_seq[i] = cbases[1]
				# don't bother otherwise

	seq['ancestral'] = anc_seq
	return seq

def calculate_gcstar(aln, locname):
	types = { 	'A': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
				'C': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
				'T': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
				'G': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'}
				}

	for sp in sps:
		if sp in aln:
			anc = aln['ancestral']
			spseq = aln[sp]

			counts = {'AT_AT': 0, 'AT_GC': 0, 'GC_GC': 0, 'GC_AT': 0}
			ancAT = anc.count('A') + anc.count('T')
			ancGC = anc.count('G') + anc.count('C') 

			for a, b in zip(anc, spseq):
				if a != 'N' and b not in ['-', 'N']:
					counts[types[a][b]] += 1

			if ancAT > 0 and ancGC > 0 and (counts['AT_GC'] + counts['GC_AT']) > 0:
				gcstar = (counts['AT_GC'] / float(ancAT)) / ((counts['AT_GC'] / float(ancAT))+(counts['GC_AT'] / float(ancGC)))
			else:
				gcstar = np.nan
			
			print('%s,%s,%s,%s,%s,%s,%.3f' % (locname, sp, ancAT, ancGC, counts['AT_GC'], counts['GC_AT'], gcstar))
						
#cpus = 4
#run_alignments(outdir)
print('locus,species,ancAT,ancGC,AT_GC,GC_AT,gcstar')
alns = glob.glob(outdir + '*aln')
for aln in alns:
	locname = re.search('(ENS\S+)\.f', aln).group(1)
	# read in alignment
	aln = read_aln(aln)
	# define ancestral sequence
	aln = get_ancestral(aln)
	# calculate gcstar per lineage
	calculate_gcstar(aln, locname)