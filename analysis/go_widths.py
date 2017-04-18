import pandas as pd
import numpy as np
import re
import argparse

parser = argparse.ArgumentParser(description='cline-fitting')
parser.add_argument('--contact', help='the contact for which to run the analysis')
args = parser.parse_args()
con = args.contact

d = pd.read_csv("mart_export.txt")

go = list(d['GO Term Accession'].unique())

c = pd.read_csv('%s.fitting2.csv' % con)
c['gene'] = [re.search('^([^_]+)', x).group(1) for x in c.locus]
c = c[c.transect == 'cline']
c = c.groupby('gene').aggregate(np.mean)

out = '%s.gene_ontology.csv' % con
o = open(out, 'w')
o.write('term,num_genes,go_mean,back_mean,sig,q2.5,q97.5\n')
for term in go:
	genes = list(d.ix[d['GO Term Accession'] == term, 'Ensembl Protein ID'].unique())
	sub1 = c[c.index.isin(genes)]
	sub2 = c[~c.index.isin(genes)]

	sim = []
	for ix in range(0, 1000):
		t = c.sample(sub1.shape[0], replace=False)
		sim.append(t.w.mean())

	go_mean = sub1.w.mean()
	q025, q975 = np.percentile(sim, [2.5, 97.5])
	per = len([x for x in sim if x <= go_mean]) / 1000.
	o.write('%s,%s,%.1f,%.1f,%.3f,%.1f,%.1f\n' % (term, sub1.shape[0], go_mean, 
		   sub2.w.mean(), per, q025, q975))
o.close()