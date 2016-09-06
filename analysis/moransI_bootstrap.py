import re
import numpy as np

contacts = ['sjo', 'nBMB', 'gillies', 'carlia']
dist = [0, 100, 200, 500, 1000, 2000, 3000, 4000, 5000]

print 'contact\tdist\tmoransI\tImin\tImax\tnum'
for contact in contacts:
	file = "/Users/sonal/thesisWork/introgression/LD/LD_%s.out" % contact

	d = {}
	f = open(file, 'r')
	for l in f:
		l = l.rstrip()
		
		t = re.split('\t', l)
		# need to change 6 to 2 for more local version of this
		if t[2] not in d:
			d[t[2]] = []
		d[t[2]].append(t)

	f.close()

	for dist1, dist2 in zip(dist, dist[1:]):
	
		nLoci = 0
		sumNumerator = 0
		numerator = []
		nComparisons = 0
		sumDenominator = 0

		for c in d:
			loci = d[c]
			if len(loci) > 1:
				nLoci += len(loci)  
			
				for array in loci:
					value = float(array[4]) ** 2
					sumDenominator += value

				for ix, array1 in enumerate(loci):
					for array2 in loci[ix + 1:]:
						# switch back to 3 for more local version
						locdist = abs(float(array1[3]) - float(array2[3]))
						if locdist > dist1 and locdist <= dist2:
								nComparisons += 1
								value = float(array1[4]) * float(array2[4])
								sumNumerator += value
								numerator.append(value)

		if nComparisons > 0:
			moranI = (nLoci * sumNumerator) / (nComparisons * sumDenominator)
			bootstrap = []
			for i in range(100):
				original = np.array(numerator)
				resample = np.floor(np.random.rand(len(original))*len(original)).astype(int)
				resampled_sum = np.sum(original[resample])
				resamp_MoransI = (nLoci * resampled_sum) / (nComparisons * sumDenominator)
				bootstrap.append(resamp_MoransI)
			print contact, "\t", dist2, "\t", moranI, "\t", np.min(bootstrap), '\t', np.max(bootstrap), '\t', nComparisons
			
	
