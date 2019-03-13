import msprime
import pyslim
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ts',type=str,help='unmutated tree sequence')
parser.add_argument('posn',type=int,help='position of the site of interest')
parser.add_argument('out',type=str)
args = parser.parse_args()

ts = pyslim.load(args.ts+'.trees')
posn = args.posn
# find the local tree
interval = (-2,-1)
for tree in ts.trees():
	interval = tree.interval
	if (interval[0] <= posn and interval[1] > posn):
		break

# !! assumes only one variant in the ts (i.e. input is before muts dropped!)
genotypes = list(ts.variants())[0].genotypes
print(genotypes)
n = len(genotypes)
codes = {}
for i in range(n):
	for j in range(i+1,n):
		codes[(i,j)] = genotypes[i] + genotypes[j]	

times = {}
memoized_samples = {}
def check_memo_samples(node):
	global tree
	global memoized_samples
	try:
		samples = memoized_samples[node]
	except:
		samples = list(tree.samples(node))
		if len(samples) == 0:
			samples = [node]
		memoized_samples[node] = samples
	return samples	

#print(tree.draw(format="unicode"))
for node in range(n,len(list(tree.nodes()))):

	children = tree.children(node)
	samples = []
	for child in children:
		samples.append(check_memo_samples(child))
	flattened_samples = [item for sublist in samples for item in sublist]
	memoized_samples[node] = flattened_samples
	t = tree.time(node)
	#print(node,children,t,samples,flattened_samples)
	for k in range(len(samples)):
		left_samples = samples[k]
		for l in range(k+1,len(samples)):
			right_samples = samples[l]
			for i in left_samples:
				for j in right_samples:
					times[tuple(sorted((i,j)))] = t	

for key in codes.keys():
	if codes[key] == 2:
		print(times[key])

