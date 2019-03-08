import msprime
import pyslim
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('out',type=str)
parser.add_argument('-q','--q',action='store_true')
parser.add_argument('-n','--n',type=int,default=2000,help='HAPLOID sample size')
parser.add_argument('-u','--u',type=float,default=1e-8,help='mut rate per bp per gen')
parser.add_argument('-r','--r',type=float,default=1e-8,help='per bp recomb rate')
parser.add_argument('-l','--l',type=int,default=2e5,help='length of locus in bp')
parser.add_argument('-c','--c',type=float,default=0.01,help='minimum MAF of standing var')
parser.add_argument('-s','--s',type=float,default=0,help='selection coefficient. only use under Mendelian Inheritance model; do not use this option with --pg!')
parser.add_argument('--pg',type=float,nargs=2,default=None,help='Use for polygenic selection. arg1 is I, the selection differential on trait; arg2 is k, the num loci. s is random (see e.g. Edge&Coop 2018)')
args = parser.parse_args()

def throw_mut_on_tree(ts):
	global args

	n = args.n
	l = args.l
	r = args.r
	q = args.q
	c = args.c
	
	# find total tree length times sequence extent
	tree_sizes = np.array([t.total_branch_length * (np.ceil(t.interval[1]) - np.ceil(t.interval[0])) for t in ts.trees()])
	tree_sizes /= sum(tree_sizes)

	# pick the tree
	tree_index = np.random.choice(ts.num_trees, size=1, p=tree_sizes)
	t = ts.first()
	for (i,t) in enumerate(ts.trees()):
	    if i == tree_index:
	    	break

	assert(t.index == tree_index)

	# pick the branch
	cpicked = -1
	while cpicked < c:
		treeloc = t.total_branch_length * np.random.uniform()
		for mut_n in t.nodes():
		    if mut_n != t.root:
		        treeloc -= t.branch_length(mut_n)
		        if treeloc <= 0:
		            cpicked = t.num_samples(mut_n)/(n)
		            print(cpicked)
		            break

	# pick the location on the sequence
	mut_base = 0.0 + np.random.randint(low=np.ceil(t.interval[0]), high=np.ceil(t.interval[1]), size=1)

	# the following assumes that there's no other mutations in the tree sequence
	assert(ts.num_sites == 0)

	# the mutation metadata
	mut_md = pyslim.MutationMetadata(mutation_type=1, selection_coeff=0.0, population=1, slim_time=1)

	tables = ts.tables
	site_id = tables.sites.add_row(position=mut_base, ancestral_state=b'')
	tables.mutations.add_row(site=site_id, node=mut_n, derived_state='1',
	        metadata=pyslim.encode_mutation([mut_md]))

	mut_ts = pyslim.load_tables(tables)

	# genotypes
	#out_slim_targets = open('%s.slim.targets'%(out),'w')
	#for i,g in enumerate(mut_ts.genotype_matrix()[0]):
	#	if g == 1:
	#		#print(i) 
	#		out_slim_targets.write('%d\n'%(i))	
	#out_slim_targets.close()
	print(mut_ts.genotype_matrix())
	print('%d / %d' %(np.sum(mut_ts.genotype_matrix()),n))
	freq = np.sum(mut_ts.genotype_matrix())/(n)

	return mut_base, freq, mut_ts 


ts = pyslim.annotate_defaults(msprime.simulate(args.n, recombination_rate = args.r, length=args.l),
                              model_type="WF", slim_generation=1)

mut_base, freq, mut_ts = throw_mut_on_tree(ts)

# save treeseq 
mut_ts.dump("%s.trees"%(args.out))
# save slim command
#out_slim = open('%s.slim.cmd'%(out),'w')
basename = args.out
#out_slim.write('./slim -d \"basename=\'%s\'\" ssv.slim'%(basename))
#out_slim.close()

s = args.s

c = args.c
w = np.sum(1/np.arange(np.ceil(args.n*c),np.floor(args.n*(1-c))+1))
if args.pg != None:
	I = args.pg[0]
	k = args.pg[1]
	beta = np.random.normal(0,np.sqrt(w/k))	
	s = beta * I
	print(I,k,beta,s)
np.savetxt(basename+'.metadata',np.array([s,freq,mut_base]))

print('./slim -d r=%.3e -d l=%d -d u=%.3e -d \"basename=\'%s\'\" -d s=%f ssv.slim'%(args.r,args.l,args.u,basename,s))




