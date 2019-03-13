import msprime
import pyslim
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ts',type=str,help='unmutated tree sequence')
parser.add_argument('out',type=str)
parser.add_argument('-u','--u',type=float,default=1e-8,help='mut rate per bp per gen')
args = parser.parse_args()

ts = pyslim.load(args.ts+'.out.trees')
mut_ts = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=args.u, keep=True))
mut_ts.dump(args.out+'.out.mut.trees')
