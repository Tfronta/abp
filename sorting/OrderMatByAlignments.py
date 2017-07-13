#!/usr/bin/env python

import argparse
import sys
ap = argparse.ArgumentParser(description="Order snv mat by read order from another file")
ap.add_argument("mat", help="snv matrix")
ap.add_argument("reads", help="read order")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

reads = open(args.reads)
matFile = open(args.mat)
outFile= open(args.out, 'w')

freqLine = matFile.readline()
readMat = {}
for line in matFile:
    v = line.split()
    readMat[v[1]] = v[0]

for line in reads:
    (r,c)= line.split()[0:2]
    
    if (r in readMat):
        sys.stdout.write(readMat[r] + "\t" +  r + "\t" + c + "\n")
