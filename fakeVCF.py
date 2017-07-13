#!/usr/bin/env python


import argparse
import numpy as np
import sys
import pickle

ap = argparse.ArgumentParser(description="Sort by haplotype")
ap.add_argument("mat", help="snv matrix file")
ap.add_argument("--vcf", help="outputvcf vcf")
#args = ap.parse_args('assembly.consensus.fragments.snv.mat.categorized')

args = ap.parse_args()
matFile = open(args.mat)


lenOfFakeVCF = 0

for line in matFile:
    line=line.split("\t")
    lenOfFakeVCF = len(line[0]) 
    break


#print(lenOfFakeVCF)

outVCF=""
for i in range(lenOfFakeVCF):
   line = "1F\t{}\tsnp{}\tC\tG\t.\tPASS\t*\tGT:GQ\t0/1/1/1:100\n".format(i+1,i+1)
   outVCF += line

#print(outVCF)


f = open(args.vcf, "w+")
f.write(outVCF)

