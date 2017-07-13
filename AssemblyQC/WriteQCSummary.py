#!/usr/bin/env python

import sys
import numpy as np
from Bio import SeqIO
if (len(sys.argv) != 7):
  print "Usage: WriteQCSummary.py assembly vcfFile covFile window maxCov"
  sys.exit(0)

assemblyFile = open(sys.argv[1])
vcfFile = open(sys.argv[2])
covFile = open(sys.argv[3])
window = int(sys.argv[4])
maxCov = int(sys.argv[5])
dirName = sys.argv[6]
vcf = [l.strip() for l in vcfFile]

allCov = [l.split() for l in covFile]

cov = [int(v[2]) for v in allCov]


i = 0
span = 0
maxSpan = 0
spanStart = 0
spanEnd = 0
last = min(window, len(cov))
total = sum(cov[0:last])
maxTotal = 0
avgCov = np.mean(cov)
medianCov = np.median(cov)
maxAvg = total/float(last)
maxIndex = 0
for i in range(0,len(cov) - window - 1):
    total -= cov[i]
    total += cov[i+window]
    avg = total / float(window)
    if (avg > maxAvg):
        maxIndex = i
        maxAvg = avg


contigs = list(SeqIO.parse(assemblyFile, "fasta"))
avgContigLength = 0
if (len(contigs) > 0):
    contigLengths = [len(c) for c in contigs]
    avgContigLength = sum(contigLengths)/len(contigLengths)


print dirName + "\t" + str(len(contigs)) + "\t" + "{:2.2f}".format(avgContigLength) + "\t" + str(len(contigs[0])) + "\t" + str(len(vcf)) + "\t" + str(maxAvg) + "\t" + str(maxIndex)  + "\t{:2.2f}\t{:2.2f}".format(avgCov, medianCov)
