#!/usr/bin/env python

import argparse
import sys
ap = argparse.ArgumentParser(description="Filter sam by read list")
ap.add_argument("sam", help="Sam file.")
ap.add_argument("readlist", help="Read list.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

readListFile = open(args.readlist)
outFile      = open(args.out, 'w')
readList = {l.strip():True for l in readListFile}
sam = open(args.sam)
for line in sam:
    if (line[0] == '@'):
        outFile.write(line)
        continue
    v = line.split()
    if v[0] in readList:
        outFile.write(line)
outFile.close()



