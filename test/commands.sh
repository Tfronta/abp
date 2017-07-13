#!/bin/bash

# this one is the one you want to use SUDIPTO, first it creates a fakevcf from the matrix file
# and then the next two scripts do the actual CC
mkdir -p withFakeVCF
cd withFakeVCF

../../fakeVCF.py ../assembly.consensus.fragments.snv.mat.categorized --vcf test.vcf

../../PairedSNVs.py ../assembly.consensus.fragments.snv.mat.categorized --maxCov 350 --minCov 50 --graph mi.gml --adj min.adj --minNShared 5 --minLRT 1.5 --vcf  test.vcf

../../MinDisagreeCluster.py --graph mi.gml --cuts mi.gml.cuts --sites mi.gml.sites --factor 2 --out mi.cuts.gml --swap 1000

cd ..

# these commands are just to make sure the fake vcf gives similar results to the real vcf 
mkdir -p original
cd original

../../PairedSNVs.py ../assembly.consensus.fragments.snv.mat.categorized --maxCov 350 --minCov 50 --graph mi.gml --adj min.adj --minNShared 5 --minLRT 1.5 --vcf assembly.consensus.nucfreq.vcf

../../MinDisagreeCluster.py --graph mi.gml --cuts mi.gml.cuts --sites mi.gml.sites --factor 2 --out mi.cuts.gml --swap 1000


cd ..


# cating the results 
echo ""
echo "original"
cat original/mi.gml.cuts
echo "fake vcf"
cat withFakeVCF/mi.gml.cuts




