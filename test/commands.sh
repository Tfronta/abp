#!/bin/bash

mkdir -p withFakeVCF
cd withFakeVCF

../../fakeVCF.py ../assembly.consensus.fragments.snv.mat.categorized --vcf test.vcf

../../PairedSNVs.py ../assembly.consensus.fragments.snv.mat.categorized --maxCov 350 --minCov 50 --graph mi.gml --adj min.adj --minNShared 5 --minLRT 1.5 --vcf  test.vcf

../../MinDisagreeCluster.py --graph mi.gml --cuts mi.gml.cuts --sites mi.gml.sites --factor 2 --out mi.cuts.gml --swap 1000


cd ..
mkdir -p original
cd original

../../PairedSNVs.py ../assembly.consensus.fragments.snv.mat.categorized --maxCov 350 --minCov 50 --graph mi.gml --adj min.adj --minNShared 5 --minLRT 1.5 --vcf assembly.consensus.nucfreq.vcf

../../MinDisagreeCluster.py --graph mi.gml --cuts mi.gml.cuts --sites mi.gml.sites --factor 2 --out mi.cuts.gml --swap 1000





cd ..
echo "original"
cat original/mi.gml.cuts
echo "fake vcf"
cat withFakeVCF/mi.gml.cuts




