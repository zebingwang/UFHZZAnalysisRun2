#!/bin/sh
setenv SCRAM_ARCH slc6_amd64_gcc491
cmsrel CMSSW_7_4_12_patch2
cd CMSSW_7_4_12_patch2/src
cmsenv
git clone https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git
git clone -b 74x-root6 https://github.com/cms-analysis/HiggsAnalysis-CombinedLi
mit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZM
atrixElement
cd ZZMatrixElement
git checkout -b from-V00-02-01 V00-02-01
cd ..
scram b -j 8
