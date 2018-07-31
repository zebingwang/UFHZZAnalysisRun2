#!/bin/sh
git cms-merge-topic rafaellopesdesa:EgammaAnalysis80_EGMSmearer_Moriond17_23Jan
cd EgammaAnalysis/ElectronTools/data
git clone https://github.com/ECALELFS/ScalesSmearings.git
cd -
git cms-merge-topic cms-met:METRecipe_8020 -u
git cms-merge-topic cms-met:METRecipe_80X_part2 -u
git clone https://github.com/VBF-HZZ/KinZfitter.git
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa
git cms-merge-topic -u perrozzi:HTXS_clean
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
cp UFHZZAnalysisRun2/Utilities/fixcrab.patch ZZMatrixElement
(cd ZZMatrixElement ; git checkout -b from-v204 v2.0.4 ; git apply fixcrab.patch ; . setup.sh -j 12)
#NOTE need to copy MELA library from peviously compiled area (the permission to access cms-higgs.cern.ch is gone, somehow) 
# e.g. scp from /raid/raid9/ahmad/CMSSW_8_0_26_patch1/src/src/ZZMatrixElement/MELA/data/slc6_amd64_gcc530/libmcfm_703.so 
# or cp lxplus.cern.ch:~mschen/public/libmcfm_703.so-80X-v2.0.4 libmcfm_703.so
scram b -j 8
