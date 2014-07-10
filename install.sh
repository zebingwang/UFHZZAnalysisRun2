#!/bin/sh

scramv1 project CMSSW CMSSW_5_3_9_patch3
cd CMSSW_5_3_9_patch3/src
eval `scramv1 runtime -sh`
git clone https://github.com/VBF-HZZ/UFHZZAnalysis8TeV.git
cd UFHZZAnalysis8TeV
git checkout -b testProd origin/testProd
cd ../
tar xvzf UFHZZAnalysis8TeV/requiredPackages.tgz
scram b -j12


