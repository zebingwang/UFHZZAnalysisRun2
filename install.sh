#!/bin/sh

scramv1 project CMSSW CMSSW_7_0_6_patch1
cd CMSSW_7_0_6_patch1/src
eval `scramv1 runtime -sh`
git clone https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git
cd UFHZZAnalysisRun2
git checkout -b csa14 origin/csa14
cd ../
scram b -j12



