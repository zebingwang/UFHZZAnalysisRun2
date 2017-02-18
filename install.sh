#!/bin/sh
git clone -b 80X https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git
git cms-merge-topic rafaellopesdesa:EgammaAnalysis80_EGMSmearer_Moriond17_23Jan
cd EgammaAnalysis/ElectronTools/data
git clone https://github.com/ECALELFS/ScalesSmearings.git
git checkout Moriond17_23Jan_v1
cd -
git cms-merge-topic cms-met:METRecipe_8020 -u
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
cp UFHZZAnalysisRun2/Utilities/fixcrab.patch ZZMatrixElement
(cd ZZMatrixElement ; git checkout -b from-v203 v2.0.4 ; git apply fixcrab.patch ; . setup.sh -j 12)
git clone https://github.com/VBF-HZZ/KinZfitter.git
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa
scram b -j 8
