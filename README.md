HZZ Analyzer for CMS Run2

------

To install:

setenv SCRAM_ARCH slc6_amd64_gcc530

cmsrel CMSSW_8_0_10

cd CMSSW_8_0_10/src

cmsenv

git cms-init

echo /PhysicsTools/PatUtils/ >> .git/info/sparse-checkout

git cms-merge-topic cms-met:metTool80X

git clone -b 80X https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git

git clone -b 74x-root6 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement

cd ZZMatrixElement ; git checkout -b from-master master

cd ..

git clone https://github.com/mhl0116/KinZfitter.git

git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V2 KaMuCa

scram b -j 8

cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_80X_cfg.py

