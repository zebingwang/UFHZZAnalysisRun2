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

git remote add -f -t smearings80X shervin86 https://github.com/shervin86/cmssw.git

git cherry-pick f3b0b0140483c336212baa035cf9a820a016a799

git cherry-pick a5aaeb7a598800ae98e88ea1a952ecd1d66aa059

git cherry-pick c7ac16dd88969510d2d6d6ea2c4702e0108bf151

git cherry-pick 054a90830c77423ee673204611522018ace69c5d

git cms-addpkg EgammaAnalysis/ElectronTools

cd EgammaAnalysis/ElectronTools/data

git clone -b ICHEP2016_approval_4fb https://github.com/ECALELFS/ScalesSmearings.git

cd -

git clone -b 80X https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git

git clone -b 74x-root6 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement

cd ZZMatrixElement ; git checkout -b from-master master

cd ..

git clone https://github.com/VBF-HZZ/KinZfitter.git

git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V3 KaMuCa

scram b -j 8

cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_80X_cfg.py

