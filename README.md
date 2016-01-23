HZZ Analyzer for CMS Run2

------


To install:

setenv SCRAM_ARCH slc6_amd64_gcc493

cmsrel CMSSW_7_6_3

cd CMSSW_7_6_3/src

cmsenv

git clone https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git

git clone -b 74x-root6 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement

cd ZZMatrixElement

git checkout -b from-V00-02-01-patch1 V00-02-01-patch1

cd ..

git clone https://github.com/tocheng/KinZfitter.git

scram b -j 8

cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_76X_cfg.py


Instructions on how to use git:

4.) Development of the codes and put your codes to git:

4.1.) Inside directory:

      CMSSW_xxxx/src/UFHZZAnalysisRun2

4.2.) Make a new branch for yourself, choose the branch name as you want, but do not be the same as others’:

4.2.1.) Create a branch of your own :

    git branch your_own_branch_name 

4.2.2.) Change the codes of the working area to be the new branch

     git checkout your_own_branch_name 

4.2.3.) Put your new branch to the to remote git repository https://github.com/

       git push origin your_own_branch_name:your_own_branch_name

4.3.) Update your changes to github after you made some modifications:

4.3.1.) Commit change to your local repository:

    git commit -m “note of the change” -a 

4.3.2.) Update the remote repository:

    git push origin

4.4.) If you want to add some new files or delete some old files:

      git add a_new_file

      git rm an_old_file

      git commit -m “note of the change add new files, and delete old files” -a

      git push origin 

New module for electron MVA ID:

SlimmedElectronMvaIDProducer

The instruction for the module can be found in 
SlimmedElectronMvaIDProducer/doc/

