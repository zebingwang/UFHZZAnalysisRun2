HZZ Analyzer for CMS Run2

------

To install:

setenv SCRAM_ARCH slc6_amd64_gcc530

cmsrel CMSSW_8_0_26_patch1

cd CMSSW_8_0_26_patch1/src

cmsenv

git cms-init

git clone -b 80X https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git

cp UFHZZAnalysisRun2/install.sh .

./install.sh

cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_80X_Moriond_cfg_1.py

cp UFHZZAnalysisRun2/Utilities/crab/* .

voms-proxy-init --valid=168:00

source /cvmfs/cms.cern.ch/crab3/crab.sh

python SubmitCrabJobs.py -t "myTask_Data" -d datasets_2016ReReco.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_80X_M1703Feb_2l_cfg.py

or similary for MC:

python SubmitCrabJobs.py -t "myTask_MC" -d datasets_Summer16_25ns_MiniAOD.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_80X_M17_4l_cfg.py

You can use manageCrabTask.py to check the status, resubmit, or kill your task. E.g. after submitting:

nohup python -u manageCrabTask.py -t resultsAna_Data_M17_Feb19 -r -l >& managedata.log &

This will start an infinite loop of running crab resubmit on all of your tasks, then sleep for 30min. You should kill the process once all of your tasks are done. Once all of your tasks are done, you should run the following command to purge your crab cache so that it doesn't fill up:

python manageCrabTask.py -t resultsAna_Data_M17_Feb19 -p


