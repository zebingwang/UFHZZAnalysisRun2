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

cp UFHZZAnalysisRun2/Utilities/SubmitCrabJobs.py .

cp UFHZZAnalysisRun2/Utilities/crabConfig_TEMPLATE.py .

cp UFHZZAnalysisRun2/Utilities/das_client.py .

cp UFHZZAnalysisRun2/Utilities/manageCrabTask.py .

cp UFHZZAnalysisRun2/Utilities/hadd.py .

cp UFHZZAnalysisRun2/Utilities/lcgDelDir.py .

cp UFHZZAnalysisRun2/Utilities/Run2016_09Jul.txt .

cp UFHZZAnalysisRun2/Utilities/datasets_Run2016.txt .

cp UFHZZAnalysisRun2/Utilities/datasets_Spring16_25ns_MiniAOD.txt .

source /cvmfs/cms.cern.ch/crab3/crab.sh

edit crabConfig_TEMPLATE.py to point to your working directory

python SubmitCrabJobs.py -t "myTask_Data_2016" -d datasets_Run2016.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_80X_cfg.py

or similary for MC:

python SubmitCrabJobs.py -t "myTask_MC_2016" -d datasets_Spring16_25ns_MiniAOD.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_80X_cfg.py

You can use manageCrabTask.py to check the status, resubmit, or kill your task.
