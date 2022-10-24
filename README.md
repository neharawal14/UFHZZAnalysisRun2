# HZZ Analyzer for CMS Run2 for CMSSW_12_4_3


SCRAM_ARCH=slc7_amd64_gcc900; export SCRAM_ARCH

cmsrel CMSSW_12_4_3

cd CMSSW_12_4_3/src/

cmsenv

git cms-init

git clone https://github.com/JHUGen/JHUGenMELA

cd JHUGenMELA/MELA/

./setup.sh -j 8

eval $(./setup.sh env) # remember to do this at every login

cd ../../

git clone -b 12_4_3 https://ferrico@github.com/ferrico/UFHZZAnalysisRun2.git

git clone -b UL20_10_6_26  https://github.com/ferrico/KinZfitter.git

scramv1 b -j 8


# ####### side note #######
# 1) Muon, electron and jets corrections have been disabled
# 2) electronID needs to be updated
# ####### side note #######




##### OLD STUFF #####
------

To install:

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_12

cd CMSSW_10_6_12/src

cmsenv

git cms-init

git clone -b UL_19 https://ferrico@github.com/ferrico/UFHZZAnalysisRun2.git

cp UFHZZAnalysisRun2/install*.sh .

./install_2.sh

cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_102X_2018_Legacy_cfg.py
cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_102X_2017_Legacy_cfg.py
cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_102X_2016_Legacy_cfg.py

cp UFHZZAnalysisRun2/Utilities/crab/* .

voms-proxy-init --valid=168:00
#probably need "voms-proxy-init -voms cms -rfc"

source /cvmfs/cms.cern.ch/crab3/crab.sh

python SubmitCrabJobs.py -t "myTask_Data" -d datasets_2016ReReco.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_80X_M1703Feb_2l_cfg.py

or similary for MC:

python SubmitCrabJobs.py -t "myTask_MC" -d datasets_Summer16_25ns_MiniAOD.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_80X_M17_4l_cfg.py

You can use manageCrabTask.py to check the status, resubmit, or kill your task. E.g. after submitting:

nohup python -u manageCrabTask.py -t resultsAna_Data_M17_Feb19 -r -l >& managedata.log &

This will start an infinite loop of running crab resubmit on all of your tasks, then sleep for 30min. You should kill the process once all of your tasks are done. Once all of your tasks are done, you should run the following command to purge your crab cache so that it doesn't fill up:

python manageCrabTask.py -t resultsAna_Data_M17_Feb19 -p

UFHZZ4LAna/python/templateMC_102X_Legacy16_4l_cfg.py
UFHZZ4LAna/python/templateMC_102X_Legacy17_4l_cfg.py
UFHZZ4LAna/python/templateMC_102X_Legacy18_4l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy16_3l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy17_3l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy18_3l_cfg.py
