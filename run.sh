#!/bin/bash
#./run.sh excl_ttbar 16006 1 QED 15 lephad
startMsg='Job started on '`date`
echo $startMsg

####### USER SETTINGS ###########
basedir=/afs/cern.ch/user/p/pchatagn/TTbarInterface/GEN_10_6_10_ProtonReco
FPMCfolder=/afs/cern.ch/user/p/pchatagn/FPMC/build/
outfolder=/eos/user/p/pchatagn/OutputDirFPMC/
CMSSWFolder=/afs/cern.ch/user/p/pchatagn/CMSSW_10_6_16/
#source /afs/cern.ch/user/m/mpitt/.zshrc
#################################
puleupDir=\'/eos/cms/store/mc/RunIISummer19ULPrePremix/Neutrino_E-10_gun/PREMIX/UL17_106X_mc2017_realistic_v6-v1/100000\'
nEvents=200

setup10() {
    echo scram p CMSSW CMSSW_10_6_10
    export SCRAM_ARCH=slc7_amd64_gcc700
    if [ -r CMSSW_10_6_10/src ] ; then 
       echo release CMSSW_10_6_10 already exists
    else
       echo No release of CMSSW
      scram p CMSSW CMSSW_10_6_10
      curl file:///afs/cern.ch/user/m/mpitt/public/exclusive_top/pythia_fragment --retry 2 --create-dirs -o CMSSW_10_6_10/src/Configuration/GenProduction/python/FSQ-RunIISummer19UL17GEN-00000-fragment.py
    fi
    cd CMSSW_10_6_10/src   
    eval `scram runtime -sh`
    scram b
    cd ../../
}

setup9UL() {
    echo scram p CMSSW CMSSW_9_4_14_UL_patch1
    if [ -r CMSSW_9_4_14_UL_patch1/src ] ; then 
       echo release CMSSW_9_4_14_UL_patch1 already exists
    else
      scram p CMSSW CMSSW_9_4_14_UL_patch1
    fi
    cd CMSSW_9_4_14_UL_patch1/src
    eval `scram runtime -sh`
    cd ../../
}

# check number of arguments
if [ "$#" -ne 6 ]; then
  echo $#
  echo "Usage: $0 NAME PROC IDX TYPINT NFLUX MODE" >&2
  echo "Example: $0 excl_ttbar 16006 1 QED 15 lephad" >&2
  exit 1
fi

echo ./run.sh $1 $2 $3 $4 $5 $6
echo
proc=${1}_${6}_${4}
mode=$6

#setup cmsenv
cd $CMSSWFolder
echo in CMSSW folder:
pwd
eval `scramv1 runtime -sh`
cd -
echo CMSENV setup done

mkdir -p ${outfolder}/MINIAOD/${proc}

procid=$2
ii=$(printf "%02d" `expr ${3} + 0`) #don't know what will be the input, 05 or 5
seed1=`expr 17673 + ${3}`
seed2=`expr 63565 + ${3}`
typeint=$4
nflux=$5 

folder=${proc}_${ii}
#create a temporary where samples will be generated
mkdir -pv $folder
cd $folder

#creates FPMC settings card
JOcard=card_${proc}_${ii}.jo
cp "$basedir/template_card.txt" $JOcard
cp "$basedir/addColorFlow.py" .
cp "$basedir/ChangeProcess.py" .
echo 'add settings to default card:'
echo "LHEFILE     '"${proc}".lhe'" >> $JOcard
echo 'IPROC       '${procid} >> $JOcard
echo "TYPINT      '"${typeint}"'" >> $JOcard
echo 'NRN1        '${seed1} >> $JOcard
echo 'NRN2        '${seed2} >> $JOcard
echo 'MAXEV       '${nEvents} >> $JOcard
echo 'NFLUX       '${nflux} >> $JOcard
echo
cat $JOcard



#simulation of aa->tt production
echo execute the code
$FPMCfolder/fpmc-lhe < $JOcard



echo fix colorflow in the output file
python addColorFlow.py ${proc}.lhe

#replace Xsec error until the bug will be fixed
sed -i 's/-306/E-06/g' ${proc}_fix.lhe

#MadSpin
echo running madspin
cp "$basedir/header_for_madspin.txt" "$basedir/runMadSpin.sh" .
bash runMadSpin.sh ${proc}_fix.lhe $mode
#bash runMadSpin_23april.sh ${proc}_fix.lhe $mode
echo End madspin


#CMSSW (set up in lxplus7)
setup10


echo "GEN-SIM starting"
cmsDriver.py step1 --filein file:events_filter.lhe  --fileout file:stepLHE.root --mc --eventcontent LHE --datatier LHE --conditions 106X_mc2017_realistic_v6 --step NONE --python_filename step0_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nEvents
cmsRun step0_cfg.py

cmsDriver.py Configuration/GenProduction/python/FSQ-RunIISummer19UL17GEN-00000-fragment.py --filein file:stepLHE.root --fileout file:stepSIM.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN,SIM --nThreads 8 --geometry DB:Extended --era Run2_2017 --python_filename step1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nEvents
#add pythia flag
sed -i "/^            'MultipartonInteractions:coreFraction=0.63', .*/a\            'ProcessLevel:all = off'," step1_cfg.py
#sed -i "/^            'SigmaTotal:mode = 0',  .*/a\            'Next:numberShowLHA = on'," step1_cfg.py
#sed -i "/^            'ColourReconnection:range=5.176',.*/a\            'Main:showAllStatistics = on'," step1_cfg.py
cmsRun step1_cfg.py

#echo EXIT using exit command
#exit

echo "DIGI-RECO starting"
cmsDriver.py step1 --filein file:stepSIM.root  --fileout file:stepDR.root --pileup_input file:dummy  --mc --eventcontent PREMIXRAW --runUnscheduled --datatier GEN-SIM-DIGI --conditions 106X_mc2017_realistic_v6 --step DIGI,DATAMIX,L1,DIGI2RAW --procModifiers premix_stage2 --nThreads 8 --geometry DB:Extended --datamix PreMix --era Run2_2017 --python_filename step2_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nEvents
#update pileup samples
sed -i "/process.mixData.input.fileNames/a import os, random\nbaseDir=$puleupDir\nfileList=[os.path.join(baseDir.replace('/eos/cms',''),f) for f in os.listdir(baseDir)]\nrandom.shuffle(fileList)\nprocess.mixData.input.fileNames = cms.untracked.vstring(fileList)" step2_cfg.py
cmsRun step2_cfg.py

echo "HLT starting"
setup9UL
cmsDriver.py step1 --filein file:stepDR.root  --fileout file:stepHLT.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions 94X_mc2017_realistic_v15 --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' --step HLT:2e34v40 --nThreads 8 --geometry DB:Extended --era Run2_2017 --python_filename stepHLT_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nEvents
cmsRun stepHLT_cfg.py

echo "AOD starting"
setup10
cmsDriver.py step1 --filein file:stepHLT.root  --fileout file:stepAOD.root --mc --eventcontent AODSIM --runUnscheduled --datatier AODSIM --conditions 106X_mc2017_realistic_v6 --step RAW2DIGI,L1Reco,RECO,RECOSIM --nThreads 8 --geometry DB:Extended --era Run2_2017 --python_filename step3_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nEvents
cmsRun step3_cfg.py

echo "MINIAOD starting"
cmsDriver.py step1 --filein file:stepAOD.root  --fileout file:miniAOD.root --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions 106X_mc2017_realistic_v6 --step PAT --nThreads 8 --geometry DB:Extended --era Run2_2017 --python_filename step4_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nEvents
cmsRun step4_cfg.py

echo done with production, ls
ls -hs

#Copy output file to $outfolder
echo "Copying to storage"
echo cp miniAOD.root  ${outfolder}/MINIAOD/${proc}/${proc}_${ii}.root
cp miniAOD.root  ${outfolder}/MINIAOD/${proc}/${proc}_${ii}.root

cd ../
#remove the temporary where samples were generated
rm -rf ${proc}_${ii}
echo $startMsg
echo job finished on `date`

