from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = "SignalDiffractiveTTbar"
config.General.workArea = "grid"
config.General.transferOutputs=True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/eos/home-p/pchatagn/CMSSW_10_6_24/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py"
config.JobType.disableAutomaticOutputCollection = False
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 1200
config.JobType.pyCfgParams = ['runOnData=False','era=era2017','ListVars=ttbar_sys']
config.JobType.inputFiles = ['jec_MC.db','jer_MC.db','muoncorr_db.txt','jecUncSources.txt','qg_db.db']

config.section_("Data")
config.Data.inputDataset = "/eos/user/p/pchatagn/OutputDirFPMC/MINIAOD/exclu_ttbar_lephad_QCD"
config.Data.ignoreLocality = False
config.Data.inputDBS = "global"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 4
config.Data.publication = False
config.Data.outLFNDirBase = "/eos/user/p/pchatagn/OutputNtuplizer/813afaa/"

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
