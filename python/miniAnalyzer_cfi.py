import FWCore.ParameterSet.Config as cms

ANALYSISVARS={
   'basic': [],
   'full':  ['generator','leptons','photons','jets','pflow','met','pps','rawmu'],
   'full_sys':  ['generator','leptons','photons','jets','pflow','met','pps','rawmu','systematics'],
   'lepton':['leptons'],
   'pps_zx':['leptons','pps','generator'],
   'pps_yx':['photons','pps','generator'],
   'pps':   ['pps'],
   'ttbar': ['leptons','jets','met','pps','generator','pflow'],
   'dilep': ['leptons','pps','generator','pflow'],
   'dilep_sys': ['leptons','pps','generator','pflow','systematics'],
   'dilep_data': ['leptons','pps','pflow'],
   'ttbar_data': ['leptons','jets','met','pps','pflow'],
   'ttbar_sys': ['leptons','jets','met','pps','generator','systematics','pflow'],
   'lowmu': ['leptons','jets','met','pps','generator','pflow'],
   'lowmu_data': ['leptons','photons','jets','met','pps','pflow'],
   'lowmu_sys': ['leptons','photons','jets','met','pps','generator','systematics','pflow'],
   'test':  ['generator','pps']
}
ANALYSISTRIGGERLISTS={
    2018:[
        'HLT_Ele32_WPTight_Gsf_v',
        'HLT_Ele35_WPTight_Gsf_v',
        'HLT_Ele38_WPTight_Gsf_v',
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
        'HLT_IsoMu24_v',
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
        'HLT_Photon50_R9Id90_HE10_IsoM_v',
        'HLT_Photon75_R9Id90_HE10_IsoM_v',
        'HLT_Photon90_R9Id90_HE10_IsoM_v',
        'HLT_Photon120_R9Id90_HE10_IsoM_v',
        'HLT_Photon165_R9Id90_HE10_IsoM_v',
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_v',
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_v',
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v',
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_v3',
        'HLT_Photon175_v',
        'HLT_Photon200_v',
        'HLT_Photon300_NoHE_v'
    ],
    2017:['HLT_Ele32_WPTight_Gsf_v',
          'HLT_Ele32_WPTight_Gsf_L1DoubleEG_v',
          'HLT_Ele35_WPTight_Gsf_v',
          'HLT_Ele38_WPTight_Gsf_v',
          'HLT_Ele40_WPTight_Gsf_v',
          'HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v',
          'HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v',
          'HLT_IsoMu24_v',
          'HLT_IsoMu24_eta2p1_v',
          'HLT_IsoMu27_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
          'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Photon150_v',
          'HLT_Photon175_v',
          'HLT_Photon200_v',
          'HLT_Photon50_R9Id90_HE10_IsoM_v',
          'HLT_Photon75_R9Id90_HE10_IsoM_v',
          'HLT_Photon90_R9Id90_HE10_IsoM_v',
          'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v',
          'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_v',
          'HLT_PFJet40_v',
          'HLT_PFJet60_v',
          'HLT_PFJet80_v',
          'HLT_PFJet140_v',
          'HLT_PFJet200_v',
          'HLT_PFJet260_v',
          'HLT_PFJet320_v',
          'HLT_PFJet400_v',
          'HLT_PFJet450_v',
          'HLT_PFJet500_v',
          'HLT_PFJet550_v',
          'HLT_PFJetFwd40_v',
          'HLT_PFJetFwd60_v',
          'HLT_PFJetFwd80_v',
          'HLT_PFJetFwd140_v',
          'HLT_PFJetFwd200_v',
          'HLT_PFJetFwd260_v',
          'HLT_PFJetFwd320_v',
          'HLT_PFJetFwd400_v',
          'HLT_PFJetFwd450_v',
          'HLT_PFJetFwd500_v',
          'HLT_Mu50_v',
          'HLT_TkMu50_v',
          'HLT_ZeroBias_v',
          'HLT_HIMu12_v',
          'HLT_HIMu15_v',
          'HLT_HIEle15_WPLoose_Gsf_v1',
          'HLT_HIEle20_WPLoose_Gsf_v1',
          'HLT_HIPhoton50_HoverELoose_v1',
          'HLT_HIPhoton60_HoverELoose_v1'
      ],
      #https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2017
      'ttbar2017':['HLT_Ele35_WPTight_Gsf_v',
          'HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v',
          'HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v',
          'HLT_IsoMu27_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
          'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v',
          'HLT_PFHT430_SixJet40_BTagCSV_p080_v ',
          'HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v',
          'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v',
          'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v'
      ],      
      'dilep2017':['HLT_Ele35_WPTight_Gsf_v',
          'HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v',
          'HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v',
          'HLT_IsoMu27_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
          'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v'
      ],  
      'lowmu2017':['HLT_HIPFJet140_v',
          'HLT_HIPFJetFwd140_v',
          'HLT_HIPhoton60_HoverELoose_v',
          'HLT_HIEle15_WPLoose_Gsf_v',
          'HLT_HIMu15_v',
          'HLT_AK4PFJet120_v',
          'HLT_Ele20_WPLoose_Gsf_v',
          'HLT_Mu17_v'
      ],          
    2016:['HLT_Ele32_eta2p1_WPTight_Gsf_v',
          'HLT_IsoMu24_v',
          'HLT_IsoMu24_eta2p1_v',
          'HLT_IsoTkMu24_v',
          'HLT_IsoTkMu24_eta2p1_v',
          'HLT_Mu50_v',
          'HLT_TkMu50_v',
          'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
          'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
          'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
          'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
          'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
          'HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
          'HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
          'HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v',
          'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
          'HLT_Photon36_v',
          'HLT_Photon36_R9Id90_HE10_IsoM_v',
          'HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF',
          'HLT_Photon50_v',
          'HLT_Photon50_R9Id90_HE10_IsoM_v',
          'HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF',
          'HLT_Photon75_v',
          'HLT_Photon75_R9Id90_HE10_IsoM_v',
          'HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF',
          'HLT_Photon90_v',
          'HLT_Photon90_R9Id90_HE10_IsoM_v',
          'HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF',
          'HLT_Photon120_v',
          'HLT_Photon120_R9Id90_HE10_IsoM_v',
          'HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_v',
          'HLT_Photon165_HE10',
          'HLT_Photon165_R9Id90_HE10_IsoM_v',
          'HLT_Photon175_v',
          'HLT_Photon250_NoHE_v',
          'HLT_Photon300_NoHE_v',
          'HLT_PFJet40_v',
          'HLT_PFJet60_v',
          'HLT_PFJet80_v',
          'HLT_PFJet140_v',
          'HLT_PFJet200_v',
          'HLT_PFJet260_v',
          'HLT_PFJet320_v',
          'HLT_PFJet400_v',
          'HLT_PFJet450_v',
          'HLT_PFJet500_v'
      ]
}

ANALYSISJETIDS={
    2018:'tightLepVeto',
    2017:'tightLepVeto',
    2016:'looseID'
}

ANALYSISRUNS={
    'era2017'   :297050,
    'era2017_B' :297050,
    'era2017_C' :299368,
    'era2017_C1':299368,
    'era2017_C2':300806,
    'era2017_D' :302031,
    'era2017_E' :303825,
    'era2017_F' :305044,
    'era2017_F1':305044,
    'era2017_F2':305178,
    'era2017_F3':305967,
    'era2017_H' :306936
}


analysis = cms.EDAnalyzer("MiniAnalyzer",
                          saveTree               = cms.bool(True),
                          ListVars               = cms.vstring(ANALYSISVARS['full']),
                          FilterType             = cms.string("full"),
                          applyFilt              = cms.bool(True),
						  runNumber              = cms.untracked.int32(1),
                          triggerBits            = cms.InputTag("TriggerResults","","HLT"),
                          prescales              = cms.InputTag("patTrigger"),
                          l1prescales            = cms.InputTag("patTrigger","l1min"),
                          triggersToUse          = cms.vstring(ANALYSISTRIGGERLISTS[2017]),
                          rho                    = cms.InputTag("fixedGridRhoFastjetAll"),
                          vertices               = cms.InputTag("offlineSlimmedPrimaryVertices"),                          
                          muons                  = cms.InputTag("slimmedMuons"),
                          RoccoR                 = cms.string("muoncorr_db.txt"),   
                          PPS_pixelEff           = cms.string("TopLJets2015/TopAnalysis/data/era2017/pixelEfficiencies_multiRP.root"),                        
                          electrons              = cms.InputTag("slimmedElectrons"),
                          photons                = cms.InputTag("slimmedPhotons"),
                          jets                   = cms.InputTag('updatedPatJetsUpdatedJECBTag'),
                          jetIdToUse             = cms.string(ANALYSISJETIDS[2017]),
                          EGIDVersion            = cms.string("2"),
                          jecUncSources          = cms.vstring("AbsoluteStat","AbsoluteScale","AbsoluteMPFBias","Fragmentation","SinglePionECAL","SinglePionHCAL","FlavorPureGluon","FlavorPureQuark","FlavorPureCharm","FlavorPureBottom","TimePtEta","RelativeJEREC1","RelativeJEREC2","RelativeJERHF","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeBal","RelativeFSR","RelativeStatFSR","RelativeStatEC","RelativeStatHF","PileUpDataMC","PileUpPtRef","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF"),
                          jecUncFile             = cms.string('jecUncSources.txt'),
                          metFilterBits          = cms.InputTag("TriggerResults","","PAT"),
                          metFiltersToUse        = cms.vstring('Flag_HBHENoiseFilter',
                                                               'Flag_HBHENoiseIsoFilter',
                                                               'Flag_EcalDeadCellTriggerPrimitiveFilter',
                                                               'Flag_goodVertices',
                                                               'Flag_eeBadScFilter',
                                                               'Flag_globalTightHalo2016Filter'), 
                          badChCandFilter        = cms.InputTag('BadChargedCandidateFilter'),
                          badPFMuonFilter        = cms.InputTag('BadPFMuonFilter'),
                          PUprotons              = cms.InputTag(''),
                          mets                   = cms.InputTag('slimmedMETsModifiedMET'),                          
                          pfCands                = cms.InputTag('packedPFCandidates'),
                          ctppsLocalTracks       = cms.InputTag('ctppsLocalTrackLiteProducer'),
                          tagRecoProtons         = cms.InputTag('ctppsProtons',"singleRP"),
                          tagMultiRecoProtons    = cms.InputTag('ctppsProtons',"multiRP"),
                          #extra info from AODs
                          tagStripHits           = cms.InputTag('totemRPRecHitProducer'),
                          tagStripPatterns       = cms.InputTag('totemRPUVPatternFinder'),
                          )
