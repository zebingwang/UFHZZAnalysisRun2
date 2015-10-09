import FWCore.ParameterSet.Config as cms

process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag='MCRUN2_74_V9::All'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(DUMMYFILELIST)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DUMMYFILENAME.root")
)

# clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
				     src = cms.InputTag("slimmedMuons"),
				     preselection = cms.string("track.isNonnull"),
				     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
				     fractionOfSharedSegments = cms.double(0.499),
				     )

# Electron MVA ID producer
process.mvaNonTrigV025nsPHYS14 = cms.EDProducer("SlimmedElectronMvaIDProducer",
                                     electronsCollection = cms.InputTag("slimmedElectrons","","PAT"),
                                     method = cms.string("BDTSimpleCat"),
                                     mvaWeightFile = cms.vstring(
        "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB1_5_oldscenario2phys14FIX_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB2_5_oldscenario2phys14FIX_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EE_5_oldscenario2phys14FIX_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB1_10_oldscenario2phys14FIX_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB2_10_oldscenario2phys14FIX_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EE_10_oldscenario2phys14FIX_BDT.weights.xml"
        ),
                                     Trig = cms.bool(False),
                                     )
     

process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("mvaNonTrigV025nsPHYS14","NonTrig"),
                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              correctedJetSrc = cms.untracked.InputTag("slimmedJets"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"), #or selectedVertices 
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(DUMMYCROSSSECTION),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              reweightForPU = cms.untracked.bool(True),
                              triggerSrc = cms.untracked.InputTag("TriggerResults","","HLT"),
                              triggerList = cms.untracked.vstring(
                                            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                            'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                            'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
                                            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
                                            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
                                            'HLT_TripleMu_12_10_5_v',
                                            'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
                                            'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
                                            'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                            'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                            'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
                                            'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
                                            'HLT_Ele32_eta2p1_WP75_Gsf_v',
                                            'HLT_Ele22_eta2p1_WP75_Gsf_v',
                                            'HLT_Ele27_WP85_Gsf_v',
                                            'HLT_Ele27_eta2p1_WP75_Gsf_v',
                                            'HLT_IsoMu20_v',
                                            'HLT_IsoTkMu20_v',
                              ),
                              verbose = cms.untracked.bool(False)              
#                              verbose = cms.untracked.bool(True)              
                             )

process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

process.p = cms.Path(#process.reCorrectedPatJets
                     process.fsrPhotonSequence*
		     process.boostedMuons*
                     process.mvaNonTrigV025nsPHYS14*
                     process.Ana
                     )
