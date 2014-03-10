import FWCore.ParameterSet.Config as cms

process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring()
myfilelist.extend( [

'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_1.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_10.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_11.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_12.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_13.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_14.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_2.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_3.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_4.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_5.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_6.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_7.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_8.root',
'file:/home/msnowball/newHiggsPatuplizer/RegressionTest/CMSSW_5_3_4/src/MissedEvents/PAT/missedEvents_53X_9.root'


]
)

process.source = cms.Source("PoolSource",fileNames = myfilelist)
#process.source.eventsToProcess = cms.untracked.VEventRange("1:912643")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Data53X.root")
                                   )

process.AnaAfterHlt = cms.EDAnalyzer('UFHZZ4LAna',
                             photonSrc    = cms.untracked.InputTag("cleanPatPhotons"),
                             electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                             muonSrc      = cms.untracked.InputTag("rochesterMuons"),
                             jetSrc       = cms.untracked.InputTag("cleanPatJets"),
                             metSrc       = cms.untracked.InputTag("patMETsPF"),
                             vertexSrc    = cms.untracked.InputTag("goodOfflinePrimaryVertices"), #or selectedVertices 
                             isMC         = cms.untracked.bool(False),
 #                            isSignal     = cms.untracked.bool(True),
 #                            mH           = cms.untracked.uint32(120),
 #                            CrossSection = cms.untracked.double(1),#0.00242802 ),
 #                            FilterEff    = cms.untracked.double(1),
 #                            weightEvents = cms.untracked.bool(True),
                             elRhoSrc     = cms.untracked.InputTag("kt6PFJets", "rho","RECO"),
                             muRhoSrc     = cms.untracked.InputTag("kt6PFJetsCentralNeutral", "rho","RECO"),
                             mZ2Low       = cms.untracked.double(4),
 #                            reweightForPU = cms.untracked.bool(False),
 #                            interactiveRun = cms.untracked.bool(True),
                             doFsrRecovery = cms.untracked.bool(True),
                             _muPtCut = cms.untracked.double(5),
                             _elecPtCut = cms.untracked.double(7)
                             )



# Trigger
process.hltHighLevel = cms.EDFilter("HLTHighLevel",
                                    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                    HLTPaths = cms.vstring( 'HLT_Mu17_TkMu8*',
                                                            'HLT_Mu17_Mu8*',
                                                            #'HLT_DoubleMu7_v*',
                                                            #'HLT_DoubleMu3_v*',
                                                            #'HLT_Ele10_LW_LR1',
                                                            #'HLT_Ele15_SW_LR1',
                                                            #'HLT_Ele15_SW_CaloEleId_L1R',
                                                            #'HLT_Ele17_SW_TightEleId_L1R_v*',
                                                            #'HLT_Ele17_SW_TighterEleIdIsol_L1R_v*',
                                                            #'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
                                                            #'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL*',
                                                            'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*',
                                                            'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*',
                                                            'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*',
                                                            'HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL*'
                                                            ),
                                    # provide list of HLT paths (or patterns) you want
                                    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                    andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true  
                                    throw = cms.bool(False)    # throw exception on unknown path names 
                                    )



process.p = cms.Path(process.hltHighLevel
                     *process.AnaAfterHlt
                     )






