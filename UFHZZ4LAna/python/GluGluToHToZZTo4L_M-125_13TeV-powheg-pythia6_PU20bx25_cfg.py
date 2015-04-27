import FWCore.ParameterSet.Config as cms

process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag='START53_V23::All'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring()
myfilelist.extend( [

       # the sync file
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3295EF7C-2070-E411-89C4-7845C4FC35DB.root'

       # all files
       #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/148E558C-946F-E411-AFA7-7845C4FC3A52.root',
       #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3295EF7C-2070-E411-89C4-7845C4FC35DB.root',
       #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/40729A64-946F-E411-B6CE-7845C4FC3A52.root',
       #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/60908275-2070-E411-94E2-7845C4FC3A58.root',
       #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/84965C0D-946F-E411-99FA-848F69FD2955.root',
       #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/A2ABB939-3670-E411-BCA8-848F69FD4553.root',
       #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/C0E81726-3670-E411-BE66-008CFA007CE0.root'

]
)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
#                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                             eventsToProcess = cms.untracked.VEventRange('1:1-1:1000')
#                             eventsToProcess = cms.untracked.VEventRange('1:66718-1:66718')
#                             eventsToProcess = cms.untracked.VEventRange('1:171190-1:171190')

                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT_testAna.root")                              
                                   )
'''
process.reCorrectedPatJets = cms.EDProducer("PatJetReCorrector",
                                            jets = cms.InputTag('selectedPatJets'),
                                            payload = cms.string('AK5PF'),
                                            rho = cms.InputTag('kt6PFJets', 'rho','RECO'),
                                            levels = cms.vstring('L1FastJet','L2Relative','L3Absolute')
                                            )
'''

process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              correctedJetSrc = cms.untracked.InputTag("slimmedJets"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"), #or selectedVertices 
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(0.00544046 ),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              reweightForPU = cms.untracked.bool(True),
                              verbose = cms.untracked.bool(False)              
#                              verbose = cms.untracked.bool(True)              
                             )

process.AnaAfterHlt = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              correctedJetSrc = cms.untracked.InputTag("slimmedJets"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"), #or selectedVertices 
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(0.00544046 ),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              reweightForPU = cms.untracked.bool(True),
#                              verbose = cms.untracked.bool(True)
                              verbose = cms.untracked.bool(False)                          
                             )


# Trigger
process.hltHighLevel = cms.EDFilter("HLTHighLevel",
                                    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                    HLTPaths = cms.vstring(
                                         'HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v1',
                                         'HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1',
                                         'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1',
                                         'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1',
                                         'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1',
                                         'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1'
                                    ),
                                    # provide list of HLT paths (or patterns) you want
                                    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                    andOr = cms.bool(True),              # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true  
                                    throw = cms.bool(False)              # throw exception on unknown path names 
                                    )


process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

##clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
				     src = cms.InputTag("slimmedMuons"),
				     preselection = cms.string("track.isNonnull"),
				     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
				     fractionOfSharedSegments = cms.double(0.499),
				     )


process.p = cms.Path(#process.reCorrectedPatJets
                     process.fsrPhotonSequence*
		     process.boostedMuons*
                     process.Ana*
                     process.hltHighLevel*
                     process.AnaAfterHlt
                     )


