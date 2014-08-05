import FWCore.ParameterSet.Config as cms

elecId = cms.untracked.string("noEID")
#elecId = cms.untracked.string("eidLoose")
#elecId = cms.untracked.string("eidTight")
#inputFile = 'file:GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT.root'
inputFile = 'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod//GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0881ABEB-2709-E411-9E42-00145EDD7581.root'
outputFile = cms.string("GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT_testAna.root")

#inputFile='https://cms-service-dqm.web.cern.ch//eos/cms/store/cmst3/user/gpetrucc/miniAOD/v1/VBF_HToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT.root'
#outputFile = cms.string("VBF_HToZZTo4L_M-125_13TeV-powheg-pythia6_PU20bx25_PAT_testAna.root")

#inputFile='root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WH_ZH_HToZZ_4LFilter_M-125_13TeV_pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000//0ABFFCEF-AA09-E411-8022-001E673970FD.root'
#outputFile=cms.string("WH_ZH_HToZZ_4LFilter_M-125_13TeV_pythia6_testAna.root");

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

inputFile

]
)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = outputFile
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
                              muonSrc      = cms.untracked.InputTag("slimmedMuons"),
                              #correctedJetSrc = cms.untracked.InputTag("reCorrectedPatJets"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"), #or selectedVertices 
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.uint32(125),
                              CrossSection = cms.untracked.double(0.00544046 ),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              reweightForPU = cms.untracked.bool(True),       
                              doVarDump    = cms.untracked.bool(True),
                              doFsrRecovery = cms.untracked.bool(True),
                              elecID       = elecId
                             )

process.AnaAfterHlt = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              muonSrc      = cms.untracked.InputTag("slimmedMuons"),
                              #correctedJetSrc = cms.untracked.InputTag("reCorrectedPatJets"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"), #or selectedVertices 
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.uint32(125),
                              CrossSection = cms.untracked.double(0.00544046 ),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              reweightForPU = cms.untracked.bool(True),            
                              doVarDump    = cms.untracked.bool(True),
                              doFsrRecovery = cms.untracked.bool(True),
                              elecID       = elecId
                             )




# Trigger
process.hltHighLevel = cms.EDFilter("HLTHighLevel",
                                    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                    HLTPaths = cms.vstring('HLT_Mu17_TkMu8*',
                                                            'HLT_Mu17_Mu8*',
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



process.p = cms.Path(#process.reCorrectedPatJets
                     process.Ana
                     *process.hltHighLevel
                     *process.AnaAfterHlt
                     )






