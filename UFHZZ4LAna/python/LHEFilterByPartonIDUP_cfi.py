import FWCore.ParameterSet.Config as cms

LHEFilterByPartonIDUP = cms.EDFilter('LheIDupFilter',
  skipLheIDup = cms.untracked.int32(5), # pdgID of particles to skip (e.g. |pdgID|=5 for b quarks)
  numLheIDup = cms.untracked.int32(2)  # required number of particles with given pdgID in the event
)
