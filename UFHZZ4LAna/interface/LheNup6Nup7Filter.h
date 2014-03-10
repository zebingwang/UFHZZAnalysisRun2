#ifndef LheNup6Nup7Filter_h
#define LheNup6Nup7Filter_h
// Filter to skip events with LHE NUP==6 || LHE NUP==7, all the other events should be selected

// system include files
#include <memory>

// other include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"


class LheNup6Nup7Filter : public edm::EDFilter 
{
  public:
  explicit LheNup6Nup7Filter(const edm::ParameterSet&);
  ~LheNup6Nup7Filter();

  private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  bool skipLheNup6Nup7_;
  int nSkippedLheNup6Nup7;
  int nScannedLhe;
  int nScanned;
	
};

#endif
