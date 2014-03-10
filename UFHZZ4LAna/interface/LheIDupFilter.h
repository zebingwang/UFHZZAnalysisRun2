#ifndef LheIDupFilter_h
#define LheIDupFilter_h
// Filter to skip events with given number of particles with given id LHE IDUP, all the other events should pass

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


class LheIDupFilter : public edm::EDFilter 
{
  public:
  explicit LheIDupFilter(const edm::ParameterSet&);
  ~LheIDupFilter();

  private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  int skipLheIDup_;
  int numLheIDup_;
  int nSkippedLheIDup;
  int nScannedLhe;
  int nScanned;
	
};

#endif
