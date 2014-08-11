// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>

class FSRPhotonProducer : public edm::EDProducer 
{
  public:
    explicit FSRPhotonProducer(const edm::ParameterSet&);
    ~FSRPhotonProducer();

  private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::EDGetTokenT<reco::CandidateView> srcCands_;
    double ptThresh_;
    bool   extractMuonFSR_;
};

FSRPhotonProducer::FSRPhotonProducer(const edm::ParameterSet& iConfig):
  srcCands_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("srcCands"))),
  ptThresh_( iConfig.getParameter<double>("ptThresh") ),
  extractMuonFSR_(iConfig.getParameter<bool>("extractMuonFSR"))
{
  produces<reco::PFCandidateCollection>(); 
}


FSRPhotonProducer::~FSRPhotonProducer()
{
}


void FSRPhotonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<reco::PFCandidateCollection> comp( new reco::PFCandidateCollection );

  edm::Handle<reco::CandidateView> cands;
  iEvent.getByToken(srcCands_, cands);

  for( reco::CandidateView::const_iterator c = cands->begin(); c != cands->end(); ++c ) 
  {
    if (c->charge()==0 && c->pdgId() == 22 && c->pt() > ptThresh_)  
    {
      comp->push_back( reco::PFCandidate(0, c->p4(), reco::PFCandidate::gamma));
      comp->back().setStatus(0);
    }

    // don't extract muon FSR photons for now.
    //if (extractMuonFSR_ && abs(c->pdgId()) == 13) 
    //{
    //  if (c->ecalEnergy() > 0) 
    //  { 
    //    // Proper version, with massless photons (but still using the muon direction)
    //    reco::Particle::PolarLorentzVector p4(c->ecalEnergy() * c->pt()/c->p(), c->eta(), c->phi(), 0.);
    //    // Improper version below, with massive photons (precise copy of Patrick's code)
    //    // reco::Particle::LorentzVector p4 = c->p4(); p4 *= c->ecalEnergy()/c->energy();
    //    if (p4.pt() > ptThresh_) 
    //    {
    //       comp->push_back( reco::PFCandidate(0, reco::Particle::LorentzVector(p4), reco::PFCandidate::gamma) );
    //    }
    //  }
    //}
  }

  iEvent.put( comp );
}


//define this as a plug-in
DEFINE_FWK_MODULE(FSRPhotonProducer);
