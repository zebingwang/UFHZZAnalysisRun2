#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "TMath.h"

class PhotonPFIsoCalculator : public edm::EDProducer {
public:
  explicit PhotonPFIsoCalculator(const edm::ParameterSet&);
  ~PhotonPFIsoCalculator();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::EDGetTokenT<reco::CandidateView> leptonLabel_;
  edm::EDGetTokenT<reco::CandidateView> pfLabel_;
  StringCutObjectSelector<reco::Candidate> endcapDefinition_;
  StringCutObjectSelector<reco::Candidate> pfSelection_;
  double deltaR_, deltaRself_, deltaZ_;
  double directional_;
  double vetoConeEndcaps_;
  bool   debug_;

};



PhotonPFIsoCalculator::PhotonPFIsoCalculator(const edm::ParameterSet& iConfig) :
  leptonLabel_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("leptonLabel"))),
  pfLabel_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("pfLabel"))),
  endcapDefinition_(iConfig.existsAs<std::string>("endcapDefinition") ? iConfig.getParameter<std::string>("endcapDefinition") : "abs(eta) > 1.479", true),
  pfSelection_(iConfig.getParameter<std::string>("pfSelection"), true),
  deltaR_(iConfig.getParameter<double>("deltaR")),
  deltaRself_(iConfig.getParameter<double>("deltaRself")),
  deltaZ_(iConfig.existsAs<double>("deltaZ") ? iConfig.getParameter<double>("deltaZ") : 0),
  directional_(iConfig.getParameter<bool>("directional")),
  vetoConeEndcaps_(iConfig.getParameter<double>("vetoConeEndcaps")),
  debug_(iConfig.getUntrackedParameter<bool>("debug",false))
{
  produces<edm::ValueMap<float> >().setBranchAlias("pfMuIso");
}

void PhotonPFIsoCalculator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  edm::Handle<reco::CandidateView> photons;
  iEvent.getByToken(leptonLabel_,photons);

  edm::Handle<reco::CandidateView> pfH;
  iEvent.getByToken(pfLabel_,pfH);

  std::vector<float> isoV;
  std::auto_ptr<edm::ValueMap<float> > isoM(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler isoF(*isoM);

  if (debug_) std::cout << "Run " << iEvent.id().run() << ", Event " << iEvent.id().event() << std::endl;

  for( reco::CandidateView::const_iterator ph = photons->begin(); ph != photons->end(); ++ph )
  {
    if (debug_) std::cout << leptonLabel_.encode() << " with pt = " << ph->pt() << ", eta = " << ph->eta() << ", phi = " << ph->phi() << std::endl;

    Double_t ptSum =0.;  
    math::XYZVector isoAngleSum;
    std::vector<math::XYZVector> coneParticles;

    for( reco::CandidateView::const_iterator pf = pfH->begin(); pf!=pfH->end(); ++pf )
    {   

      double dr = deltaR(pf->p4(), pf->p4()) ;
      if (dr >= deltaR_) continue;
        
      if (!pfSelection_(*pf)) continue; 

      if (debug_) 
        std::cout << "   pfCandidate of pdgId " << pf->pdgId() 
          << ", pt = " << pf->pt() << ", dr = " << dr << ", dz = " << (pf->vz() - ph->vz()) 
          << " is in cone... " << std::endl;

      if (deltaZ_ > 0 && fabs(pf->vz() - ph->vz()) > deltaZ_) continue;

      if (deltaR(pf->p4(), ph->p4()) < deltaRself_) continue;

      // dR Veto for Gamma: no-one in EB, dR > 0.08 in EE
      if (vetoConeEndcaps_ > 0 && endcapDefinition_(*ph) && dr < vetoConeEndcaps_) continue;

      if (debug_) std::cout << "          ...and passes all vetos, so it's added to the sum." << std::endl;
      // scalar sum
      ptSum += pf->pt();

      // directional sum
      math::XYZVector transverse( pf->eta() - ph->eta()
              , reco::deltaPhi(pf->phi(), ph->phi())
              , 0);
      transverse *= pf->pt() / transverse.rho();
      if (transverse.rho() > 0) 
      {
          isoAngleSum += transverse;
          coneParticles.push_back(transverse);
      }
    }
    if (directional_) 
    {
      double directionalPT = 0;
      for (unsigned int iPtcl = 0; iPtcl < coneParticles.size(); ++iPtcl)
      {
        directionalPT += pow(TMath::ACos( coneParticles[iPtcl].Dot(isoAngleSum) / coneParticles[iPtcl].rho() / isoAngleSum.rho() ),2) * coneParticles[iPtcl].rho();
      }
      isoV.push_back(directionalPT);
    } 
    else isoV.push_back(ptSum);
  }

  isoF.insert(photons,isoV.begin(),isoV.end());

  isoF.fill();
  iEvent.put(isoM);

}

PhotonPFIsoCalculator::~PhotonPFIsoCalculator() { }
void PhotonPFIsoCalculator::beginJob() { }
void PhotonPFIsoCalculator::endJob() { }
DEFINE_FWK_MODULE(PhotonPFIsoCalculator);
