// -*- C++ -*-
//
// Package:    UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer
// Class:      KalmanMuonCalibrationsProducer

// system include files
#include <memory>
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// Kalman Muon Corrections 
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

//
// class declaration
//

class KalmanMuonCalibrationsProducer : public edm::EDProducer {
   public:
      explicit KalmanMuonCalibrationsProducer(const edm::ParameterSet&);
      ~KalmanMuonCalibrationsProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      // Kalman Muon Calibrator 
      KalmanMuonCalibrator *kalmanMuonCalibrator;

      edm::EDGetToken muonsCollection_;
      bool isMC;
      bool isSync;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
KalmanMuonCalibrationsProducer::KalmanMuonCalibrationsProducer(const edm::ParameterSet& iConfig)
{

   muonsCollection_ = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonsCollection"));
   isMC = iConfig.getParameter<bool>("isMC");
   isMC = iConfig.getParameter<bool>("isSync");

   if (isMC) {
       kalmanMuonCalibrator = new KalmanMuonCalibrator("MC_76X_13TeV");
   } else {
       kalmanMuonCalibrator = new KalmanMuonCalibrator("DATA_76X_13TeV");
   }

   produces<std::vector<pat::Muon> >();      

}


KalmanMuonCalibrationsProducer::~KalmanMuonCalibrationsProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
KalmanMuonCalibrationsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   edm::Handle<vector<pat::Muon>> muonsCollection;
   iEvent.getByToken(muonsCollection_ , muonsCollection);
   
   // output muons
   std::vector<pat::Muon> * patMuons = new std::vector<pat::Muon>();
   
   // input muons
   const vector<pat::Muon> * theMuons = muonsCollection.product();
   unsigned int nbMuon =  theMuons->size();

   TLorentzVector p4;
   for(unsigned i = 0 ; i < nbMuon; i++){

       pat::Muon mu = theMuons->at(i); 
       
       double oldpt=mu.pt();
       double newpt=oldpt;
       double oldpterr=mu.muonBestTrack()->ptError();
       double newpterr=oldpterr;

       if (!isMC) {
           newpt = kalmanMuonCalibrator->getCorrectedPt(oldpt,mu.eta(),mu.phi(),mu.charge());
           newpterr = kalmanMuonCalibrator->getCorrectedError(oldpt,mu.eta(),oldpterr/oldpt);
       } else {
           double unsmearednewpt = kalmanMuonCalibrator->getCorrectedPt(oldpt,mu.eta(),mu.phi(),mu.charge());
           newpt = kalmanMuonCalibrator->smear(unsmearednewpt,mu.eta());
           newpterr = kalmanMuonCalibrator->getCorrectedErrorAfterSmearing(oldpt,mu.eta(),oldpterr/oldpt);
       }
       
       mu.addUserFloat("correctedPtError",newpterr);
       patMuons->push_back(mu);    
       
       p4.SetPtEtaPhiM(newpt, mu.eta(), mu.phi(), mu.mass());
       patMuons->back().setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mu.mass()));

       //std::cout<<"muon old pt: "<<mu.pt()<<" +/- "<<oldpterr<<" new pt: "<<newpt<<" +/- "<<newpterr<<std::endl;

   }
   
   // add the muons to the event output
   std::auto_ptr<std::vector<pat::Muon> > ptr(patMuons);
   iEvent.put(ptr);
   
}

// ------------ method called once each job just before starting event loop  ------------
void 
KalmanMuonCalibrationsProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KalmanMuonCalibrationsProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
KalmanMuonCalibrationsProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
KalmanMuonCalibrationsProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
KalmanMuonCalibrationsProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
KalmanMuonCalibrationsProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
KalmanMuonCalibrationsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KalmanMuonCalibrationsProducer);
