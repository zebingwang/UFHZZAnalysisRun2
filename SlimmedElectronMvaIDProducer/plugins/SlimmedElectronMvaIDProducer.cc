// -*- C++ -*-
//
// Package:    UFHZZAnalysisRun2/SlimmedElectronMvaIDProducer
// Class:      SlimmedElectronMvaIDProducer
// 
/**\class SlimmedElectronMvaIDProducer SlimmedElectronMvaIDProducer.cc UFHZZAnalysisRun2/SlimmedElectronMvaIDProducer/plugins/SlimmedElectronMvaIDProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tongguang Cheng
//         Created:  Fri, 03 Apr 2015 13:52:09 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

// electron MvaID related
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

//needed if create ID from reco/aod
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//
// class declaration
//

class SlimmedElectronMvaIDProducer : public edm::EDProducer {
   public:
      explicit SlimmedElectronMvaIDProducer(const edm::ParameterSet&);
      ~SlimmedElectronMvaIDProducer();

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

     edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
     edm::EDGetToken electronsToken_;
     edm::InputTag electronsCollection_;

     string method_;
     vector<string> mvaWeightFiles_;
     bool Trig_;
    

     EGammaMvaEleEstimatorCSA14* mvaID_;

     std::string idname; 

     GreaterByPt<pat::Electron> pTComparator_;

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
SlimmedElectronMvaIDProducer::SlimmedElectronMvaIDProducer(const edm::ParameterSet& iConfig):
   mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap")))
{

   electronsCollection_ = iConfig.getParameter<edm::InputTag>("electronsCollection");
   electronsToken_ = mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronsCollection"));
   Trig_ = iConfig.getParameter<bool>("Trig");

    if(Trig_){ idname = "Trig";}    
    if(!Trig_){ idname = "NonTrig";    }

   //register your products

   produces<edm::ValueMap<float> >();
   produces<std::vector<pat::Electron> >(idname);      

}


SlimmedElectronMvaIDProducer::~SlimmedElectronMvaIDProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SlimmedElectronMvaIDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   edm::Handle<vector<pat::Electron>> electronsCollection;
   iEvent.getByLabel(electronsCollection_ , electronsCollection);

   edm::Handle<edm::View<reco::GsfElectron> > gsfelectrons;
   iEvent.getByToken(electronsToken_,gsfelectrons);
        
   // electron mva values            
   edm::Handle<edm::ValueMap<float> > mvaValues;
   iEvent.getByToken(mvaValuesMapToken_,mvaValues);
   
   // output electrons
   std::vector<pat::Electron> * patElectrons = new std::vector<pat::Electron>();
   
   // output valuemap
   std::auto_ptr<edm::ValueMap<float> > ID(new edm::ValueMap<float>() );

   std::vector<float> values;

   // input electrons
   const vector<pat::Electron> * theElectrons = electronsCollection.product();
   unsigned int nbElectron =  theElectrons->size();
   values.reserve(theElectrons->size());    

   for(unsigned i = 0 ; i < nbElectron; i++){

        const auto gsf = gsfelectrons->ptrAt((size_t)i);

        float idvalue = (*mvaValues)[gsf]; 

        pat::Electron anElectron = theElectrons->at(i); 

        std::vector<pat::Electron::IdPair> ids;
        pat::Electron::IdPair id;
        //std::pair<std::string,float> id;
        id.first  =  idname;    
        id.second =  idvalue;
        
        ids.push_back(id);     
        anElectron.setElectronIDs(ids);
        
        patElectrons->push_back(anElectron);    
        
        values.push_back( idvalue );
    }
    
    // add the value map to the input electron collection
    edm::ValueMap<float>::Filler filler(*ID);
    filler.insert(electronsCollection, values.begin(), values.end() );
    
    filler.fill();
    iEvent.put(ID);
    
    // sort electrons in pt
    //std::sort(patElectrons->begin(), patElectrons->end(), pTComparator_);
    
    // add the electrons to the event output
    std::auto_ptr<std::vector<pat::Electron> > ptr(patElectrons);
    iEvent.put(ptr,idname);
    
}

// ------------ method called once each job just before starting event loop  ------------
void 
SlimmedElectronMvaIDProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimmedElectronMvaIDProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
SlimmedElectronMvaIDProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SlimmedElectronMvaIDProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SlimmedElectronMvaIDProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SlimmedElectronMvaIDProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SlimmedElectronMvaIDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SlimmedElectronMvaIDProducer);
