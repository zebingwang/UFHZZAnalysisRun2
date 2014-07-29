#ifndef HZZ4LPHOTONTREE_H
#define HZZ4LPHOTONTREE_H

//system includes
#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <iomanip>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"

 // user include files 
 #include "FWCore/Framework/interface/Frameworkfwd.h"
 #include "FWCore/Framework/interface/EDAnalyzer.h"
 #include "FWCore/Framework/interface/Event.h"
 #include "FWCore/Framework/interface/MakerMacros.h"

 #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
 #include "DataFormats/PatCandidates/interface/PackedCandidate.h" // for miniAOD

 #include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/ServiceRegistry/interface/Service.h"
 #include "CommonTools/UtilAlgos/interface/TFileService.h"
 #include "TROOT.h"
 #include "TH1.h"
 #include "TH2.h"
 #include "TTree.h"
 #include "TMath.h"
 #include "TString.h"
 #include "TLorentzVector.h"
 #include "FWCore/Framework/interface/LuminosityBlock.h"
 #include "FWCore/Framework/interface/Frameworkfwd.h"
 #include "FWCore/Framework/interface/EDAnalyzer.h"
 #include "FWCore/Framework/interface/Event.h"
 #include "FWCore/Framework/interface/MakerMacros.h"
 #include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/ServiceRegistry/interface/Service.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "Math/VectorUtil.h"
 #include "DataFormats/Common/interface/MergeableCounter.h"
 #include "DataFormats/Math/interface/Point3D.h"
 #include "DataFormats/Common/interface/RefToBase.h"
 #include "DataFormats/Math/interface/deltaR.h"
 #include "DataFormats/Math/interface/deltaPhi.h"

 // PAT
 #include "DataFormats/PatCandidates/interface/Electron.h"
 #include "DataFormats/PatCandidates/interface/Photon.h"
 #include "DataFormats/PatCandidates/interface/Muon.h"
 #include "DataFormats/PatCandidates/interface/Tau.h"
 #include "DataFormats/PatCandidates/interface/Jet.h"
 #include "DataFormats/PatCandidates/interface/MET.h"
 #include "DataFormats/PatCandidates/interface/TriggerEvent.h"
 #include "DataFormats/Provenance/interface/Timestamp.h"
 #include "DataFormats/Candidate/interface/Candidate.h"
 #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
 #include "UFHZZAnalysis8TeV/UFHZZ4LAna/interface/HZZ4LHelper.h"

class HZZ4LPhotonTree
{

 public:

  HZZ4LPhotonTree(TString treeName,edm::Service<TFileService> fs);
  ~HZZ4LPhotonTree();


  void fillPhotonDumpTree(std::vector<pat::PackedCandidate> photons,std::vector<double> deltaRVec, const edm::Event& iEvent, const reco::Vertex *&vertex); // miniAOD

 private:
  
  double Event, Run, LumiSect;
  double pT, eta, phi, pX, pY, pZ;
  //double relIso, isoNH, isoCH, isoCHPU, isoPhot; // miniAOD
  double energy, deltaR;

  HZZ4LHelper myHelper;
  
  TTree *photonTree;
    

};

#endif

#ifndef HZZ4LPHOTONTREE_CC
#define HZZ4LPHOTONTREE_CC


HZZ4LPhotonTree::HZZ4LPhotonTree(TString treeName,edm::Service<TFileService> fs)
{
  TFileDirectory photons = fs->mkdir("photons");
  photonTree = photons.make<TTree>(treeName,treeName);
  photonTree->Branch("Event",&Event,"Event/D");
  photonTree->Branch("Run",&Run,"Run/D");
  photonTree->Branch("Lumi",&LumiSect,"Lumi/D");
  photonTree->Branch("pT",&pT,"pT/D");
  photonTree->Branch("eta",&eta,"eta/D");
  photonTree->Branch("phi",&phi,"phi/D");
  photonTree->Branch("pX",&pX,"pX/D");
  photonTree->Branch("pY",&pY,"pY/D");
  photonTree->Branch("pZ",&pZ,"pZ/D");
  //photonTree->Branch("relIso",&relIso,"relIso/D"); // miniAOD
  //photonTree->Branch("isoNH",&isoNH,"isoNH/D"); // miniAOD
  //photonTree->Branch("isoCH",&isoCH,"isoCH/D"); // miniAOD
  //photonTree->Branch("isoCHPU",&isoCHPU,"isoCHPU/D"); // miniAOD
  //photonTree->Branch("isoPhot",&isoPhot,"isoPhot/D"); // miniAOD
  photonTree->Branch("energy",&energy,"energy/D");
  photonTree->Branch("deltaR",&deltaR,"deltaR/D");

 
}


HZZ4LPhotonTree::~HZZ4LPhotonTree()
{

}


void HZZ4LPhotonTree::fillPhotonDumpTree(std::vector<pat::PackedCandidate> photons,std::vector<double> deltaRVec, // miniAOD
					 const edm::Event& iEvent, const reco::Vertex *&vertex)
{

  using namespace std;

  
  Run = iEvent.id().run();
  Event = iEvent.id().event();
  LumiSect = iEvent.id().luminosityBlock();

  for(unsigned int i = 0; i < photons.size(); i++)
    {
      pT = photons[i].pt();
      eta = photons[i].eta();
      phi = photons[i].phi();
      pX = photons[i].px();
      pY = photons[i].py();
      pZ = photons[i].pz();
      energy = photons[i].energy();
      //isoNH = photons[i].userFloat("fsrPhotonPFIsoNHad03"); // miniAOD
      //isoCH = photons[i].userFloat("fsrPhotonPFIsoChHad03pt02") ; // miniAOD
      //isoCHPU = photons[i].userFloat("fsrPhotonPFIsoChHadPU03pt02") ; // miniAOD
      //isoPhot = photons[i].userFloat("fsrPhotonPFIsoPhoton03") ; // miniAOD
      deltaR = deltaRVec[i];
      //relIso = (photons[i].userFloat("fsrPhotonPFIsoChHad03pt02") + // miniAOD
      //          photons[i].userFloat("fsrPhotonPFIsoNHad03") + // miniAOD
      //          photons[i].userFloat("fsrPhotonPFIsoPhoton03") + // miniAOD
      //          photons[i].userFloat("fsrPhotonPFIsoChHadPU03pt02"))/photons[i].pt();  // miniAOD

      photonTree->Fill();

    }
      
}



#endif
