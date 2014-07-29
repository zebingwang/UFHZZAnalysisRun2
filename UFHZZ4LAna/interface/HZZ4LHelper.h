#ifndef HZZ4LHELPER_H
#define HZZ4LHELPER_H

#define PI 3.14159

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
#include "TLorentzRotation.h"

// user include files 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
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
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//Boost
#include "CommonTools/CandUtils/interface/CenterOfMassBooster.h"
#include "CommonTools/CandUtils/interface/Booster.h"



class HZZ4LHelper
{

 public:

  HZZ4LHelper();
  ~HZZ4LHelper();
  
  std::vector< pat::Muon > goodMuons2012_Iso(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut, double rho, double isocut, const reco::Vertex *&vertex);
  std::vector< pat::Electron > goodElectrons2012_Iso(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, double rho, double isocut, std::string elecID, const reco::Vertex *&vertex);
  
  std::vector< pat::Muon > goodMuons2012_Iso(std::vector< pat::Muon > Muons, double muPtCut, double rho, double isocut, const reco::Vertex *&vertex);
  std::vector< pat::Electron > goodElectrons2012_Iso(std::vector< pat::Electron > Electrons, double elecPtCut, double rho, double isocut, std::string elecID, const reco::Vertex *&vertex);

  std::vector<pat::Muon> allMuons(edm::Handle<edm::View<pat::Muon> > Muons);
  std::vector<pat::Electron> allElectrons(edm::Handle<edm::View<pat::Electron> > Electrons);
  std::vector<pat::Muon> goodBaselineMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut);
  std::vector<pat::Electron> goodBaselineElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID);
  std::vector<pat::Muon> goodMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut);
  std::vector<pat::Electron> goodElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID);
  std::vector<pat::Muon> goodLooseMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut);
  std::vector<pat::Electron> goodLooseElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID);
  std::vector<pat::Muon> goodTightMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut);
  std::vector<pat::Electron> goodTightElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID);
  std::vector<pat::Muon> goodMuonsAllCuts(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut);
  std::vector<pat::Electron> goodElectronsAllCuts(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID);

  std::vector<pat::Muon> goodMuons2012(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut, double muIsoCut,double Rho, std::string isoType, std::string PUCorrType,const reco::Vertex *&vertex);
  std::vector<pat::Electron> goodElectrons2012(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID, double elecIsoCut,double Rho, std::string isoType, std::string PUCorrType, const reco::Vertex *&vertex);

  std::vector<pat::Muon> goodMuons2012_noIso(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut, const reco::Vertex *&vertex);
  std::vector<pat::Electron> goodElectrons2012_noIso(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID, const reco::Vertex *&vertex);

  std::vector<pat::Muon> goodMuons2012_noIso(std::vector<pat::Muon> Muons, double muPtCut, const reco::Vertex *&vertex);
  std::vector<pat::Electron> goodElectrons2012_noIso(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID, const reco::Vertex *&vertex,const edm::Event& iEvent);

  std::vector<pat::Muon> goodLooseMuons2012(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut);
  std::vector<pat::Electron> goodLooseElectrons2012(edm::Handle<edm::View<pat::Electron> > Electrons, double elPtCut);
  

  std::vector<pat::Muon> goodLooseMuons(std::vector<pat::Muon> Muons, double muPtCut);
  std::vector<pat::Electron> goodLooseElectrons(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID);
  std::vector<pat::Muon> goodTightMuons(std::vector<pat::Muon>  Muons, double muPtCut);
  std::vector<pat::Electron> goodTightElectrons(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID);



  void cleanOverlappingLeptons(std::vector<pat::Muon> &Muons, std::vector<pat::Electron> &Electrons,const reco::Vertex *&vertex);
  void m2lCutAt4GeV(std::vector<pat::Muon> &Muons, std::vector<pat::Electron> &Electrons);
  void m2lCutAt4GeV(std::vector<pat::Muon> inputMuons, std::vector<pat::Electron> inputElectrons, std::vector<pat::Muon> &outputMuons, std::vector<pat::Electron> &outputElectrons);
  bool passedM2lCut(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double lowM2lCut, double highM2lCut);
  bool passedM2lCut_OS(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double lowM2lCut, double highM2lCut);
  bool passedM2lAndPtCut(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double lowM2lCut, double highM2lCut);
  double minMass2l(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, bool &sameFlavor, bool &sameSign);
  double maxMass2l(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons );

  double minMass3l(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons);
  double M3lSoftestLep(TLorentzVector l1, TLorentzVector l2, TLorentzVector l3, TLorentzVector l4);

  double largestLepMomentum(TLorentzVector l1, TLorentzVector l2, TLorentzVector l3, TLorentzVector l4);

  double getSumOf2LeastIso(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr);
  double getLeastIso1(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr);
  double getLeastIso1(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr, double &pt);
  double getLeastIso2(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr);

  int passedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut, double muonRhoCorr, double elecRhoCorr);
  int passedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut, TString isoType, TString muPUCorr, TString elPUCorr, double Rho);
  int passedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut, TString isoType, TString muPUCorr, TString elPUCorr, double muonRho, double elecRho, double &leastIso1, double &leastIso2 );
  
    
  int passedEarlyIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr, double earlyIsoCut);
  int passedInvertedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr,
			      double chosenIsoCut);
  int passedInvertedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut,
			      TString isoType, TString muPUCorr, TString elPUCorr, double Rho);


  bool passedPtandEarlyIso(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double leadingPtCut, double subleadingPtCut,
			   double earlyIsoCut, double muonRhoCorr, double elecRhoCorr );
  bool passedPtandEarlyIso(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double leadingPtCut, double subleadingPtCut,
			   double earlyIsoCut, TString isoType, TString muPUCorr, TString elPUCorr, double Rho );
  bool passedPtIsoAndM2l(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double leadingPtCut, double subleadingPtCut,
			 double earlyIsoCut, TString isoType, TString muPUCorr, TString elPUCorr, double Rho );

  double getSIP3D(pat::Muon muon);
  double getSIP3D(pat::Electron electron);
  TString getStringFromDouble(double number);

  double getCorrWeightForSChan(double m4l);
  bool isPFMuon(pat::Muon muon,edm::Handle<edm::View<reco::PFCandidate> > pf);
  double pfIso(pat::Muon muon, double Rho);
  double pfIso(pat::Electron elec, double Rho);
  double pfIsoFSR(pat::Muon muon, double Rho, double photPt);
  double pfIsoFSR(pat::Electron elec, double Rho, double photPt);


  enum MuonEffectiveAreaType {
    kMuTrkIso03, 
    kMuEcalIso03, 
    kMuHcalIso03, 
    kMuTrkIso05, 
    kMuEcalIso05, 
    kMuHcalIso05, 
    kMuChargedIso03, 
    kMuGammaIso03, 
    kMuNeutralHadronIso03, 
    kMuGammaAndNeutralHadronIso03,
    kMuGammaIso03Tight, 
    kMuNeutralHadronIso03Tight, 
    kMuGammaAndNeutralHadronIso03Tight,
    kMuChargedIso04, 
    kMuGammaIso04, 
    kMuNeutralHadronIso04, 
    kMuGammaAndNeutralHadronIso04,
    kMuGammaIso04Tight, 
    kMuNeutralHadronIso04Tight, 
    kMuGammaAndNeutralHadronIso04Tight,
    kMuGammaIsoDR0p0To0p1,
    kMuGammaIsoDR0p1To0p2,
    kMuGammaIsoDR0p2To0p3,
    kMuGammaIsoDR0p3To0p4,
    kMuGammaIsoDR0p4To0p5,
    kMuNeutralHadronIsoDR0p0To0p1,
    kMuNeutralHadronIsoDR0p1To0p2,
    kMuNeutralHadronIsoDR0p2To0p3,
    kMuNeutralHadronIsoDR0p3To0p4,
    kMuNeutralHadronIsoDR0p4To0p5,
    kMuGammaIso05,
    kMuNeutralIso05
  };
  
  enum MuonEffectiveAreaTarget {
    kMuNoCorr,
    kMuEAData2011,
    kMuEASummer11MC,
    kMuEAFall11MC,
    kMuEAData2012
  };

  enum ElecEffectiveAreaType {
    kEGammaNeutralHadIso04,
    kENeutralHadIso04,
    kEGammaNeutralIso04,
    kEGammaIso04,
    kENeutralHadIso03,
    kEGammaNeutralHadIso03,
    kEGammaIso03    
  };

  enum ElecEffectiveAreaTarget {
    kElNoCorr,
    kElEAData2012,
    kElEAData2011,
    kElEAFall11MC
  };

  double MuonEffArea(MuonEffectiveAreaType type, double SCEta, MuonEffectiveAreaTarget EffectiveAreaTarget);
  double ElecEffArea(ElecEffectiveAreaType type, double SCEta, ElecEffectiveAreaTarget EffectiveAreaTarget);

  MuonEffectiveAreaType muEAtype;
  MuonEffectiveAreaTarget muEAtarget;
  ElecEffectiveAreaType elEAtype;
  ElecEffectiveAreaTarget elEAtarget;

  /*
  void getMuonEffAreaType(MuonEffectiveAreaType type);
  void getMuonEffAreaTarget(MuonEffectiveAreaTarget target);
  void getElecEffAreaType(ElecEffectiveAreaType type);
  void getElecEffAreaTarget(ElecEffectiveAreaTarget target);
  */

};


#endif


#ifndef HZZ4LHELPER_CC
#define HZZ4LHELPER_CC

HZZ4LHelper::HZZ4LHelper()
{

  muEAtype = kMuGammaAndNeutralHadronIso04;
  muEAtarget = kMuEAData2012;
  elEAtype = kEGammaNeutralHadIso04;
  elEAtarget = kElEAData2012;


  //declarations

}


HZZ4LHelper::~HZZ4LHelper()
{

  //destructor ---do nothing

}



std::vector<pat::Muon> HZZ4LHelper::allMuons(edm::Handle<edm::View<pat::Muon> > Muons)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  
  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      bestMuons.push_back(*mu);
    }
  
  return bestMuons;
  
}





//Returns vector of good electrons for analysis
std::vector<pat::Electron> HZZ4LHelper::allElectrons(edm::Handle<edm::View<pat::Electron> > Electrons)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      bestElectrons.push_back(*elec);
    }
  
  
  return bestElectrons;
  
  
}



//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodBaselineMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  //double muChi2Cut = 10.0;
  //int muMuonHits = 0;
  //int muNumMatches = 1;
  int muTrackerHits = 10;
  /**************************************/
  
  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( abs(getSIP3D(*mu)) < 100 )
	{
	  if( (double)mu->trackIso()/mu->pt() < 0.7 )
	    { 
	      if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isGlobalMuon() == 1 && mu->track()->numberOfValidHits() > muTrackerHits )
		{
		  bestMuons.push_back(*mu);
		}
	    }
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodBaselineElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  //int missingHitsCut = 2; // for miniAOD
  /**************************************/
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(getSIP3D(*elec)) < 100 )
	{
	  if( (double)elec->dr03TkSumPt()/elec->pt() < 0.7 )
	    {
	      if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
		{
		  //if( elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < missingHitsCut ) // for miniAOD
		    //{
		    if (elec->electronID(elecID) == 1 || elec->electronID(elecID) == 3 || elec->electronID(elecID) == 5 || elec->electronID(elecID) == 7 ||
		       elec->electronID(elecID) == 9 || elec->electronID(elecID) == 11 || elec->electronID(elecID) == 13 || elec->electronID(elecID) == 15)
		      {
			bestElectrons.push_back(*elec);
		      }
		    //}
		}
	    }
	}
    }
  
  
  return bestElectrons;
  
  
}


//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  double muChi2Cut = 10.0;
  int muMuonHits = 0;
  int muNumMatches = 1;
  //int muTrackerHits = 10;
  /**************************************/
  
  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isGlobalMuon() == 1 && mu->isTrackerMuon() == 1 )
	{
	  if( mu->normChi2() < muChi2Cut && mu->globalTrack()->hitPattern().numberOfValidMuonHits() > muMuonHits)
	    {
	      if( mu->numberOfMatches() > muNumMatches )
	  	{
		  bestMuons.push_back(*mu);
	  	}
	    }
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  // int missingHitsCuts = 2; // for miniAOD
  /**************************************/
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
	{
	  //if( elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < missingHitsCuts ) // for miniAOD
	    //{
	      //if (elec->electronID(elecID) == 1 || elec->electronID(elecID) == 3 || elec->electronID(elecID) == 5 || elec->electronID(elecID) == 7 ||
	      //  elec->electronID(elecID) == 9 || elec->electronID(elecID) == 11 || elec->electronID(elecID) == 13 || elec->electronID(elecID) == 15)
	      if(elec->electronID(elecID) == 5 || elec->electronID(elecID) == 7 || elec->electronID(elecID) == 13 || elec->electronID(elecID) == 15)
		{
		  bestElectrons.push_back(*elec);
		}
	    //}
	}
    }
  
  
  return bestElectrons;
  
  
}




//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodLooseMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  int muTrackerHits = 10;
  /**************************************/
  
  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( abs(getSIP3D(*mu)) < 100 )
	{
	  if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isGlobalMuon() == 1 && mu->track()->numberOfValidHits() > muTrackerHits )
	    {
	      bestMuons.push_back(*mu);
	    }
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodLooseElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  // int missingHitsCuts = 2; // for miniAOD
  /**************************************/
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(getSIP3D(*elec)) < 100 )
	{
	  if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
	    {
	      //if( elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < missingHitsCuts ) // for miniAOD
		//{
		  if (elec->electronID(elecID) == 1 || elec->electronID(elecID) == 3 || elec->electronID(elecID) == 5 || elec->electronID(elecID) == 7 ||
		      elec->electronID(elecID) == 9 || elec->electronID(elecID) == 11 || elec->electronID(elecID) == 13 || elec->electronID(elecID) == 15)
		    {
		      bestElectrons.push_back(*elec);
		    }
		//}
	    }
	}
    }
  
  return bestElectrons;
  
  
}



//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodTightMuons(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  double muChi2Cut = 10.0;
  int muMuonHits = 0;
  int muNumMatches = 1;
  int muTrackerHits = 10;
  /**************************************/
  
  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( abs(getSIP3D(*mu)) < 100 )
	{
	  if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isGlobalMuon() == 1 && mu->isTrackerMuon() == 1 )
	    {
	      if( mu->normChi2() < muChi2Cut && mu->globalTrack()->hitPattern().numberOfValidMuonHits() > muMuonHits &&  mu->track()->numberOfValidHits() > muTrackerHits)
		{
		  if( mu->numberOfMatches() > muNumMatches )
		    {
		      bestMuons.push_back(*mu);
		    }
		}
	    }
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodTightElectrons(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  // int missingHitsCuts = 2; // for miniAOD
  /**************************************/
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(getSIP3D(*elec)) < 100 )
	{
	  if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
	    {
	      //if( elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < missingHitsCuts ) // for miniAOD
		//{
		  if(elec->electronID(elecID) == 5 || elec->electronID(elecID) == 7 || elec->electronID(elecID) == 13 || elec->electronID(elecID) == 15)
		    {
		      bestElectrons.push_back(*elec);
		    }
		//}
	    }
	}
    }
  
  
  return bestElectrons;
  
  
}


///////////////////////////////////////////////////////////


//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodMuonsAllCuts(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  int muTrackerHits = 10;
  double sipCut = 4.0;
  /**************************************/
  
  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( abs(getSIP3D(*mu)) < sipCut )
	{
	  if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isGlobalMuon() == 1 && mu->track()->numberOfValidHits() > muTrackerHits )
	    {
	      bestMuons.push_back(*mu);
	    }
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodElectronsAllCuts(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  //int missingHitsCuts = 2; // for miniAOD
  double sipCut = 4.0;
  /**************************************/
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(getSIP3D(*elec)) < sipCut )
	{
	  if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
	    {
	      //if( elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < missingHitsCuts ) // for miniAOD
		//{
		  if (elec->electronID(elecID) == 1 || elec->electronID(elecID) == 3 || elec->electronID(elecID) == 5 || elec->electronID(elecID) == 7 ||
		      elec->electronID(elecID) == 9 || elec->electronID(elecID) == 11 || elec->electronID(elecID) == 13 || elec->electronID(elecID) == 15)
		    {
		      bestElectrons.push_back(*elec);
		    }
		//}
	    }
	}
    }
  
  return bestElectrons;
  
  
}


//////////////////////////////////////////////////////////////

/// 2012
//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodMuons2012(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut, double muIsoCut, double Rho, std::string isoType, std::string PUCorrType, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  //int muTrackerHits = 10;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/
  double PUCorr = 0;


  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isPFMuon() == 1 && (mu->isGlobalMuon() || mu->isTrackerMuon() ) )
	{
	  if( abs(getSIP3D(*mu)) < sipCut )
	    {
	      if( fabs(mu->innerTrack()->dxy(vertex->position())) < dxyCut )
		{ 
		  if( fabs(mu->innerTrack()->dz(vertex->position())) < dzCut )
		    {
		      
		      double iso = 999;
		      if( isoType == "PF" && PUCorrType == "dB"){
			iso = (mu->chargedHadronIso()+max(mu->photonIso()+mu->neutralHadronIso()-0.5*mu->userIsolation("PfPUChargedHadronIso"),0.0))/mu->pt();
			
		      }
		      else if( isoType == "PF" && PUCorrType == "PFEffAreaRho"){
			PUCorr = Rho*MuonEffArea(muEAtype,mu->eta(),muEAtarget);
			iso = (mu->chargedHadronIso()+max(mu->photonIso()+mu->neutralHadronIso()-PUCorr,0.0))/mu->pt();
			
		      }
		      else{ 
			string err = "HZZ4LHelper::goodMuons2012 --- unknown isolation type\n";
			err+="Possibilites are:\n";
			err+="\t PF  with PFEffAreaRho or dB\n";
			throw invalid_argument( err );
		      }
		      if(iso < muIsoCut)
			{
			  bestMuons.push_back(*mu);
			}
		    }
		}
	    }
	}
    }


  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodElectrons2012(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID, double elecIsoCut,double Rho, std::string isoType, std::string PUCorrType, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  int missingHitsCuts = 2;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/
  double PUCorr = 0;

  /*
  ////////Emanuele////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = -0.12;
  double bdtCutLoPt_14_25 = 0.39;
  double bdtCutHiPt_0_08 = 0.49;
  double bdtCutHiPt_08_14 = 0.046;
  double bdtCutHiPt_14_25 = 0.69;
  */

  /*
  ////////Si+Duncan////////
  double bdtCutLoPt_0_08 = 0.369;
  double bdtCutLoPt_08_14 = -0.025;
  double bdtCutLoPt_14_25 = 0.531;
  double bdtCutHiPt_0_08 = 0.735;
  double bdtCutHiPt_08_14 = 0.467;
  double bdtCutHiPt_14_25 = 0.795;
  */
  
  ////////Christophe////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = 0.004;
  double bdtCutLoPt_14_25 = 0.295;
  double bdtCutHiPt_0_08 = 0.5;
  double bdtCutHiPt_08_14 = 0.12;
  double bdtCutHiPt_14_25 = 0.6;
  
  ///MVA ID = mvaNonTrigV0
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(getSIP3D(*elec)) < sipCut && fabs(elec->gsfTrack()->dxy(vertex->position())) < dxyCut && fabs(elec->gsfTrack()->dz(vertex->position())) < dzCut )
	{
	  if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
	    { 
	      double iso = 999;
              if( isoType == "PF" && PUCorrType == "dB"){
		iso = (elec->chargedHadronIso()+max(elec->photonIso()+elec->neutralHadronIso()-0.5*elec->userIsolation("PfPUChargedHadronIso"),0.0))/elec->pt();
	      }
	      else if( isoType == "PF" && PUCorrType == "PFEffAreaRho"){
		PUCorr = Rho*ElecEffArea(elEAtype,elec->superCluster()->eta(),elEAtarget);
		iso = (elec->chargedHadronIso()+max(elec->photonIso()+elec->neutralHadronIso()-PUCorr,0.0))/elec->pt();
	      }
              else{
                string err = "HZZ4LHelper::goodElectrons2012 --- unknown isolation type\n";
                err+="Possibilites are:\n";
                err+="\t PF  with PFEffAreaRho or dB\n";
                throw invalid_argument( err );
              }
	      
	      if(iso < elecIsoCut)
		{
		  double cutVal = 1000;
		  if( elec->pt() > 0 && elec->pt() < 10.0 )
		    {
		      if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutLoPt_0_08;}
		      if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutLoPt_08_14;}
		      if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutLoPt_14_25;}
		    }
		  if( elec->pt() >= 10.0 )
		    {
		      if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutHiPt_0_08;}
		      if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutHiPt_08_14;}
		      if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutHiPt_14_25;}
		    }
		  if(cutVal == 1000){ cout << "Problem in HZZ4LHelper::goodElectrons2012!" << endl; continue;}
		  //int misHits = elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
		  int misHits = 0; // for miniAOD
                  if(elec->electronID(elecID) > cutVal && misHits < missingHitsCuts) 
		    {
		      bestElectrons.push_back(*elec);
		    }
		}
	    }
	}
    }
  
  return bestElectrons;
  
  
}


/////////////////////////////////////////////////////

std::vector<pat::Muon> HZZ4LHelper::goodMuons2012_noIso(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/


  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isPFMuon() == 1 && (mu->isGlobalMuon() || mu->isTrackerMuon() ) )
	{
	  if( abs(getSIP3D(*mu)) < sipCut )
	    {
	      if( fabs(mu->innerTrack()->dxy(vertex->position())) < dxyCut )
		{ 
		  if( fabs(mu->innerTrack()->dz(vertex->position())) < dzCut )
		    {
		      bestMuons.push_back(*mu);
		    }
		}
	    }
	}
    }


  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodElectrons2012_noIso(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, std::string elecID, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  int missingHitsCuts = 2;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/

  /*
  ////////Emanuele////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = -0.12;
  double bdtCutLoPt_14_25 = 0.39;
  double bdtCutHiPt_0_08 = 0.49;
  double bdtCutHiPt_08_14 = 0.046;
  double bdtCutHiPt_14_25 = 0.69;
  */

  /*
  ////////Si+Duncan////////
  double bdtCutLoPt_0_08 = 0.369;
  double bdtCutLoPt_08_14 = -0.025;
  double bdtCutLoPt_14_25 = 0.531;
  double bdtCutHiPt_0_08 = 0.735;
  double bdtCutHiPt_08_14 = 0.467;
  double bdtCutHiPt_14_25 = 0.795;
  */
  
  ////////Moriond////////
  /*
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = 0.004;
  double bdtCutLoPt_14_25 = 0.295;
  double bdtCutHiPt_0_08 = 0.5;
  double bdtCutHiPt_08_14 = 0.12;
  double bdtCutHiPt_14_25 = 0.6;
  */
  
  //Legacy
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = 0.004;
  double bdtCutLoPt_14_25 = 0.295;
  double bdtCutHiPt_0_08 = -0.34;
  double bdtCutHiPt_08_14 = -0.65;
  double bdtCutHiPt_14_25 = 0.6;
  
  ///MVA ID = mvaNonTrigV0
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(getSIP3D(*elec)) < sipCut && fabs(elec->gsfTrack()->dxy(vertex->position())) < dxyCut && fabs(elec->gsfTrack()->dz(vertex->position())) < dzCut )
	{
	  if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
	    { 
	      //cout << "mva: " << elec->electronID(elecID) << endl;
	      double cutVal = 1000;
	      if( elec->pt() > 0 && elec->pt() < 10.0 )
		{
		  if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutLoPt_0_08;}
		  if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutLoPt_08_14;}
		  if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutLoPt_14_25;}
		}
	      if( elec->pt() >= 10.0 )
		{
		  if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutHiPt_0_08;}
		  if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutHiPt_08_14;}
		  if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutHiPt_14_25;}
		}
	      if(cutVal == 1000){ cout << "ele->superCluster()->eta() > 2.5" << endl; continue;}
              //int misHits = elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); //for miniAOD
              int misHits = 0;
	      if(elec->electronID(elecID) > cutVal && misHits < missingHitsCuts)
		{
		  bestElectrons.push_back(*elec);
		}
	    }
	}
    }

  
  return bestElectrons;
  
  
}




std::vector<pat::Muon> HZZ4LHelper::goodMuons2012_noIso(std::vector<pat::Muon> Muons, double muPtCut, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/


  for(unsigned int i = 0; i < Muons.size(); i++)
    {

      //cout << Muons[i].pt() << "  " << Muons[i].eta() << "  " << Muons[i].isPFMuon() << endl;

      if( Muons[i].pt() > muPtCut && abs(Muons[i].eta()) < muEtaCut && Muons[i].isPFMuon() == 1 && (Muons[i].isGlobalMuon() || Muons[i].isTrackerMuon() ) )
	{
	  if( abs(getSIP3D(Muons[i])) < sipCut )
	    {
	      if( fabs(Muons[i].innerTrack()->dxy(vertex->position())) < dxyCut )
		{ 
		  if( fabs(Muons[i].innerTrack()->dz(vertex->position())) < dzCut )
		    {
		      bestMuons.push_back(Muons[i]);
		    }
		}
	    }
	}
    }


  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodElectrons2012_noIso(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID, const reco::Vertex *&vertex,const edm::Event& iEvent)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  int missingHitsCuts = 2;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/

  /*
  ////////Emanuele////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = -0.12;
  double bdtCutLoPt_14_25 = 0.39;
  double bdtCutHiPt_0_08 = 0.49;
  double bdtCutHiPt_08_14 = 0.046;
  double bdtCutHiPt_14_25 = 0.69;
  */

  /*
  ////////Si+Duncan////////
  double bdtCutLoPt_0_08 = 0.369;
  double bdtCutLoPt_08_14 = -0.025;
  double bdtCutLoPt_14_25 = 0.531;
  double bdtCutHiPt_0_08 = 0.735;
  double bdtCutHiPt_08_14 = 0.467;
  double bdtCutHiPt_14_25 = 0.795;
  */
  
  ////////Christophe////////
  /*
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = 0.004;
  double bdtCutLoPt_14_25 = 0.295;
  double bdtCutHiPt_0_08 = 0.5;
  double bdtCutHiPt_08_14 = 0.12;
  double bdtCutHiPt_14_25 = 0.6;
  */

  //Legacy
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = 0.004;
  double bdtCutLoPt_14_25 = 0.295;
  double bdtCutHiPt_0_08 = -0.34;
  double bdtCutHiPt_08_14 = -0.65;
  double bdtCutHiPt_14_25 = 0.6;
  
  ///MVA ID = mvaNonTrigV0
  
  for(unsigned int i = 0; i < Electrons.size(); i++)
    {
      if( abs(getSIP3D(Electrons[i])) < sipCut && fabs(Electrons[i].gsfTrack()->dxy(vertex->position())) < dxyCut && fabs(Electrons[i].gsfTrack()->dz(vertex->position())) < dzCut )
	{
	  if( Electrons[i].pt() > elecPtCut && abs(Electrons[i].eta()) < elecEtaCut)
	    { 
              //cout << "mva: " << Electrons[i].electronID(elecID) << endl;

	      double cutVal = 1000;
	      if( Electrons[i].pt() > 0 && Electrons[i].pt() < 10.0 )
		{
		  if(fabs(Electrons[i].superCluster()->eta()) > 0 && fabs(Electrons[i].superCluster()->eta()) < 0.8 ){ cutVal = bdtCutLoPt_0_08;}
		  if(fabs(Electrons[i].superCluster()->eta()) >= 0.8 && fabs(Electrons[i].superCluster()->eta()) < 1.479 ){ cutVal = bdtCutLoPt_08_14;}
		  if(fabs(Electrons[i].superCluster()->eta()) >= 1.479){ cutVal = bdtCutLoPt_14_25;}
		}
	      if( Electrons[i].pt() >= 10.0 )
		{
		  if(fabs(Electrons[i].superCluster()->eta()) > 0 && fabs(Electrons[i].superCluster()->eta()) < 0.8 ){ cutVal = bdtCutHiPt_0_08;}
		  if(fabs(Electrons[i].superCluster()->eta()) >= 0.8 && fabs(Electrons[i].superCluster()->eta()) < 1.479 ){ cutVal = bdtCutHiPt_08_14;}
		  if(fabs(Electrons[i].superCluster()->eta()) >= 1.479){ cutVal = bdtCutHiPt_14_25;}
		}
	      //if(cutVal == 1000){ cout << "Electrons[i].superCluster()->eta() > 2.5, event " << iEvent.id().event() 
	      //<< "pt " << Electrons[i].pt() << " eta " << Electrons[i].eta() << " scEta " << Electrons[i].superCluster()->eta() << endl; continue;}
	      //int misHits = Electrons[i].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniADO
              int misHits = 0; // for miniAOD
	      if(Electrons[i].electronID(elecID) > cutVal && misHits < missingHitsCuts)
		{
		  bestElectrons.push_back(Electrons[i]);
		}
	    }
	}
    }

  
  return bestElectrons;
  
  
}








//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodLooseMuons2012(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( (mu->isGlobalMuon() || mu->isTrackerMuon() || mu->isPFMuon()) && fabs(mu->eta()) < 2.4 && mu->pt() > muPtCut)
	{
	  //cout << "dB: " << mu->userIsolation("PfPUChargedHadronIso") << endl;
	  //double PUCorrVal = 0.5*mu->userIsolation("PfPUChargedHadronIso");
	  //double tmpIso = (mu->chargedHadronIso()+max(mu->photonIso()+mu->neutralHadronIso()-PUCorrVal,0.0))/mu->pt();
	  //cout << tmpIso << endl;

	  bestMuons.push_back(*mu);
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodLooseElectrons2012(edm::Handle<edm::View<pat::Electron> > Electrons, double elPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(elec->eta()) < 2.5 && elec->pt() > elPtCut)
	{
	  bestElectrons.push_back(*elec);
	}
    }
  
  return bestElectrons;
  
  
}



//////////////////////////////////////////////////////////////


//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodLooseMuons(std::vector<pat::Muon> Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  int muTrackerHits = 10;
  /**************************************/
  
  for(unsigned int i = 0; i < Muons.size(); i++)
    {
      if( abs(getSIP3D(Muons[i])) < 100 )
	{
	  if( Muons[i].pt() > muPtCut && abs(Muons[i].eta()) < muEtaCut && Muons[i].isGlobalMuon() == 1 && Muons[i].track()->numberOfValidHits() > muTrackerHits )
	    {
	      bestMuons.push_back(Muons[i]);
	    }
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodLooseElectrons(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  //int missingHitsCuts = 2; // for miniAOD
  /**************************************/
  
  for(unsigned int i = 0; i < Electrons.size(); i++ )
    {
      if( abs(getSIP3D(Electrons[i])) < 100 )
	{
	  if( Electrons[i].pt() > elecPtCut && abs(Electrons[i].eta()) < elecEtaCut)
	    {
	      //if( Electrons[i].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < missingHitsCuts ) // for miniAOD
		//{
		  if (Electrons[i].electronID(elecID) == 1 || Electrons[i].electronID(elecID) == 3 || Electrons[i].electronID(elecID) == 5 || Electrons[i].electronID(elecID) == 7 ||
		      Electrons[i].electronID(elecID) == 9 || Electrons[i].electronID(elecID) == 11 || Electrons[i].electronID(elecID) == 13 || Electrons[i].electronID(elecID) == 15)
		    {
		      bestElectrons.push_back(Electrons[i]);
		    }
		//}
	    }
	}
    }
  
  return bestElectrons;
  
  
}



//Returns vector of good muons for analysis
std::vector<pat::Muon> HZZ4LHelper::goodTightMuons(std::vector<pat::Muon> Muons, double muPtCut)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Muon> bestMuons;
  
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  double muChi2Cut = 10.0;
  int muMuonHits = 0;
  int muNumMatches = 1;
  int muTrackerHits = 10;
  /**************************************/
  
  for(unsigned int i = 0; i < Muons.size(); i++ )
    {
      if( abs(getSIP3D(Muons[i])) < 100 )
	{
	  if( Muons[i].pt() > muPtCut && abs(Muons[i].eta()) < muEtaCut && Muons[i].isGlobalMuon() == 1 && Muons[i].isTrackerMuon() == 1 )
	    {
	      if( Muons[i].normChi2() < muChi2Cut && Muons[i].globalTrack()->hitPattern().numberOfValidMuonHits() > muMuonHits &&  Muons[i].track()->numberOfValidHits() > muTrackerHits)
		{
		  if( Muons[i].numberOfMatches() > muNumMatches )
		    {
		      bestMuons.push_back(Muons[i]);
		    }
		}
	    }
	}
    }
  
  return bestMuons;
  
}


std::vector<pat::Electron> HZZ4LHelper::goodTightElectrons(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector<pat::Electron> bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  //int missingHitsCuts = 2; // for miniAOD
  /**************************************/
  
  for(unsigned int i = 0; i < Electrons.size(); i++ )
    {
      if( abs(getSIP3D(Electrons[i])) < 100 )
	{
	  if( Electrons[i].pt() > elecPtCut && abs(Electrons[i].eta()) < elecEtaCut)
	    {
	      //if( Electrons[i].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < missingHitsCuts )
		//{
		  if(Electrons[i].electronID(elecID) == 5 || Electrons[i].electronID(elecID) == 7 || Electrons[i].electronID(elecID) == 13 || Electrons[i].electronID(elecID) == 15)
		    {
		      bestElectrons.push_back(Electrons[i]);
		    }
		//}
	    }
	}
    }
  
  
  return bestElectrons;
  
  
}


/////////////////////////////////////////

void 
HZZ4LHelper::cleanOverlappingLeptons(std::vector<pat::Muon> &Muons, std::vector<pat::Electron> &Electrons,const reco::Vertex *&vertex)
{
  
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  double dxyCut = 0.2;
  double dzCut = 0.5;
  double tmpDeltR = 100;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      if( (Muons[i].isPFMuon() == 1 || Muons[i].isGlobalMuon() == 1) && fabs(Muons[i].muonBestTrack()->dxy(vertex->position())) < dxyCut  // for miniAOD
	  && fabs(Muons[i].muonBestTrack()->dz(vertex->position())) < dzCut && Muons[i].pt() > 5)
	{
	  for( unsigned int j = 0; j < Electrons.size(); j++ )
	    {
	      tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Electrons[j].eta(),Electrons[j].phi());
	      if( tmpDeltR < 0.05 )
		{
		  Electrons.erase(Electrons.begin()+j);
		}
	    }
	}
    }
      

}


void 
HZZ4LHelper::m2lCutAt4GeV(std::vector<pat::Muon> &Muons, std::vector<pat::Electron> &Electrons)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  
  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = i; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      if( (Muons[i].p4() + Muons[j].p4()).M() < 4.0 )
		{
		  Muons.erase(Muons.begin()+j);
		  Muons.erase(Muons.begin()+i);
		  if(Muons.size() < 2){ break;} //Empty pointer protection
		}
	      
	    }
	}
      
    }

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = i; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      if( (Electrons[i].p4() + Electrons[j].p4()).M() < 4.0 )
		{
		  Electrons.erase(Electrons.begin()+j);
		  Electrons.erase(Electrons.begin()+i);
		  if(Electrons.size() < 2 ){ break;} //Empty pointer protection
                }
	      
	    }
	}
    }


  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( (Muons[i].p4() + Electrons[j].p4()).M() < 4.0 )
	    {
	      Muons.erase(Muons.begin()+i);
	      Electrons.erase(Electrons.begin()+j);
	      if(Electrons.size() < 2 || Muons.size() < 2){break;}
	    }
	  
	}
    }
  

}


void 
HZZ4LHelper::m2lCutAt4GeV(std::vector<pat::Muon> inputMuons, std::vector<pat::Electron> inputElectrons, std::vector<pat::Muon> &outputMuons, std::vector<pat::Electron> &outputElectrons)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  vector<unsigned int> ignoredMuons;
  vector<unsigned int> ignoredElectrons;
  

  for( unsigned int i = 0; i < inputMuons.size(); i++ )
    {
      for( unsigned int j = i; j < inputMuons.size(); j++ )
        {
          if( i != j )
            {
	      if( (inputMuons[i].p4() + inputMuons[j].p4()).M() < 4.0 )
		{
		  ignoredMuons.push_back(i);
		  ignoredMuons.push_back(j);
		}
	    }
	}
      
    }

  for( unsigned int i = 0; i < inputElectrons.size(); i++ )
    {
      for( unsigned int j = i; j < inputElectrons.size(); j++ )
	{
	  if( i != j )
	    {
	      if( (inputElectrons[i].p4() + inputElectrons[j].p4()).M() < 4.0 )
		{
		  ignoredElectrons.push_back(i);
		  ignoredElectrons.push_back(j);
		}
	      
	    }
	}
    }





  for( unsigned int i = 0; i < inputMuons.size(); i++ )
    {
      for( unsigned int j = 0; j < inputElectrons.size(); j++ )
	{
	  if( (inputMuons[i].p4() + inputElectrons[j].p4()).M() < 4.0 )
	    {
	      ignoredMuons.push_back(i);
	      ignoredElectrons.push_back(j);
	    }
	  
	}
    }
  
  bool passedMu = true;
  for( unsigned int i = 0; i < inputMuons.size(); i++)
    {
      passedMu = true;
      for( unsigned int j = 0; j < ignoredMuons.size(); j++)
	{
	  if( i == ignoredMuons[j]){ passedMu = false;}
	}
      if(passedMu){outputMuons.push_back(inputMuons[i]);}
    }


  bool passedEl = true;
  for( unsigned int i = 0; i < inputElectrons.size(); i++)
    {
      passedEl = true;
      for( unsigned int j = 0; j < ignoredElectrons.size(); j++)
	{
	  if( i == ignoredElectrons[j]){ passedEl = false;}
	}
      if(passedEl){outputElectrons.push_back(inputElectrons[i]);}
    }



}



bool 
HZZ4LHelper::passedM2lCut(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double lowM2lCut, double highM2lCut)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  bool passed = false;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = i; j < Muons.size(); j++ )
        {
	  if( (Muons[i].p4() + Muons[j].p4()).M() > lowM2lCut)
	    {
	      passed = true;
	    }
	}
    }
	    
  
  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = i; j < Electrons.size(); j++ )
	{
	  if( (Electrons[i].p4() + Electrons[j].p4()).M() > lowM2lCut)
	    {
	      passed = true;
	    }
	}
    }


  return passed;

}

bool 
HZZ4LHelper::passedM2lCut_OS(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double lowM2lCut, double highM2lCut)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  bool passed = false;
  int counter = 0;
  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = i; j < Muons.size(); j++ )
        {
	  if(Muons[i].charge() * Muons[j].charge() == -1)
	    {
	      if( (Muons[i].p4() + Muons[j].p4()).M() > lowM2lCut)
		{
		  counter++;
		}
	    }
	}
    }
      
  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = i; j < Electrons.size(); j++ )
	{
          if(Electrons[i].charge() * Electrons[j].charge() == -1)
	    {
	      if( (Electrons[i].p4() + Electrons[j].p4()).M() > lowM2lCut)
		{
		  counter++;
		}
	    }
	}
    }


  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
        {
	  if(Muons[i].charge() * Electrons[j].charge() == -1)
	    {
	      if( (Muons[i].p4() + Electrons[j].p4()).M() > lowM2lCut)
		{
		  counter++;
		}
	    }
	}
    }




  if( counter >= 4 )passed = true;

  return passed;

}


bool 
HZZ4LHelper::passedM2lAndPtCut(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double lowM2lCut, double highM2lCut)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  bool passed = false;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = i; j < Muons.size(); j++ )
        {
          if( i != j  && (Muons[i].charge()*Muons[j].charge() == -1))
            {
	      if( (Muons[i].p4() + Muons[j].p4()).M() > lowM2lCut && (Muons[i].p4() + Muons[j].p4()).M() < highM2lCut )
		{
		  if( (Muons[i].pt() > 20 && Muons[j].pt() > 10) || (Muons[i].pt() > 10 && Muons[j].pt() > 20) )
		    {
		      passed = true;
		    }
		  
		}
	    }
	}
    }
  
  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = i; j < Electrons.size(); j++ )
	{
	  if( i != j  && (Electrons[i].charge()*Electrons[j].charge() == -1))
	    {
	      if( (Electrons[i].p4() + Electrons[j].p4()).M() > lowM2lCut && (Electrons[i].p4() + Electrons[j].p4()).M() < highM2lCut  )
		{
		  if( (Electrons[i].pt() > 20 && Electrons[j].pt() > 10) || (Electrons[i].pt() > 10 && Electrons[j].pt() > 20) )
		    {
		      passed = true;
		    }
		}
	    }
	}
    }

  return passed;

}







double 
HZZ4LHelper::minMass2l(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, bool &sameFlavor, bool &sameSign )
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  double min2l = 1000;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = i; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      if( (Muons[i].p4() + Muons[j].p4()).M() < min2l )
		{
		  min2l = (Muons[i].p4() + Muons[j].p4()).M();
		  if(Muons[i].charge()*Muons[j].charge() == 1){sameSign = true;}
		  else{ sameSign = false;}
		  sameFlavor = true;
		}
	    }
	}
    }
	    

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = i; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      if( (Electrons[i].p4() + Electrons[j].p4()).M() < min2l )
		{
		  min2l = (Electrons[i].p4() + Electrons[j].p4()).M();
		  if(Electrons[i].charge()*Electrons[j].charge() == 1){sameSign = true;}
		  else{ sameSign = false;}
		  sameFlavor = true;
                }
	      
	    }
	}
    }


  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      if( (Muons[i].p4() + Electrons[j].p4()).M() < min2l )
		{
                  min2l = (Muons[i].p4() + Electrons[j].p4()).M();
		  if(Muons[i].charge()*Electrons[j].charge() == 1){sameSign = true;}
		  else{ sameSign = false;}
		  sameFlavor = false;

                }
	      
	    }
	}
    }



  return min2l;


}


double 
HZZ4LHelper::maxMass2l(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons )
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  double max2l = 0;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = i; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      if( (Muons[i].p4() + Muons[j].p4()).M() > max2l )
		{
		  max2l = (Muons[i].p4() + Muons[j].p4()).M();
		}
	    }
	}
    }
	    

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = i; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      if( (Electrons[i].p4() + Electrons[j].p4()).M() > max2l )
		{
		  max2l = (Electrons[i].p4() + Electrons[j].p4()).M();
		}
	    }
	}
    }


  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      if( (Muons[i].p4() + Electrons[j].p4()).M() > max2l )
		{
                  max2l = (Muons[i].p4() + Electrons[j].p4()).M();
                }
	    }
	}
    }



  return max2l;


}


double 
HZZ4LHelper::minMass3l(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons )
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  vector<math::XYZTLorentzVector> fourVecs; 
  
  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      fourVecs.push_back(Muons[i].p4());
    }
	    

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      fourVecs.push_back(Electrons[i].p4());
    }

  double tmpMinM3l = 10000;
  double minM3l = 10000;
  for( unsigned int j = 0; j < fourVecs.size(); j++ )
    {
      for( unsigned int k = 0; k < fourVecs.size(); k++ )
	{
	  for( unsigned int l = 0; l < fourVecs.size(); l++ )
	    {
	      if( l != k && l != j && j != k )
		{
		  tmpMinM3l = (fourVecs[j] + fourVecs[k] + fourVecs[l]).M();
		  
		  if( tmpMinM3l < minM3l ){ minM3l = tmpMinM3l; }
		}
	    }
	}
    }


  return minM3l;


}


double 
HZZ4LHelper::largestLepMomentum(TLorentzVector l1, TLorentzVector l2, TLorentzVector l3, TLorentzVector l4)
{
  using namespace std;
  
  //Define the c.o.m frame
  TLorentzVector v4_Z = l1 + l2 + l3 + l4;
  
  // beta parameter of Lorentz rotation: pay attention to the negetive sign!!!
  double zx = -v4_Z.Px()/v4_Z.E(); 
  double zy = -v4_Z.Py()/v4_Z.E(); 
  double zz = -v4_Z.Pz()/v4_Z.E();
  
  TLorentzRotation COM(zx,zy,zz);
  
  //Lorentz tranformation on leptons
  TLorentzVector v4_LepPrime1 = l1.Transform(COM); TLorentzVector v4_LepPrime2 = l2.Transform(COM);
  TLorentzVector v4_LepPrime3 = l3.Transform(COM); TLorentzVector v4_LepPrime4 = l4.Transform(COM);
  
  // find max lep momentum in Z rest frame 
  double p1 =  v4_LepPrime1.P(); 
  double p2 =  v4_LepPrime2.P(); 
  double p3 =  v4_LepPrime3.P(); 
  double p4 =  v4_LepPrime4.P();
  
  double P_com[4] = {p1,p2,p3,p4};
  
  double p_max = 0.0;  //int max_lep = -1;
  for (int p = 0; p<4; p++)
    {
      if(P_com[p]>=p_max) { p_max = P_com[p]; /*max_lep = p;*/ }
    }
  

  return p_max;
  
}

double
HZZ4LHelper::M3lSoftestLep(TLorentzVector l1, TLorentzVector l2, TLorentzVector l3, TLorentzVector l4)
{
  
  using namespace std;
  
  //Define the c.o.m frame
  TLorentzVector v4_Z = l1 + l2 + l3 + l4;
  
  // beta parameter of Lorentz rotation: pay attention to the negetive sign!!!
  double zx = -v4_Z.Px()/v4_Z.E(); 
  double zy = -v4_Z.Py()/v4_Z.E(); 
  double zz = -v4_Z.Pz()/v4_Z.E();
  
  TLorentzRotation COM(zx,zy,zz);
  
  //Lorentz tranformation on leptons
  TLorentzVector v4_LepPrime1 = l1.Transform(COM); TLorentzVector v4_LepPrime2 = l2.Transform(COM);
  TLorentzVector v4_LepPrime3 = l3.Transform(COM); TLorentzVector v4_LepPrime4 = l4.Transform(COM);
  
  // find max lep momentum in Z rest frame 
  double p1 =  v4_LepPrime1.P(); 
  double p2 =  v4_LepPrime2.P(); 
  double p3 =  v4_LepPrime3.P(); 
  double p4 =  v4_LepPrime4.P();
  
  double P_com[4] = {p1,p2,p3,p4};

  TLorentzVector p_trans[4] = {v4_LepPrime1,v4_LepPrime2,v4_LepPrime3,v4_LepPrime4};
  
  double p_max = 0.0;  int max_lep = -1;
  for (int p = 0; p<4; p++)
    {
      if(P_com[p]>=p_max) { p_max = P_com[p]; max_lep = p; }
    }

  TLorentzVector vecc;
  for( int q = 0; q < 4; q++ )
    {
      if( max_lep != q )
	{
	  vecc += p_trans[q];
	}
    }
     
  double softest = vecc.M();
  
  return softest;

}



double
HZZ4LHelper::getSumOf2LeastIso(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  double tmpDeltR = 100;
  double tmpIso = 1000;
  vector<double> iso_;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = i; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Muons[j].eta(),Muons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpIso = ( (Muons[i].trackIso() - Muons[j].pt()) + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0) )/Muons[i].pt();
		}
	      else{
		tmpIso = (Muons[i].trackIso() + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0))/Muons[i].pt();
	      }
	    }
	}
      iso_.push_back(tmpIso);
    }
	    

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = i; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      tmpDeltR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Electrons[j].eta(),Electrons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpIso =  ( (Electrons[i].dr03TkSumPt()-Electrons[j].pt()) + max( (Electrons[i].dr03EcalRecHitSumEt()-Electrons[j].et())
			      + Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt(); 
		}
	      else{
		tmpIso = (Electrons[i].dr03TkSumPt() + max(Electrons[i].dr03EcalRecHitSumEt() + Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt();
	      }
	    }
	}
      iso_.push_back(tmpIso);
    }


  double leastIso1 = 0;
  double leastIso2 = 0;

  unsigned int taken = 100;
  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( iso_[i] > leastIso1 ){ leastIso1 = iso_[i]; taken = i; }
    }

  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( i != taken )
	{
	  if( iso_[i] > leastIso2 ){ leastIso2 = iso_[i]; }
	}
    }

  double leastIso = leastIso1 + leastIso2;

  return leastIso;

}



double 
HZZ4LHelper::getLeastIso1(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  double tmpIso = 1000;
  vector<double> iso_;

  double leastIso1 = 0;
  double leastIso2 = 0;


  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      tmpIso = ( Muons[i].trackIso()  + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0) )/Muons[i].pt();
      iso_.push_back(tmpIso);
    }
	    

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      tmpIso =  ( Electrons[i].dr03TkSumPt() + max( Electrons[i].dr03EcalRecHitSumEt() + Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt(); 
      iso_.push_back(tmpIso);
    }

  unsigned int taken = 100;
  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( iso_[i] > leastIso1 ){ leastIso1 = iso_[i]; taken = i; }
    }
  
  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( i != taken )
	{
	  if( iso_[i] > leastIso2 ){ leastIso2 = iso_[i]; }
	}
    }
  
  
  return leastIso1;

}

double 
HZZ4LHelper::getLeastIso1(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr, double &pT)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  double tmpIso = 1000;
  vector<double> iso_;
  
  double leastIso1 = 0;
  

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      tmpIso = ( Muons[i].trackIso()  + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0) )/Muons[i].pt();
      if( tmpIso > leastIso1 ){ leastIso1 = tmpIso; pT = Muons[i].pt(); }
     
    }
	    

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      tmpIso =  ( Electrons[i].dr03TkSumPt() + max( Electrons[i].dr03EcalRecHitSumEt() + Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt(); 
      if( tmpIso > leastIso1 ){leastIso1 = tmpIso; pT = Electrons[i].pt();}
      
    }
  
  return leastIso1;

}


double 
HZZ4LHelper::getLeastIso2(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  double tmpIso = 1000;
  vector<double> iso_;

  double leastIso1 = 0;
  double leastIso2 = 0;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      tmpIso = ( Muons[i].trackIso()  + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0) )/Muons[i].pt();
      iso_.push_back(tmpIso);
    }
	    

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      tmpIso =  ( Electrons[i].dr03TkSumPt() + max( Electrons[i].dr03EcalRecHitSumEt() + Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt(); 
      iso_.push_back(tmpIso);
    }

  unsigned int taken = 100;
  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( iso_[i] > leastIso1 ){ leastIso1 = iso_[i]; taken = i; }
    }
  
  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( i != taken )
	{
	  if( iso_[i] > leastIso2 ){ leastIso2 = iso_[i]; }
	}
    }
  
  
  return leastIso2;

}

int 
HZZ4LHelper::passedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut, double muonRhoCorr, double elecRhoCorr )
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  int tmpIsoCounter = 0;
  double tmpDeltR = 100;
  double tmpIso = 1000;

  double tmpEtCorr = 0;
  double tmpPtCorr = 0;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Muons[j].eta(),Muons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Muons[j].pt();
		  //tmpEtCorr += Muons[j].et();
		}
	    }
	}

      for( unsigned int k = 0; k < Electrons.size(); k++ )
	{
	  tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(), Electrons[k].eta(),Electrons[k].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Electrons[k].pt();
	      tmpEtCorr += Electrons[k].et();
	    }
	}
    

      tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) + Muons[i].hcalIso() - muonRhoCorr,0.0) )/Muons[i].pt();

      if( tmpIso < chosenIsoCut ){ tmpIsoCounter++;}

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }
	    

  tmpPtCorr = 0;
  tmpEtCorr = 0;

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      tmpDeltR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Electrons[j].eta(),Electrons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Electrons[j].pt();
		  tmpEtCorr += Electrons[j].et();

		}
	    }
	}

      for( unsigned int k = 0; k < Muons.size(); k++ )
        {
	  tmpDeltR = deltaR(Muons[k].eta(),Muons[k].phi(),Electrons[i].eta(), Electrons[i].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Muons[k].pt();
	      //tmpEtCorr += Muons[k].et();
	    }
	}
    

      tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								+ Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt(); 

      if( tmpIso < chosenIsoCut ){ tmpIsoCounter++;}

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }


  return tmpIsoCounter;

}


int 
HZZ4LHelper::passedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut, TString isoType, TString muPUCorr, TString elPUCorr, double Rho)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  int tmpIsoCounter = 0;
  double tmpDeltR = 100;
  double tmpIso = 1000;

  double tmpEtCorr = 0;
  double tmpPtCorr = 0;

  double PUCorrVal = 0;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Muons[j].eta(),Muons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Muons[j].pt();
		  /////tmpEtCorr += Muons[j].et();
		}
	    }
	}

      for( unsigned int k = 0; k < Electrons.size(); k++ )
	{
	  tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(), Electrons[k].eta(),Electrons[k].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Electrons[k].pt();
	      tmpEtCorr += Electrons[k].et();
	    }
	}
      

      if(muPUCorr == "dB" && isoType == "PF"){
	PUCorrVal = 0.5*Muons[i].userIsolation("PfPUChargedHadronIso");
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "PF" ){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "DET"){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) + 
							    Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "PF" ){ 
	PUCorrVal = Rho*MuonEffArea(muEAtype,Muons[i].eta(),muEAtarget);
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "DET"){
	double effA = 0;
        if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
        if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
        PUCorrVal = Rho*effA;
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) +
							    Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "PFEffAreaRho" && isoType == "PF"){
	PUCorrVal = Rho*MuonEffArea(muEAtype,Muons[i].eta(),muEAtarget);
        tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRhoBL" && isoType == "DET"){
	double effA_Ecal = 0, effA_Hcal = 0;
	if( abs(Muons[i].eta()) < 1.479){effA_Ecal = 0.074; effA_Hcal = 0.022;}
        if( abs(Muons[i].eta()) > 1.479){effA_Ecal = 0.045; effA_Hcal = 0.030;}
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr-effA_Ecal*Rho ) +
							    Muons[i].hcalIso()-effA_Hcal*Rho,0.0) )/Muons[i].pt();
      }
      else{
	string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	err+="Possibilites are:\n";
	err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	throw invalid_argument( err );
      }
      
      if( tmpIso < chosenIsoCut ){ tmpIsoCounter++;}
      
      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }
	    

  tmpPtCorr = 0;
  tmpEtCorr = 0;

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      tmpDeltR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Electrons[j].eta(),Electrons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Electrons[j].pt();
		  tmpEtCorr += Electrons[j].et();

		}
	    }
	}

      for( unsigned int k = 0; k < Muons.size(); k++ )
        {
	  tmpDeltR = deltaR(Muons[k].eta(),Muons[k].phi(),Electrons[i].eta(), Electrons[i].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Muons[k].pt();
	      //tmpEtCorr += Muons[k].et();
	    }
	}
    

      if(elPUCorr == "dB" && isoType == "PF"){
	PUCorrVal = 0.5*Electrons[i].userIsolation("PfPUChargedHadronIso");
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "rho" && isoType == "PF"){ 
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "rho" && isoType == "DET"){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								  + Electrons[i].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRho" && isoType == "PF"){ 
	double effA = 0;
	if( Electrons[i].isEB()){effA = 0.085;}//are there numbers for electrons?
	if( Electrons[i].isEE()){effA = 0.060;}
	PUCorrVal = Rho*effA;
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRho" && isoType == "DET"){ 
	double effA = 0;
	if( Electrons[i].isEB()){effA = 0.085;}//are there numbers for electrons?
	if( Electrons[i].isEE()){effA = 0.060;}
	PUCorrVal = Rho*effA;
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								  + Electrons[i].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "PFEffAreaRho" && isoType == "PF"){
	PUCorrVal = Rho*ElecEffArea(elEAtype,Electrons[i].eta(),elEAtarget); 
        tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRhoBL" && isoType == "DET"){double effA_Ecal = 0, effA_Hcal = 0;
	if( Electrons[i].isEB()){effA_Ecal = 0.101; effA_Hcal = 0.021;}
        if( Electrons[i].isEE()){effA_Ecal = 0.046; effA_Hcal = 0.040;}
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr-effA_Ecal*Rho)
								  + Electrons[i].dr03HcalTowerSumEt()-effA_Hcal*Rho, 0.0))/Electrons[i].pt();
      }
      else{
	string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	err+="Possibilites are:\n";
	err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	throw invalid_argument( err );
      }
      
      
      if( tmpIso < chosenIsoCut ){ tmpIsoCounter++;}

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }


  return tmpIsoCounter;

}


int 
HZZ4LHelper::passedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut, TString isoType, TString muPUCorr, TString elPUCorr, double muonRho,double elecRho,double &leastIso1, double &leastIso2 )
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  int tmpIsoCounter = 0;
  double tmpDeltR = 100;
  double tmpIso = 1000;
  vector<double> iso_;

  leastIso1 = 0;
  leastIso2 = 0;

  double tmpEtCorr = 0;
  double tmpPtCorr = 0;

  double PUCorrVal = 0;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Muons[j].eta(),Muons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Muons[j].pt();
		  
		}
	    }
	}

      for( unsigned int k = 0; k < Electrons.size(); k++ )
	{
	  tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(), Electrons[k].eta(),Electrons[k].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Electrons[k].pt();
	      tmpEtCorr += Electrons[k].et();
	    }
	}
      
      if(muPUCorr == "dB" && isoType == "PF"){
	PUCorrVal = 0.5*Muons[i].userIsolation("PfPUChargedHadronIso");
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "PF" ){
	PUCorrVal = muonRho*PI*0.3*0.3;
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "DET"){
	PUCorrVal = muonRho*PI*0.3*0.3;
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) + 
							    Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "PF" ){ 
	double effA = 0;
	if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
	if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
	PUCorrVal = muonRho*effA;
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "DET"){
	double effA = 0;
        if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
        if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
        PUCorrVal = muonRho*effA;
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) +
							    Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "PFEffAreaRho" && isoType == "PF"){
	PUCorrVal = muonRho*MuonEffArea(muEAtype,Muons[i].eta(),muEAtarget);
        tmpIso = (Muons[i].chargedHadronIso()+ max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRhoBL" && isoType == "DET"){
	double effA_Ecal = 0, effA_Hcal = 0;
	if( abs(Muons[i].eta()) < 1.479){effA_Ecal = 0.074; effA_Hcal = 0.022;}
        if( abs(Muons[i].eta()) > 1.479){effA_Ecal = 0.045; effA_Hcal = 0.030;}
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr-effA_Ecal*muonRho ) +
							    Muons[i].hcalIso()-effA_Hcal*muonRho,0.0) )/Muons[i].pt();
      }
      else{
	string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	err+="Possibilites are:\n";
	err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	throw invalid_argument( err );
      }
      
      if( tmpIso < chosenIsoCut ){ tmpIsoCounter++;}
      iso_.push_back(tmpIso);

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }
	    

  tmpPtCorr = 0;
  tmpEtCorr = 0;

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      tmpDeltR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Electrons[j].eta(),Electrons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Electrons[j].pt();
		  tmpEtCorr += Electrons[j].et();

		}
	    }
	}

      for( unsigned int k = 0; k < Muons.size(); k++ )
        {
	  tmpDeltR = deltaR(Muons[k].eta(),Muons[k].phi(),Electrons[i].eta(), Electrons[i].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Muons[k].pt();
	      
	    }
	}
    

      if(elPUCorr == "dB" && isoType == "PF"){
	PUCorrVal = 0.5*Electrons[i].userIsolation("PfPUChargedHadronIso");
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "rho" && isoType == "PF"){ 
	PUCorrVal = elecRho*PI*0.3*0.3;
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "rho" && isoType == "DET"){
	PUCorrVal = elecRho*PI*0.3*0.3;
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								  + Electrons[i].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRho" && isoType == "PF"){ 
	double effA = 0;
	if( Electrons[i].isEB()){effA = 0.085;}//are there numbers for electrons?
	if( Electrons[i].isEE()){effA = 0.060;}
	PUCorrVal = elecRho*effA;
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRho" && isoType == "DET"){ 
	double effA = 0;
	if( Electrons[i].isEB()){effA = 0.085;}//are there numbers for electrons?
	if( Electrons[i].isEE()){effA = 0.060;}
	PUCorrVal = elecRho*effA;
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								  + Electrons[i].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "PFEffAreaRho" && isoType == "PF"){
	PUCorrVal = elecRho*ElecEffArea(elEAtype,Electrons[i].eta(),elEAtarget);
        tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRhoBL" && isoType == "DET"){double effA_Ecal = 0, effA_Hcal = 0;
	if( Electrons[i].isEB()){effA_Ecal = 0.101; effA_Hcal = 0.021;}
        if( Electrons[i].isEE()){effA_Ecal = 0.046; effA_Hcal = 0.040;}
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr-effA_Ecal*elecRho)
								  + Electrons[i].dr03HcalTowerSumEt()-effA_Hcal*elecRho, 0.0))/Electrons[i].pt();
      }
      else{
	string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	err+="Possibilites are:\n";
	err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	throw invalid_argument( err );
      }
      
      
      if( tmpIso < chosenIsoCut ){ tmpIsoCounter++;}
      iso_.push_back(tmpIso);
      
      tmpPtCorr = 0;
      tmpEtCorr = 0;
      
    }
  
  
  
  unsigned int taken = 100;
  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( iso_[i] > leastIso1 ){ leastIso1 = iso_[i]; taken = i; }
    }

  for( unsigned int i = 0; i < iso_.size(); i++)
    {
      if( i != taken )
        {
          if( iso_[i] > leastIso2 ){ leastIso2 = iso_[i]; }
	}
    }
    
  
  return tmpIsoCounter;

}



int 
HZZ4LHelper::passedEarlyIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr, double earlyIsoCut)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  int tmpIsoCounter = 0;
  double tmpIso = 1000;
 
  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      tmpIso = (Muons[i].trackIso() + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0))/Muons[i].pt();
      if( tmpIso < earlyIsoCut ){ tmpIsoCounter++;}
    }

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      tmpIso = (Electrons[i].dr03TkSumPt() + max(Electrons[i].dr03EcalRecHitSumEt() + Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt();
      if( tmpIso < earlyIsoCut ){ tmpIsoCounter++;}
    }

  return tmpIsoCounter;

}



int 
HZZ4LHelper::passedInvertedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double muonRhoCorr, double elecRhoCorr, double chosenIsoCut)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  int tmpIsoCounter = 0;
  double tmpDeltR = 100;
  double tmpIso = 1000;
  vector<double> iso_;

  double tmpEtCorr = 0;
  double tmpPtCorr = 0;

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Muons[j].eta(),Muons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Muons[j].pt();
		  //tmpEtCorr += Muons[j].et();
		}
	    }
	}

      for( unsigned int k = 0; k < Electrons.size(); k++ )
	{
	  tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(), Electrons[k].eta(),Electrons[k].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Electrons[k].pt();
	      tmpEtCorr += Electrons[k].et();
	    }
	}
    

      tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) + Muons[i].hcalIso() - muonRhoCorr,0.0) )/Muons[i].pt();
      if( tmpIso > chosenIsoCut ){ tmpIsoCounter++;}
      iso_.push_back(tmpIso);

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }
	    

  tmpPtCorr = 0;
  tmpEtCorr = 0;

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      tmpDeltR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Electrons[j].eta(),Electrons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Electrons[j].pt();
		  tmpEtCorr += Electrons[j].et();

		}
	    }
	}

      for( unsigned int k = 0; k < Muons.size(); k++ )
        {
	  tmpDeltR = deltaR(Muons[k].eta(),Muons[k].phi(),Electrons[i].eta(), Electrons[i].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Muons[k].pt();
	      //tmpEtCorr += Muons[k].et();
	    }
	}
    

      tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								+ Electrons[i].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[i].pt(); 
      if( tmpIso > chosenIsoCut ){ tmpIsoCounter++;}
      iso_.push_back(tmpIso);

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }

  return tmpIsoCounter;

}


int 
HZZ4LHelper::passedInvertedIsolation(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double chosenIsoCut,TString isoType, TString muPUCorr, TString elPUCorr, double Rho)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  int tmpIsoCounter = 0;
  double tmpDeltR = 100;
  double tmpIso = 1000;
  vector<double> iso_;

  double tmpEtCorr = 0;
  double tmpPtCorr = 0;

  double PUCorrVal = 0;


  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      for( unsigned int j = 0; j < Muons.size(); j++ )
        {
          if( i != j )
            {
	      tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Muons[j].eta(),Muons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Muons[j].pt();
		  //tmpEtCorr += Muons[j].et();
		}
	    }
	}

      for( unsigned int k = 0; k < Electrons.size(); k++ )
	{
	  tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(), Electrons[k].eta(),Electrons[k].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Electrons[k].pt();
	      tmpEtCorr += Electrons[k].et();
	    }
	}
    
      if(muPUCorr == "dB" && isoType == "PF"){
	PUCorrVal = 0.5*Muons[i].userIsolation("PfPUChargedHadronIso");
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "PF" ){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "DET"){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) + 
							    Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "PF" ){ 
	double effA = 0;
	if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
	if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
	PUCorrVal = Rho*effA;
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "DET"){
	double effA = 0;
        if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
        if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
        PUCorrVal = Rho*effA;
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr) +
							    Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "PFEffAreaRho" && isoType == "PF"){
	PUCorrVal = Rho*MuonEffArea(muEAtype,Muons[i].eta(),muEAtarget);
        tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRhoBL" && isoType == "DET"){
	double effA_Ecal = 0, effA_Hcal = 0;
	if( abs(Muons[i].eta()) < 1.479){effA_Ecal = 0.074; effA_Hcal = 0.022;}
        if( abs(Muons[i].eta()) > 1.479){effA_Ecal = 0.045; effA_Hcal = 0.030;}
	tmpIso = ( (Muons[i].trackIso() - tmpPtCorr) + max( (Muons[i].ecalIso()-tmpEtCorr-effA_Ecal*Rho ) +
							    Muons[i].hcalIso()-effA_Hcal*Rho,0.0) )/Muons[i].pt();
      }
      else{
	string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	err+="Possibilites are:\n";
	err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	throw invalid_argument( err );
      }

      if( tmpIso > chosenIsoCut ){ tmpIsoCounter++;}
      iso_.push_back(tmpIso);

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }
	    

  tmpPtCorr = 0;
  tmpEtCorr = 0;

  for( unsigned int i = 0; i < Electrons.size(); i++ )
    {
      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  if( i != j )
	    {
	      tmpDeltR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Electrons[j].eta(),Electrons[j].phi());
	      if( tmpDeltR < 0.3 )
		{
		  tmpPtCorr += Electrons[j].pt();
		  tmpEtCorr += Electrons[j].et();

		}
	    }
	}

      for( unsigned int k = 0; k < Muons.size(); k++ )
        {
	  tmpDeltR = deltaR(Muons[k].eta(),Muons[k].phi(),Electrons[i].eta(), Electrons[i].phi());
	  if( tmpDeltR < 0.3 )
	    {
	      tmpPtCorr += Muons[k].pt();
	      //tmpEtCorr += Muons[k].et();
	    }
	}
    
      if(elPUCorr == "dB" && isoType == "PF"){
	PUCorrVal = 0.5*Electrons[i].userIsolation("PfPUChargedHadronIso");
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "rho" && isoType == "PF"){ 
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "rho" && isoType == "DET"){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								  + Electrons[i].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRho" && isoType == "PF"){ 
	double effA = 0;
	if( Electrons[i].isEB()){effA = 0.085;}//are there numbers for electrons?
	if( Electrons[i].isEE()){effA = 0.060;}
	PUCorrVal = Rho*effA;
	tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRho" && isoType == "DET"){ 
	double effA = 0;
	if( Electrons[i].isEB()){effA = 0.085;}//are there numbers for electrons?
	if( Electrons[i].isEE()){effA = 0.060;}
	PUCorrVal = Rho*effA;
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr)
								  + Electrons[i].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "PFEffAreaRho" && isoType == "PF"){
	PUCorrVal = Rho*ElecEffArea(elEAtype,Electrons[i].eta(),elEAtarget);
        tmpIso = (Electrons[i].chargedHadronIso()+max(Electrons[i].photonIso()+Electrons[i].neutralHadronIso()-PUCorrVal,0.0))/Electrons[i].pt();
      }
      else if(elPUCorr == "EffAreaRhoBL" && isoType == "DET"){double effA_Ecal = 0, effA_Hcal = 0;
	if( Electrons[i].isEB()){effA_Ecal = 0.101; effA_Hcal = 0.021;}
        if( Electrons[i].isEE()){effA_Ecal = 0.046; effA_Hcal = 0.040;}
	tmpIso =  ( (Electrons[i].dr03TkSumPt()-tmpPtCorr) + max( (Electrons[i].dr03EcalRecHitSumEt()-tmpEtCorr-effA_Ecal*Rho)
								  + Electrons[i].dr03HcalTowerSumEt()-effA_Hcal*Rho, 0.0))/Electrons[i].pt();
      }
      else{
	string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	err+="Possibilites are:\n";
	err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	throw invalid_argument( err );
      }

      if( tmpIso > chosenIsoCut ){ tmpIsoCounter++;}
      iso_.push_back(tmpIso);

      tmpPtCorr = 0;
      tmpEtCorr = 0;

    }

  return tmpIsoCounter;

}






bool 
HZZ4LHelper::passedPtandEarlyIso(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double leadingPtCut, double subleadingPtCut,
				      double earlyIsoCut, double muonRhoCorr, double elecRhoCorr )
{

  using namespace std;

  bool passed = false;
  int passedLeading = 0;
  int passedSubLeading = 0;
  
  unsigned int pickedMu = 999;
  unsigned int pickedE =  999;


  ///////////////////////Leading//////////////////////

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      if( Muons[i].pt() > leadingPtCut )
	{
	  if( (Muons[i].trackIso() + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0))/Muons[i].pt() < earlyIsoCut ) 
	    {
	      passedLeading++;
	      pickedMu = i;
	      break;
	    }
	}
    }

  for( unsigned int j = 0; j < Electrons.size(); j++ )
    {
      if( Electrons[j].pt() > leadingPtCut )
	{
	  if( (Electrons[j].dr03TkSumPt() + max(Electrons[j].dr03EcalRecHitSumEt() + Electrons[j].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[j].pt() < earlyIsoCut )
	    {
	      passedLeading++;
	      pickedE = j;
	      break;
	    }
	}
    }

  ///////////////////////SubLeading//////////////////////

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      if( Muons[i].pt() > subleadingPtCut )
	{
	  if( (Muons[i].trackIso() + max(Muons[i].ecalIso() + Muons[i].hcalIso() - muonRhoCorr,0.0))/Muons[i].pt() < earlyIsoCut ) 
	    {
	      if( i != pickedMu )
		{
		  passedSubLeading++;
		}
	    }
	}
    }

  for( unsigned int j = 0; j < Electrons.size(); j++ )
    {
      if( Electrons[j].pt() > subleadingPtCut )
	{
	  if( (Electrons[j].dr03TkSumPt() + max(Electrons[j].dr03EcalRecHitSumEt() + Electrons[j].dr03HcalTowerSumEt() - elecRhoCorr, 0.0))/Electrons[j].pt() < earlyIsoCut )
	    {
	      if( j != pickedE )
		{
		  passedSubLeading++;
		}
	    }
	}
    }

  if( (passedLeading > 0 && passedSubLeading > 0)  or (pickedE<999 and pickedMu<999)){ passed = true; }

  return passed;

}

bool 
HZZ4LHelper::passedPtandEarlyIso(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons, double leadingPtCut, double subleadingPtCut, double earlyIsoCut, TString isoType, TString muPUCorr, TString elPUCorr, double Rho )
{

  using namespace std;

  bool passed = false;
  int passedLeading = 0;
  int passedSubLeading = 0;
  
  unsigned int pickedMu = 999;
  unsigned int pickedE =  999;

  double tmpIso = 1000;
  double PUCorrVal = 0;

  ///////////////////////Leading//////////////////////


  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      if(muPUCorr == "dB" && isoType == "PF"){
	PUCorrVal = 0.5*Muons[i].userIsolation("PfPUChargedHadronIso");
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "PF" ){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "rho" && isoType == "DET"){
	PUCorrVal = Rho*PI*0.3*0.3;
	tmpIso = ( Muons[i].trackIso() + max( Muons[i].ecalIso() + Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "PF" ){ 
	double effA = 0;
	if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
	if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
	PUCorrVal = Rho*effA;
	tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRho" && isoType == "DET"){
	double effA = 0;
        if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
        if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
        PUCorrVal = Rho*effA;
	tmpIso = ( Muons[i].trackIso()  + max( Muons[i].ecalIso() + Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
      }
      else if(muPUCorr == "PFEffAreaRho" && isoType == "PF"){
	PUCorrVal = Rho*MuonEffArea(muEAtype,Muons[i].eta(),muEAtarget);
        tmpIso = (Muons[i].chargedHadronIso()+ max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
      }
      else if(muPUCorr == "EffAreaRhoBL" && isoType == "DET"){
	double effA_Ecal = 0, effA_Hcal = 0;
	if( abs(Muons[i].eta()) < 1.479){effA_Ecal = 0.074; effA_Hcal = 0.022;}
        if( abs(Muons[i].eta()) > 1.479){effA_Ecal = 0.045; effA_Hcal = 0.030;}
	tmpIso = ( Muons[i].trackIso()  + max( (Muons[i].ecalIso()-effA_Ecal*Rho ) + Muons[i].hcalIso()-effA_Hcal*Rho,0.0) )/Muons[i].pt();
      }
      else{
	string err = "HZZ4LHelper::passedPtandEarlyIso --- unknown isolation type\n";
	err+="Possibilites are:\n";
	err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	throw invalid_argument( err );
      }

      if( Muons[i].pt() > leadingPtCut )
	{
	  if( tmpIso < earlyIsoCut ) 
	    {
	      passedLeading++;
	      pickedMu = i;
	      break;
	    }
	      }
	}

      for( unsigned int j = 0; j < Electrons.size(); j++ )
	{
	  
	  if(elPUCorr == "dB" && isoType == "PF"){
	    PUCorrVal = 0.5*Electrons[j].userIsolation("PfPUChargedHadronIso");
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "rho" && isoType == "PF"){ 
	    PUCorrVal = Rho*PI*0.3*0.3;
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "rho" && isoType == "DET"){
	    PUCorrVal = Rho*PI*0.3*0.3;
	    tmpIso =  ( Electrons[j].dr03TkSumPt() + max( Electrons[j].dr03EcalRecHitSumEt()
								      + Electrons[j].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "EffAreaRho" && isoType == "PF"){ 
	    double effA = 0;
	    if( Electrons[j].isEB()){effA = 0.085;}//are there numbers for electrons?
	    if( Electrons[j].isEE()){effA = 0.060;}
	    PUCorrVal = Rho*effA;
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "EffAreaRho" && isoType == "DET"){ 
	    double effA = 0;
	    if( Electrons[j].isEB()){effA = 0.085;}//are there numbers for electrons?
	    if( Electrons[j].isEE()){effA = 0.060;}
	    PUCorrVal = Rho*effA;
	    tmpIso =  ( Electrons[j].dr03TkSumPt() + max( Electrons[j].dr03EcalRecHitSumEt() + Electrons[j].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "PFEffAreaRho" && isoType == "PF"){
	    PUCorrVal = Rho*ElecEffArea(elEAtype,Electrons[j].eta(),elEAtarget);
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "EffAreaRhoBL" && isoType == "DET"){double effA_Ecal = 0, effA_Hcal = 0;
	    if( Electrons[j].isEB()){effA_Ecal = 0.101; effA_Hcal = 0.021;}
	    if( Electrons[j].isEE()){effA_Ecal = 0.046; effA_Hcal = 0.040;}
	    tmpIso =  ( Electrons[j].dr03TkSumPt() + max( (Electrons[j].dr03EcalRecHitSumEt()-effA_Ecal*Rho)
							  + Electrons[j].dr03HcalTowerSumEt()-effA_Hcal*Rho, 0.0))/Electrons[j].pt();
	  }
	  else{
	    string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	    err+="Possibilites are:\n";
	    err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	    err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	    throw invalid_argument( err );
	  }
      
	  if( Electrons[j].pt() > leadingPtCut )
	    {
	      if( tmpIso < earlyIsoCut )
		{
		  passedLeading++;
		  pickedE = j;
		  break;
		}
	    }
	}

  ///////////////////////SubLeading//////////////////////

  for( unsigned int i = 0; i < Muons.size(); i++ )
    {
      if( Muons[i].pt() > subleadingPtCut )
	{
	  if(muPUCorr == "dB" && isoType == "PF"){
	    PUCorrVal = 0.5*Muons[i].userIsolation("PfPUChargedHadronIso");
	    tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
	  }
	  else if(muPUCorr == "rho" && isoType == "PF" ){
	    PUCorrVal = Rho*PI*0.3*0.3;
	    tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
	  }
	  else if(muPUCorr == "rho" && isoType == "DET"){
	    PUCorrVal = Rho*PI*0.3*0.3;
	    tmpIso = ( Muons[i].trackIso() + max( Muons[i].ecalIso() + Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
	  }
	  else if(muPUCorr == "EffAreaRho" && isoType == "PF" ){ 
	    double effA = 0;
	    if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
	    if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
	    PUCorrVal = Rho*effA;
	    tmpIso = (Muons[i].chargedHadronIso()+max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
	  }
	  else if(muPUCorr == "EffAreaRho" && isoType == "DET"){
	    double effA = 0;
	    if( abs(Muons[i].eta()) < 1.479){effA = 0.085;}
	    if( abs(Muons[i].eta()) > 1.479){effA = 0.060;}
	    PUCorrVal = Rho*effA;
	    tmpIso = ( Muons[i].trackIso()  + max( Muons[i].ecalIso() + Muons[i].hcalIso() - PUCorrVal,0.0) )/Muons[i].pt();
	  }
	  else if(muPUCorr == "PFEffAreaRho" && isoType == "PF"){
	    PUCorrVal = Rho*MuonEffArea(muEAtype,Muons[i].eta(),muEAtarget);
	    tmpIso = (Muons[i].chargedHadronIso()+ max(Muons[i].photonIso()+Muons[i].neutralHadronIso()-PUCorrVal,0.0))/Muons[i].pt();
	  }
	  else if(muPUCorr == "EffAreaRhoBL" && isoType == "DET"){
	    double effA_Ecal = 0, effA_Hcal = 0;
	    if( abs(Muons[i].eta()) < 1.479){effA_Ecal = 0.074; effA_Hcal = 0.022;}
	    if( abs(Muons[i].eta()) > 1.479){effA_Ecal = 0.045; effA_Hcal = 0.030;}
	    tmpIso = ( Muons[i].trackIso()  + max( (Muons[i].ecalIso()-effA_Ecal*Rho ) + Muons[i].hcalIso()-effA_Hcal*Rho,0.0) )/Muons[i].pt();
	  }
	  else{
	    string err = "HZZ4LHelper::passedPtandEarlyIso --- unknown isolation type\n";
	    err+="Possibilites are:\n";
	    err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	    err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	    throw invalid_argument( err );
	  }

      
	  if( tmpIso < earlyIsoCut ) 
	    {
	      if( i != pickedMu )
		{
		  passedSubLeading++;
		}
	    }
	}
    }
  
  for( unsigned int j = 0; j < Electrons.size(); j++ )
    {
	  if(elPUCorr == "dB" && isoType == "PF"){
	    PUCorrVal = 0.5*Electrons[j].userIsolation("PfPUChargedHadronIso");
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "rho" && isoType == "PF"){ 
	    PUCorrVal = Rho*PI*0.3*0.3;
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "rho" && isoType == "DET"){
	    PUCorrVal = Rho*PI*0.3*0.3;
	    tmpIso =  ( Electrons[j].dr03TkSumPt() + max( Electrons[j].dr03EcalRecHitSumEt()
							  + Electrons[j].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "EffAreaRho" && isoType == "PF"){ 
	    double effA = 0;
	    if( Electrons[j].isEB()){effA = 0.085;}//are there numbers for electrons?
	    if( Electrons[j].isEE()){effA = 0.060;}
	    PUCorrVal = Rho*effA;
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "EffAreaRho" && isoType == "DET"){ 
	    double effA = 0;
	    if( Electrons[j].isEB()){effA = 0.085;}//are there numbers for electrons?
	    if( Electrons[j].isEE()){effA = 0.060;}
	    PUCorrVal = Rho*effA;
	    tmpIso =  ( Electrons[j].dr03TkSumPt() + max( Electrons[j].dr03EcalRecHitSumEt() + Electrons[j].dr03HcalTowerSumEt() - PUCorrVal, 0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "PFEffAreaRho" && isoType == "PF"){
	    PUCorrVal = Rho*ElecEffArea(elEAtype,Electrons[j].eta(),elEAtarget);
	    tmpIso = (Electrons[j].chargedHadronIso()+max(Electrons[j].photonIso()+Electrons[j].neutralHadronIso()-PUCorrVal,0.0))/Electrons[j].pt();
	  }
	  else if(elPUCorr == "EffAreaRhoBL" && isoType == "DET"){double effA_Ecal = 0, effA_Hcal = 0;
	    if( Electrons[j].isEB()){effA_Ecal = 0.101; effA_Hcal = 0.021;}
	    if( Electrons[j].isEE()){effA_Ecal = 0.046; effA_Hcal = 0.040;}
	    tmpIso =  ( Electrons[j].dr03TkSumPt() + max( (Electrons[j].dr03EcalRecHitSumEt()-effA_Ecal*Rho)
							  + Electrons[j].dr03HcalTowerSumEt()-effA_Hcal*Rho, 0.0))/Electrons[j].pt();
	  }
	  else{
	    string err = "HZZ4LHelper::passedIsolation --- unknown isolation type\n";
	    err+="Possibilites are:\n";
	    err+="\t DET with rho, EffAreaRho, or EffAreaRhoBL\n";
	    err+="\t PF  with rho, EffAreaRho, PFEffAreaRho, or dB\n"; 
	    throw invalid_argument( err );
	  }
      

	  if( Electrons[j].pt() > subleadingPtCut )
	    {
	      if( tmpIso < earlyIsoCut )
		{
		  if( j != pickedE )
		    {
		      passedSubLeading++;
		    }
		}
	    }
    }
  
  if( (passedLeading > 0 && passedSubLeading > 0) || (pickedE<999 and pickedMu<999) ){ passed = true; }
  
  return passed;

}



double HZZ4LHelper::getSIP3D(pat::Muon muon)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  double ip = fabs(muon.dB(pat::Muon::PV3D));
  double ipError = muon.edB(pat::Muon::PV3D);

  double sip = ip/ipError;

  return sip;

}

double HZZ4LHelper::getSIP3D(pat::Electron electron)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  double ip = fabs(electron.dB(pat::Electron::PV3D));
  double ipError = electron.edB(pat::Electron::PV3D);

  double sip = ip/ipError;

  return sip;

}




TString HZZ4LHelper::getStringFromDouble(double number)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  char _string[192];
  sprintf(_string,"%5.1f",number);
  TString Number_string = _string;

  return Number_string;

}


double HZZ4LHelper::getCorrWeightForSChan(double m4l)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;

  if( m4l > 100 ){return 1.0;}

  double MZ0 = 91.2;
  double Const = 0.8;
  double A = 3.95; 
  double B = 2.6;

  double func = Const + A*exp( -pow(m4l-MZ0,2)/pow(B,2));

  return func;

}




bool HZZ4LHelper::isPFMuon(pat::Muon muon,edm::Handle<edm::View<reco::PFCandidate> > pf)
{

  //using namespace edm;
  using namespace pat;
  using namespace std;


  bool ispf = false;
  bool matchByReference = true;
  
    edm::View<reco::PFCandidate>::const_iterator pfit, pfbegin = pf->begin(), pfend = pf->end();
    /// Now loop
    //for (size_t i = 0, n = muons->size(); i < n; ++i) {
      // read the edm reference to the muon
    //RefToBase<pat::Muon> muonRef = muon->refAt(i);
    edm::Ptr<reco::Candidate> muonRef1 = muon.originalObjectRef();

    //Reco::MuonRef muonRef1 = muon.muonRef();
    //edm::Ptr<reco::Candidate> = muon.
    /// Loop on PF
      for (pfit = pfbegin; pfit != pfend; ++pfit) {
	// Get muon ref, skip those that don't have one
	reco::MuonRef pfRef = pfit->muonRef(); 

	if (pfRef.isNull()){/*cout << "No PFCandRef\n";*/ continue;}
	//if (!pfCut_(*pfit)){cout << "No PFCandRef\n"; continue;}

	// Perform the matching
	if (matchByReference) {
	  if (pfRef.id() != muonRef1.id()) /*cout << "No PFCandRef\n";*/continue; //throw cms::Exception("Configuration") 
	    // << "Cannot match by reference the muons from " << muons_.encode() << " (id " << muonRef.id() << ")" 
	    //<< " with the ones referenced by " << pf_.encode() << " (id " << pfRef.id() << ")\n";
	  if (pfRef.key() == muonRef1.key()) {
	    cout << "Found PFCandRef\n";
	    ispf = true;
	    //out->push_back(muonRef);
	    break;
	  }
	} else if (std::abs(muonRef1->eta() - pfRef->eta()) < 1e-5 && 
		   std::abs(muonRef1->phi() - pfRef->phi()) < 1e-5 &&
		   std::abs(muonRef1->pt()  - pfRef->pt() ) < 1e-5) {
	  ispf = true;
	  cout << "Found PFCandRef\n";
	  //out->push_back(muonRef);
	  break;
	}
      }
      //}

  return ispf;

}


double HZZ4LHelper::pfIso(pat::Muon muon, double Rho)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  //double PUCorr = Rho*MuonEffArea(muEAtype,muon.eta(),muEAtarget);
  double PUCorr = 0.5*muon.userIsolation("PfPUChargedHadronIso");
  double iso = (muon.chargedHadronIso()+max(muon.photonIso()+muon.neutralHadronIso()-PUCorr,0.0))/muon.pt();

  
  return iso;
}


double HZZ4LHelper::pfIso(pat::Electron elec, double Rho)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  double PUCorr = Rho*ElecEffArea(elEAtype,elec.superCluster()->eta(),elEAtarget);
  double iso = (elec.chargedHadronIso()+max(elec.photonIso()+elec.neutralHadronIso()-PUCorr,0.0))/elec.pt();
  
  return iso;
}

double HZZ4LHelper::pfIsoFSR(pat::Muon muon, double Rho, double photPt)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  //double PUCorr = Rho*MuonEffArea(muEAtype,muon.eta(),muEAtarget);
  double PUCorr= 0.5*muon.userIsolation("PfPUChargedHadronIso");
  double iso = (muon.chargedHadronIso()+max(muon.photonIso()-photPt+muon.neutralHadronIso()-PUCorr,0.0))/muon.pt();
  
  return iso;
}


double HZZ4LHelper::pfIsoFSR(pat::Electron elec, double Rho, double photPt)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  double PUCorr = Rho*ElecEffArea(elEAtype,elec.superCluster()->eta(),elEAtarget);
  double iso = (elec.chargedHadronIso()+max(elec.photonIso()-photPt+elec.neutralHadronIso()-PUCorr,0.0))/elec.pt();
  
  return iso;
}


double HZZ4LHelper::MuonEffArea(MuonEffectiveAreaType type, double SCEta, MuonEffectiveAreaTarget EffectiveAreaTarget)
{

  double EffectiveArea = 0;
  
  if (EffectiveAreaTarget == kMuNoCorr) {
    return 0.0;
  }
  
  //2012 Data Effective Areas
  else if (EffectiveAreaTarget == kMuEAData2012) {
    if (type == kMuGammaIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.50419;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.30582;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.19765;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.28723;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.52529;
      if (fabs(SCEta) >= 2.3 )                        EffectiveArea = 0.48818;
    }
    if (type == kMuNeutralHadronIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.16580;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.25904;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.24695;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.22021;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.34045;
      if (fabs(SCEta) >= 2.3 )                        EffectiveArea = 0.21592;
    }
    if (type == kMuGammaAndNeutralHadronIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.674;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.565;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.442;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.515;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.821;
      if (fabs(SCEta) >= 2.3 )                        EffectiveArea = 0.660;
    }
    if (type == kMuGammaAndNeutralHadronIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.382;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.317;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.242;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.326;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.462;
      if (fabs(SCEta) >= 2.3 )                        EffectiveArea = 0.372;
    }
    if (type == kMuGammaAndNeutralHadronIso04Tight){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.340;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.310;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.315;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.415;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.658;
      if (fabs(SCEta) >= 2.3 )                        EffectiveArea = 0.405;
    }
    if (type == kMuGammaAndNeutralHadronIso03Tight){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.207;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.183;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.177;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.271;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.348;
      if (fabs(SCEta) >= 2.3 )                        EffectiveArea = 0.246;
    }
  }
  //2011 Data Effective Areas
  else if (EffectiveAreaTarget == kMuEAData2011) {
      
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3 ) EffectiveArea = 0.005;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.011;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.3 ) EffectiveArea = 0.011;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.023;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.014;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.021;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.023;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.028;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.032;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.028;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.033;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.052;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.001;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.007;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.014;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.038;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.033;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.045;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.051;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.114;
    }
    /// BEGIN FROM SLIDE 11 OF  https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=188494
    /// NOTE: to be used with the rho from ALL pf candidates within |eta|<2.5
    if (type == kMuGammaIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.049;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.5 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 1.5 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.022;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.034;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 2.3 )                      EffectiveArea = 0.048;
    }
    if (type == kMuGammaIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.085;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.5 ) EffectiveArea = 0.052;
      if (fabs(SCEta) >= 1.5 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.070;
      if (fabs(SCEta) >= 2.3 )                      EffectiveArea = 0.081;
    }
    if (type == kMuNeutralHadronIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.5 ) EffectiveArea = 0.039;
      if (fabs(SCEta) >= 1.5 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.044;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.047;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 2.3 )                      EffectiveArea = 0.065;
    }
    if (type == kMuNeutralHadronIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.046;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.5 ) EffectiveArea = 0.067;
      if (fabs(SCEta) >= 1.5 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.074;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.083;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.095;
      if (fabs(SCEta) >= 2.3 )                      EffectiveArea = 0.105;
    }
    if (type == kMuGammaAndNeutralHadronIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.076;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.5 ) EffectiveArea = 0.070;
      if (fabs(SCEta) >= 1.5 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.067;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.082;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.097;
      if (fabs(SCEta) >= 2.3 )                      EffectiveArea = 0.115;
    }
    if (type == kMuGammaAndNeutralHadronIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.132;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.5 ) EffectiveArea = 0.120;
      if (fabs(SCEta) >= 1.5 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.114;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.139;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.168;
      if (fabs(SCEta) >= 2.3 )                      EffectiveArea = 0.189;
    }
    if (type == kMuGammaIso05){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.05317;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.03502;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.03689;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.05221;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.06668;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.0744;
    }
    if (type == kMuNeutralIso05) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.06408;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.07557;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.08864;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.11492;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.13784;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.18745;
    }
  }
       
  //Summer11 MC Effective Areas
  else if (EffectiveAreaTarget == kMuEASummer11MC) {
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.006;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.015;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.023;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.036;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.033;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.041;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.062;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.052;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.066;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.093;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.000;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.003;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.005;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.013;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.025;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.044;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.030;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.038;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.048;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.118;
    }
  } 
  //Fall11 MC Effective Areas
  else if (EffectiveAreaTarget == kMuEAFall11MC) {
    if (type == kMuGammaIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.003;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.011;
    }
    if (type == kMuGammaIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.008;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.006;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.019;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.024;
    }
    if (type == kMuGammaIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.020;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.012;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.022;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.027;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.034;
    }
    if (type == kMuGammaIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.042;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.033;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.022;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.059;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.068;
    }
    if (type == kMuGammaIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.060;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.043;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.036;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.055;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.092;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.115;
    }
    if (type == kMuNeutralHadronIsoDR0p0To0p1) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.002;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.004;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.010;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.014;
    }
    if (type == kMuNeutralHadronIsoDR0p1To0p2) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.005;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.007;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.009;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.015;
      if (fabs(SCEta) >= 2.3  )                       EffectiveArea = 0.017;
    }
    if (type == kMuNeutralHadronIsoDR0p2To0p3) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.009;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.015;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.016;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.018;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.022;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.026;
    }
    if (type == kMuNeutralHadronIsoDR0p3To0p4) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.013;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.021;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.032;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.042;
    }
    if (type == kMuNeutralHadronIsoDR0p4To0p5) {
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.017;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.026;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.035;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.046;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.063;
      if (fabs(SCEta) >= 2.3  ) EffectiveArea = 0.135;
    }
  
    return EffectiveArea;
    
  }


  return EffectiveArea;


}





double HZZ4LHelper::ElecEffArea(ElecEffectiveAreaType type, double SCEta, ElecEffectiveAreaTarget EffectiveAreaTarget)
{

  //From: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection
  double EffectiveArea = 0;

  if (EffectiveAreaTarget == kElNoCorr) {
    return 0.0;
  }
  
  if(EffectiveAreaTarget == kElEAData2012)
    {
      
      if (type == kEGammaNeutralHadIso04){
	if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.19;
	if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.25;
	if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.12;
	if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.21;
	if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.27; 
	if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.44;
	if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.52;
      }
	
      return EffectiveArea;
    }

  if(EffectiveAreaTarget == kElEAData2011){

    if (type == kENeutralHadIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.044;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.065;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.068;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.057;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.058; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.061;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.11;
    }
    if (type == kEGammaIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.14;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.13;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.079;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.13;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.15; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.16;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.18;
    }
    if (type == kEGammaNeutralHadIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.18;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.20;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.15;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.19;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.21; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.22;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.29;
    }

    if (type == kENeutralHadIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.024;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.037;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.023;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.023; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.021;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.021;
    }
    if (type == kEGammaIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.081;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.084;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.048;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.089;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.092; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.097;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.11;
    }
    if (type == kEGammaNeutralHadIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.10;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.12;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.085;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.11;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.12; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.12;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.13;
    }

    return EffectiveArea;
  }

  if(EffectiveAreaTarget == kElEAFall11MC){

    if (type == kENeutralHadIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.041;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.068;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.075;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.068;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.072; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.077;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.14;
    }
    if (type == kEGammaIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.14;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.14;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.084;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.15;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.20; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.22;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.26;
    }
    if (type == kEGammaNeutralHadIso04){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.18;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.21;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.16;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.22;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.27; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.30;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.41;
    }

    if (type == kENeutralHadIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.022;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.039;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.040;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.028;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.027; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.024;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.030;
    }
    if (type == kEGammaIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.084;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.090;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.049;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.099;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.122; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.132;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.155;
    }
    if (type == kEGammaNeutralHadIso03){
      if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )   EffectiveArea = 0.11;
      if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.13;
      if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.089;
      if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.13;
      if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )   EffectiveArea = 0.15; 
      if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )   EffectiveArea = 0.16;
      if (fabs(SCEta) >= 2.4 )                        EffectiveArea = 0.19;
    }

    return EffectiveArea;
  }
  
  return EffectiveArea;


}

/*
void HZZ4LHelper::getMuonEffAreaType(MuoneEffectiveAreaType type){type = muEAtype;}
void HZZ4LHelper::getMuonEffAreaTarget(MuonEffectiveAreaTarget target){target = muEAtarget;}
void HZZ4LHelper::getElecEffAreaType(ElecEffectiveAreaType type){type = elEAtype;}
void HZZ4LHelper::getElecEffAreaTarget(ElecEffectiveAreaTarget target){target = elEAtarget;}
*/



/////////////////////////////////////////////////////

std::vector< pat::Muon > HZZ4LHelper::goodMuons2012_Iso(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut, double rho, double isocut, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector< pat::Muon > bestMuons;
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/


  for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu)
    {
      if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isPFMuon() == 1 && (mu->isGlobalMuon() || mu->isTrackerMuon() ) )
	{
	  if( abs(getSIP3D(*mu)) < sipCut )
	    {
	      //if( fabs(mu->innerTrack()->dxy(vertex->position())) < dxyCut ) // for miniAOD
	      if( fabs(mu->muonBestTrack()->dxy(vertex->position())) < dxyCut ) // for miniAOD
		{ 
		  //if( fabs(mu->innerTrack()->dz(vertex->position())) < dzCut ) // for miniAOD
		  if( fabs(mu->muonBestTrack()->dz(vertex->position())) < dzCut ) // for miniAOD
		    {
                      if(pfIso(*mu,rho) < isocut)
                       {
		         bestMuons.push_back(*mu);
                       }
		    }
		}
	    }
	}
    }


  return bestMuons;
  
}


std::vector< pat::Electron > HZZ4LHelper::goodElectrons2012_Iso(edm::Handle<edm::View<pat::Electron> > Electrons, double elecPtCut, double rho, double isocut, std::string elecID, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;
  
  vector< pat::Electron > bestElectrons;
  
  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  int missingHitsCuts = 2;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/

  /*
  ////////Emanuele////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = -0.12;
  double bdtCutLoPt_14_25 = 0.39;
  double bdtCutHiPt_0_08 = 0.49;
  double bdtCutHiPt_08_14 = 0.046;
  double bdtCutHiPt_14_25 = 0.69;
  */

  /*
  ////////Si+Duncan////////
  double bdtCutLoPt_0_08 = 0.369;
  double bdtCutLoPt_08_14 = -0.025;
  double bdtCutLoPt_14_25 = 0.531;
  double bdtCutHiPt_0_08 = 0.735;
  double bdtCutHiPt_08_14 = 0.467;
  double bdtCutHiPt_14_25 = 0.795;
  */
  
  ////////Christophe////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = 0.004;
  double bdtCutLoPt_14_25 = 0.295;
  double bdtCutHiPt_0_08 = 0.5;
  double bdtCutHiPt_08_14 = 0.12;
  double bdtCutHiPt_14_25 = 0.6;
  
  ///MVA ID = mvaNonTrigV0
  
  for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec)
    {
      if( abs(getSIP3D(*elec)) < sipCut && fabs(elec->gsfTrack()->dxy(vertex->position())) < dxyCut && fabs(elec->gsfTrack()->dz(vertex->position())) < dzCut )
	{
	  if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
	    { 
	      double cutVal = 1000;
	      if( elec->pt() > 0 && elec->pt() < 10.0 )
		{
		  if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutLoPt_0_08;}
		  if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutLoPt_08_14;}
		  if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutLoPt_14_25;}
		}
	      if( elec->pt() >= 10.0 )
		{
		  if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutHiPt_0_08;}
		  if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutHiPt_08_14;}
		  if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutHiPt_14_25;}
		}
	      //if(cutVal == 1000){ cout << "ele->superCluster()->eta() > 2.5" << endl; continue;}
	      //int misHits = elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
	      int misHits = 0; // for miniAOD
	      if(elec->electronID(elecID) > cutVal && misHits < missingHitsCuts)
		{
                  if(pfIso(*elec,rho) < isocut)
                    {
		     bestElectrons.push_back(*elec);
                    }	  
		}
	    }
	}
    }

  
  return bestElectrons;
  
  
}



///////////////////////////////

std::vector< pat::Muon > HZZ4LHelper::goodMuons2012_Iso(std::vector< pat::Muon > Muons, double muPtCut, double rho, double isocut, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  vector<pat::Muon > bestMuons;
  /********** M U O N  C U T S **********/
  double muEtaCut = 2.4;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/


  for(unsigned int i=0; i < Muons.size(); i++)
    {
      pat::Muon *mu = &(Muons[i]);
      if( mu->pt() > muPtCut && abs(mu->eta()) < muEtaCut && mu->isPFMuon() == 1 && (mu->isGlobalMuon() || mu->isTrackerMuon() ) )
        {
          if( abs(getSIP3D(*mu)) < sipCut )
            {
              if( fabs(mu->innerTrack()->dxy(vertex->position())) < dxyCut )
                {
                  if( fabs(mu->innerTrack()->dz(vertex->position())) < dzCut )
                    {
                      if(pfIso(*mu,rho) < isocut)
                       {
                         bestMuons.push_back(*mu);
                       }
                    }
                }
            }
        }
    }


  return bestMuons;
}


std::vector< pat::Electron > HZZ4LHelper::goodElectrons2012_Iso(std::vector< pat::Electron > Electrons, double elecPtCut, double rho, double isocut, std::string elecID, const reco::Vertex *&vertex)
{
  //using namespace edm;
  using namespace pat;
  using namespace std;

  vector< pat::Electron > bestElectrons;

  /****** E L E C T R O N  C U T S ******/
  double elecEtaCut = 2.5;
  int missingHitsCuts = 2;
  double sipCut = 4.0;
  double dxyCut = 0.5;
  double dzCut = 1;
  /**************************************/

  /*
  ////////Emanuele////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = -0.12;
  double bdtCutLoPt_14_25 = 0.39;
  double bdtCutHiPt_0_08 = 0.49;
  double bdtCutHiPt_08_14 = 0.046;
  double bdtCutHiPt_14_25 = 0.69;
  */
  /*
  ////////Si+Duncan////////
  double bdtCutLoPt_0_08 = 0.369;
  double bdtCutLoPt_08_14 = -0.025;
  double bdtCutLoPt_14_25 = 0.531;
  double bdtCutHiPt_0_08 = 0.735;
  double bdtCutHiPt_08_14 = 0.467;
  double bdtCutHiPt_14_25 = 0.795;
  */

  ////////Christophe////////
  double bdtCutLoPt_0_08 = 0.47;
  double bdtCutLoPt_08_14 = 0.004;
  double bdtCutLoPt_14_25 = 0.295;
  double bdtCutHiPt_0_08 = 0.5;
  double bdtCutHiPt_08_14 = 0.12;
  double bdtCutHiPt_14_25 = 0.6;

  ///MVA ID = mvaNonTrigV0

  for(unsigned int i = 0; i < Electrons.size(); i++)
    {
      pat::Electron *elec = &(Electrons[i]);

    if( abs(getSIP3D(*elec)) < sipCut && fabs(elec->gsfTrack()->dxy(vertex->position())) < dxyCut && fabs(elec->gsfTrack()->dz(vertex->position())) < dzCut )
        {
          if( elec->pt() > elecPtCut && abs(elec->eta()) < elecEtaCut)
            {
              double cutVal = 1000;
              if( elec->pt() > 0 && elec->pt() < 10.0 )
                {
                  if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutLoPt_0_08;}
                  if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutLoPt_08_14;}
                  if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutLoPt_14_25;}
                }
              if( elec->pt() >= 10.0 )
                {
                  if(fabs(elec->superCluster()->eta()) > 0 && fabs(elec->superCluster()->eta()) < 0.8 ){ cutVal = bdtCutHiPt_0_08;}
                  if(fabs(elec->superCluster()->eta()) >= 0.8 && fabs(elec->superCluster()->eta()) < 1.479 ){ cutVal = bdtCutHiPt_08_14;}
                  if(fabs(elec->superCluster()->eta()) >= 1.479 && fabs(elec->superCluster()->eta()) < elecEtaCut ){ cutVal = bdtCutHiPt_14_25;}
                }
              //if(cutVal == 1000){ cout << "ele->superCluster()->eta() > 2.5" << endl; continue;}
              //int misHits = elec->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
              int misHits = 0; // for miniAOD
              if(elec->electronID(elecID) > cutVal && misHits < missingHitsCuts)
                {
                  if(pfIso(*elec,rho) < isocut)
                    {
                     bestElectrons.push_back(*elec);
                    }
                }
            }
        }
    }


  return bestElectrons;

}




#endif
