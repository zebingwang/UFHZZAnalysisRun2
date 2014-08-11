// -*- C++ -*-
//
// Package:    UFHZZ4LAna
// Class:      UFHZZ4LAna
// 
/*
  class UFHZZ4LAna UFHZZ4LAna.cc UFHZZAnalysis8TeV/UFHZZ4LAna/src/UFHZZ4LAna.cc
  
  Description: UF HZZ4L Analysis Analyzer. Works in CMSSW 53X
  
  Implementation: Full analysis step for HZZ4L in 2012.
  
  Last updated: 05.06.2013 --- Matt Snowball
  
*/
//
// Original Author:  Matthew Snowball, snowball@phys.ufl.edu
//
//

// system include files
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

#define PI 3.14159

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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"

// PAT
#include "DataFormats/PatCandidates/interface/PFParticle.h"
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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"  
//#include "CMGTools/External/interface/PileupJetIdAlgo.h"  // removed for miniAOD
//#include "CMGTools/External/interface/PileupJetIdentifier.h"  // removed for miniAOD
// Reco
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"


//Angles
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LAngles.h"
//#include "ZZMatrixElement/MELA/interface/Mela.h" // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/PseudoMELA.h"  // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/SpinOneEvenMELA.h"  // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/SpinOneOddMELA.h"  // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/SpinTwoMinimalMELA.h"  // removed for miniAOD

//MEKD
//#include "ZZMatrixElement/MEKD/interface/MEKD.h"  // removed for miniAOD
//#include "ZZMatrixElement/MEKD/interface/MEKD_MG.h"  // removed for miniAOD
//#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"  // removed for miniAOD

//Helper
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LHelper.h"
//Muons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonAna.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonTree.h"
//Electrons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LElectronTree.h"
//Photons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPhotonTree.h"
//Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJetTree.h"
//Final Leps
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LFinalLepTree.h"
//Sip
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LSipAna.h"
//PU
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPileUp.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//Iso
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LIsoEff.h"
//Mass Error
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMassErr.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPerLepResolution.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LDiLepResolution.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LResolution.h"

//Signal Efficiency
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LSigEff.h"
//Z->4L
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LZto4LAna.h"
//GEN
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LGENAna.h"
//VBF Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJets.h"

//
// class declaration
//

class UFHZZ4LAna : public edm::EDAnalyzer {
public:
  explicit UFHZZ4LAna(const edm::ParameterSet&);
  ~UFHZZ4LAna();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, const edm::EventSetup& iSetup);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup);
  
  void findHiggsCandidate(std::vector<pat::Muon> &candMuons, std::vector<pat::Electron> &candElectrons, 
                  std::vector<pat::PFParticle> fsrPhotons,
                  std::vector<double> deltaRVec,
                  std::vector<pat::Muon> &selectedMuons, std::vector<pat::Electron> &selectedElectrons,
                  std::vector<pat::PFParticle > &selectedFsrPhotons, const edm::Event& iEvent);

  double getMinDeltaR(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons);
  void plotMinDeltaRemu(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons);
   
  void bookStepPlots();
  void fillStepPlots();
  void bookResolutionHistograms();
  void fillResolutionHistograms(edm::Handle<edm::View<pat::Muon> > muons);
  bool findZ(std::vector<pat::PFParticle> photons, std::vector<double> deltaRVec, 
             pat::Muon &muon1, pat::Muon &muon2,int taken1, int &taken2,
             int &assocMuon,math::XYZTLorentzVector &ZVec, math::XYZTLorentzVector &photVec, bool &foundPhoton);
  
  bool findZ(std::vector<pat::PFParticle> photons, std::vector<double> deltaRVec, 
             pat::Electron &electron1, pat::Electron &electron2,int taken1, int &taken2, 
             int &assocElec, math::XYZTLorentzVector &ZVec, math::XYZTLorentzVector &photVec, bool &foundPhoton);
  
  //--------------------------------------------    
  void findHiggsCandidate_MixFlavour(std::vector<pat::Muon> &candMuons, std::vector<pat::Electron> &candElectrons,
                             std::vector<pat::Muon> &selectedMuons, std::vector<pat::Electron> &selectedElectrons, 
                             bool noChargeReq = false);
  int RecoFourMixEvent;
  bool mixedFlavorCharge;
  //--------------------------------------------
  

  //MELA
  HZZ4LAngles angles;
  // Mela *mela;   // removed for miniAOD


  //PseudoMELA pseudoMela;
  //SpinOneEvenMELA spin1EvMela;
  //SpinOneOddMELA spin1OddMela;
  //SpinTwoMinimalMELA spin2Mela;

  double melaLD, mela_Sig, mela_Bkg;
  float melaLD_Pt, mela_Sig_Pt, mela_Bkg_Pt;
  float melaLD_Y, mela_Sig_Y, mela_Bkg_Y;
  float melaLD_PtY, mela_Sig_PtY, mela_Bkg_PtY;
  float pseudoMelaLD, pseudoMela_SM, pseudoMela_PS;
  float spin1EvMelaLD, spin1EvMela_SM, spin1EvMela_S1E;
  float spin1OddMelaLD, spin1OddMela_SM, spin1OddMela_S1O;
  float spin2MelaLD, spin2Mela_SM, spin2Mela_S2;

  double pdfSigM4l, pdfSigM4l_noFSR, pdfBkgM4l, pdfBkgM4l_noFSR;
  double pdfSigM4l_ScaleUp, pdfSigM4l_ScaleUp_noFSR, pdfBkgM4l_ScaleUp, pdfBkgM4l_ScaleUp_noFSR;
  double pdfSigM4l_ScaleDown, pdfSigM4l_ScaleDown_noFSR, pdfBkgM4l_ScaleDown, pdfBkgM4l_ScaleDown_noFSR;
  double pdfSigM4l_ResUp, pdfSigM4l_ResUp_noFSR, pdfBkgM4l_ResUp, pdfBkgM4l_ResUp_noFSR;
  double pdfSigM4l_ResDown, pdfSigM4l_ResDown_noFSR, pdfBkgM4l_ResDown, pdfBkgM4l_ResDown_noFSR;
 
  float p0plus_melaNorm;   // higgs, analytic distribution, normalized as for normal MELA distribution     
  float p0plus_mela;   // higgs, analytic distribution          
  float p0minus_mela;  // pseudoscalar, analytic distribution 
  float p0plus_VAJHU;  // higgs, vector algebra, JHUgen
  float p0minus_VAJHU; // pseudoscalar, vector algebra, JHUgen
  float p0plus_VAMCFM; // higgs, vector algebra, MCFM
  float p1_mela;  // zprime, analytic distribution not yet functional
  float p1_VAJHU; // zprime, vector algebra, JHUgen, not yet functional
  float p2_mela; // graviton, analytic distribution 
  float p2_VAJHU; // graviton, vector algebra, JHUgen,
  float mela_bkg_analytic;  // background,  analytic distribution 
  float mela_bkg_VAMCFM; // background, vector algebra, MCFM
  float mela_bkg_ggzz_VAMCFM; // background, vector algebra, MCFM
  float mela_bkg_VAMCFMNorm; // background, vector algebra, MCFM, Normalized 
  float mela_p0_pt; // multiplicative probability for signal pt
  float mela_p0_y; // multiplicative probability for signal y
  float mela_bkg_pt; // multiplicative probability for bkg pt
  float mela_bkg_y; // multiplicative probability for bkg y
  float p0plus_m4l;
  float bkg_m4l;

  //MEKD
  // MEMs *MEMsnoPDFs_noFSR; // removed for miniAOD
  // MEMs *MEMsnoPDFs;  // removed for miniAOD
  
  double interferenceWeight;
  double JHUKD_H_qqZZ_noPDF, JHU_ME_H, MCFM_ME_qqZZ;
  double JHUKD_H_h0M_noPDF,  JHU_ME_h0M;
  double JHUKD_H_h0P_noPDF,  JHU_ME_h0P;
  double JHUKD_H_h1M_noPDF,  JHU_ME_h1M;
  double JHUKD_H_h1P_noPDF,  JHU_ME_h1P;
  double JHUKD_H_ggh2P_noPDF, JHU_ME_ggh2P;
  double JHUKD_H_qqh2P_noPDF, JHU_ME_qqh2P;
  
  double JHUKD_H_h2hP_noPDF, JHU_ME_h2hP;
  double JHUKD_H_h2hM_noPDF, JHU_ME_h2hM;
  double JHUKD_H_h2bP_noPDF, JHU_ME_h2bP;
  double JHUKD_H_h2P_prodInd_noPDF, JHU_ME_h2P_prodInd;
  double JHUKD_H_h1P_prodInd_noPDF, JHU_ME_h1P_prodInd;
  double JHUKD_H_h1M_prodInd_noPDF, JHU_ME_h1M_prodInd;

  double MEKD_noPDF_noFSR, MEKD_ME_H_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR;
  double MEKD_h0M_ZZ_noPDF_noFSR, MEKD_ME_h0M_noPDF_noFSR;
  double MEKD_h0P_ZZ_noPDF_noFSR, MEKD_ME_h0P_noPDF_noFSR;
  double MEKD_h1M_ZZ_noPDF_noFSR, MEKD_ME_h1M_noPDF_noFSR;
  double MEKD_h1P_ZZ_noPDF_noFSR, MEKD_ME_h1P_noPDF_noFSR;
  double MEKD_qqh2P_ZZ_noPDF_noFSR,MEKD_ME_qqh2P_noPDF_noFSR;
  double MEKD_ggh2P_ZZ_noPDF_noFSR,MEKD_ME_ggh2P_noPDF_noFSR;

  double MEKD_h2hP_ZZ_noPDF_noFSR, MEKD_ME_h2hP_noPDF_noFSR;
  double MEKD_h2hM_ZZ_noPDF_noFSR, MEKD_ME_h2hM_noPDF_noFSR;
  double MEKD_h2bP_ZZ_noPDF_noFSR, MEKD_ME_h2bP_noPDF_noFSR;
  double MEKD_h2P_ZZ_prodInd_noPDF_noFSR, MEKD_ME_h2P_prodInd_noPDF_noFSR;
  double MEKD_h1P_ZZ_prodInd_noPDF_noFSR, MEKD_ME_h1P_prodInd_noPDF_noFSR;
  double MEKD_h1M_ZZ_prodInd_noPDF_noFSR, MEKD_ME_h1M_prodInd_noPDF_noFSR;
  
  double JHUKD_H_qqZZ_noPDF_noFSR,  JHU_ME_H_noPDF_noFSR,   JHU_ME_ZZ_noPDF_noFSR;
  double JHUKD_h0M_ZZ_noPDF_noFSR,  JHU_ME_h0M_noPDF_noFSR;
  double JHUKD_h0P_ZZ_noPDF_noFSR,  JHU_ME_h0P_noPDF_noFSR;
  double JHUKD_h1M_ZZ_noPDF_noFSR,  JHU_ME_h1M_noPDF_noFSR;
  double JHUKD_h1P_ZZ_noPDF_noFSR,  JHU_ME_h1P_noPDF_noFSR;
  double JHUKD_qqh2P_ZZ_noPDF_noFSR,JHU_ME_qqh2P_noPDF_noFSR;
  double JHUKD_ggh2P_ZZ_noPDF_noFSR,JHU_ME_ggh2P_noPDF_noFSR;
  
  double MEKD_wPDF, MEKD_ME_H_wPDF, MEKD_ME_ZZ_wPDF;
  double MEKD_noPDF, MEKD_ME_H_noPDF, MEKD_ME_ZZ_noPDF;
  double MEKD_wPDF_noFSR, MEKD_ME_H_wPDF_noFSR, MEKD_ME_ZZ_wPDF_noFSR;
 
  //Helper Class
  HZZ4LHelper helper;
  
  //Muon Class
  //HZZ4LMuonAna muonAna;

  //Dump Class
  HZZ4LMuonTree *muonDump;
  HZZ4LElectronTree *electronDump;
  HZZ4LPhotonTree *photonDump;
  HZZ4LFinalLepTree *finalLepDump;
  HZZ4LJetTree *jetDump;
  
  //Sip Class
  HZZ4LSipAna sipAna;

  //Iso Class
  HZZ4LIsoEff *isoEff;
  //HZZ4LIsoEff *isoEff_Z4l;

  //Mass Err
  HZZ4LMassErr massErr;

  //Sig Eff
  HZZ4LSigEff *sigEff_4,*sigEff_6,*sigEff_8,*sigEff_9,*sigEff_10,*sigEff_11,*sigEff_12;

  //Zto4L
  HZZ4LZto4LAna Zto4LAna;
  HZZ4LZto4LAna *Zto4LAnaOP;

  //GEN
  HZZ4LGENAna genAna;
  
  //VBF
  HZZ4LJets jetHelper;

  //PU Reweighting
  edm::LumiReWeighting *lumiWeight;
  HZZ4LPileUp pileUp;

  //Saved Events Trees
  TTree *passedEventsTree_All;
  void bookPassedEventTree(TString treeName, TTree *tree);
  void setTreeVariables( const edm::Event&, const edm::EventSetup&, 
                         std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, 
                         std::vector<pat::Jet> selectedVBFJets, std::vector<pat::Jet> correctedVBFJets);
  void setGENVariables(std::vector<reco::GenParticle> Higgs, 
                       std::vector<reco::GenParticle> Zs, 
                       std::vector<reco::GenParticle> leptonsS1, std::vector<reco::GenParticle> leptonsS3);
  void setGENMatchedVariables(std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons);

  //Variables
  bool notDuplicateEvent;
  ULong64_t Run, Event, LumiSect;
  double mass4l, mass4e, mass4mu, mass2e2mu, pT4l, massZ1, massZ2;
  double eta4l, phi4l;
  int idL1, idL2, idL3, idL4;
  double mvaL1, mvaL2, mvaL3, mvaL4;
  double SipL1, SipL2, SipL3, SipL4;
  double IPL1, IPL2, IPL3, IPL4;
  double dIPL1,dIPL2,dIPL3,dIPL4;
  double pTL1, pTL2, pTL3, pTL4;
  double pXL1, pXL2, pXL3, pXL4;
  double pYL1, pYL2, pYL3, pYL4;
  double pZL1, pZL2, pZL3, pZL4;
  double EL1, EL2, EL3, EL4;
  double pTL1FSR, pTL2FSR, pTL3FSR, pTL4FSR;
  double pXL1FSR, pXL2FSR, pXL3FSR, pXL4FSR;
  double pYL1FSR, pYL2FSR, pYL3FSR, pYL4FSR;
  double pZL1FSR, pZL2FSR, pZL3FSR, pZL4FSR;
  double EL1FSR, EL2FSR, EL3FSR, EL4FSR;
  double EZ1, EZ2, pTZ1, pTZ2, pXZ1, pXZ2;
  double pYZ1, pYZ2, pZZ1, pZZ2;
  int chargeL1, chargeL2, chargeL3, chargeL4;
  double etaL1, etaL2, etaL3, etaL4;
  double phiL1, phiL2, phiL3, phiL4;
  double isoTrackL1, isoTrackL2, isoTrackL3, isoTrackL4;
  double isoEcalL1, isoEcalL2, isoEcalL3, isoEcalL4;
  double isoHcalL1, isoHcalL2, isoHcalL3, isoHcalL4;
  double isoNHL1, isoNHL2, isoNHL3, isoNHL4;
  double isoCHL1, isoCHL2, isoCHL3, isoCHL4;
  double isoPhotL1, isoPhotL2, isoPhotL3, isoPhotL4;
  double RelIsoL1, RelIsoL2, RelIsoL3, RelIsoL4;
  double RelIsoUCL1, RelIsoUCL2, RelIsoUCL3, RelIsoUCL4;
  int missingHitsL1, missingHitsL2, missingHitsL3, missingHitsL4;
  bool ecalSeedL1, ecalSeedL2, ecalSeedL3, ecalSeedL4;
  double muRho, elRho;
  double worstIso, worstSIP;
  float cosTheta1, cosTheta2, Phi, cosThetaStar, phiStar1, phiStar2, phiStar12, Phi1, Phi2;
  int nJets, nVtx, nPhotons;
  double metVal, rawRho;
  double m4lErrUCSD, m4lErrUCSDCorr, m4lErrUF, m4lErrUFCorr, m4lErrUFADCorr;
  
  bool FSR_Z1, FSR_Z2;
  double FSRPhot1_Pt, FSRPhot2_Pt;
  double FSRPhot1_eta, FSRPhot2_eta;
  double FSRPhot1_phi, FSRPhot2_phi;
  double FSRPhot1_Px, FSRPhot2_Px;
  double FSRPhot1_Py, FSRPhot2_Py;
  double FSRPhot1_Pz, FSRPhot2_Pz;
  double FSRPhot1_E, FSRPhot2_E;


  //GEN
  int GENidL1_S1, GENidL2_S1, GENidL3_S1, GENidL4_S1;
  double GENpTL1_S1, GENpTL2_S1, GENpTL3_S1, GENpTL4_S1;
  double GENpXL1_S1, GENpXL2_S1, GENpXL3_S1, GENpXL4_S1;
  double GENpYL1_S1, GENpYL2_S1, GENpYL3_S1, GENpYL4_S1;
  double GENpZL1_S1, GENpZL2_S1, GENpZL3_S1, GENpZL4_S1;
  double GENEL1_S1, GENEL2_S1, GENEL3_S1, GENEL4_S1;
  double GENEZ1, GENEZ2, GENpTZ1, GENpTZ2;
  double GENpXZ1, GENpXZ2, GENpYZ1, GENpYZ2, GENpZZ1, GENpZZ2;
  int GENchargeL1_S1, GENchargeL2_S1, GENchargeL3_S1, GENchargeL4_S1;
  double GENetaL1_S1, GENetaL2_S1, GENetaL3_S1, GENetaL4_S1;
  double GENphiL1_S1, GENphiL2_S1, GENphiL3_S1, GENphiL4_S1;

  int GENidL1_S3, GENidL2_S3, GENidL3_S3, GENidL4_S3;
  double GENpTL1_S3, GENpTL2_S3, GENpTL3_S3, GENpTL4_S3;
  double GENpXL1_S3, GENpXL2_S3, GENpXL3_S3, GENpXL4_S3;
  double GENpYL1_S3, GENpYL2_S3, GENpYL3_S3, GENpYL4_S3;
  double GENpZL1_S3, GENpZL2_S3, GENpZL3_S3, GENpZL4_S3;
  double GENEL1_S3, GENEL2_S3, GENEL3_S3, GENEL4_S3;
  int GENchargeL1_S3, GENchargeL2_S3, GENchargeL3_S3, GENchargeL4_S3;
  double GENetaL1_S3, GENetaL2_S3, GENetaL3_S3, GENetaL4_S3;
  double GENphiL1_S3, GENphiL2_S3, GENphiL3_S3, GENphiL4_S3;

  double GENMH, GENM4L, GENMZ1, GENMZ2;

  int idL1_GENMatched, idL2_GENMatched, idL3_GENMatched, idL4_GENMatched;
  double pTL1_GENMatched, pTL2_GENMatched, pTL3_GENMatched, pTL4_GENMatched;
  double pXL1_GENMatched, pXL2_GENMatched, pXL3_GENMatched, pXL4_GENMatched;
  double pYL1_GENMatched, pYL2_GENMatched, pYL3_GENMatched, pYL4_GENMatched;
  double pZL1_GENMatched, pZL2_GENMatched, pZL3_GENMatched, pZL4_GENMatched;
  double EL1_GENMatched, EL2_GENMatched, EL3_GENMatched, EL4_GENMatched;
  double EZ1_GENMatched, EZ2_GENMatched, pTZ1_GENMatched, pTZ2_GENMatched;
  double pXZ1_GENMatched, pXZ2_GENMatched, pYZ1_GENMatched, pYZ2_GENMatched;
  double pZZ1_GENMatched, pZZ2_GENMatched;
  int chargeL1_GENMatched, chargeL2_GENMatched, chargeL3_GENMatched, chargeL4_GENMatched;
  double etaL1_GENMatched, etaL2_GENMatched, etaL3_GENMatched, etaL4_GENMatched;
  double phiL1_GENMatched, phiL2_GENMatched, phiL3_GENMatched, phiL4_GENMatched;
  double m4l_GENMatched;
  
  //Calculate Angles
  TLorentzVector HP4;
  TLorentzVector Z1P4, L11P4, L12P4;
  TLorentzVector Z2P4, L21P4, L22P4;
  TLorentzVector L11P4_noFSR, L12P4_noFSR, L21P4_noFSR, L22P4_noFSR;
  
  math::XYZTLorentzVector Lep1,Lep2,Lep3,Lep4;

  double theta12, theta13, theta14;
  
  //All leptons 
  int extraLep_id[10];
  double extraLep_pT[10], extraLep_iso[10], extraLep_e[10];
  double extraLep_pX[10], extraLep_pY[10], extraLep_pZ[10];
  double extraLep_eta[10], extraLep_phi[10], extraLep_sip[10];
  double extraLep_chIso[10], extraLep_nhIso[10], extraLep_phIso[10];
  
  
  // ----------member data ---------------------------
  std::map<std::string,TH1F*> histContainer_;
  std::map<std::string,TH2F*> histContainer2D_;
 
  // mass err
  std::map<TString,TH1F*> hContainer_;
  std::map<TString,TH2F*> hContainer2D_;
  std::map<TString,TH3F*> hContainer3D_;
  //////////////////////////////////
 
  const double Zmass;
  
  //MUST HAVE THESE FOR HIGGS CANDIDATE SELECTION 
  bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;
  bool eventPassedPtAndIsoCuts, foundHiggsCandidate, twoLooseIsoLeptons;
  bool passedM4lCut, isIsolated, passedSIP3D, passedPtCuts;
  double mZ1, mZ2, m4l, m4lNoFSR;
  math::XYZTLorentzVector Z1Vec, Z2Vec, HiggsCandVec,HiggsCandVecNoFSR;
  bool twoLep_ID,fourLep_Cleaned;
  bool passedFullSelection, passedZ4lSelection, passedQCDcut;

 
  //Input tags
  edm::InputTag photonSrc_;
  edm::InputTag elecSrc_;
  edm::InputTag muonSrc_;
  //edm::InputTag correctedJetSrc_; // removed for miniAOD
  edm::InputTag jetSrc_;
  edm::InputTag metSrc_;
  edm::InputTag vertexSrc_;
  edm::InputTag muRhoSrc_;
  edm::InputTag elRhoSrc_;

  //Taken from config
  double mZ1Low, mZ2Low;
  double mZ1High, mZ2High;
  double m4lLowCut;
  std::string elecID;
  bool isMC, isSignal;
  unsigned int mH;
  double CrossSection, FilterEff, Lumi;
  bool weightEvents;
  double isoCut, earlyIsoCut, sip3dCut;
  double leadingPtCut, subleadingPtCut;
  
  //Counters
  int nEventsTotal,nEventsSkimmed;
  double nEvAfterSkim, nEvPassedHlt, nEvPassedPtCut;
  double nEvWith2IsoLep, nEvWith4lep, nEvBeforeZCuts;
  double nEvAfterM4lCut, nEvAfterZ1Cut, nEvAfterZ2Cut, nEvAfterZZCut;
  double nEvAfterIso, nEvAfterIso_comb, nEvAfterSIP, nEvAfterMelaCut;
  double nEvAfterZ1Formed, nEvAfterZ1Formed_2mu, nEvAfterZ1Formed_2e;
  double nEvAfterZ2Formed, nEvAfterZ2Formed_4mu, nEvAfterZ2Formed_4e, nEvAfterZ2Formed_2e2mu;
  double nEv4GoodLep, nEv4GoodLep_4mu, nEv4GoodLep_4e, nEv4GoodLep_2e2mu;
  double nEvAfterVBFJet1, nEvAfterVBFJet1_4mu, nEvAfterVBFJet1_4e, nEvAfterVBFJet1_2e2mu;
  double nEvAfterVBFJet2, nEvAfterVBFJet2_4mu, nEvAfterVBFJet2_4e, nEvAfterVBFJet2_2e2mu;
  double nEvAfterVBFJetCuts, nEvAfterVBFJetCuts_4mu, nEvAfterVBFJetCuts_4e, nEvAfterVBFJetCuts_2e2mu;

  double nEvAfterMelaCut_4e, nEvAfterMelaCut_4mu, nEvAfterMelaCut_2e2mu;

  double nEvAfterPsMelaCut_4e, nEvAfterPsMelaCut_4mu, nEvAfterPsMelaCut_2e2mu, nEvAfterPsMelaCut;
  double nEvAfterGrMelaCut_4e, nEvAfterGrMelaCut_4mu, nEvAfterGrMelaCut_2e2mu, nEvAfterGrMelaCut;

  double nEvBeforeZCuts_4mu, nEvAfterM4lCut_4mu, nEvAfterZ1Cut_2mu;
  double nEvAfterZ2Cut_4mu, nEvAfterIso_4mu, nEvAfterSIP_4mu;
  
  double nEvBeforeZCuts_4e, nEvAfterM4lCut_4e, nEvAfterZ1Cut_2e;
  double nEvAfterZ2Cut_4e, nEvAfterIso_4e, nEvAfterSIP_4e;
  
  double nEvAfterZZCut_4e, nEvAfterZZCut_4mu, nEvAfterZZCut_2e2mu;

  double nEvBeforeZCuts_2e2mu, nEvAfterM4lCut_2e2mu, nEvAfterZ1Cut_2e2mu;
  double nEvAfterZ2Cut_2e2mu, nEvAfterIso_2e2mu, nEvAfterSIP_2e2mu;
  
  double nEvAfterId, nEvAfterId_4mu, nEvAfterId_4e, nEvAfterId_2e2mu;
  double nEvAfterCleaning, nEvAfterCleaning_4mu, nEvAfterCleaning_4e, nEvAfterCleaning_2e2mu;
  double nEvPassedPtCut_4mu, nEvPassedPtCut_4e, nEvPassedPtCut_2e2mu;
  double nEvPassedPtCut2_4mu, nEvPassedPtCut2_4e, nEvPassedPtCut2_2e2mu,nEvPassedPtCut2;
  double nEvAfterZ4lCut_4mu, nEvAfterZ4lCut_4e, nEvAfterZ4lCut_2e2mu, nEvAfterZ4lCut;
  double gen4mu,gen4e,gen2e2mu;
  double gen4muPseudo, gen4ePseudo,gen2e2muPseudo; 

  double counterMuZ1ph, counterElZ1ph;
  double counterMuZ2ph, counterElZ2ph;

  double nEvWith1FSRZ, nEvWith1FSRZ_4e, nEvWith1FSRZ_4mu, nEvWith1FSRZ_2e2mu;
  double nEvWith2FSRZ, nEvWith2FSRZ_4e, nEvWith2FSRZ_4mu, nEvWith2FSRZ_2e2mu;

  double nEvFSRPtLt4, nEvFSRPtGt4dR0p5MatchFsrISO, nEvFSRPtGt4dR0p07Match; 
  double nEvFSRPtLt4Zmm, nEvFSRPtGt4dR0p5MatchFsrISOZmm, nEvFSRPtGt4dR0p07MatchZmm; 
  double nEvFSRPtLt4Zee, nEvFSRPtGt4dR0p5MatchFsrISOZee, nEvFSRPtGt4dR0p07MatchZee; 

  // register to the TFileService
  edm::Service<TFileService> fs;

  //Event Weights
  double eventWeight, scaleWeight, MC_weight;
  
  //Rho Correction
  const double deltaRmu, deltaRe;
  double Rho, muonRho, elecRho;
  
  double leastIso, highestSip;
  
  //data duplicate event
  std::vector<ULong64_t> runVec, lumiVec, eventVec;
  
  //pt Cuts
  double _elecPtCut, _muPtCut;

  //ID's only
  bool tightIdsOnly;
  bool looseIdsOnly;
  

  //PU Reweighting
  bool reweightForPU;
  double npv, BX;
  double PU_weight, PT_weight;
  bool interactiveRun;
  std::string PUVersion;

  int finalState;

  double minM3l, Z4lmaxP, minDeltR, m3l_soft;

  double minMass2Lep, maxMass2Lep;
  double thetaPhoton,thetaPhoton_deg, thetaPhotonZ, thetaPhotonZ_deg;
  double theta12_deg, theta13_deg, theta14_deg;

  int counter_4mu, counter_4e, counter_2e2mu, counter_2mu2e;
  int LMcounter_4mu, LMcounter_4e, LMcounter_2e2mu, LMcounter_2mu2e;
  int Z4lcounter_4mu, Z4lcounter_4e, Z4lcounter_2e2mu, Z4lcounter_2mu2e;

  double thetaZ2, etaZ2;
 
  bool doVarDump, doBlinding, doFsrRecovery;


  //Resolution
  bool bStudyResolution;
  bool bStudyDiLeptonResolution;
  bool bStudyFourLeptonResolution;
  double massErrorUCSD, massErrorUCSDCorr, massErrorUF, massErrorUFCorr, massErrorUFADCorr;
  string eventType;
  HZZ4LPerLepResolution * PerLepReso;
  HZZ4LDiLepResolution * DiLepReso;
  HZZ4LResolution * FourLepReso;

  //VBF
  bool VBFJet1, VBFJet2;
  double pTVBFJet1, pTVBFJet2;
  double etaVBFJet1, etaVBFJet2;
  double phiVBFJet1, phiVBFJet2;
  double massVBFJet1, massVBFJet2;
  double VBFDiJetMass, VBFDeltaEta;
  double FisherDiscrim;

  TString tmpEvent;

};


UFHZZ4LAna::UFHZZ4LAna(const edm::ParameterSet& iConfig) :
  mixedFlavorCharge(iConfig.getUntrackedParameter<bool>("mixedFlavorCharge",false)),
  histContainer_(),
  histContainer2D_(),
  Zmass(91.1876),
  photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
  elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
  muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
  // correctedJetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("correctedJetSrc" )), // removed for miniAOD
  jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" )),
  metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc" )),
  vertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")),
  muRhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc")),
  elRhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc")),
  mZ1Low(iConfig.getUntrackedParameter<double>("mZ1Low",40)),
  mZ2Low(iConfig.getUntrackedParameter<double>("mZ2Low",0)),
  mZ1High(iConfig.getUntrackedParameter<double>("mZ1High",120)),
  mZ2High(iConfig.getUntrackedParameter<double>("mZ2High",120)),
  m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",100)),
  elecID(iConfig.getUntrackedParameter<std::string>("elecID","eidTight")),
  isMC(iConfig.getUntrackedParameter<bool>("isMC",true)),
  isSignal(iConfig.getUntrackedParameter<bool>("isSignal",false)),
  mH(iConfig.getUntrackedParameter<unsigned int>("mH",0)),
  CrossSection(iConfig.getUntrackedParameter<double>("CrossSection",1.0)),
  FilterEff(iConfig.getUntrackedParameter<double>("FilterEff",1.0)),
  Lumi(iConfig.getUntrackedParameter<double>("Lumi",1.0)),
  weightEvents(iConfig.getUntrackedParameter<bool>("weightEvents",false)),
  isoCut(iConfig.getUntrackedParameter<double>("isoCut",0.4)),
  earlyIsoCut(iConfig.getUntrackedParameter<double>("earlyIsoCut",0.4)),
  sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
  leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
  subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
  deltaRmu(iConfig.getUntrackedParameter<double>("deltaRmu",0.3)),
  deltaRe(iConfig.getUntrackedParameter<double>("deltaRe",0.3)),
  _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",5)),
  _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",3)),
  tightIdsOnly(iConfig.getUntrackedParameter<bool>("tightIdsOnly",false)),
  looseIdsOnly(iConfig.getUntrackedParameter<bool>("looseIdsOnly",true)),
  reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
  interactiveRun(iConfig.getUntrackedParameter<bool>("interactiveRun",false)),
  PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Legacy53X")),
  doVarDump(iConfig.getUntrackedParameter<bool>("doVarDump",false)),
  doBlinding(iConfig.getUntrackedParameter<bool>("doBlinding",false)),
  doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",false)),
  bStudyResolution(iConfig.getUntrackedParameter<bool>("bStudyResolution",false)),
  bStudyDiLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyDiLeptonResolution",false)),
  bStudyFourLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyFourLeptonResolution",false))
{
  
  if(!isMC){reweightForPU = false;}
  
  // book histograms:
  //ex: histContainer_["photons"]=fs->make<TH1F>("photons", "photon multiplicity",   10, 0,  10);
  
  histContainer_["extraParticles"]=fs->make<TH1F>("extraParticles", "Extra Leptons after Iso; N Leptons; N Events", 15, 0, 15);

  histContainer_["NEVENTS"]=fs->make<TH1F>("nEvents","nEvents in Sample",2,0,2);
  histContainer_["NVTX"]=fs->make<TH1F>("nVtx","Number of Vertices",61,-0.5,60.5);
  histContainer_["NVTX_RW"]=fs->make<TH1F>("nVtx_ReWeighted","Number of Vertices",61,-0.5,60.5);
  histContainer_["NINTERACT"]=fs->make<TH1F>("nInteractions","Number of True Interactions",61,-0.5,60.5);
  histContainer_["NINTERACT_RW"]=fs->make<TH1F>("nInteraction_ReWeighted","Number of True Interactions",61,-0.5,60.5);
  histContainer_["PUweights"]=fs->make<TH1F>("PUweights","PUweights",240,0,80);

  histContainer_["minMass2l"]=fs->make<TH1F>("minMass2l","Minimum Mass of 2 Leptons",1000,0,100);
  histContainer_["minMass2l_SS_SF"]=fs->make<TH1F>("minMass2l_SS_SF","Minimum Mass of 2 Leptons",1000,0,100);
  histContainer_["minMass2l_OS_SF"]=fs->make<TH1F>("minMass2l_OS_SF","Minimum Mass of 2 Leptons",1000,0,100);
  histContainer_["minMass2l_SS_OF"]=fs->make<TH1F>("minMass2l_SS_OF","Minimum Mass of 2 Leptons",1000,0,100);
  histContainer_["minMass2l_OS_OF"]=fs->make<TH1F>("minMass2l_OS_OF","Minimum Mass of 2 Leptons",1000,0,100);


  passedEventsTree_All = new TTree("passedEvents","passedEvents");
  sigEff_4 = new HZZ4LSigEff("SigEff_4");
  sigEff_6 = new HZZ4LSigEff("SigEff_6");
  sigEff_8 = new HZZ4LSigEff("SigEff_8");
  sigEff_9 = new HZZ4LSigEff("SigEff_9");
  sigEff_10 = new HZZ4LSigEff("SigEff_10");
  sigEff_11 = new HZZ4LSigEff("SigEff_11");
  sigEff_12 = new HZZ4LSigEff("SigEff_12");
  

  
  //MUST HAVE THESE FOR HIGGS CANDIDATE SELECTION
  RecoFourMuEvent = false;
  RecoFourEEvent = false;
  RecoTwoETwoMuEvent = false;
  RecoTwoMuTwoEEvent = false;
  eventPassedPtAndIsoCuts = false;
  foundHiggsCandidate = false;
  twoLooseIsoLeptons = false;
  passedM4lCut = false;
  isIsolated = false;
  passedSIP3D = false;
  passedPtCuts = false;
  twoLep_ID = false;
  fourLep_Cleaned = false;
  mZ1 = 0;
  mZ2 = 0;
  m4l = 0;
  passedFullSelection = false;
  passedZ4lSelection = false;
  passedQCDcut = false;

  eta4l = 0;  phi4l = 0;
  mvaL1 = 0; mvaL2 = 0; mvaL3 = 0; mvaL4 = 0;
  missingHitsL1 = -1; missingHitsL2 = -1; missingHitsL3 = -1; missingHitsL4 = -1;
  ecalSeedL1 = false; ecalSeedL2 = false; ecalSeedL3 = false; ecalSeedL4 = false;
  pTL1 = 0;  pTL2 = 0;  pTL3 = 0;  pTL4 = 0;
  pXL1 = 0;  pXL2 = 0;  pXL3 = 0;  pXL4 = 0;
  pYL1 = 0;  pYL2 = 0;  pYL3 = 0;  pYL4 = 0;
  pZL1 = 0;  pZL2 = 0;  pZL3 = 0;  pZL4 = 0;
  EL1 = 0;   EL2 = 0;   EL3 = 0;   EL4 = 0;
  pTL1FSR = 0;  pTL2FSR = 0;  pTL3FSR = 0;  pTL4FSR = 0;
  pXL1FSR = 0;  pXL2FSR = 0;  pXL3FSR = 0;  pXL4FSR = 0;
  pYL1FSR = 0;  pYL2FSR = 0;  pYL3FSR = 0;  pYL4FSR = 0;
  pZL1FSR = 0;  pZL2FSR = 0;  pZL3FSR = 0;  pZL4FSR = 0;
  EL1FSR = 0;   EL2FSR = 0;   EL3FSR = 0;   EL4FSR = 0;
  EZ1 = 0;   EZ2 = 0;  
  pTZ1 = 0;  pTZ2 = 0;  
  pXZ1 = 0;  pXZ2 = 0;
  pYZ1 = 0;  pYZ2 = 0;
  pZZ1 = 0;  pZZ2 = 0;
  chargeL1 = -999;  chargeL2 = -999;  chargeL3 = -999;  chargeL4 = -999;
  etaL1 = 0;  etaL2 = 0;  etaL3 = 0;  etaL4 = 0;
  phiL1 = 0;  phiL2 = 0;  phiL3 = 0;  phiL4 = 0;
  m4lErrUCSD = 0; m4lErrUCSDCorr = 0; m4lErrUF = 0; m4lErrUFCorr = 0; m4lErrUFADCorr = 0;
  isoNHL1 = 0; isoNHL2 = 0; isoNHL3 = 0; isoNHL4 = 0;
  isoCHL1 = 0; isoCHL2 = 0; isoCHL3 = 0; isoCHL4 = 0;
  isoPhotL1 = 0; isoPhotL2 = 0; isoPhotL3 = 0; isoPhotL4 = 0;
  FSR_Z1 = false; FSR_Z2 = false;
  FSRPhot1_Pt = -999, FSRPhot2_Pt = -999;
  FSRPhot1_eta = -999, FSRPhot2_eta = -999;
  FSRPhot1_phi = -999, FSRPhot2_phi = -999;
  FSRPhot1_Px = -999; FSRPhot2_Px = -999;
  FSRPhot1_Py = -999; FSRPhot2_Py = -999;
  FSRPhot1_Pz = -999; FSRPhot2_Pz = -999;
  FSRPhot1_E = -999; FSRPhot2_E = -999;

  //Calculated Angles
  cosTheta1 = 0; cosTheta2 = 0; Phi = 0; cosThetaStar = 0;
  phiStar1 = 0;  phiStar2 = 0; phiStar12 = 0; Phi1 = 0; Phi2 = 0;

 
  //Total event counter
  nEventsTotal = 0; nEventsSkimmed = 0; nEvAfterSkim = 0; nEvPassedHlt = 0;
  nEvPassedPtCut = 0; nEvWith2IsoLep = 0; nEvWith4lep = 0; nEvBeforeZCuts = 0;
  nEvAfterM4lCut = 0; nEvAfterZ1Cut = 0; nEvAfterZ2Cut = 0; nEvAfterIso = 0;

  nEvAfterZZCut = 0; nEvAfterZZCut_4e = 0; nEvAfterZZCut_4mu = 0; nEvAfterZZCut_2e2mu = 0;
  nEvAfterMelaCut = 0; nEvAfterMelaCut_4e = 0; nEvAfterMelaCut_4mu = 0; nEvAfterMelaCut_2e2mu = 0;

  nEvAfterPsMelaCut_4e = 0; nEvAfterPsMelaCut_4mu = 0; nEvAfterPsMelaCut_2e2mu = 0; nEvAfterPsMelaCut = 0;
  nEvAfterGrMelaCut_4e = 0; nEvAfterGrMelaCut_4mu = 0; nEvAfterGrMelaCut_2e2mu = 0; nEvAfterGrMelaCut = 0;

  nEvAfterIso_comb = 0; nEvAfterSIP = 0;
  
  nEvAfterId = 0; nEvAfterId_4mu = 0; nEvAfterId_4e = 0; nEvAfterId_2e2mu = 0;
  nEvAfterCleaning = 0; nEvAfterCleaning_4mu = 0; nEvAfterCleaning_4e = 0; nEvAfterCleaning_2e2mu = 0;
  nEvPassedPtCut_4mu = 0; nEvPassedPtCut_4e = 0; nEvPassedPtCut_2e2mu = 0;
  nEvPassedPtCut2_4mu = 0; nEvPassedPtCut2_4e = 0; nEvPassedPtCut2_2e2mu = 0; nEvPassedPtCut2 = 0;
  nEvAfterZ4lCut_4mu = 0; nEvAfterZ4lCut_4e = 0; nEvAfterZ4lCut_2e2mu = 0; nEvAfterZ4lCut = 0;

  nEvAfterZ1Formed = 0; nEvAfterZ1Formed_2mu = 0; nEvAfterZ1Formed_2e = 0;
  nEvAfterZ2Formed = 0; nEvAfterZ2Formed_4mu = 0; nEvAfterZ2Formed_4e = 0; nEvAfterZ2Formed_2e2mu = 0;
  nEv4GoodLep = 0; nEv4GoodLep_4mu = 0; nEv4GoodLep_4e = 0; nEv4GoodLep_2e2mu = 0;
  nEvAfterVBFJet1 = 0; nEvAfterVBFJet1_4mu = 0; nEvAfterVBFJet1_4e = 0; nEvAfterVBFJet1_2e2mu = 0;
  nEvAfterVBFJet2 = 0; nEvAfterVBFJet2_4mu = 0; nEvAfterVBFJet2_4e = 0; nEvAfterVBFJet2_2e2mu = 0;
  nEvAfterVBFJetCuts = 0; nEvAfterVBFJetCuts_4mu = 0; nEvAfterVBFJetCuts_4e = 0; nEvAfterVBFJetCuts_2e2mu = 0;

  counterMuZ1ph = 0; counterElZ1ph = 0;
  counterMuZ2ph = 0; counterElZ2ph = 0;

  nEvWith1FSRZ = 0; nEvWith1FSRZ_4e = 0; nEvWith1FSRZ_4mu = 0; nEvWith1FSRZ_2e2mu = 0;
  nEvWith2FSRZ = 0; nEvWith2FSRZ_4e = 0; nEvWith2FSRZ_4mu = 0; nEvWith2FSRZ_2e2mu = 0;

  nEvFSRPtLt4 = 0; nEvFSRPtGt4dR0p5MatchFsrISO = 0; nEvFSRPtGt4dR0p07Match = 0;
  nEvFSRPtLt4Zmm = 0; nEvFSRPtGt4dR0p5MatchFsrISOZmm = 0; nEvFSRPtGt4dR0p07MatchZmm = 0;
  nEvFSRPtLt4Zee = 0; nEvFSRPtGt4dR0p5MatchFsrISOZee = 0; nEvFSRPtGt4dR0p07MatchZee = 0;

  nEvBeforeZCuts_4mu = 0;
  nEvAfterM4lCut_4mu = 0;
  nEvAfterZ1Cut_2mu = 0;
  nEvAfterZ2Cut_4mu = 0;
  nEvAfterIso_4mu = 0;
  nEvAfterSIP_4mu = 0;
  
  
  nEvBeforeZCuts_4e = 0;
  nEvAfterM4lCut_4e = 0;
  nEvAfterZ1Cut_2e = 0;
  nEvAfterZ2Cut_4e = 0;
  nEvAfterIso_4e = 0;
  nEvAfterSIP_4e = 0;
  
  
  nEvBeforeZCuts_2e2mu = 0;
  nEvAfterM4lCut_2e2mu = 0;
  nEvAfterZ1Cut_2e2mu = 0;
  nEvAfterZ2Cut_2e2mu = 0;
  nEvAfterIso_2e2mu = 0;
  nEvAfterSIP_2e2mu = 0;


  nEvAfterZZCut = 0;
  nEvAfterZZCut_4e = 0;
  nEvAfterZZCut_4mu = 0;
  nEvAfterZZCut_2e2mu = 0;

  nEvAfterMelaCut = 0;
  nEvAfterMelaCut_4e = 0;
  nEvAfterMelaCut_4mu = 0;
  nEvAfterMelaCut_2e2mu = 0;
  nEvAfterPsMelaCut_4e = 0;
  nEvAfterPsMelaCut_4mu = 0;
  nEvAfterPsMelaCut_2e2mu = 0;
  nEvAfterPsMelaCut = 0;
  nEvAfterGrMelaCut_4e = 0;
  nEvAfterGrMelaCut_4mu = 0;
  nEvAfterGrMelaCut_2e2mu = 0;
  nEvAfterGrMelaCut = 0;

  gen4mu = 0;gen4e = 0;gen2e2mu = 0;
  gen4muPseudo = 0;gen4ePseudo = 0; gen2e2muPseudo = 0; 

  
  //Event weight... =1 if nothing is specified in cfg
  eventWeight = 1.0;
  MC_weight = 1.0;
  if(weightEvents && isMC){scaleWeight = CrossSection*FilterEff;}
  else{scaleWeight = 1.0;}
  
  Rho = 0;
  muonRho = 0;
  elecRho = 0;
  
  leastIso = 0;
  highestSip = 0;
  

  //if(interactiveRun)
  //{ lumiWeight = new edm::LumiReWeighting( "hists/mcFlat10_Fall11.root", "hists/Data_PU_Fall11.root", "pileup_mc", "pileup" );}
  //else
  //{ lumiWeight = new edm::LumiReWeighting( "UFHZZAnalysis8TeV/UFHZZ4LAna/hists/mcFlat10_Fall11.root", 
  //                                         "UFHZZAnalysis8TeV/UFHZZ4LAna/hists/Data_PU_Fall11.root",
  //					     "pileup_mc", "pileup" ); } 

  npv = -1;
  BX = 0;
  PU_weight = 1;
  PT_weight = 1;


  isoEff = new HZZ4LIsoEff("");

  Zto4LAnaOP = new HZZ4LZto4LAna("_offPeak");
  
  //mela = new Mela(false,8);
  //MEKDwithPDFs = new MEKD(8,"CTEQ6L");
  //MEKDnoPDFs = new MEKD(8,"");
  //MEKDnoPDFs_noFSR = new MEKD(8,"");
  //MEMsnoPDFs_noFSR = new MEMs(8); //removed for miniAOD
  //MEMsnoPDFs = new MEMs(8);//removed for miniAOD
  

  muonDump = new HZZ4LMuonTree("muonDumpTree",fs);
  electronDump = new HZZ4LElectronTree("electronDumpTree",fs);
  photonDump = new HZZ4LPhotonTree("photonDumpTree",fs);
  finalLepDump = new HZZ4LFinalLepTree("finalLepDumpTree",fs);
  jetDump = new HZZ4LJetTree("jetDumpTree",fs);
 
  counter_4mu = 0;    counter_4e = 0;    counter_2e2mu = 0;    counter_2mu2e = 0;
  LMcounter_4mu = 0;  LMcounter_4e = 0;  LMcounter_2e2mu = 0;  LMcounter_2mu2e = 0;
  Z4lcounter_4mu = 0; Z4lcounter_4e = 0; Z4lcounter_2e2mu = 0; Z4lcounter_2mu2e = 0;

  if(bStudyResolution) 
  {
    PerLepReso = new HZZ4LPerLepResolution();
    if(bStudyDiLeptonResolution) {  DiLepReso = new HZZ4LDiLepResolution(); }
    if(bStudyFourLeptonResolution)  { FourLepReso = new HZZ4LResolution(); }
  }

 
}


UFHZZ4LAna::~UFHZZ4LAna()
{
  
  //destructor --- don't do anything here
  
}



// ------------ method called for each event  ------------
void
UFHZZ4LAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace pat;
  //using namespace MEMNames; // removed for miniAOD

  // ======= Get Collections ======= //

  // vertex collection
  edm::Handle<reco::VertexCollection> vertex;
  iEvent.getByLabel(vertexSrc_,vertex);
  const reco::Vertex *PV = 0;
  if(!vertex->empty() && vertex->size() > 0) PV = &(vertex->at(0));
  
  // PU collection
  if(isMC && reweightForPU)
  {
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
      
    std::vector<PileupSummaryInfo>::const_iterator PVI;
     
    npv = -1;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) 
    {
      BX = PVI->getBunchCrossing();
      if(BX == 0) 
      { 
        //npv = PVI->getPU_NumInteractions();//old
        npv = PVI->getTrueNumInteractions(); 
        continue;
      }
    }
  }
  histContainer_["NINTERACT"]->Fill(npv);
  //if(isMC && reweightForPU){PU_weight = lumiWeight->weight( npv );}
  if(isMC && reweightForPU)
  {PU_weight = pileUp.getPUWeight(npv,PUVersion);}
  else
  { PU_weight = 1.0;}
  //cout << "PU WEIGHT: " <<  PU_weight << endl;
  histContainer_["NINTERACT_RW"]->Fill(npv,PU_weight);

  // photon collection 
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByLabel(photonSrc_,photons);
  
  // electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);
  
  // muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);
  
  // met collection 
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_,mets);
    
  // beamspot collection
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByLabel("offlineBeamSpot",beamspot);
  
  // Rho Correction
  edm::Handle<double> eventRhoMu;
  iEvent.getByLabel(muRhoSrc_,eventRhoMu);
  muonRho = *eventRhoMu;

  edm::Handle<double> eventRhoE;
  iEvent.getByLabel(elRhoSrc_,eventRhoE);
  elecRho = *eventRhoE;

  // Particle Flow Cands
  edm::Handle<edm::View<pat::PackedCandidate> > pfCands; // modified for miniAOD
  iEvent.getByLabel("packedPFCandidates",pfCands); // modified for miniAOD

  edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;  //remove for miniAOD
  iEvent.getByLabel("boostedFsrPhotons",photonsForFsr); //remove for miniAOD
  
  // GEN collection
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("prunedGenParticles", genParticles); // modified for miniAOD

  //VBF 
  //edm::Handle<edm::View<pat::Jet> > correctedJets; //removed for miniAOD
  //iEvent.getByLabel(correctedJetSrc_,correctedJets); //removed for miniAOD

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);
  
  //edm::Handle<edm::ValueMap<float> > puJetIdMva; // removed for miniAOD
  //iEvent.getByLabel("puJetMva","full53xDiscriminant",puJetIdMva); // removed for miniAOD
  //iEvent.getByLabel("puJetMva","fullDiscriminant",puJetIdMva);

  //edm::Handle<edm::ValueMap<int> > puJetIdFlag; // removed for miniAOD
  //iEvent.getByLabel("puJetMva","full53xId",puJetIdFlag);  // removed for miniAOD
  //iEvent.getByLabel("puJetMva","fullId",puJetIdFlag);  
  

  // ============= Set Variables ============= //
  
  //MUST HAVE THESE FOR HIGGS CANDIDATE SELECTION 
  RecoFourMuEvent = false;
  RecoFourEEvent = false;
  RecoTwoETwoMuEvent = false;
  RecoTwoMuTwoEEvent = false;
  eventPassedPtAndIsoCuts = false;
  foundHiggsCandidate = false;
  twoLooseIsoLeptons = false;
  passedM4lCut = false;
  isIsolated = false;
  passedSIP3D = false;
  passedPtCuts = false;
  twoLep_ID = false;
  fourLep_Cleaned = false;
  passedFullSelection = false;
  passedZ4lSelection = false;
  passedQCDcut = false;
  mZ1 = 0;
  mZ2 = 0;
  m4l = 0;
  m4lNoFSR = 0;
  FSR_Z1 = false; FSR_Z2 = false;
  FSRPhot1_Pt = -999; FSRPhot2_Pt = -999;
  FSRPhot1_eta = -999; FSRPhot2_eta = -999;
  FSRPhot1_phi = -999; FSRPhot2_phi = -999;
  FSRPhot1_Px = -999; FSRPhot2_Px = -999;
  FSRPhot1_Py = -999; FSRPhot2_Py = -999;
  FSRPhot1_Pz = -999; FSRPhot2_Pz = -999;
  FSRPhot1_E = -999; FSRPhot2_E = -999;

  VBFJet1 = false; VBFJet2 = false;
  pTVBFJet1 = -999; pTVBFJet2 = -999;
  etaVBFJet1 = -999; etaVBFJet2 = -999;
  phiVBFJet1 = -999; phiVBFJet2 = -999;
  massVBFJet1 = -999; massVBFJet2 = -999;
  VBFDiJetMass = -999; VBFDeltaEta = -999;
  FisherDiscrim = -999;

  GENidL1_S1 = -999; GENidL2_S1 = -999; GENidL3_S1 = -999; GENidL4_S1 = -999;
  GENpTL1_S1 = -999; GENpTL2_S1 = -999; GENpTL3_S1 = -999; GENpTL4_S1 = -999;
  GENpXL1_S1 = -999; GENpXL2_S1 = -999; GENpXL3_S1 = -999; GENpXL4_S1 = -999;
  GENpYL1_S1 = -999; GENpYL2_S1 = -999; GENpYL3_S1 = -999; GENpYL4_S1 = -999;
  GENpZL1_S1 = -999; GENpZL2_S1 = -999; GENpZL3_S1 = -999; GENpZL4_S1 = -999;
  GENEL1_S1 = -999; GENEL2_S1 = -999; GENEL3_S1 = -999; GENEL4_S1 = -999;
  GENEZ1 = -999; GENEZ2 = -999;
  GENpTZ1 = -999; GENpTZ2 = -999;
  GENpXZ1 = -999; GENpXZ2 = -999;
  GENpYZ1 = -999; GENpYZ2 = -999;
  GENpZZ1 = -999; GENpZZ2 = -999;
  GENchargeL1_S1 = -999; GENchargeL2_S1 = -999; GENchargeL3_S1 = -999; GENchargeL4_S1 = -999;
  GENetaL1_S1 = -999; GENetaL2_S1 = -999; GENetaL3_S1 = -999; GENetaL4_S1 = -999;
  GENphiL1_S1 = -999; GENphiL2_S1 = -999; GENphiL3_S1 = -999; GENphiL4_S1 = -999;

  GENidL1_S3 = -999; GENidL2_S3 = -999; GENidL3_S3 = -999; GENidL4_S3 = -999;
  GENpTL1_S3 = -999; GENpTL2_S3 = -999; GENpTL3_S3 = -999; GENpTL4_S3 = -999;
  GENpXL1_S3 = -999; GENpXL2_S3 = -999; GENpXL3_S3 = -999; GENpXL4_S3 = -999;
  GENpYL1_S3 = -999; GENpYL2_S3 = -999; GENpYL3_S3 = -999; GENpYL4_S3 = -999;
  GENpZL1_S3 = -999; GENpZL2_S3 = -999; GENpZL3_S3 = -999; GENpZL4_S3 = -999;
  GENEL1_S3 = -999; GENEL2_S3 = -999; GENEL3_S3 = -999; GENEL4_S3 = -999;
  GENchargeL1_S3 = -999; GENchargeL2_S3 = -999; GENchargeL3_S3 = -999; GENchargeL4_S3 = -999;
  GENetaL1_S3 = -999; GENetaL2_S3 = -999; GENetaL3_S3 = -999; GENetaL4_S3 = -999;
  GENphiL1_S3 = -999; GENphiL2_S3 = -999; GENphiL3_S3 = -999; GENphiL4_S3 = -999;

  GENMH = -999; GENMZ1 = -999; GENMZ2 = -999; GENM4L = -999;
  
  idL1_GENMatched = -999; idL2_GENMatched = -999; idL3_GENMatched = -999; idL4_GENMatched = -999;
  pTL1_GENMatched = -999; pTL2_GENMatched = -999; pTL3_GENMatched = -999; pTL4_GENMatched = -999;
  pXL1_GENMatched = -999; pXL2_GENMatched = -999; pXL3_GENMatched = -999; pXL4_GENMatched = -999;
  pYL1_GENMatched = -999; pYL2_GENMatched = -999; pYL3_GENMatched = -999; pYL4_GENMatched = -999;
  pZL1_GENMatched = -999; pZL2_GENMatched = -999; pZL3_GENMatched = -999; pZL4_GENMatched = -999;
  EL1_GENMatched = -999; EL2_GENMatched = -999; EL3_GENMatched = -999; EL4_GENMatched = -999;
  EZ1_GENMatched = -999; EZ2_GENMatched = -999;
  pTZ1_GENMatched = -999; pTZ2_GENMatched = -999;
  pXZ1_GENMatched = -999; pXZ2_GENMatched = -999;
  pYZ1_GENMatched = -999; pYZ2_GENMatched = -999;
  pZZ1_GENMatched = -999; pZZ2_GENMatched = -999;
  chargeL1_GENMatched = -999; chargeL2_GENMatched = -999; chargeL3_GENMatched = -999; chargeL4_GENMatched = -999;
  etaL1_GENMatched = -999; etaL2_GENMatched = -999; etaL3_GENMatched = -999; etaL4_GENMatched = -999;
  phiL1_GENMatched = -999; phiL2_GENMatched = -999; phiL3_GENMatched = -999; phiL4_GENMatched = -999;
  m4l_GENMatched = -999;

  melaLD = -1; mela_Sig = -1; mela_Bkg = -1;
  melaLD_Pt = -1; mela_Sig_Pt = -1; mela_Bkg_Pt = -1;
  melaLD_Y = -1; mela_Sig_Y = -1; mela_Bkg_Y = -1;
  melaLD_PtY = -1; mela_Sig_PtY = -1; mela_Bkg_PtY = -1;

  pdfSigM4l = -1; pdfSigM4l_noFSR = -1; pdfBkgM4l = -1; pdfBkgM4l_noFSR = -1;
  pdfSigM4l_ScaleUp = -1; pdfSigM4l_ScaleUp_noFSR = -1; pdfBkgM4l_ScaleUp = -1; pdfBkgM4l_ScaleUp_noFSR = -1;
  pdfSigM4l_ScaleDown = -1; pdfSigM4l_ScaleDown_noFSR = -1; pdfBkgM4l_ScaleDown = -1; pdfBkgM4l_ScaleDown_noFSR = -1;
  pdfSigM4l_ResUp = -1; pdfSigM4l_ResUp_noFSR = -1; pdfBkgM4l_ResUp = -1; pdfBkgM4l_ResUp_noFSR = -1;
  pdfSigM4l_ResDown = -1; pdfSigM4l_ResDown_noFSR = -1; pdfBkgM4l_ResDown = -1; pdfBkgM4l_ResDown_noFSR = -1;

  pseudoMelaLD = -1; pseudoMela_SM = -1; pseudoMela_PS = -1;
  spin1EvMelaLD = -1; spin1EvMela_SM = -1; spin1EvMela_S1E = -1;
  spin1OddMelaLD = -1; spin1OddMela_SM = -1; spin1OddMela_S1O = -1;
  spin2MelaLD = -1; spin2Mela_SM = -1; spin2Mela_S2 = -1;

  p0plus_melaNorm = -999;
  p0plus_mela = -999;   
  p0minus_mela = -999;  
  p0plus_VAJHU = -999;  
  p0minus_VAJHU = -999; 
  p0plus_VAMCFM = -999; 
  p1_mela = -999;  
  p1_VAJHU = -999; 
  p2_mela = -999; 
  p2_VAJHU = -999; 
  mela_bkg_analytic = -999;  
  mela_bkg_VAMCFM = -999;
  mela_bkg_ggzz_VAMCFM = -999;
  mela_bkg_VAMCFMNorm = -999; 
  mela_p0_pt = -999; 
  mela_p0_y = -999; 
  mela_bkg_pt = -999;
  mela_bkg_y = -999;
  p0plus_m4l = -999;
  bkg_m4l = -999;

  interferenceWeight = 1;
  JHUKD_H_qqZZ_noPDF = -999; JHU_ME_H = -999; MCFM_ME_qqZZ = -999;
  JHUKD_H_h0M_noPDF = -999; JHU_ME_h0M = -999;
  JHUKD_H_h0P_noPDF = -999; JHU_ME_h0P = -999;
  JHUKD_H_h1M_noPDF = -999; JHU_ME_h1M = -999;
  JHUKD_H_h1P_noPDF = -999; JHU_ME_h1P = -999;
  JHUKD_H_ggh2P_noPDF = -999; JHU_ME_ggh2P = -999;
  JHUKD_H_qqh2P_noPDF = -999; JHU_ME_qqh2P = -999;

  JHUKD_H_h2hP_noPDF = -999; JHU_ME_h2hP = -999;
  JHUKD_H_h2hM_noPDF = -999; JHU_ME_h2hM = -999;
  JHUKD_H_h2bP_noPDF = -999; JHU_ME_h2bP = -999;
  JHUKD_H_h2P_prodInd_noPDF = -999; JHU_ME_h2P_prodInd = -999;
  JHUKD_H_h1P_prodInd_noPDF = -999; JHU_ME_h1P_prodInd = -999;
  JHUKD_H_h1M_prodInd_noPDF = -999; JHU_ME_h1M_prodInd = -999;
  
  MEKD_noPDF_noFSR = -999; MEKD_ME_H_noPDF_noFSR = -999; MEKD_ME_ZZ_noPDF_noFSR = -999;
  MEKD_h0M_ZZ_noPDF_noFSR = -999; MEKD_ME_h0M_noPDF_noFSR = -999;
  MEKD_h0P_ZZ_noPDF_noFSR = -999; MEKD_ME_h0P_noPDF_noFSR = -999;
  MEKD_h1M_ZZ_noPDF_noFSR = -999; MEKD_ME_h1M_noPDF_noFSR = -999;
  MEKD_h1P_ZZ_noPDF_noFSR = -999; MEKD_ME_h1P_noPDF_noFSR = -999;
  MEKD_qqh2P_ZZ_noPDF_noFSR = -999;MEKD_ME_qqh2P_noPDF_noFSR = -999;
  MEKD_ggh2P_ZZ_noPDF_noFSR = -999;MEKD_ME_ggh2P_noPDF_noFSR = -999;

  MEKD_h2hP_ZZ_noPDF_noFSR = -999; MEKD_ME_h2hP_noPDF_noFSR = -999;
  MEKD_h2hM_ZZ_noPDF_noFSR = -999; MEKD_ME_h2hM_noPDF_noFSR = -999;
  MEKD_h2bP_ZZ_noPDF_noFSR = -999; MEKD_ME_h2bP_noPDF_noFSR = -999;
  MEKD_h2P_ZZ_prodInd_noPDF_noFSR = -999; MEKD_ME_h2P_prodInd_noPDF_noFSR = -999;
  MEKD_h1P_ZZ_prodInd_noPDF_noFSR = -999; MEKD_ME_h1P_prodInd_noPDF_noFSR = -999;
  MEKD_h1M_ZZ_prodInd_noPDF_noFSR = -999; MEKD_ME_h1M_prodInd_noPDF_noFSR = -999;

  JHUKD_H_qqZZ_noPDF_noFSR = -999;  JHU_ME_H_noPDF_noFSR = -999;   JHU_ME_ZZ_noPDF_noFSR = -999;
  JHUKD_h0M_ZZ_noPDF_noFSR = -999;  JHU_ME_h0M_noPDF_noFSR = -999;
  JHUKD_h0P_ZZ_noPDF_noFSR = -999;  JHU_ME_h0P_noPDF_noFSR = -999;
  JHUKD_h1M_ZZ_noPDF_noFSR = -999;  JHU_ME_h1M_noPDF_noFSR = -999;
  JHUKD_h1P_ZZ_noPDF_noFSR = -999;  JHU_ME_h1P_noPDF_noFSR = -999;
  JHUKD_qqh2P_ZZ_noPDF_noFSR = -999;JHU_ME_qqh2P_noPDF_noFSR = -999;
  JHUKD_ggh2P_ZZ_noPDF_noFSR = -999;JHU_ME_ggh2P_noPDF_noFSR = -999;
    
  MEKD_wPDF = -999; MEKD_ME_H_wPDF = -999; MEKD_ME_ZZ_wPDF = -999;
  MEKD_noPDF = -999; MEKD_ME_H_noPDF = -999; MEKD_ME_ZZ_noPDF = -999;
  MEKD_wPDF_noFSR = -999; MEKD_ME_H_wPDF_noFSR = -999; MEKD_ME_ZZ_wPDF_noFSR = -999;
  

  minDeltR = 1000;
  leastIso = 0;
  highestSip = 0;
  theta12 = 1000;
  theta13 = 1000;
  theta14 = 1000;
  theta12_deg = 1000;
  theta13_deg = 1000;
  theta14_deg = 1000;
  thetaPhoton = 1000;
  thetaPhoton_deg = 1000;
  thetaPhotonZ = 1000;
  thetaPhotonZ_deg = 1000;
  thetaZ2 = 1000;
  etaZ2 = 1000;

  minMass2Lep = 0;
  maxMass2Lep = 0;
  
  massErrorUCSD = 0; massErrorUCSDCorr = 0; massErrorUF = 0; massErrorUFCorr = 0; massErrorUFADCorr = 0;
  eventType = "";

  minM3l = 1000;
  Z4lmaxP = 0;
  minDeltR = 1000;
  finalState = -1;
  m3l_soft = 0;
  // ====================== Do Analysis ======================== //

  //Event weighting
  if( isSignal ){ eventWeight = PU_weight*PT_weight; }
  else if ( !isMC ){eventWeight = 1.0; }
  else{ eventWeight = PU_weight;}

  if(isMC && isSignal)
  {
    nEventsTotal += eventWeight;
    sigEff_4->advanceSigDenCounters(genParticles,eventType,eventWeight);
    sigEff_6->advanceSigDenCounters(genParticles,eventType,eventWeight);
    sigEff_8->advanceSigDenCounters(genParticles,eventType,eventWeight);
    sigEff_9->advanceSigDenCounters(genParticles,eventType,eventWeight);	 
    sigEff_10->advanceSigDenCounters(genParticles,eventType,eventWeight);	
    sigEff_11->advanceSigDenCounters(genParticles,eventType,eventWeight);
    sigEff_12->advanceSigDenCounters(genParticles,eventType,eventWeight);
  }



  //Check for duplicate events in data
  notDuplicateEvent = true;
  if(!isMC)
  { 
    ULong64_t runId = iEvent.id().run();
    ULong64_t eventId = iEvent.id().event();
    ULong64_t lumiId = iEvent.id().luminosityBlock();    
    for (unsigned int n = 0; n < runVec.size(); n++) 
    {
      if(runId == runVec[n] && lumiId == lumiVec[n] && eventId == eventVec[n]){notDuplicateEvent = false;}
    } 
    if (notDuplicateEvent) 
    {
      runVec.push_back(runId);
      lumiVec.push_back(lumiId);
      eventVec.push_back(eventId);
    }
  }


  if( notDuplicateEvent && !vertex->empty())
  {
    
    nEvPassedHlt += eventWeight;

    //N Vertex 
    nVtx = vertex->size();
    histContainer_["NVTX"]->Fill(nVtx);
    histContainer_["NVTX_RW"]->Fill(nVtx,PU_weight);
    
    vector<pat::Muon> AllMuons;
    vector<pat::Electron> AllElectrons;
  
    AllMuons = helper.goodLooseMuons2012(muons,_muPtCut);
    AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut);

    //bool skimStep1 = helper.passedPtandEarlyIso(AllMuons, AllElectrons, leadingPtCut, subleadingPtCut, 10000, muonRho, elecRho);
    //bool skimStep2 = helper.passedM2lCut(AllMuons, AllElectrons, 40, 120);

    // no more skim
    //    if(skimStep1 && skimStep2)
      
    vector<pat::Muon> recoMuons;
    vector<pat::Electron> recoElectrons;
    helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);

    recoMuons = helper.goodMuons2012_noIso(AllMuons,_muPtCut,PV);
    recoElectrons = helper.goodElectrons2012_noIso(AllElectrons,_elecPtCut,elecID,PV,iEvent);

    if( (recoMuons.size() + recoElectrons.size()) >= 2 ){twoLep_ID = true;}

    if(doVarDump)
    {
      muonDump->fillMuonDumpTree(AllMuons,iEvent,muonRho,PV);
      electronDump->fillElectronDumpTree(AllElectrons,iEvent,elecRho,PV,elecID);
    }	      

    if(twoLep_ID)
    {
      // Mass Resolution Study
      if(bStudyResolution)
      {
        vector< pat::Muon > recoIsoMuons;
        vector< pat::Electron > recoIsoElectrons;
	    
        recoIsoMuons = helper.goodMuons2012_Iso(AllMuons,_muPtCut, muonRho, isoCut, PV);
        recoIsoElectrons = helper.goodElectrons2012_Iso(AllElectrons,_elecPtCut, elecRho, isoCut, elecID,PV);
	    
        PerLepReso->fillHistograms(hContainer_, hContainer2D_, hContainer3D_, recoIsoElectrons, recoIsoMuons, eventWeight, !isMC);
        if(bStudyDiLeptonResolution) { DiLepReso->fillHistograms(hContainer_, hContainer2D_, recoIsoElectrons, recoIsoMuons, eventWeight, !isMC); }
      }

	
      nEvAfterId += eventWeight;
      if(isSignal) sigEff_4->advanceSigNumCounters_ID(eventType,eventWeight);

      //Plot Min Mass(2l)
      bool sameFlavorMin = false, sameSignMin = false;
      minMass2Lep = helper.minMass2l(recoMuons,recoElectrons, sameFlavorMin, sameSignMin);
      histContainer_["minMass2l"]->Fill(minMass2Lep,eventWeight);    
      if(sameFlavorMin)
      {
        if(sameSignMin){histContainer_["minMass2l_SS_SF"]->Fill(minMass2Lep,eventWeight);}
        else{histContainer_["minMass2l_OS_SF"]->Fill(minMass2Lep,eventWeight);}
      }
      else
      {
        if(sameSignMin){histContainer_["minMass2l_SS_OF"]->Fill(minMass2Lep,eventWeight);}
        else{histContainer_["minMass2l_OS_OF"]->Fill(minMass2Lep,eventWeight);}
      }
	
      nPhotons = 0;
      std::vector<pat::PFParticle> fsrPhotons; std::vector<double> deltaRVec;
      if(doFsrRecovery)
      {
        for(edm::View<pat::PFParticle>::const_iterator phot=photonsForFsr->begin(); phot!=photonsForFsr->end(); ++phot)
        {
          bool matched = false; double chosenDeltaRPh = 999;
          for(int i=0; i<(int)recoElectrons.size(); i++)
          {
            double tmpDeltaREPh = deltaR(recoElectrons[i].eta(), recoElectrons[i].phi(), phot->eta(),phot->phi());
            double fsrDeltaPhi = fabs(deltaPhi(phot->phi(),recoElectrons[i].phi()));
            double fsrDeltaEta = fabs(phot->eta()-recoElectrons[i].eta());
            if( tmpDeltaREPh < 0.15){matched = true;}	
            if( fsrDeltaPhi < 2 && fsrDeltaEta < 0.05 ){matched = true;}
            if( tmpDeltaREPh < chosenDeltaRPh ){chosenDeltaRPh = tmpDeltaREPh;}
          }
          for(int i=0; i<(int)recoMuons.size(); i++)
          {
            double tmpDeltaRMPh = deltaR(recoMuons[i].eta(), recoMuons[i].phi(), phot->eta(),phot->phi());
            if( tmpDeltaRMPh < chosenDeltaRPh ){chosenDeltaRPh = tmpDeltaRMPh;}
          }

          if( (phot->pt() > 2 && chosenDeltaRPh < 0.07) || (phot->pt() > 4 && chosenDeltaRPh < 0.5) )
          {
            if(!matched && fabs(phot->eta()) < 2.4){ fsrPhotons.push_back(*phot); deltaRVec.push_back(chosenDeltaRPh); nPhotons++;}
          }
        }
      }	   

      if(doVarDump) photonDump->fillPhotonDumpTree(fsrPhotons,deltaRVec,iEvent,PV);



      // =========== BEGIN Extra Event Info =========== //

      int lepCounter = 0;
      for(unsigned int i = 0; i < recoMuons.size(); i++)
      {
        extraLep_id[lepCounter] = recoMuons[i].pdgId();
        extraLep_pT[lepCounter] = recoMuons[i].pt();
        extraLep_iso[lepCounter] = helper.pfIso(recoMuons[i],muonRho);
        extraLep_e[lepCounter] = recoMuons[i].energy();
        extraLep_pX[lepCounter] = recoMuons[i].px();
        extraLep_pY[lepCounter] = recoMuons[i].py();
        extraLep_pZ[lepCounter] = recoMuons[i].pz();
        extraLep_eta[lepCounter] = recoMuons[i].eta();
        extraLep_phi[lepCounter] = recoMuons[i].phi();
        extraLep_pX[lepCounter] = recoMuons[i].px();
        extraLep_pY[lepCounter] = recoMuons[i].py();
        extraLep_pZ[lepCounter] = recoMuons[i].pz();
        extraLep_chIso[lepCounter] = recoMuons[i].chargedHadronIso();
        extraLep_nhIso[lepCounter] = recoMuons[i].neutralHadronIso();
        extraLep_phIso[lepCounter] = recoMuons[i].photonIso();
        extraLep_sip[lepCounter] = helper.getSIP3D(recoMuons[i]);
        lepCounter++; 
      }
      
      for(unsigned int i = 0; i < recoElectrons.size(); i++)
      {
        extraLep_id[lepCounter] = recoElectrons[i].pdgId();
        extraLep_pT[lepCounter] = recoElectrons[i].pt();
        extraLep_iso[lepCounter] = helper.pfIso(recoElectrons[i],elecRho);
        extraLep_e[lepCounter] = recoElectrons[i].energy();
        extraLep_pX[lepCounter] = recoElectrons[i].px();
        extraLep_pY[lepCounter] = recoElectrons[i].py();
        extraLep_pZ[lepCounter] = recoElectrons[i].pz();
        extraLep_eta[lepCounter] = recoElectrons[i].eta();
        extraLep_phi[lepCounter] = recoElectrons[i].phi();
        extraLep_pX[lepCounter] = recoElectrons[i].px();
        extraLep_pY[lepCounter] = recoElectrons[i].py();
        extraLep_pZ[lepCounter] = recoElectrons[i].pz();
        extraLep_chIso[lepCounter] = recoElectrons[i].chargedHadronIso();
        extraLep_nhIso[lepCounter] = recoElectrons[i].neutralHadronIso();
        extraLep_phIso[lepCounter] = recoElectrons[i].photonIso();
        extraLep_sip[lepCounter] = helper.getSIP3D(recoElectrons[i]);
        lepCounter++;
      }
      
      for(unsigned int i = (recoElectrons.size() + recoMuons.size()); i < 10; i++)
      {
        extraLep_id[lepCounter] = 0;
        extraLep_pT[lepCounter] = 0;
        extraLep_iso[lepCounter] = 0;
        extraLep_e[lepCounter] = 0;
        extraLep_pX[lepCounter] = 0;
        extraLep_pY[lepCounter] = 0;
        extraLep_pZ[lepCounter] = 0;
        extraLep_eta[lepCounter] = 0;
        extraLep_phi[lepCounter] = 0;
        extraLep_pX[lepCounter] = 0;
        extraLep_pY[lepCounter] = 0;
        extraLep_pZ[lepCounter] = 0;
        extraLep_chIso[lepCounter] = 0;
        extraLep_nhIso[lepCounter] = 0;
        extraLep_phIso[lepCounter] = 0;
        extraLep_sip[lepCounter] = 0;
        lepCounter++;
      }
      

      nJets = 0;
      //nJets
      for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet)
      {
        if(jet->pt()>10)
        {
          bool matchesElec = false;
          bool matchesMuon = false;
          for(unsigned int i = 0; i < recoMuons.size(); i++)
          {
            if( jet->eta() >= (recoMuons[i].eta() - 0.1) && jet->eta() <= (recoMuons[i].eta() + 0.1) )
            {
              if( jet->phi() >= (recoMuons[i].phi() - 0.1) && jet->phi() <= (recoMuons[i].phi() + 0.1) )
              {
                matchesMuon = true;
              }
            }
          }

          for(unsigned int i = 0; i < recoElectrons.size(); i++)
          {
            if( jet->eta() >= (recoElectrons[i].eta() - 0.1) && jet->eta() <= (recoElectrons[i].eta() + 0.1) )
            {
              if( jet->phi() >= (recoElectrons[i].phi() - 0.1) && jet->phi() <= (recoElectrons[i].phi() + 0.1) )
              {
                matchesElec = true;
              }
            }   
          }
		    
          if( matchesElec == false && matchesMuon == false ){ nJets++;}

        } 
      } 

      //MET
      metVal= mets->empty() ? 0 : (*mets)[0].et();
	            
      // =========== END Extra Event Info  =========== //                  

      vector<pat::Muon> selectedMuons;
      vector<pat::Electron> selectedElectrons;

      RecoFourMixEvent = 0;
      //if(mixedFlavorCharge){findHiggsCandidate_MixFlavour(recoMuons,recoElectrons,selectedMuons,selectedElectrons, true);}
      std::vector<pat::PFParticle> selectedFsrPhotons; 
      findHiggsCandidate(recoMuons,recoElectrons,fsrPhotons,deltaRVec,selectedMuons,selectedElectrons,selectedFsrPhotons,iEvent);
	  
      if( foundHiggsCandidate )
      {	    
        //VBF Jets
        vector<pat::Jet> selectedVBFJets, correctedVBFJets;
        double tempDeltaR = -999;
        //if(doVarDump) jetDump->fillJetDumpTree(jets,correctedJets,iEvent); // for miniAOD
        if(doVarDump) jetDump->fillJetDumpTree(jets,jets,iEvent); // for miniAOD

        for(unsigned int i = 0; i < jets->size(); ++i) 
        {
          const pat::Jet & patjet = jets->at(i);
          //const pat::Jet & correctedJet = correctedJets->at(i); // for miniAOD
          const pat::Jet & correctedJet = jets->at(i); // for miniAOD
          //float mva   = (*puJetIdMva)[jets->refAt(i)];
          //int  idflag = (*puJetIdFlag)[jets->refAt(i)]; // for miniAOD
          //cout << "jet " << i << " pt " << patjet.pt() << " eta " << patjet.eta() << " PU JetID MVA " << mva << endl; 
          //if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose ) )  // for miniAOD
          if (patjet.userFloat("pileupJetId:fullDiscriminant")) // for miniAOD
          {
            //PF ID
            bool looseIDPass = false; 
            if(jetHelper.patjetID(correctedJet) == 1) looseIDPass = true;
            if(looseIDPass)
            {
              bool isDeltaR = true;
              for(unsigned int phIndex = 0; phIndex < selectedFsrPhotons.size(); phIndex++)
              {
                tempDeltaR = deltaR(patjet.eta(),patjet.phi(),selectedFsrPhotons[phIndex].eta(),selectedFsrPhotons[phIndex].phi());
                if (tempDeltaR < 0.5) isDeltaR = false;
              }
              if(correctedJet.pt() > 30 && fabs(patjet.eta()) < 4.7 && isDeltaR)
              {
                selectedVBFJets.push_back(patjet);
                correctedVBFJets.push_back(correctedJet);
              }
            }
          }
        }
      
        thetaZ2 = Z2Vec.Theta();
        etaZ2 = Z2Vec.Eta(); 

        passedPtCuts = helper.passedPtandEarlyIso(selectedMuons, selectedElectrons, leadingPtCut, subleadingPtCut, 1000, "PF", "dB","PFEffAreaRho", muonRho);

        if( passedPtCuts )
        {
          if(isSignal)
          {
            if(RecoFourMuEvent) sigEff_4->advanceSigNumCounters_PT2010(eventType,"reco4mu",eventWeight);
            if(RecoFourEEvent) sigEff_4->advanceSigNumCounters_PT2010(eventType,"reco4e",eventWeight);
            if(RecoTwoMuTwoEEvent || RecoTwoETwoMuEvent) sigEff_4->advanceSigNumCounters_PT2010(eventType,"reco2e2mu",eventWeight);
          }

          nEvPassedPtCut2 += eventWeight;		  
          if( RecoFourMuEvent ) nEvPassedPtCut2_4mu += eventWeight;
          if( RecoFourEEvent ) nEvPassedPtCut2_4e += eventWeight;
          if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ) nEvPassedPtCut2_2e2mu += eventWeight;
  
          // ========= m2l > 4 cut ========= //
          fourLep_Cleaned = helper.passedM2lCut_OS(selectedMuons,selectedElectrons,Z2Vec.M()< 12 ?0:4 ,10000);
          passedQCDcut = helper.passedM2lCut_OS(selectedMuons,selectedElectrons,4,10000);
          // ================================ //
          if(fourLep_Cleaned)
          {
            if(isSignal)
            {
              if(RecoFourMuEvent){sigEff_4->advanceSigNumCounters_4GeV(eventType,"reco4mu",eventWeight);}
              if(RecoFourEEvent){sigEff_4->advanceSigNumCounters_4GeV(eventType,"reco4e",eventWeight);}
              if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_4->advanceSigNumCounters_4GeV(eventType,"reco2e2mu",eventWeight);}
            }

            nEvAfterCleaning += eventWeight;
            if( RecoFourMuEvent ) nEvAfterCleaning_4mu += eventWeight;
            if( RecoFourEEvent ) nEvAfterCleaning_4e += eventWeight;
            if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ) nEvAfterCleaning_2e2mu += eventWeight;
	
            //Min Mass 3L
            minM3l = helper.minMass3l(selectedMuons,selectedElectrons);

            //M4L Error
            massErrorUCSD = massErr.getMassResolution(selectedElectrons, selectedMuons, selectedFsrPhotons);
            massErrorUCSDCorr = massErr.getMassResolutionCorr(selectedElectrons, selectedMuons, selectedFsrPhotons, true, !isMC);
            if(RecoFourMuEvent)
            {
              massErrorUF = massErr.calc4muErr(selectedMuons,selectedFsrPhotons,false,!isMC);
              massErrorUFCorr = massErr.calc4muErr(selectedMuons,selectedFsrPhotons,true,!isMC);
            }
            if(RecoFourEEvent)
            {
              massErrorUF = massErr.calc4eErr(selectedElectrons,selectedFsrPhotons,false,!isMC);
              massErrorUFCorr = massErr.calc4eErr(selectedElectrons,selectedFsrPhotons,true,!isMC);                      
            }
            if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent)
            { 
              massErrorUF = massErr.calc2e2muErr(selectedElectrons,selectedMuons,selectedFsrPhotons,false,!isMC);
              massErrorUFCorr = massErr.calc2e2muErr(selectedElectrons,selectedMuons,selectedFsrPhotons,true,!isMC);               
            }
	
            //MinDeltaR
            minDeltR = getMinDeltaR(selectedMuons, selectedElectrons);
		
            if( m4l > m4lLowCut )
            {
              passedM4lCut = true;
              if(isSignal)
              {
                if( RecoFourMuEvent ){sigEff_4->advanceSigNumCounters_M4L(eventType, "reco4mu",eventWeight);}
                if( RecoFourEEvent  ){sigEff_4->advanceSigNumCounters_M4L(eventType, "reco4e",eventWeight);}
                if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){sigEff_4->advanceSigNumCounters_M4L(eventType, "reco2e2mu",eventWeight); }
              }
            }
		
            // ========= BEGIN SIP ======== //
            for(unsigned int i = 0; i < selectedMuons.size(); i++)
            {  
              if( helper.getSIP3D(selectedMuons[i]) > highestSip ){ highestSip = helper.getSIP3D(selectedMuons[i]); }
            }
            for(unsigned int i = 0; i < selectedElectrons.size(); i++)
            {
              if( helper.getSIP3D(selectedElectrons[i]) > highestSip ){ highestSip = helper.getSIP3D(selectedElectrons[i]); }
            }
            // ========= END SIP Cuts ======== //            
		
            // ========= BEGIN ISO Cuts ========= //
            double leastIso1, leastIso2;
            helper.passedIsolation(selectedMuons,selectedElectrons,isoCut,"PF","dB","PFEffAreaRho",muonRho,elecRho,leastIso1,leastIso2);
            worstIso = leastIso1;	    
            // ========= END ISO Cuts ======== //
  
            if(RecoFourMuEvent) 
            {
              HiggsCandVecNoFSR = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedMuons[2].p4() + selectedMuons[3].p4();
            }
            if(RecoFourEEvent)
            {
              HiggsCandVecNoFSR = selectedElectrons[0].p4() + selectedElectrons[1].p4() + selectedElectrons[2].p4() + selectedElectrons[3].p4();
            }
            if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent)
            {
              HiggsCandVecNoFSR = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedElectrons[0].p4() + selectedElectrons[1].p4();
            }
	      
            //Set All the Variables for Saved Trees --- must be done after variables are available
            setTreeVariables(iEvent, iSetup, selectedMuons, selectedElectrons, selectedVBFJets, correctedVBFJets);

            //Calculate Angles
            HP4.SetPxPyPzE(HiggsCandVec.Px(), HiggsCandVec.Py(), HiggsCandVec.Pz(), HiggsCandVec.E() );
            Z1P4.SetPxPyPzE(pXZ1,pYZ1,pZZ1,EZ1);
            Z2P4.SetPxPyPzE(pXZ2,pYZ2,pZZ2,EZ2);

            int tmpIdL1,tmpIdL2,tmpIdL3,tmpIdL4;
            if(chargeL1 < 0){ L11P4.SetPxPyPzE(Lep1.Px(),Lep1.Py(),Lep1.Pz(),Lep1.E()); tmpIdL1 = idL1;}
            else{ L11P4.SetPxPyPzE(Lep2.Px(),Lep2.Py(),Lep2.Pz(),Lep2.E()); tmpIdL1 = idL2;}
            if(chargeL2 > 0){ L12P4.SetPxPyPzE(Lep2.Px(),Lep2.Py(),Lep2.Pz(),Lep2.E()); tmpIdL2 = idL2;}
            else{ L12P4.SetPxPyPzE(Lep1.Px(),Lep1.Py(),Lep1.Pz(),Lep1.E()); tmpIdL2 = idL1;}
            if(chargeL3 < 0){ L21P4.SetPxPyPzE(Lep3.Px(),Lep3.Py(),Lep3.Pz(),Lep3.E()); tmpIdL3 = idL3;}
            else{ L21P4.SetPxPyPzE(Lep4.Px(),Lep4.Py(),Lep4.Pz(),Lep4.E()); tmpIdL3 = idL4;}
            if(chargeL4 > 0){ L22P4.SetPxPyPzE(Lep4.Px(),Lep4.Py(),Lep4.Pz(),Lep4.E()); tmpIdL4 = idL4;}
            else{ L22P4.SetPxPyPzE(Lep3.Px(),Lep3.Py(),Lep3.Pz(),Lep3.E()); tmpIdL4 = idL3;}

            //MELA
            //mela->computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4,cosThetaStar,cosTheta1,cosTheta2,Phi,phiStar1,melaLD,
            //      mela_Sig, mela_Bkg,false /*withPt*/,false/*withY*/);
            //mela->computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4,cosThetaStar,cosTheta1,cosTheta2,Phi,phiStar1,melaLD_Pt,
            //      mela_Sig_Pt, mela_Bkg_Pt,true/*withPt*/,false/*withY*/);
            //mela->computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4,cosThetaStar,cosTheta1,cosTheta2,Phi,phiStar1,melaLD_Y,
            //      mela_Sig_Y, mela_Bkg_Y,false/*withPt*/,true/*withY*/);
            //mela->computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4,cosThetaStar,cosTheta1,cosTheta2,Phi,phiStar1,melaLD_PtY,
            //	    mela_Sig_PtY, mela_Bkg_PtY,true/*withPt*/,true/*withY*/);
            /*
            int tmpFlavor;
            if(RecoFourEEvent) tmpFlavor = 1;
            else if(RecoFourMuEvent) tmpFlavor = 2;
            else tmpFlavor = 3;
	      
            float tmpHPt = HP4.Pt();
            float tmpHM  = HP4.M();
            float tmpMassZ1 = massZ1;
            float tmpMassZ2 = massZ2;
            float tmpHY = HP4.Eta();

            //doesnt work with MEMs version
            mela->computeP(tmpHM, tmpMassZ1, tmpMassZ2,cosThetaStar,cosTheta1,cosTheta2,Phi,phiStar1,
                          //signal probabilities
                          p0plus_melaNorm,   // higgs, analytic distribution, normalized as for normal MELA distribution     
                          p0plus_mela,   // higgs, analytic distribution          
                          p0minus_mela,  // pseudoscalar, analytic distribution 
                          p0plus_VAJHU,  // higgs, vector algebra, JHUgen
                          p0minus_VAJHU, // pseudoscalar, vector algebra, JHUgen
                          p0plus_VAMCFM,// higgs, vector algebra, MCFM
                          p1_mela,  // zprime, analytic distribution 
                          p1_VAJHU, // zprime, vector algebra, JHUgen,
                          p2_mela , // graviton, analytic distribution 
                          p2_VAJHU, // graviton, vector algebra, JHUgen,
                          //backgrounds
                          mela_bkg_analytic,  // background,  analytic distribution 
                          mela_bkg_VAMCFM, // background, vector algebra, MCFM
                          mela_bkg_ggzz_VAMCFM, // background, vector algebra, MCFM for ggZZ
                          mela_bkg_VAMCFMNorm, // background, vector algebra, MCFM, Normalized 
                          //pt/rapidity
                          mela_p0_pt, // multiplicative probability for signal pt
                          mela_p0_y, // multiplicative probability for signal y
                          mela_bkg_pt, // multiplicative probability for bkg pt
                          mela_bkg_y, // multiplicative probability for bkg y
                          //supermela
                          p0plus_m4l,  // signal m4l probability as in datacards
                          bkg_m4l,     // backgroun m4l probability as in datacards
                          //optional input parameters
                          tmpHPt,
                          tmpHY,
                          tmpFlavor // 1:4e, 2:4mu, 3:2e2mu (for interference effects)
                          );
	      
	      
            pseudoMela.computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4, pseudoMelaLD, pseudoMela_SM, pseudoMela_PS);
            spin1EvMela.computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4, spin1EvMelaLD, spin1EvMela_SM, spin1EvMela_S1E);
            spin1OddMela.computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4, spin1OddMelaLD, spin1OddMela_SM, spin1OddMela_S1O);
            spin2Mela.computeKD(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4, spin2MelaLD, spin2Mela_SM, spin2Mela_S2);

            //MEKD
            MEKDnoPDFs->computeMEs(L11P4, tmpIdL1, L12P4, tmpIdL2, L21P4, tmpIdL3, L22P4, tmpIdL4);
            MEKDnoPDFs->computeKD( (string)"SMHiggs", (string)"ZZ", MEKD_noPDF, MEKD_ME_H_noPDF, MEKD_ME_ZZ_noPDF);
            MEKDnoPDFs->computeKD( (string)"Higgs0M", (string)"ZZ", MEKD_PSZZ_noPDF, MEKD_PSZZ_ME_PS_noPDF, MEKD_PSZZ_ME_ZZ_noPDF);
            MEKDnoPDFs->computeKD( (string)"Graviton2PM", (string)"ZZ", MEKD_S2ZZ_noPDF, MEKD_S2ZZ_ME_S2_noPDF, MEKD_S2ZZ_ME_ZZ_noPDF);
*/

            if(chargeL1 < 0){ L11P4_noFSR.SetPxPyPzE(pXL1,pYL1,pZL1,EL1); tmpIdL1 = idL1;}
            else{ L11P4_noFSR.SetPxPyPzE(pXL2,pYL2,pZL2,EL2); tmpIdL1 = idL2;}
            if(chargeL2 > 0){ L12P4_noFSR.SetPxPyPzE(pXL2,pYL2,pZL2,EL2); tmpIdL2 = idL2;}
            else{ L12P4_noFSR.SetPxPyPzE(pXL1,pYL1,pZL1,EL1); tmpIdL2 = idL1;}
            if(chargeL3 < 0){ L21P4_noFSR.SetPxPyPzE(pXL3,pYL3,pZL3,EL3); tmpIdL3 = idL3;}
            else{ L21P4_noFSR.SetPxPyPzE(pXL4,pYL4,pZL4,EL4); tmpIdL3 = idL4;}
            if(chargeL4 > 0){ L22P4_noFSR.SetPxPyPzE(pXL4,pYL4,pZL4,EL4); tmpIdL4 = idL4;}
            else{ L22P4_noFSR.SetPxPyPzE(pXL3,pYL3,pZL3,EL3); tmpIdL4 = idL3;}

            //MEKDnoPDFs_noFSR->computeMEs(L11P4_noFSR, tmpIdL1, L12P4_noFSR, tmpIdL2, L21P4_noFSR, tmpIdL3, L22P4_noFSR, tmpIdL4);
            //MEKDnoPDFs_noFSR->computeKD( (string)"SMHiggs", (string)"ZZ", MEKD_noPDF_noFSR, MEKD_ME_H_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            //MEKDnoPDFs_noFSR->computeKD( (string)"Higgs0M", (string)"ZZ", MEKD_PSZZ_noPDF_noFSR, MEKD_PSZZ_ME_PS_noPDF_noFSR, MEKD_PSZZ_ME_ZZ_noPDF_noFSR);
            //MEKDnoPDFs_noFSR->computeKD( (string)"Graviton2PM", (string)"ZZ", MEKD_S2ZZ_noPDF_noFSR, MEKD_S2ZZ_ME_S2_noPDF_noFSR, MEKD_S2ZZ_ME_ZZ_noPDF_noFSR);
  
            vector<TLorentzVector> P4s_noFSR, P4s;
            vector<int> tmpIDs;
            P4s.push_back(L11P4); P4s_noFSR.push_back(L11P4_noFSR);
            P4s.push_back(L12P4); P4s_noFSR.push_back(L12P4_noFSR);
            P4s.push_back(L21P4); P4s_noFSR.push_back(L21P4_noFSR);
            P4s.push_back(L22P4); P4s_noFSR.push_back(L22P4_noFSR);
            tmpIDs.push_back(tmpIdL1);
            tmpIDs.push_back(tmpIdL2);
            tmpIDs.push_back(tmpIdL3);
            tmpIDs.push_back(tmpIdL4);

            /*  removed for miniAOD , 70 lines below
            MEMsnoPDFs->computeMEs(P4s,tmpIDs);
            interferenceWeight = MEMsnoPDFs->getMELAWeight();
            MEMsnoPDFs->computePm4l(P4s,tmpIDs, MEMNames::kNone, pdfSigM4l, pdfBkgM4l);
            MEMsnoPDFs->computePm4l(P4s,tmpIDs, MEMNames::kScaleUp, pdfSigM4l_ScaleUp, pdfBkgM4l_ScaleUp);
            MEMsnoPDFs->computePm4l(P4s,tmpIDs, MEMNames::kResolUp, pdfSigM4l_ResUp, pdfBkgM4l_ResUp);
            MEMsnoPDFs->computePm4l(P4s,tmpIDs, MEMNames::kScaleDown, pdfSigM4l_ScaleDown, pdfBkgM4l_ScaleDown);
            MEMsnoPDFs->computePm4l(P4s,tmpIDs, MEMNames::kResolDown, pdfSigM4l_ResDown, pdfBkgM4l_ResDown);

            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, kqqZZ,  kMCFM, &MEMs::probRatio, JHUKD_H_qqZZ_noPDF,  JHU_ME_H, MCFM_ME_qqZZ);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k0minus, kJHUGen, &MEMs::probRatio, JHUKD_H_h0M_noPDF,   JHU_ME_H, JHU_ME_h0M);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k0hplus,       kJHUGen, &MEMs::probRatio, JHUKD_H_h0P_noPDF,   JHU_ME_H, JHU_ME_h0P);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k1minus,       kJHUGen, &MEMs::probRatio, JHUKD_H_h1M_noPDF,   JHU_ME_H, JHU_ME_h1M);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k1plus,        kJHUGen, &MEMs::probRatio, JHUKD_H_h1P_noPDF,   JHU_ME_H, JHU_ME_h1P);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k2mplus_gg,    kJHUGen, &MEMs::probRatio, JHUKD_H_ggh2P_noPDF, JHU_ME_H, JHU_ME_ggh2P);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k2mplus_qqbar, kJHUGen, &MEMs::probRatio, JHUKD_H_qqh2P_noPDF, JHU_ME_H, JHU_ME_qqh2P);
            MEMsnoPDFs->computeKD(kSMHiggs, kMELA_HCP, kqqZZ,  kMELA_HCP, &MEMs::probRatio, melaLD,  mela_Sig, mela_Bkg);
	    
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k2hplus,  kJHUGen, &MEMs::probRatio, JHUKD_H_h2hP_noPDF, JHU_ME_H, JHU_ME_h2hP);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k2hminus, kJHUGen, &MEMs::probRatio, JHUKD_H_h2hM_noPDF, JHU_ME_H, JHU_ME_h2hM);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k2bplus,  kJHUGen, &MEMs::probRatio, JHUKD_H_h2bP_noPDF, JHU_ME_H, JHU_ME_h2bP);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k2mplus_prodIndep, kJHUGen, &MEMs::probRatio, JHUKD_H_h2P_prodInd_noPDF,
                                  JHU_ME_H, JHU_ME_h2P_prodInd);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k1plus_prodIndep,  kJHUGen, &MEMs::probRatio, JHUKD_H_h1P_prodInd_noPDF,
                                  JHU_ME_H, JHU_ME_h1P_prodInd);
            MEMsnoPDFs->computeKD(kSMHiggs, kJHUGen, k1minus_prodIndep, kJHUGen, &MEMs::probRatio, JHUKD_H_h1M_prodInd_noPDF,
                                  JHU_ME_H, JHU_ME_h1M_prodInd);

            MEMsnoPDFs_noFSR->computeMEs(P4s_noFSR,tmpIDs);
            MEMsnoPDFs_noFSR->computePm4l(P4s_noFSR,tmpIDs, MEMNames::kNone, pdfSigM4l_noFSR, pdfBkgM4l_noFSR);
            MEMsnoPDFs_noFSR->computePm4l(P4s_noFSR,tmpIDs, MEMNames::kScaleUp, pdfSigM4l_ScaleUp_noFSR, pdfBkgM4l_ScaleUp_noFSR);
            MEMsnoPDFs_noFSR->computePm4l(P4s_noFSR,tmpIDs, MEMNames::kResolUp, pdfSigM4l_ResUp_noFSR, pdfBkgM4l_ResUp_noFSR);
            MEMsnoPDFs_noFSR->computePm4l(P4s_noFSR,tmpIDs, MEMNames::kScaleDown, pdfSigM4l_ScaleDown_noFSR, pdfBkgM4l_ScaleDown_noFSR);
            MEMsnoPDFs_noFSR->computePm4l(P4s_noFSR,tmpIDs, MEMNames::kResolDown, pdfSigM4l_ResDown_noFSR, pdfBkgM4l_ResDown_noFSR);

            MEMsnoPDFs_noFSR->computeKD(kSMHiggs, kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_noPDF_noFSR, MEKD_ME_H_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k0minus, kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_h0M_ZZ_noPDF_noFSR, MEKD_ME_h0M_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k0hplus, kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_h0P_ZZ_noPDF_noFSR, MEKD_ME_h0P_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k1minus,kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_h1M_ZZ_noPDF_noFSR, MEKD_ME_h1M_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k1plus,kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_h1P_ZZ_noPDF_noFSR, MEKD_ME_h1P_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k2mplus_qqbar,kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_qqh2P_ZZ_noPDF_noFSR,MEKD_ME_qqh2P_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k2mplus_gg,kMEKD,kqqZZ,kMEKD,&MEMs::probRatio,MEKD_ggh2P_ZZ_noPDF_noFSR,MEKD_ME_ggh2P_noPDF_noFSR,MEKD_ME_ZZ_noPDF_noFSR);


            MEMsnoPDFs_noFSR->computeKD(k2hplus, kMEKD, kqqZZ,  kMEKD, &MEMs::probRatio, MEKD_h2hP_ZZ_noPDF_noFSR, 
                                        MEKD_ME_h2hP_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k2hminus, kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_h2hM_ZZ_noPDF_noFSR, 
                                        MEKD_ME_h2hM_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k2bplus, kMEKD, kqqZZ,  kMEKD, &MEMs::probRatio, MEKD_h2bP_ZZ_noPDF_noFSR, 
                                        MEKD_ME_h2bP_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k2mplus_prodIndep, kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_h2P_ZZ_prodInd_noPDF_noFSR,
                                        MEKD_ME_h2P_prodInd_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k1plus_prodIndep, kMEKD, kqqZZ,  kMEKD, &MEMs::probRatio, MEKD_h1P_ZZ_prodInd_noPDF_noFSR,
                                        MEKD_ME_h1P_prodInd_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k1minus_prodIndep, kMEKD, kqqZZ, kMEKD, &MEMs::probRatio, MEKD_h1M_ZZ_prodInd_noPDF_noFSR,
                                        MEKD_ME_h1M_prodInd_noPDF_noFSR, MEKD_ME_ZZ_noPDF_noFSR);

            //cout << "               " << MEMsnoPDFs_noFSR->computeME(k2mplus_prodIndep, kMEKD, P4s_noFSR,tmpIDs, MEKD_ME_h2P_prodInd_noPDF_noFSR) << endl;


            // JHUGen/MCMF cross-check
            MEMsnoPDFs_noFSR->computeKD(kSMHiggs,     kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, JHUKD_H_qqZZ_noPDF_noFSR,  JHU_ME_H_noPDF_noFSR,   JHU_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k0minus,      kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, JHUKD_h0M_ZZ_noPDF_noFSR,  JHU_ME_h0M_noPDF_noFSR, JHU_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k0hplus,      kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, JHUKD_h0P_ZZ_noPDF_noFSR,  JHU_ME_h0P_noPDF_noFSR, JHU_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k1minus,      kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, JHUKD_h1M_ZZ_noPDF_noFSR,  JHU_ME_h1M_noPDF_noFSR, JHU_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k1plus,       kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, JHUKD_h1P_ZZ_noPDF_noFSR,  JHU_ME_h1P_noPDF_noFSR, JHU_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k2mplus_qqbar,kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, JHUKD_qqh2P_ZZ_noPDF_noFSR,JHU_ME_qqh2P_noPDF_noFSR,JHU_ME_ZZ_noPDF_noFSR);
            MEMsnoPDFs_noFSR->computeKD(k2mplus_gg,   kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, JHUKD_ggh2P_ZZ_noPDF_noFSR,JHU_ME_ggh2P_noPDF_noFSR,JHU_ME_ZZ_noPDF_noFSR);
            */   

            Z4lmaxP = helper.largestLepMomentum(L11P4,L12P4,L21P4,L22P4);
            vector<double> thetas = angles.angleBetweenLep(L11P4,L12P4,L21P4,L22P4);
            theta12 = thetas[0]; theta13 = thetas[1]; theta14 = thetas[2];
            theta12_deg = acos(theta12)*180/PI; theta13_deg = acos(theta13)*180/PI;
            theta14_deg = acos(theta14)*180/PI;
            thetaPhoton = angles.minAngleOfPhoton(L11P4,L12P4,L21P4,L22P4);
            thetaPhoton_deg = acos(thetaPhoton)*180/PI;
            thetaPhotonZ = angles.angleOfPhotonZframe(L11P4,L12P4,L21P4,L22P4,chargeL3);
            thetaPhotonZ_deg = acos(thetaPhotonZ)*180/PI;
            maxMass2Lep = helper.maxMass2l(recoMuons,recoElectrons);
	      
            //Fill Passed Events Tree
            if(!doBlinding || isMC)
            {
              double globalPtMuCut = 5, globalPtElCut = 7;
              bool passedPtSelection = true;
              if(Z2Vec.M() > 0)
              {
                if( abs(idL1) == 13 && pTL1 < globalPtMuCut ) passedPtSelection = false;
                if( abs(idL2) == 13 && pTL2 < globalPtMuCut ) passedPtSelection = false;
                if( abs(idL3) == 13 && pTL3 < globalPtMuCut ) passedPtSelection = false;
                if( abs(idL4) == 13 && pTL4 < globalPtMuCut ) passedPtSelection = false;
                if( abs(idL1) == 11 && pTL1 < globalPtElCut ) passedPtSelection = false;
                if( abs(idL2) == 11 && pTL2 < globalPtElCut ) passedPtSelection = false;
                if( abs(idL3) == 11 && pTL3 < globalPtElCut ) passedPtSelection = false;
                if( abs(idL4) == 11 && pTL4 < globalPtElCut ) passedPtSelection = false;
                if(passedPtSelection)
                {
                  passedZ4lSelection = true;
                  if(Z2Vec.M() > 12) passedFullSelection = true;
                }
              }
              if(isMC)
              {
                std::vector<reco::GenParticle> Higgs,Zs,leptonsS1,leptonsS3;
                genAna.fillGenEvent(genParticles,Higgs,Zs,leptonsS1,leptonsS3);
                setGENVariables(Higgs,Zs,leptonsS1,leptonsS3);
                //setGENMatchedVariables(selectedMuons,selectedElectrons); // for miniAOD
              }
              passedEventsTree_All->Fill();
            }

            if(passedM4lCut && isSignal)
            {
              if(Z2Vec.M() > 4)
              {
                if(RecoFourMuEvent){sigEff_4->advanceSigNumCounters_FINAL(eventType,"reco4mu",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoFourEEvent){sigEff_4->advanceSigNumCounters_FINAL(eventType,"reco4e",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_4->advanceSigNumCounters_FINAL(eventType,"reco2e2mu",eventWeight,selectedMuons,selectedElectrons);}
              }
  
              if(Z2Vec.M() > 6)
              {
                if(RecoFourMuEvent){sigEff_6->advanceSigNumCounters_FINAL(eventType,"reco4mu",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoFourEEvent){sigEff_6->advanceSigNumCounters_FINAL(eventType,"reco4e",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_6->advanceSigNumCounters_FINAL(eventType,"reco2e2mu",eventWeight,selectedMuons,selectedElectrons);}
              }
              if(Z2Vec.M() > 8)
              {
                if(RecoFourMuEvent){sigEff_8->advanceSigNumCounters_FINAL(eventType,"reco4mu",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoFourEEvent){sigEff_8->advanceSigNumCounters_FINAL(eventType,"reco4e",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_8->advanceSigNumCounters_FINAL(eventType,"reco2e2mu",eventWeight,selectedMuons,selectedElectrons);}
              }

              if(Z2Vec.M() > 9)
              {
                if(RecoFourMuEvent){sigEff_9->advanceSigNumCounters_FINAL(eventType,"reco4mu",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoFourEEvent){sigEff_9->advanceSigNumCounters_FINAL(eventType,"reco4e",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_9->advanceSigNumCounters_FINAL(eventType,"reco2e2mu",eventWeight,selectedMuons,selectedElectrons);}
              }

              if(Z2Vec.M() > 10)
              {
                if(RecoFourMuEvent){sigEff_10->advanceSigNumCounters_FINAL(eventType,"reco4mu",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoFourEEvent){sigEff_10->advanceSigNumCounters_FINAL(eventType,"reco4e",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_10->advanceSigNumCounters_FINAL(eventType,"reco2e2mu",eventWeight,selectedMuons,selectedElectrons);}
              }

              if(Z2Vec.M() > 11)
              {
                if(RecoFourMuEvent){sigEff_11->advanceSigNumCounters_FINAL(eventType,"reco4mu",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoFourEEvent){sigEff_11->advanceSigNumCounters_FINAL(eventType,"reco4e",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_11->advanceSigNumCounters_FINAL(eventType,"reco2e2mu",eventWeight,selectedMuons,selectedElectrons);}
              }
 
              if(Z2Vec.M() > 12)
              {
                if(RecoFourMuEvent){sigEff_12->advanceSigNumCounters_FINAL(eventType,"reco4mu",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoFourEEvent){sigEff_12->advanceSigNumCounters_FINAL(eventType,"reco4e",eventWeight,selectedMuons,selectedElectrons);}
                if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){sigEff_12->advanceSigNumCounters_FINAL(eventType,"reco2e2mu",eventWeight,selectedMuons,selectedElectrons);}
              }
            }
	      
            //if( m4l > 100 && m4l < 120 ){Zto4LAnaOP->fillZto4LHistograms(histContainer_, selectedMuons, selectedElectrons,eventWeight);}
            if( m4l < 100 )
            {
              if( m4l > 80 && m4l < 100 )
              {
                if(RecoFourMuEvent)Z4lcounter_4mu += eventWeight;
                if(RecoFourEEvent) Z4lcounter_4e += eventWeight;
                if(RecoTwoMuTwoEEvent) Z4lcounter_2mu2e += eventWeight;
                if(RecoTwoETwoMuEvent) Z4lcounter_2e2mu += eventWeight;
                m3l_soft = helper.M3lSoftestLep(L11P4,L12P4,L21P4,L22P4);
                //Zto4LAna.fillZto4LHistograms(histContainer_, selectedMuons, selectedElectrons,eventWeight);
              }
            }
     
            if( m4l > 70 && Z2Vec.M() > 12)
            {
              nEvAfterZ4lCut += eventWeight;
              if( RecoFourMuEvent ){nEvAfterZ4lCut_4mu += eventWeight; }
              if( RecoFourEEvent ){nEvAfterZ4lCut_4e += eventWeight; }
              if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){nEvAfterZ4lCut_2e2mu += eventWeight; }
            }

            if( Z2Vec.M() > 12 )
            {
              if( passedM4lCut)
              {
                if(doVarDump) finalLepDump->fillFinalLepDumpTree(selectedMuons,selectedElectrons,iEvent,muonRho,elecRho,PV,elecID);      
                sipAna.advanceSipCounters(highestSip,eventWeight);
                //Step Plots
                if(isMC || !doBlinding)
                {
                  nEvAfterM4lCut += eventWeight;
                  if( RecoFourMuEvent ){nEvAfterM4lCut_4mu += eventWeight;}
                  if( RecoFourEEvent  ){nEvAfterM4lCut_4e += eventWeight; }
                  if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterM4lCut_2e2mu += eventWeight; }
                  if(VBFJet1 || VBFJet2 || (VBFJet1 && VBFJet2))
                  {
                    nEvAfterVBFJet1 += eventWeight;
                    if( RecoFourMuEvent ){nEvAfterVBFJet1_4mu += eventWeight;}
                    if( RecoFourEEvent  ){nEvAfterVBFJet1_4e += eventWeight; }
                    if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterVBFJet1_2e2mu += eventWeight; }
 		      
                    if(VBFJet2)
                    {
                      nEvAfterVBFJet2 += eventWeight;
                      if( RecoFourMuEvent ){nEvAfterVBFJet2_4mu += eventWeight;}
                      if( RecoFourEEvent  ){nEvAfterVBFJet2_4e += eventWeight; }
                      if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterVBFJet2_2e2mu += eventWeight; }
	  			
                      //cout << "VBF: " << VBFDiJetMass << "   " << VBFDeltaEta << endl;
  
                      if(FisherDiscrim > 0.4)
                      {      
                        nEvAfterVBFJetCuts += eventWeight;
                        if( RecoFourMuEvent ){nEvAfterVBFJetCuts_4mu += eventWeight;}
                        if( RecoFourEEvent  ){nEvAfterVBFJetCuts_4e += eventWeight; }
                        if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterVBFJetCuts_2e2mu += eventWeight; }
                      }
                    }
                  }  

                  if(JHUKD_H_h0M_noPDF > 0.3)
                  {
                    nEvAfterPsMelaCut += eventWeight;
                    if( RecoFourMuEvent ){nEvAfterPsMelaCut_4mu += eventWeight;}
                    if( RecoFourEEvent  ){nEvAfterPsMelaCut_4e += eventWeight; }
                    if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterPsMelaCut_2e2mu += eventWeight; }
                  }
                  if(JHUKD_H_ggh2P_noPDF > 0.15)
                  {
                    nEvAfterGrMelaCut += eventWeight;
                    if( RecoFourMuEvent ){nEvAfterGrMelaCut_4mu += eventWeight;}
                    if( RecoFourEEvent  ){nEvAfterGrMelaCut_4e += eventWeight; }
                    if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterGrMelaCut_2e2mu += eventWeight; }
                  }
  			  
                  if(JHUKD_H_qqZZ_noPDF > 0.1)
                  //if(melaLD > 0.1)
                  {
                    nEvAfterMelaCut += eventWeight;
                    if( RecoFourMuEvent ){nEvAfterMelaCut_4mu += eventWeight;}
                    if( RecoFourEEvent  ){nEvAfterMelaCut_4e += eventWeight; }
                    if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterMelaCut_2e2mu += eventWeight; }
                    if(FSR_Z1 || FSR_Z2)
                    {
                      if( (FSR_Z1 && !FSR_Z2) || (FSR_Z2 && !FSR_Z1) )
                      {
                        nEvWith1FSRZ += eventWeight;
                        if(RecoFourMuEvent) nEvWith1FSRZ_4mu += eventWeight;
                        if(RecoFourEEvent) nEvWith1FSRZ_4e += eventWeight;
                        if(RecoTwoETwoMuEvent) nEvWith1FSRZ_2e2mu += eventWeight;
                        if(RecoTwoMuTwoEEvent) nEvWith1FSRZ_2e2mu += eventWeight;
                      }
                      if(FSR_Z1 && FSR_Z2)
                      {
                        nEvWith2FSRZ += eventWeight;
                        if(RecoFourMuEvent) nEvWith2FSRZ_4mu += eventWeight;
                        if(RecoFourEEvent) nEvWith2FSRZ_4e += eventWeight;
                        if(RecoTwoETwoMuEvent) nEvWith2FSRZ_2e2mu += eventWeight;
                        if(RecoTwoMuTwoEEvent) nEvWith2FSRZ_2e2mu += eventWeight;
                      }
                    } 
                  }
 			  
                  if((Z1Vec.M() > 60 && Z1Vec.M() < 120) && (Z2Vec.M() > 60 && Z2Vec.M() < 120))
                  {
                    nEvAfterZZCut += eventWeight;
                    if( RecoFourMuEvent ){nEvAfterZZCut_4mu += eventWeight;}
                    if( RecoFourEEvent  ){nEvAfterZZCut_4e += eventWeight; }
                    if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterZZCut_2e2mu += eventWeight; }
                  }
                }
	      
                //Iso Efficiency
                if(isMC) isoEff->advanceIsoCounters(m4l, leastIso1, leastIso2, isSignal,eventWeight);
 		      
                if( m4l >180 )
                {
                  if( RecoFourMuEvent ){counter_4mu += eventWeight;}
                  if( RecoFourEEvent ){counter_4e += eventWeight;}
                  if( RecoTwoETwoMuEvent ){counter_2e2mu += eventWeight;}
                  if( RecoTwoMuTwoEEvent ){counter_2mu2e += eventWeight;}
                }
                if( (m4l > 100 && m4l <= 180) && (((!isMC && !doBlinding) || (m4l < 300 && m4l > 140) || m4l < 110) || isMC) )
                {
                  if( RecoFourMuEvent ){LMcounter_4mu += eventWeight;}
                  if( RecoFourEEvent ){LMcounter_4e += eventWeight;}
                  if( RecoTwoETwoMuEvent ){LMcounter_2e2mu += eventWeight;}
                  if( RecoTwoMuTwoEEvent ){LMcounter_2mu2e += eventWeight;}
                }
	 		    

                //Extra Particles
                int nPart = recoMuons.size() + recoElectrons.size() - 4;
                histContainer_["extraParticles"]->Fill(nPart,eventWeight);
  		      			    
                //Isolationg Efficiency
                nEvAfterIso += eventWeight;
                //Signal Efficiency
                if( RecoFourMuEvent ){nEvAfterIso_4mu += eventWeight; }
                if( RecoFourEEvent  ){nEvAfterIso_4e += eventWeight; }
                if( RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent ){nEvAfterIso_2e2mu += eventWeight;}
 		 
                if(bStudyResolution && bStudyFourLeptonResolution)
                {
                  FourLepReso->fillHistograms(hContainer_, hContainer2D_, selectedElectrons, selectedMuons, selectedFsrPhotons, eventWeight,!isMC);
                }
              }//passedM4lCut		  
            }//Z2 > 12      
          }//fourLep_Cleaned	    
        }//passedPtCuts	  
      }//if HC
    }//if ID   
  }//notDuplicate

  

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void UFHZZ4LAna::beginJob()
{
  //using namespace edm;
  using namespace std;
  using namespace pat;

  bookStepPlots();

  bookPassedEventTree("passedEvents", passedEventsTree_All);

  bookResolutionHistograms();
  //muonAna.bookMuonHistograms(fs,histContainer_);
  sipAna.bookSipHistograms(fs,histContainer_);
  //sigEff->bookSigEffHistograms(fs);
  sigEff_4->makeSigEffTree();
  sigEff_6->makeSigEffTree();
  sigEff_8->makeSigEffTree();
  sigEff_9->makeSigEffTree();
  sigEff_10->makeSigEffTree();
  sigEff_11->makeSigEffTree();
  sigEff_12->makeSigEffTree();
  isoEff->bookIsoHistograms(fs);
  Zto4LAna.bookZto4LHistograms(fs,histContainer_);
  //Zto4LAnaOP->bookZto4LHistograms(fs,histContainer_);

  if(bStudyResolution)
  {
    PerLepReso->bookHistograms(fs, hContainer_, hContainer2D_, hContainer3D_);
    if(bStudyDiLeptonResolution) { DiLepReso->bookHistograms(fs, hContainer_, hContainer2D_); }
    if(bStudyFourLeptonResolution) { FourLepReso->bookHistograms(fs, hContainer_, hContainer2D_);}
  }


}

// ------------ method called once each job just after ending the event loop  ------------
void  UFHZZ4LAna::endJob() 
{
  //using namespace edm;
  using namespace std;
  using namespace pat;

  gen4mu = sigEff_4->getNGen4mu();
  gen4e = sigEff_4->getNGen4e();
  gen2e2mu = sigEff_4->getNGen2e2mu();
   
  gen4muPseudo = sigEff_4->getNGen4muPseudo();
  gen4ePseudo = sigEff_4->getNGen4ePseudo();
  gen2e2muPseudo = sigEff_4->getNGen2e2muPseudo();


  cout << "XS:               " << CrossSection << endl
       << "Filter:           " << FilterEff << endl
       << "nEventsAfterSkim: " << nEvAfterSkim << endl
       << "nEvPassedHLT:     " << nEvPassedHlt << endl
       << "Weight:           " << scaleWeight << endl;

  cout << endl << endl;
  cout << "  4L/ 4e/ 4mu/ 2e2mu  " << endl;
  cout << "nEvTotalReco          " << endl
       << nEventsTotal << endl 
       << "nEvTotalGen           " << endl
       << gen4mu+gen4e+gen2e2mu << "/ " << gen4e << "/ " << gen4mu << "/ " << gen2e2mu << endl
       << "nEvAccepted:          " << endl
       << gen4muPseudo+gen4ePseudo+gen2e2muPseudo << "/ " << gen4ePseudo << "/ " << gen4muPseudo << "/ " << gen2e2muPseudo << endl
       << "nEvAfterId:           " << endl
       << nEvAfterId << endl
       << "nEvAfterZ1Formed:     " << endl 
       << nEvAfterZ1Formed << "/ " << nEvAfterZ1Formed_2e << "/ " << nEvAfterZ1Formed_2mu << endl
       << "nEvAfterZ1Cut:        " << endl
       << nEvAfterZ1Cut << "/ " << nEvAfterZ1Cut_2e << "/ " << nEvAfterZ1Cut_2mu << endl
       << counterElZ1ph << "/ " << counterMuZ1ph << endl
       << "nEv4GoodLep:          " << endl 
       << nEv4GoodLep << "/ " << nEv4GoodLep_4e << "/ " << nEv4GoodLep_4mu << "/ " << nEv4GoodLep_2e2mu << endl
       << "nEvAfterZ2Formed:     " << endl
       << nEvAfterZ2Formed << "/ " << nEvAfterZ2Formed_4e << "/ " << nEvAfterZ2Formed_4mu << "/ " << nEvAfterZ2Formed_2e2mu << endl
       << "nEvAfterZ2:           " << endl
       << nEvAfterZ2Cut << "/ " << nEvAfterZ2Cut_4e << "/ " << nEvAfterZ2Cut_4mu << "/ " << nEvAfterZ2Cut_2e2mu << endl
       << counterElZ2ph << "/ " << counterMuZ2ph << endl
       << "nEvPassedPtCut:       " << endl
       << nEvPassedPtCut2 << "/ " << nEvPassedPtCut2_4e << "/ " << nEvPassedPtCut2_4mu << "/ " << nEvPassedPtCut2_2e2mu << endl
       << "nEvAfter4GeV:         " << endl
       << nEvAfterCleaning << "/ " << nEvAfterCleaning_4e << "/ " << nEvAfterCleaning_4mu << "/ " << nEvAfterCleaning_2e2mu << endl
       << "nEvAfterZ4lCut:       " << endl
       << nEvAfterZ4lCut << "/ " << nEvAfterZ4lCut_4e << "/ " << nEvAfterZ4lCut_4mu << "/ " << nEvAfterZ4lCut_2e2mu << endl
       << "nEvAfterM4lCut:       " << endl
       << nEvAfterM4lCut << "/ " << nEvAfterM4lCut_4e << "/ " << nEvAfterM4lCut_4mu << "/ " << nEvAfterM4lCut_2e2mu << endl
       << "nEvAfterVBFJet1:       " << endl
       << nEvAfterVBFJet1 << "/ " << nEvAfterVBFJet1_4e << "/ " << nEvAfterVBFJet1_4mu << "/ " << nEvAfterVBFJet1_2e2mu << endl
       << "nEvAfterVBFJet2:       " << endl
       << nEvAfterVBFJet2 << "/ " << nEvAfterVBFJet2_4e << "/ " << nEvAfterVBFJet2_4mu << "/ " << nEvAfterVBFJet2_2e2mu << endl
       << "nEvAfterVBFJetCuts:       " << endl
       << nEvAfterVBFJetCuts << "/ " << nEvAfterVBFJetCuts_4e << "/ " << nEvAfterVBFJetCuts_4mu << "/ " << nEvAfterVBFJetCuts_2e2mu << endl
       << "nEvAfterZZCut:       " << endl
       << nEvAfterZZCut << "/ " << nEvAfterZZCut_4e << "/ " << nEvAfterZZCut_4mu << "/ " << nEvAfterZZCut_2e2mu << endl
       << "nEvAfterMelaCut:       " << endl
       << nEvAfterMelaCut << "/ " << nEvAfterMelaCut_4e << "/ " << nEvAfterMelaCut_4mu << "/ " << nEvAfterMelaCut_2e2mu << endl
       << "nEvAfterPseudoMelaCut:       " << endl
       << nEvAfterPsMelaCut << "/ " << nEvAfterPsMelaCut_4e << "/ " << nEvAfterPsMelaCut_4mu << "/ " << nEvAfterPsMelaCut_2e2mu << endl
       << "nEvAfterGraviMelaCut:       " << endl
       << nEvAfterGrMelaCut << "/ " << nEvAfterGrMelaCut_4e << "/ " << nEvAfterGrMelaCut_4mu << "/ " << nEvAfterGrMelaCut_2e2mu << endl
       << "nEvWith1FSR: " << endl
       << nEvWith1FSRZ << "/ " << nEvWith1FSRZ_4e << "/ " << nEvWith1FSRZ_4mu << "/ " << nEvWith1FSRZ_2e2mu << endl
       << "nEvWith2FSR: " << endl
       << nEvWith2FSRZ<< "/ " << nEvWith2FSRZ_4e << "/ " << nEvWith2FSRZ_4mu << "/ " << nEvWith2FSRZ_2e2mu <<endl;



  // If events should be weighted
  // Be sure to add your weighted histograms here
  if( weightEvents )
  {
    if(histContainer_["extraParticles"]->GetEntries() > 0){histContainer_["extraParticles"]->Scale(scaleWeight);}
    if(histContainer_["minMass2l"]->GetEntries() > 0){histContainer_["minMass2l"]->Scale(scaleWeight);}
    if(histContainer_["minMass2l_SS_SF"]->GetEntries() > 0){histContainer_["minMass2l_SS_SF"]->Scale(scaleWeight);}
    if(histContainer_["minMass2l_OS_SF"]->GetEntries() > 0){histContainer_["minMass2l_OS_SF"]->Scale(scaleWeight);}
    if(histContainer_["minMass2l_SS_OF"]->GetEntries() > 0){histContainer_["minMass2l_SS_OF"]->Scale(scaleWeight);}
    if(histContainer_["minMass2l_OS_OF"]->GetEntries() > 0){histContainer_["minMass2l_OS_OF"]->Scale(scaleWeight);}
  }//weightEvents


  histContainer_["NEVENTS"]->SetBinContent(1,nEventsTotal);
  histContainer_["NEVENTS"]->GetXaxis()->SetBinLabel(1,"N Events in Sample");

  fillStepPlots();
   
  sipAna.plotSipHistograms(histContainer_);
  //sigEff->plotSigEffHistograms();
  isoEff->plotIsoHistograms(scaleWeight);
  //Zto4LAna.plotZto4LHistograms(histContainer_ , weightEvents,scaleWeight);
  //Zto4LAnaOP->plotZto4LHistograms(histContainer_ , weightEvents,scaleWeight);

  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@ END @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

}

void UFHZZ4LAna::beginRun(edm::Run const&, const edm::EventSetup& iSetup)
{
  massErr.init(iSetup);
}


// ------------ method called when ending the processing of a run  ------------
void  UFHZZ4LAna::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void UFHZZ4LAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void UFHZZ4LAna::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup)
{
  using namespace std;

  // Keep track of all the events run over
  edm::Handle<edm::MergeableCounter> numEventsCounter;
  lumiSeg.getByLabel("nEventsTotal", numEventsCounter);
  //if( isMC )
  //{
  //  if(numEventsCounter.isValid())
  //  {
  //    nEventsTotal += numEventsCounter->value;
  //  }
  //}

  edm::Handle<edm::MergeableCounter> numEventsSkimmedCounter;
  lumiSeg.getByLabel("nEvents2LSkim", numEventsSkimmedCounter);
  if(numEventsSkimmedCounter.isValid())
  {
    nEventsSkimmed += numEventsSkimmedCounter->value;
  }
 
  if(!isSignal)
  {
    edm::Handle<edm::MergeableCounter> numEventsAfterSkimCounter;
    lumiSeg.getByLabel("nEventsTriLep", numEventsAfterSkimCounter);

    if(numEventsAfterSkimCounter.isValid())
    {
      nEvAfterSkim += numEventsAfterSkimCounter->value;
    }
  }

}



// ============================ UF Functions ============================= //
//Find Z1,Z2, and Higgs candidate
//Pass good leptons for analysis as candMuons and candElectrons
//Pass empty vectors of pat leptons as selectedMuons and selectedElectrons
// these will be filled in the function and then useable for more analysis.
void UFHZZ4LAna::findHiggsCandidate(std::vector<pat::Muon> &candMuons, std::vector<pat::Electron> &candElectrons,
                           std::vector<pat::PFParticle> fsrPhotons, std::vector<double> deltaRVec,  // for miniAOD
                           std::vector<pat::Muon> &selectedMuons, std::vector<pat::Electron> &selectedElectrons, 
                           std::vector<pat::PFParticle> &selectedFsrPhotons,const edm::Event& iEvent )  // for miniAOD
{
  using namespace pat;
  using namespace std;

  bool Z1isMuons = false;
  bool Z1isElectrons = false;

  bool Z2isMuons = false;
  bool Z2isElectrons = false;

  double dm = 0;
  double ZmassDiff = 1000;
  const double Zmass = 91.1876;

  int nCandMuons = candMuons.size();
  int nCandElectrons = candElectrons.size();
  int takenMu_1 = 1000, takenMu_2 = 1000;
  int takenE_1 = 1000, takenE_2 = 1000;
  int takenZ1_1 = 1000, takenZ1_2 = 1000;
  int takenZ2_1 = 1000, takenZ2_2 = 1000;

  int takenMuTmp1 = 1000, takenMuTmp2 = 1000;
  int takenETmp1 = 1000, takenETmp2 = 1000;
  int takenPhotEl1 = 999, takenPhotEl2 = 999;
  int takenPhotMu1 = 999, takenPhotMu2 = 999;
  int takenPhotonZ1 = 999, takenPhotonZ2 = 999;
  int associatedMuPh1 = 999, associatedMuPh2 = 999;
  int associatedElPh1 = 999, associatedElPh2 = 999;
  int tmpAssociatedPh = 999;
  math::XYZTLorentzVector photVecZ1, photVecZ2, tmpPhotVec, tmpZVec;
  bool foundPhotZMu1 = false, foundPhotZEl1 = false;
  bool foundPhotZMu2 = false, foundPhotZEl2 = false;
  bool foundPhot = false;

  // Select Z1 by 2 leptons with same flavor opposite charge
  // with m(2l) closest to mZ
  for( int i = 0; i < nCandMuons; i++ )
  {
    for( int j = i; j < nCandMuons; j++ )
    {
      if( candMuons[i].charge() * candMuons[j].charge() == -1 )
      {
        bool foundZ1 = findZ(fsrPhotons, deltaRVec, candMuons[i], candMuons[j],999, 
                       takenPhotMu1, tmpAssociatedPh, tmpZVec, tmpPhotVec, foundPhot);
        dm = abs(Zmass-tmpZVec.M());
        if(dm < ZmassDiff && foundZ1)
        {
          Z1Vec = tmpZVec;
          ZmassDiff = dm;
          Z1isMuons = true;
          takenMu_1 = i; takenMu_2 = j;
          photVecZ1 = tmpPhotVec;
          foundPhotZMu1 = foundPhot;
          if(tmpAssociatedPh == 1) associatedMuPh1 = i;
          if(tmpAssociatedPh == 2) associatedMuPh1 = j;
          takenPhotonZ1 = takenPhotMu1;
          FSRPhot1_Px = photVecZ1.Px();
          FSRPhot1_Py = photVecZ1.Py();
          FSRPhot1_Pz = photVecZ1.Pz();
          FSRPhot1_E = photVecZ1.E();
        } // if (dm < ZmassDiff ...
      } // if( candMuons[i].charge() ...
    } // for (int j=i; ...
  } // for (int i=0; ...


  foundPhot = false;
  tmpAssociatedPh = 999;
  for( int i = 0; i < nCandElectrons; i++ )
  {
    for( int j = i; j < nCandElectrons; j++ )
    {
      if( candElectrons[i].charge() * candElectrons[j].charge() == -1 )
      {
        bool foundZ1 = findZ(fsrPhotons, deltaRVec, candElectrons[i], candElectrons[j],999, 
                             takenPhotEl1,tmpAssociatedPh, tmpZVec, tmpPhotVec,foundPhot);
        dm = abs(Zmass-tmpZVec.M());
        if(dm < ZmassDiff && foundZ1)
        {
          Z1Vec = tmpZVec;
          ZmassDiff = dm;
          Z1isMuons = false;
          Z1isElectrons = true;
          takenE_1 = i; takenE_2 = j;
          photVecZ1 = tmpPhotVec;
          foundPhotZEl1 = foundPhot;
          if(tmpAssociatedPh == 1) associatedElPh1 = i;
          if(tmpAssociatedPh == 2) associatedElPh1 = j;
          takenPhotonZ1 = takenPhotEl1;
          FSRPhot1_Px = photVecZ1.Px();
          FSRPhot1_Py = photVecZ1.Py();
          FSRPhot1_Pz = photVecZ1.Pz();
          FSRPhot1_E = photVecZ1.E();		   
        } // if (dm < ZmassDiff ..
      } // if( candElectrons[i]...
    }// for (int j=i; ..
  }//for( int i=0; ..
   
  // Keep track of whether Z1 is Muons or Electrons
  // Assign a tmp variable for pT comparisons
  if( Z1isMuons )
  {
    takenZ1_1   = takenMu_1; takenZ1_2   = takenMu_2;
    takenMuTmp1 = takenZ1_1; takenMuTmp2 = takenZ1_2;
    if( foundPhotZMu1 )
    {
      FSR_Z1 = true; 
      counterMuZ1ph++;
      if( associatedMuPh1 == takenMuTmp1 ) Lep1 = candMuons[takenMuTmp1].p4() + photVecZ1;
      else Lep1 = candMuons[takenMuTmp1].p4();
      if( associatedMuPh1 == takenMuTmp2 ) Lep2 = candMuons[takenMuTmp2].p4() + photVecZ1;
      else Lep2 = candMuons[takenMuTmp2].p4();
      FSRPhot1_Pt = photVecZ1.Pt();
      FSRPhot1_eta = photVecZ1.Eta();
      FSRPhot1_phi = photVecZ1.Phi();
      selectedFsrPhotons.push_back(fsrPhotons[takenPhotonZ1]);
    }
    else
    {
      Lep1 = candMuons[takenMuTmp1].p4();
      Lep2 = candMuons[takenMuTmp2].p4();
      FSRPhot1_Pt = -999;
      FSRPhot1_eta = -999;
      FSRPhot1_phi = -999;
    }
  }
  if( Z1isElectrons)
  {
    takenZ1_1  = takenE_1;  takenZ1_2  = takenE_2;
    takenETmp1 = takenZ1_1; takenETmp2 = takenZ1_2;
    if( foundPhotZEl1 )
    {
      FSR_Z1 = true;
      counterElZ1ph++;
      if( associatedElPh1 == takenETmp1 ) Lep1 = candElectrons[takenETmp1].p4() + photVecZ1;
      else Lep1 = candElectrons[takenETmp1].p4();
      if( associatedElPh1 == takenETmp2 ) Lep2 = candElectrons[takenETmp2].p4() + photVecZ1;
      else Lep2 = candElectrons[takenETmp2].p4();
      FSRPhot1_Pt = photVecZ1.Pt();
      FSRPhot1_eta = photVecZ1.Eta();
      FSRPhot1_phi = photVecZ1.Phi();
      selectedFsrPhotons.push_back(fsrPhotons[takenPhotonZ1]);
    }
    else
    {
      Lep1 = candElectrons[takenETmp1].p4();
      Lep2 = candElectrons[takenETmp2].p4();
      FSRPhot1_Pt = -999;
      FSRPhot1_eta = -999;
      FSRPhot1_phi = -999;
    }
  }

  if( Z1isMuons || Z1isElectrons)
  {
    nEvAfterZ1Formed += eventWeight;
    if(Z1isMuons)nEvAfterZ1Formed_2mu += eventWeight;
    if(Z1isElectrons)nEvAfterZ1Formed_2e += eventWeight;
    if((candMuons.size() + candElectrons.size()) >= 4)
    {
      nEv4GoodLep += eventWeight;
      if(candMuons.size() >= 2 && candElectrons.size() >= 2) nEv4GoodLep_2e2mu += eventWeight;
      if(candMuons.size() >= 4 ) nEv4GoodLep_4mu += eventWeight;
      if(candElectrons.size() >= 4) nEv4GoodLep_4e += eventWeight;
    }
  }

  if( Z1Vec.M() > mZ1Low && Z1Vec.M() < mZ1High && (Z1isMuons || Z1isElectrons) )
  {
    if(isSignal) sigEff_4->advanceSigNumCounters_MZ1(eventType, eventWeight);
    nEvAfterZ1Cut += eventWeight;
    if(Z1isMuons) nEvAfterZ1Cut_2mu += eventWeight;
    if(Z1isElectrons) nEvAfterZ1Cut_2e += eventWeight;
  }

  /////////////////////Z2////////////////////////////
  double sumPtZ2 = 0, sumPtZ2_tmp = 0;
  tmpAssociatedPh = 999;
  // Select Z1 by 2 leptons with same flavor opposite charge
  // with m(2l) closest to mZ
  for( int i = 0; i < nCandMuons; i++ )
  {
    for( int j = i; j < nCandMuons; j++ )
    {
      if( candMuons[i].charge() * candMuons[j].charge() == -1 )
      {
        if( i != takenMuTmp1 && i != takenMuTmp2 && j != takenMuTmp1 && j != takenMuTmp2 )
        {
          double tmpDeltaR1=999,tmpDeltaR2=999,tmpDeltaR3=999,tmpDeltaR4=999;
          if(Z1isMuons)
          {
            tmpDeltaR1 = deltaR(candMuons[i].eta(),candMuons[i].phi(),candMuons[takenMuTmp1].eta(),candMuons[takenMuTmp1].phi());
            tmpDeltaR2 = deltaR(candMuons[j].eta(),candMuons[j].phi(),candMuons[takenMuTmp1].eta(),candMuons[takenMuTmp1].phi());
            tmpDeltaR3 = deltaR(candMuons[i].eta(),candMuons[i].phi(),candMuons[takenMuTmp2].eta(),candMuons[takenMuTmp2].phi());
            tmpDeltaR4 = deltaR(candMuons[j].eta(),candMuons[j].phi(),candMuons[takenMuTmp2].eta(),candMuons[takenMuTmp2].phi());
          }
          if(tmpDeltaR1 < 0.02) continue;
          if(tmpDeltaR2 < 0.02) continue;
          if(tmpDeltaR3 < 0.02) continue;
          if(tmpDeltaR4 < 0.02) continue;

          bool foundZ2 = findZ(fsrPhotons, deltaRVec, candMuons[i], candMuons[j],takenPhotonZ1, 
                               takenPhotMu2, tmpAssociatedPh, tmpZVec, tmpPhotVec, foundPhot);
          sumPtZ2_tmp = candMuons[i].pt() + candMuons[j].pt();
		   
          if(sumPtZ2_tmp > sumPtZ2 && foundZ2)
          {
            Z2Vec = tmpZVec;
            Z2isMuons = true;
            if(candMuons[i].pt() > candMuons[j].pt()){takenZ2_1 = i; takenZ2_2 = j;}
            else{takenZ2_1 = j; takenZ2_2 = i;}
            foundPhotZMu2 = foundPhot;
            if(tmpAssociatedPh == 1) associatedMuPh2 = i;
            if(tmpAssociatedPh == 2) associatedMuPh2 = j;
            takenPhotonZ2 = takenPhotMu2;
            sumPtZ2 = sumPtZ2_tmp;
            photVecZ2 = tmpPhotVec;
            FSRPhot2_Px = photVecZ2.Px();
            FSRPhot2_Py = photVecZ2.Py();
            FSRPhot2_Pz = photVecZ2.Pz();
            FSRPhot2_E = photVecZ2.E();
          }
        }
      }
    }
  }

  foundPhot = false;
  for( int i = 0; i < nCandElectrons; i++ )
  {
    for( int j = i; j < nCandElectrons; j++ )
    {
      if( candElectrons[i].charge() * candElectrons[j].charge() == -1 )
      {
        if( i != takenETmp1 && i != takenETmp2 && j != takenETmp1 && j != takenETmp2 )
        {
          bool foundZ2 = findZ(fsrPhotons, deltaRVec, candElectrons[i], candElectrons[j],takenPhotonZ1,
                               takenPhotEl2, tmpAssociatedPh, tmpZVec, tmpPhotVec,foundPhot);
          sumPtZ2_tmp = candElectrons[i].pt() + candElectrons[j].pt();
          if(sumPtZ2_tmp > sumPtZ2 && foundZ2)
          {
            Z2Vec = tmpZVec;
            Z2isMuons = false;
            Z2isElectrons = true;
            if(candElectrons[i].pt() > candElectrons[j].pt()){takenZ2_1 = i; takenZ2_2 = j;}
            else{takenZ2_1 = j; takenZ2_2 = i;}
            photVecZ2 = tmpPhotVec;
            foundPhotZEl2 = foundPhot;
            if(tmpAssociatedPh == 1) associatedElPh2 = i;
            if(tmpAssociatedPh == 2) associatedElPh2 = j;
            takenPhotonZ2 = takenPhotEl2;
            sumPtZ2 = sumPtZ2_tmp;
            photVecZ2 = tmpPhotVec;
            FSRPhot2_Px = photVecZ2.Px();
            FSRPhot2_Py = photVecZ2.Py();
            FSRPhot2_Pz = photVecZ2.Pz();
            FSRPhot2_E = photVecZ2.E();
          }
        }
      }
    }
  }

  if( (Z2isMuons || Z2isElectrons) && (Z1Vec.M() > mZ1Low && Z1Vec.M() < mZ1High) )
  {
    nEvAfterZ2Formed += eventWeight;
    if(Z1isMuons && Z2isMuons) nEvAfterZ2Formed_4mu += eventWeight;
    if(Z1isElectrons && Z2isElectrons) nEvAfterZ2Formed_4e += eventWeight;
    if((Z1isMuons && Z2isElectrons) || (Z1isElectrons && Z2isMuons)) nEvAfterZ2Formed_2e2mu += eventWeight;
  }
   
  if(Z2isMuons)
  {
    if( foundPhotZMu2 )
    {
      FSR_Z2 = true;
      counterMuZ2ph++;
      if( associatedMuPh2 == takenZ2_1 ) Lep3 = candMuons[takenZ2_1].p4() + photVecZ2;
      else Lep3 = candMuons[takenZ2_1].p4();
      if( associatedMuPh2 == takenZ2_2 ) Lep4 = candMuons[takenZ2_2].p4() + photVecZ2;
      else Lep4 = candMuons[takenZ2_2].p4();
      FSRPhot2_Pt = photVecZ2.Pt();
      FSRPhot2_eta = photVecZ2.Eta();
      FSRPhot2_phi = photVecZ2.Phi();
      selectedFsrPhotons.push_back(fsrPhotons[takenPhotonZ2]);
    }
    else
    {
      Lep3 = candMuons[takenZ2_1].p4();
      Lep4 = candMuons[takenZ2_2].p4();
      FSRPhot2_Pt = -1;
      FSRPhot2_eta = -1;
      FSRPhot2_phi = -1;
    }
  }
  if(Z2isElectrons)
  {
    if( foundPhotZEl2 )
    {
      FSR_Z2 = true;
      counterElZ2ph++;
      if( associatedElPh2 == takenZ2_1 ) Lep3 = candElectrons[takenZ2_1].p4() + photVecZ2;
      else Lep3 = candElectrons[takenZ2_1].p4();
      if( associatedElPh2 == takenZ2_2 ) Lep4 = candElectrons[takenZ2_2].p4() + photVecZ2;
      else Lep4 = candElectrons[takenZ2_2].p4();
      FSRPhot2_Pt = photVecZ2.Pt();
      FSRPhot2_eta = photVecZ2.Eta();
      FSRPhot2_phi = photVecZ2.Phi();
      selectedFsrPhotons.push_back(fsrPhotons[takenPhotonZ2]);
    }
    else
    {
      Lep3 = candElectrons[takenZ2_1].p4();
      Lep4 = candElectrons[takenZ2_2].p4();
      FSRPhot2_Pt = -1;
      FSRPhot2_eta = -1;
      FSRPhot2_phi = -1;
    }
  }

  //Determine whether a Higgs candidate was formed
  if( Z1isMuons == true && Z2isMuons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_4mu += eventWeight;

    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High )
    {
      if( mZ2 > mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco4mu",eventWeight);
        nEvAfterZ2Cut += eventWeight;
        nEvAfterZ2Cut_4mu += eventWeight;
      }
    }

    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoFourMuEvent = true;
      selectedMuons.push_back(candMuons[takenZ1_1]);
      selectedMuons.push_back(candMuons[takenZ1_2]);
      selectedMuons.push_back(candMuons[takenZ2_1]);
      selectedMuons.push_back(candMuons[takenZ2_2]);
    }
  }
  else if( Z1isMuons == true && Z2isElectrons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_2e2mu += eventWeight;
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High)
    {
      if( mZ2 >mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco2e2mu",eventWeight);
        nEvAfterZ2Cut += eventWeight;
        nEvAfterZ2Cut_2e2mu += eventWeight;
      }
    }

    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoTwoMuTwoEEvent = true;
      selectedMuons.push_back(candMuons[takenZ1_1]);
      selectedMuons.push_back(candMuons[takenZ1_2]);
      selectedElectrons.push_back(candElectrons[takenZ2_1]);
      selectedElectrons.push_back(candElectrons[takenZ2_2]);
    }
  }
  else if( Z1isElectrons == true && Z2isMuons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_2e2mu += eventWeight;
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High )
    {
      if( mZ2 >mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco2e2mu",eventWeight);
        nEvAfterZ2Cut += eventWeight;
        nEvAfterZ2Cut_2e2mu += eventWeight;
      }
    }
    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoTwoETwoMuEvent = true;
      selectedMuons.push_back(candMuons[takenZ2_1]);
      selectedMuons.push_back(candMuons[takenZ2_2]);
      selectedElectrons.push_back(candElectrons[takenZ1_1]);
      selectedElectrons.push_back(candElectrons[takenZ1_2]);
    }
  }
  else if( Z1isElectrons == true && Z2isElectrons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_4e += eventWeight;
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High)
    {
      if( mZ2 >mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco4e",eventWeight);
        nEvAfterZ2Cut += eventWeight;
        nEvAfterZ2Cut_4e += eventWeight;
      }
    }
    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoFourEEvent = true;
      selectedElectrons.push_back(candElectrons[takenZ1_1]);
      selectedElectrons.push_back(candElectrons[takenZ1_2]);
      selectedElectrons.push_back(candElectrons[takenZ2_1]);
      selectedElectrons.push_back(candElectrons[takenZ2_2]);
    }
  }

  // If a Higgs candidate is formed, save its information
  if( foundHiggsCandidate )
  {
    //Impose tight lepton requirements if mZ is far off mass peak
    if( RecoFourMuEvent && !looseIdsOnly )
    {
      if( mZ1 > (Zmass+20) || mZ1 < (Zmass-20) )
      {
        vector<pat::Muon> muonsZ1;
        muonsZ1.push_back(selectedMuons[0]);
        muonsZ1.push_back(selectedMuons[1]);
        vector<pat::Muon> passedMuonsZ1 = helper.goodTightMuons(muonsZ1,_muPtCut);
        if( passedMuonsZ1.size() != muonsZ1.size() ){ foundHiggsCandidate = false; }
      }
      if( mZ2 > (Zmass+20) || mZ2 < (Zmass-20) )
      {
        vector<pat::Muon> muonsZ2;
        muonsZ2.push_back(selectedMuons[2]);
        muonsZ2.push_back(selectedMuons[3]);
        vector<pat::Muon> passedMuonsZ2 = helper.goodTightMuons(muonsZ2,_muPtCut);
        if( passedMuonsZ2.size() != muonsZ2.size() ){ foundHiggsCandidate = false; }
      }
    }
    if( RecoFourEEvent && !looseIdsOnly)
    {
      if( mZ1 > (Zmass+20) || mZ1 < (Zmass-20) )
      {
        vector<pat::Electron> elecsZ1;
        elecsZ1.push_back(selectedElectrons[0]);
        elecsZ1.push_back(selectedElectrons[1]);
        vector<pat::Electron> passedElecsZ1 = helper.goodTightElectrons(elecsZ1,_elecPtCut,elecID);
        if( passedElecsZ1.size() != elecsZ1.size() ){ foundHiggsCandidate = false; }
      }
      if( mZ2 > (Zmass+20) || mZ2 < (Zmass-20) )
      {
        vector<pat::Electron> elecsZ2;
        elecsZ2.push_back(selectedElectrons[2]);
        elecsZ2.push_back(selectedElectrons[3]);
        vector<pat::Electron> passedElecsZ2 = helper.goodTightElectrons(elecsZ2,_elecPtCut,elecID);
        if( passedElecsZ2.size() != elecsZ2.size() ){ foundHiggsCandidate = false; }
      }
    }
    if( RecoTwoETwoMuEvent && !looseIdsOnly)
    {
      if( mZ1 > (Zmass+20) || mZ1 < (Zmass-20) )
      {
        vector<pat::Electron> elecsZ1;
        elecsZ1.push_back(selectedElectrons[0]);
        elecsZ1.push_back(selectedElectrons[1]);
        vector<pat::Electron> passedElecsZ1 = helper.goodTightElectrons(elecsZ1,_elecPtCut,elecID);
        if( passedElecsZ1.size() != elecsZ1.size() ){ foundHiggsCandidate = false; }
      }
      if( mZ2 > (Zmass+20) || mZ2 < (Zmass-20) )
      {
        vector<pat::Muon> muonsZ2;
        muonsZ2.push_back(selectedMuons[0]);
        muonsZ2.push_back(selectedMuons[1]);
        vector<pat::Muon> passedMuonsZ2 = helper.goodTightMuons(muonsZ2,_muPtCut);
        if( passedMuonsZ2.size() != muonsZ2.size() ){ foundHiggsCandidate = false; }
      }
    }
    if( RecoTwoMuTwoEEvent && !looseIdsOnly)
    {
      if( mZ1 > (Zmass+20) || mZ1 < (Zmass-20) )
      {
        vector<pat::Muon> muonsZ1;
        muonsZ1.push_back(selectedMuons[0]);
        muonsZ1.push_back(selectedMuons[1]);
        vector<pat::Muon> passedMuonsZ1 = helper.goodTightMuons(muonsZ1,_muPtCut);
        if( passedMuonsZ1.size() != muonsZ1.size() ){ foundHiggsCandidate = false; }
      }
      if( mZ2 > (Zmass+20) || mZ2 < (Zmass-20) )
      {
        vector<pat::Electron> elecsZ2;
        elecsZ2.push_back(selectedElectrons[0]);
        elecsZ2.push_back(selectedElectrons[1]);
        vector<pat::Electron> passedElecsZ2 = helper.goodTightElectrons(elecsZ2,_elecPtCut,elecID);
        if( passedElecsZ2.size() != elecsZ2.size() ){ foundHiggsCandidate = false; }
      }
    }
    if( foundHiggsCandidate ) 
    {
      HiggsCandVec = Z1Vec + Z2Vec;
      m4l = HiggsCandVec.M();
    }
  }

}


//------------------------------------------------------------
// findHiggsCandidate without flavour requirements
//------------------------------------------------------------
void UFHZZ4LAna::findHiggsCandidate_MixFlavour(std::vector<pat::Muon> &candMuons, std::vector<pat::Electron> &candElectrons,
                                          std::vector<pat::Muon> &selectedMuons, std::vector<pat::Electron> &selectedElectrons, bool noChargeReq)
{
  using namespace pat;
  using namespace std;
    
  bool Z1isMuons = false;
  bool Z1isElectrons = false;
  bool Z1isEMu = false;
    
  bool Z2isMuons = false;
  bool Z2isElectrons = false;
  bool Z2isEMu = false;
  bool Z2isMuE = false;
    
  double dm = 0;
  double ZmassDiff = 1000;
  const double Zmass = 91.1876;
    
  int nCandMuons = candMuons.size();
  int nCandElectrons = candElectrons.size();
  int takenMu_1 = 1000, takenMu_2 = 1000;
  int takenE_1 = 1000, takenE_2 = 1000;
  int takenZ1_1 = 1000, takenZ1_2 = 1000;
  int takenZ2_1 = 1000, takenZ2_2 = 1000;
    
  int takenMuTmp1 = 1000, takenMuTmp2 = 1000;
  int takenETmp1 = 1000, takenETmp2 = 1000;
    
  // Select Z1 by 2 leptons with opposite charge (unless noChargeReq=TRUE)
  // with m(2l) closest to mZ
  for( int i = 0; i < nCandMuons; i++ )
  {
    for( int j = i; j < nCandMuons; j++ )
    {
      if(( candMuons[i].charge() * candMuons[j].charge() == -1 ) || noChargeReq)
      {
        double MuinvMass_Zm = (candMuons[i].p4()+candMuons[j].p4()).M();
        dm = abs(Zmass-MuinvMass_Zm);
        if(dm < ZmassDiff)
        {
          Z1Vec = candMuons[i].p4() + candMuons[j].p4();
          ZmassDiff = dm;
          Z1isMuons = true;
          takenMu_1 = i; takenMu_2 = j;
        }
      }
    }
  }
    
  for( int i = 0; i < nCandElectrons; i++ )
  {
    for( int j = i; j < nCandElectrons; j++ )
    {
      if(( candElectrons[i].charge() * candElectrons[j].charge() == -1 ) || noChargeReq)
      {
        double EinvMass_Zm = (candElectrons[i].p4()+candElectrons[j].p4()).M();
        dm = abs(Zmass-EinvMass_Zm);
        if(dm < ZmassDiff)
        {
          Z1Vec = candElectrons[i].p4() + candElectrons[j].p4();
          ZmassDiff = dm;
          Z1isMuons = false;
          Z1isElectrons = true;
          takenE_1 = i; takenE_2 = j;
        }
      }
    }
  }
    
  for( int i = 0; i < nCandMuons; i++ )
  {
    for( int j = i; j < nCandElectrons; j++ )
    {
      if(( candMuons[i].charge() * candElectrons[j].charge() == -1 ) || noChargeReq)
      {
        double EMuinvMass_Zm = (candMuons[i].p4()+candElectrons[j].p4()).M();
        dm = abs(Zmass-EMuinvMass_Zm);
        if(dm < ZmassDiff)
        {
          Z1Vec = candMuons[i].p4() + candElectrons[j].p4();
          ZmassDiff = dm;
          Z1isMuons = false;
          Z1isElectrons = false;
          Z1isEMu = true;
          takenMu_1 = i; takenE_1 = j;
        }
      }
    }
  }
    
  // Keep track of whether Z1 is 2xMuons, 2xElectrons or combination Elec-Muon
  // Assign a tmp variable for pT comparisons
  if( Z1isMuons )
  {
    takenZ1_1   = takenMu_1; takenZ1_2   = takenMu_2;
    takenMuTmp1 = takenZ1_1; takenMuTmp2 = takenZ1_2;
  }
  if( Z1isElectrons)
  {
    takenZ1_1  = takenE_1;  takenZ1_2  = takenE_2;
    takenETmp1 = takenZ1_1; takenETmp2 = takenZ1_2;
  }
  if( Z1isEMu)
  {
    takenZ1_1  = takenMu_1;  takenZ1_2  = takenE_1;
    takenMuTmp1 = takenZ1_1; takenETmp1 = takenZ1_2;
  }
   
  // Variables for finding Z2
  double biggestPt_Muons = 0;
  double biggestPt_Electrons = 0;
  int takenPt_Mu = 1000;
  int takenPt_E = 1000;
  
  int takenPt_Index1 = 1000;
  double biggestPt2 = 0;
  double biggestPt3 = 0;
  int takenPt_Index2 = 1000;
  
  bool highestPtisMuon = false;
  bool highestPtisElectron = false;
   
  // Find highest pT muon and electron that remains
  for( int i = 0; i < nCandMuons; i++)
  {
    if( i != takenMuTmp1 && i != takenMuTmp2 )
    {
      if( candMuons[i].pt() > biggestPt_Muons )
      {
        takenPt_Mu = i;
        biggestPt_Muons = candMuons[i].pt();
      }
    }
  }
    
  for( int j = 0; j < nCandElectrons; j++)
  {
    if( j != takenETmp1 && j != takenETmp2 )
    {
      if( candElectrons[j].pt() > biggestPt_Electrons )
      {
        takenPt_E = j;
        biggestPt_Electrons = candElectrons[j].pt();
      }
    }
  }
    
    
  // Compare highest pT muon and electron
  if( biggestPt_Muons > biggestPt_Electrons )
  {
    highestPtisMuon = true;
    takenPt_Index1 = takenPt_Mu;
  }
  else if( biggestPt_Electrons > biggestPt_Muons )
  {
    highestPtisElectron = true;
    takenPt_Index1 = takenPt_E;
  }
   
  //If the highest pT lepton left is mu 
  if( highestPtisMuon )
  {
    for( int i = 0; i < nCandMuons; i++)
    {
      if( i != takenPt_Index1 && i != takenMuTmp1 && i != takenMuTmp2 )
      {
        if(( candMuons[takenPt_Index1].charge() * candMuons[i].charge() == -1 ) || noChargeReq)
        {
          if( candMuons[i].pt() > biggestPt2 ){ takenPt_Index2 = i; biggestPt2 = candMuons[i].pt();}
        }
      }
    }
    if( biggestPt2 > 0 && ((candMuons[takenPt_Index1].charge() * candMuons[takenPt_Index2].charge() == -1) || noChargeReq) )
    {
      Z2Vec = candMuons[takenPt_Index1].p4() + candMuons[takenPt_Index2].p4();
      takenZ2_1 = takenPt_Index1; takenZ2_2 = takenPt_Index2;
      Z2isMuons = true;
    }
    for( int i = 0; i < nCandElectrons; i++)
    {
      if( i != takenETmp1 && i != takenETmp2 )
      {
        if(( candMuons[takenPt_Index1].charge() * candElectrons[i].charge() == -1 ) || noChargeReq)
        {
          if( candElectrons[i].pt() > biggestPt3 ){ takenPt_Index2 = i; biggestPt3 = candElectrons[i].pt();}
        }
      }
    }
    if( biggestPt3 > biggestPt2 && ((candMuons[takenPt_Index1].charge() * candElectrons[takenPt_Index2].charge() == -1) || noChargeReq) )
    {
      Z2Vec = candMuons[takenPt_Index1].p4() + candElectrons[takenPt_Index2].p4();
      takenZ2_1 = takenPt_Index1; takenZ2_2 = takenPt_Index2;
      Z2isMuons = false;
      Z2isMuE = true;
    }
  }
    
  //If the highest pT lepton left is e
  if( highestPtisElectron )
  {
    for( int i = 0; i < nCandElectrons; i++)
    {
      if( i != takenPt_Index1 && i != takenETmp1 && i != takenETmp2 )
      {
        if(( candElectrons[takenPt_Index1].charge() * candElectrons[i].charge() == -1 ) || noChargeReq)
        {
          if( candElectrons[i].pt() > biggestPt2 ){ takenPt_Index2 = i; biggestPt2 = candElectrons[i].pt();}
        }
      }
    }
    if( biggestPt2 > 0 && ((candElectrons[takenPt_Index1].charge() * candElectrons[takenPt_Index2].charge() == -1) || noChargeReq) )
    {
      Z2Vec = candElectrons[takenPt_Index1].p4() + candElectrons[takenPt_Index2].p4();
      takenZ2_1 = takenPt_Index1; takenZ2_2 = takenPt_Index2;
      Z2isElectrons = true;
    }
    for( int i = 0; i < nCandMuons; i++)
    {
      if( i != takenMuTmp1 && i != takenMuTmp2 )
      {
        if(( candElectrons[takenPt_Index1].charge() * candMuons[i].charge() == -1 ) || noChargeReq)
        {
          if( candMuons[i].pt() > biggestPt3 ){ takenPt_Index2 = i; biggestPt3 = candMuons[i].pt();}
        }
      }
    }
    if( biggestPt3 > biggestPt2 && ((candElectrons[takenPt_Index1].charge() * candMuons[takenPt_Index2].charge() == -1) || noChargeReq) )
    {
      Z2Vec = candElectrons[takenPt_Index1].p4() + candMuons[takenPt_Index2].p4();
      takenZ2_1 = takenPt_Index1; takenZ2_2 = takenPt_Index2;
      Z2isElectrons = false;
      Z2isEMu = true;
    }
  }
  
  //Determine whether a Higgs candidate was formed
  if( Z1isMuons == true && Z2isMuons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_4mu += eventWeight;
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High )
    {
      //nEvAfterZ1Cut += eventWeight;
      //nEvAfterZ1Cut_4mu += eventWeight;
      if( mZ2 > mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco4mu",eventWeight);
        nEvAfterZ2Cut += eventWeight;
        nEvAfterZ2Cut_4mu += eventWeight;
      }
    }
        
    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoFourMuEvent = true;
      selectedMuons.push_back(candMuons[takenZ1_1]);
      selectedMuons.push_back(candMuons[takenZ1_2]);
      selectedMuons.push_back(candMuons[takenZ2_1]);
      selectedMuons.push_back(candMuons[takenZ2_2]);
    }
  }
  else if( Z1isMuons == true && Z2isElectrons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_2e2mu += eventWeight;
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High)
    {
      //nEvAfterZ1Cut += eventWeight;
      //nEvAfterZ1Cut_2e2mu += eventWeight;
      if( mZ2 >mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco2e2mu",eventWeight);
        nEvAfterZ2Cut += eventWeight;
        nEvAfterZ2Cut_2e2mu += eventWeight;
      }
    }
    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoTwoMuTwoEEvent = true;
      selectedMuons.push_back(candMuons[takenZ1_1]);
      selectedMuons.push_back(candMuons[takenZ1_2]);
      selectedElectrons.push_back(candElectrons[takenZ2_1]);
      selectedElectrons.push_back(candElectrons[takenZ2_2]);
    }
  }    
  else if( Z1isElectrons == true && Z2isMuons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_2e2mu += eventWeight;
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High )
    {
      //nEvAfterZ1Cut += eventWeight;
      //nEvAfterZ1Cut_2e2mu += eventWeight;
      if( mZ2 >mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco2e2mu",eventWeight);
        nEvAfterZ2Cut += eventWeight;
        nEvAfterZ2Cut_2e2mu += eventWeight;
      }
    }
    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoTwoETwoMuEvent = true;
      selectedMuons.push_back(candMuons[takenZ2_1]);
      selectedMuons.push_back(candMuons[takenZ2_2]);
      selectedElectrons.push_back(candElectrons[takenZ1_1]);
      selectedElectrons.push_back(candElectrons[takenZ1_2]);
    }
  }
  else if( Z1isElectrons == true && Z2isElectrons == true )
  {
    nEvBeforeZCuts += eventWeight;
    nEvBeforeZCuts_4e += eventWeight;
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( mZ1 > mZ1Low && mZ1 < mZ1High)
    {
      //nEvAfterZ1Cut += eventWeight;
      //nEvAfterZ1Cut_4e += eventWeight;
      if( mZ2 >mZ2Low && mZ2 < mZ2High)
      {
        if(isSignal) sigEff_4->advanceSigNumCounters_MZ2(eventType, "reco4e",eventWeight);
        nEvAfterZ2Cut += eventWeight;  
        nEvAfterZ2Cut_4e += eventWeight;
      }
    }
    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      RecoFourEEvent = true;
      selectedElectrons.push_back(candElectrons[takenZ1_1]);
      selectedElectrons.push_back(candElectrons[takenZ1_2]);
      selectedElectrons.push_back(candElectrons[takenZ2_1]);
      selectedElectrons.push_back(candElectrons[takenZ2_2]);
    }
  }
  else if( (Z1isEMu && (Z2isElectrons || Z2isMuons || Z2isMuE || Z2isEMu)) || ((Z1isElectrons || Z1isMuons || Z1isEMu) && (Z2isEMu || Z2isMuE)) )
  {
    mZ1 = Z1Vec.M();
    mZ2 = Z2Vec.M();
    if( (mZ1 > mZ1Low && mZ1 < mZ1High) && (mZ2 > mZ2Low && mZ2 < mZ2High) )
    {
      foundHiggsCandidate = true;
      if( Z1isEMu && Z2isEMu )
      {
        selectedMuons.push_back(candMuons[takenZ1_1]);
        selectedElectrons.push_back(candElectrons[takenZ1_2]);
        selectedElectrons.push_back(candElectrons[takenZ2_1]);
        selectedMuons.push_back(candMuons[takenZ2_2]);
        RecoFourMixEvent = 1;
      }
      if( Z1isEMu && Z2isMuE )
      {
        selectedMuons.push_back(candMuons[takenZ1_1]);
        selectedElectrons.push_back(candElectrons[takenZ1_2]);
        selectedMuons.push_back(candMuons[takenZ2_1]);
        selectedElectrons.push_back(candElectrons[takenZ2_2]);
        RecoFourMixEvent = 2;
      }
      if( Z1isEMu && Z2isElectrons )
      {
        selectedMuons.push_back(candMuons[takenZ1_1]);
        selectedElectrons.push_back(candElectrons[takenZ1_2]);
        selectedElectrons.push_back(candElectrons[takenZ2_1]);
        selectedElectrons.push_back(candElectrons[takenZ2_2]);
        RecoFourMixEvent = 3;
      }
      if( Z1isEMu && Z2isMuons )
      {
        selectedMuons.push_back(candMuons[takenZ1_1]);
        selectedElectrons.push_back(candElectrons[takenZ1_2]);
        selectedMuons.push_back(candMuons[takenZ2_1]);
        selectedMuons.push_back(candMuons[takenZ2_2]);
        RecoFourMixEvent = 4;
      }
      if( Z1isElectrons && Z2isEMu )
      {
        selectedElectrons.push_back(candElectrons[takenZ1_1]);
        selectedElectrons.push_back(candElectrons[takenZ1_2]);
        selectedElectrons.push_back(candElectrons[takenZ2_1]);
        selectedMuons.push_back(candMuons[takenZ2_2]);
        RecoFourMixEvent = 5;
      }
      if( Z1isElectrons && Z2isMuE )
      {
        selectedElectrons.push_back(candElectrons[takenZ1_1]);
        selectedElectrons.push_back(candElectrons[takenZ1_2]);
        selectedMuons.push_back(candMuons[takenZ2_1]);
        selectedElectrons.push_back(candElectrons[takenZ2_2]);
        RecoFourMixEvent = 6;
      }
      if( Z1isMuons && Z2isEMu )
      {
        selectedMuons.push_back(candMuons[takenZ1_1]);
        selectedMuons.push_back(candMuons[takenZ1_2]);
        selectedElectrons.push_back(candElectrons[takenZ2_1]);
        selectedMuons.push_back(candMuons[takenZ2_2]);
        RecoFourMixEvent = 7;
      }
      if( Z1isMuons && Z2isMuE )
      {
        selectedMuons.push_back(candMuons[takenZ1_1]);
        selectedMuons.push_back(candMuons[takenZ1_2]);
        selectedMuons.push_back(candMuons[takenZ2_1]);
        selectedElectrons.push_back(candElectrons[takenZ2_2]);
        RecoFourMixEvent = 8;
      }
    }
  }
    
  // If a Higgs candidate is formed, save its information
  if( foundHiggsCandidate )
  {
    HiggsCandVec = Z1Vec + Z2Vec;
    m4l = HiggsCandVec.M();
  }
}

double UFHZZ4LAna::getMinDeltaR(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons)
{

  using namespace pat;
  using namespace std;

  double minDeltaR = 1000;
  vector<double> deltaR_, deltR_emu;
  double tmpDeltaR = 1000;
  
  for( unsigned int i = 0; i < Muons.size(); i++ )
  {
    for( unsigned int j = i; j < Muons.size(); j++ )
    {
      if( i != j )
      {
        tmpDeltaR = deltaR(Muons[i].eta(),Muons[i].phi(),Muons[j].eta(),Muons[j].phi());
        deltaR_.push_back(tmpDeltaR);
      }
    }
  }
      
  for( unsigned int i = 0; i < Electrons.size(); i++ )
  {
    for( unsigned int j = i; j < Electrons.size(); j++ )
    {
      if( i != j )
      {
        tmpDeltaR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Electrons[j].eta(),Electrons[j].phi());
        deltaR_.push_back(tmpDeltaR);
      }
    }
  }

  for( unsigned int i = 0; i < Electrons.size(); i++ )
  {
    for( unsigned int j = 0; j < Muons.size(); j++ )
    {
      tmpDeltaR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Muons[j].eta(),Muons[j].phi());
      deltaR_.push_back(tmpDeltaR);
    }
  }

  for( unsigned int i = 0; i < deltaR_.size(); i++)
  {
    if( deltaR_[i] < minDeltaR ){ minDeltaR = deltaR_[i]; }
  }

  return minDeltaR;
}


void UFHZZ4LAna::plotMinDeltaRemu(std::vector<pat::Muon> Muons, std::vector<pat::Electron> Electrons)
{
  using namespace pat;
  using namespace std;

  double minDeltaR = 1000;
  vector<double> deltaR_emu;
  double tmpDeltaR = 1000;

  for( unsigned int i = 0; i < Electrons.size(); i++ )
  {
    for( unsigned int j = 0; j < Muons.size(); j++ )
    {
      tmpDeltaR = deltaR(Electrons[i].eta(), Electrons[i].phi(), Muons[j].eta(),Muons[j].phi());
      deltaR_emu.push_back(tmpDeltaR);
    }
  }

  for( unsigned int i = 0; i < deltaR_emu.size(); i++)
  {
    if( deltaR_emu[i] < minDeltaR ){ minDeltaR = deltaR_emu[i]; }
  }

  //histContainer_["minDeltaRemu"]->Fill(minDeltaR,eventWeight);

}


void UFHZZ4LAna::bookStepPlots()
{
  int nSteps = 9;
  histContainer_["stepPlot_w"]=fs->make<TH1F>("stepPlotWeighted", "Events After Cuts in H->ZZ->4l; ;N Events", nSteps, 0, nSteps);
  histContainer_["stepPlot_w_4mu"]=fs->make<TH1F>("stepPlotWeighted_4mu", "Events After Cuts in H->ZZ->4mu; ;N Events", nSteps, 0, nSteps);
  histContainer_["stepPlot_w_4e"]=fs->make<TH1F>("stepPlotWeighted_4e", "Events After Cuts in H->ZZ->4e; ;N Events", nSteps, 0, nSteps);
  histContainer_["stepPlot_w_2e2mu"]=fs->make<TH1F>("stepPlotWeighted_2e2mu", "Events After Cuts in H->ZZ->2e2mu; ;N Events", nSteps, 0, nSteps);
}


void UFHZZ4LAna::fillStepPlots()
{
  TString skim  = "Skim";
  TString trig = "HLT";
  TString ID = "2L ID";
  TString Z1Cut = "40 < mZ1 < 120";
  TString Z2Cut = "12 < mZ2 < 120";
  TString twoLepPt = "2L pT(20,10)";
  TString cleaning = "All pairs m2l > 4";
  TString m4l70cut = "m(4l) > 70";
  TString m4l100cut = "m(4l) > 100";

  if(nEvAfterSkim == 0){ nEvAfterSkim = nEvPassedHlt;}
  if( isMC ){ nEvAfterSkim = nEvPassedHlt;}
  if(!isMC && doBlinding)
  {
    nEvAfterZ4lCut = -1;
    nEvAfterM4lCut = -1;
    nEvAfterZ4lCut_4mu = -1;
    nEvAfterM4lCut_4mu = -1;
    nEvAfterZ4lCut_4e = -1;
    nEvAfterM4lCut_4e = -1;
    nEvAfterZ4lCut_2e2mu = -1;
    nEvAfterM4lCut_2e2mu = -1;
  }

  //4L
  histContainer_["stepPlot_w"]->SetBinContent(1,nEvAfterSkim*scaleWeight);
  histContainer_["stepPlot_w"]->SetBinContent(2,nEvPassedHlt*scaleWeight);  
  histContainer_["stepPlot_w"]->SetBinContent(3,nEvAfterId*scaleWeight);
  histContainer_["stepPlot_w"]->SetBinContent(4,nEvAfterZ1Cut*scaleWeight);
  histContainer_["stepPlot_w"]->SetBinContent(5,nEvAfterZ2Cut*scaleWeight);
  histContainer_["stepPlot_w"]->SetBinContent(6,nEvPassedPtCut2*scaleWeight);
  histContainer_["stepPlot_w"]->SetBinContent(7,nEvAfterCleaning*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w"]->SetBinContent(8,nEvAfterZ4lCut*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w"]->SetBinContent(9,nEvAfterM4lCut*scaleWeight);

  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(1,skim);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(2,trig);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(3,ID);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(4,Z1Cut);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(5,Z2Cut);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(6,twoLepPt);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(7,cleaning);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(8,m4l70cut);
  histContainer_["stepPlot_w"]->GetXaxis()->SetBinLabel(9,m4l100cut);

  //4mu
  histContainer_["stepPlot_w_4mu"]->SetBinContent(1,nEvAfterSkim*scaleWeight);
  histContainer_["stepPlot_w_4mu"]->SetBinContent(2,nEvPassedHlt*scaleWeight);  
  histContainer_["stepPlot_w_4mu"]->SetBinContent(3,nEvAfterId*scaleWeight);
  histContainer_["stepPlot_w_4mu"]->SetBinContent(4,nEvAfterZ1Cut*scaleWeight);
  histContainer_["stepPlot_w_4mu"]->SetBinContent(5,nEvAfterZ2Cut_4mu*scaleWeight);
  histContainer_["stepPlot_w_4mu"]->SetBinContent(6,nEvPassedPtCut2_4mu*scaleWeight);
  histContainer_["stepPlot_w_4mu"]->SetBinContent(7,nEvAfterCleaning_4mu*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w_4mu"]->SetBinContent(8,nEvAfterZ4lCut_4mu*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w_4mu"]->SetBinContent(9,nEvAfterM4lCut_4mu*scaleWeight);

  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(1,skim);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(2,trig);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(3,ID);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(4,Z1Cut);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(5,Z2Cut);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(6,twoLepPt);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(7,cleaning);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(8,m4l70cut);
  histContainer_["stepPlot_w_4mu"]->GetXaxis()->SetBinLabel(9,m4l100cut);

  //4e
  histContainer_["stepPlot_w_4e"]->SetBinContent(1,nEvAfterSkim*scaleWeight);
  histContainer_["stepPlot_w_4e"]->SetBinContent(2,nEvPassedHlt*scaleWeight);  
  histContainer_["stepPlot_w_4e"]->SetBinContent(3,nEvAfterId*scaleWeight);
  histContainer_["stepPlot_w_4e"]->SetBinContent(4,nEvAfterZ1Cut*scaleWeight);
  histContainer_["stepPlot_w_4e"]->SetBinContent(5,nEvAfterZ2Cut_4e*scaleWeight);
  histContainer_["stepPlot_w_4e"]->SetBinContent(6,nEvPassedPtCut2_4e*scaleWeight);
  histContainer_["stepPlot_w_4e"]->SetBinContent(7,nEvAfterCleaning_4e*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w_4e"]->SetBinContent(8,nEvAfterZ4lCut_4e*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w_4e"]->SetBinContent(9,nEvAfterM4lCut_4e*scaleWeight);

  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(1,skim);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(2,trig);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(3,ID);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(4,Z1Cut);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(5,Z2Cut);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(6,twoLepPt);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(7,cleaning);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(8,m4l70cut);
  histContainer_["stepPlot_w_4e"]->GetXaxis()->SetBinLabel(9,m4l100cut);

  //2e2mu
  histContainer_["stepPlot_w_2e2mu"]->SetBinContent(1,nEvAfterSkim*scaleWeight);
  histContainer_["stepPlot_w_2e2mu"]->SetBinContent(2,nEvPassedHlt*scaleWeight);  
  histContainer_["stepPlot_w_2e2mu"]->SetBinContent(3,nEvAfterId*scaleWeight);
  histContainer_["stepPlot_w_2e2mu"]->SetBinContent(4,nEvAfterZ1Cut*scaleWeight);
  histContainer_["stepPlot_w_2e2mu"]->SetBinContent(5,nEvAfterZ2Cut_2e2mu*scaleWeight);
  histContainer_["stepPlot_w_2e2mu"]->SetBinContent(6,nEvPassedPtCut2_2e2mu*scaleWeight);
  histContainer_["stepPlot_w_2e2mu"]->SetBinContent(7,nEvAfterCleaning_2e2mu*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w_2e2mu"]->SetBinContent(8,nEvAfterZ4lCut_2e2mu*scaleWeight);
  if(!doBlinding || isMC)histContainer_["stepPlot_w_2e2mu"]->SetBinContent(9,nEvAfterM4lCut_2e2mu*scaleWeight);

  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(1,skim);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(2,trig);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(3,ID);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(4,Z1Cut);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(5,Z2Cut);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(6,twoLepPt);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(7,cleaning);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(8,m4l70cut);
  histContainer_["stepPlot_w_2e2mu"]->GetXaxis()->SetBinLabel(9,m4l100cut);

}


void UFHZZ4LAna::bookPassedEventTree(TString treeName, TTree *tree)
{     
  tree->Branch("Run",&Run,"Run/l");
  tree->Branch("Event",&Event,"Event/l");
  tree->Branch("LumiSect",&LumiSect,"LumiSect/l");
  tree->Branch("EventType",&tmpEvent);
  tree->Branch("mass4l",&mass4l,"mass4l/D");
  tree->Branch("ZZMass",&mass4l,"ZZMass/D");
  tree->Branch("mass4lNoFSR",&m4lNoFSR,"mass4lNoFSR/D");
  tree->Branch("mass4mu",&mass4mu,"mass4mu/D");
  tree->Branch("mass4e",&mass4e,"mass4e/D");
  tree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/D");
  tree->Branch("minMass3l",&minM3l,"minMass3l/D");
  tree->Branch("M3l_softestLep",&m3l_soft,"M3l_softestLep/D");
  tree->Branch("Z4l_maxP",&Z4lmaxP,"Z4l_maxP/D");
  tree->Branch("minDeltaR",&minDeltR,"minDeltaR/D");
  tree->Branch("Z4l_thetaPhoton",&thetaPhoton,"Z4l_thetaPhoton/D");
  tree->Branch("Z4l_thetaPhoton_deg",&thetaPhoton_deg,"Z4l_thetaPhoton_deg/D");
  tree->Branch("Z4l_thetaPhotonZ",&thetaPhotonZ,"Z4l_thetaPhotonZ/D");
  tree->Branch("Z4l_thetaPhotonZ_deg",&thetaPhotonZ_deg,"Z4l_thetaPhotonZ_deg/D");
  tree->Branch("Z4l_theta12",&theta12,"Z4l_theta12/D");
  tree->Branch("Z4l_theta13",&theta13,"Z4l_theta13/D");
  tree->Branch("Z4l_theta14",&theta14,"Z4l_theta14/D");
  tree->Branch("Z4l_theta12_deg",&theta12_deg,"Z4l_theta12_deg/D");
  tree->Branch("Z4l_theta13_deg",&theta13_deg,"Z4l_theta13_deg/D");
  tree->Branch("Z4l_theta14_deg",&theta14_deg,"Z4l_theta14_deg/D");
  tree->Branch("Z4l_minMass2l",&minMass2Lep,"Z4l_minMass2l/D");
  tree->Branch("Z4l_maxMass2l",&maxMass2Lep,"Z4l_maxMass2l/D");
  tree->Branch("thetaZ2",&thetaZ2,"thetaZ2/D");
  tree->Branch("etaZ2",&etaZ2,"etaZ2/D");
  tree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
  tree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
  tree->Branch("passedQCDcut",&passedQCDcut,"passedQCDcut/O");
  tree->Branch("finalState",&finalState,"finalState/I");
  tree->Branch("mixedFinalState",&RecoFourMixEvent,"mixedFinalState/I");
  tree->Branch("scaleWeight",&scaleWeight,"scaleWeight/D");
  tree->Branch("eventWeight",&eventWeight,"eventWeight/D");
  tree->Branch("MC_weight",&MC_weight,"MC_weight/D");
  tree->Branch("pT4l",&pT4l,"pT4l/D");
  tree->Branch("massZ1",&massZ1,"massZ1/D");
  tree->Branch("massZ2",&massZ2,"massZ2/D");  
  tree->Branch("Z1Mass",&massZ1,"Z1Mass/D");
  tree->Branch("Z2Mass",&massZ2,"Z2Mass/D");
  tree->Branch("eta4l",&eta4l,"eta4l/D");
  tree->Branch("phi4l",&phi4l,"phi4l/D");
  tree->Branch("idL1",&idL1,"idL1/I");
  tree->Branch("idL2",&idL2,"idL2/I");
  tree->Branch("idL3",&idL3,"idL3/I");
  tree->Branch("idL4",&idL4,"idL4/I");
  tree->Branch("missingHitsL1",&missingHitsL1,"missingHitsL1/I");
  tree->Branch("missingHitsL2",&missingHitsL2,"missingHitsL2/I");
  tree->Branch("missingHitsL3",&missingHitsL3,"missingHitsL3/I");
  tree->Branch("missingHitsL4",&missingHitsL4,"missingHitsL4/I");
  tree->Branch("ecalSeedL1",&ecalSeedL1,"ecalSeedL1/O");
  tree->Branch("ecalSeedL2",&ecalSeedL2,"ecalSeedL2/O");
  tree->Branch("ecalSeedL3",&ecalSeedL3,"ecalSeedL3/O");
  tree->Branch("ecalSeedL4",&ecalSeedL4,"ecalSeedL4/O");
  tree->Branch("EL1",&EL1,"EL1/D");
  tree->Branch("EL2",&EL2,"EL2/D");
  tree->Branch("EL3",&EL3,"EL3/D");
  tree->Branch("EL4",&EL4,"EL4/D");
  tree->Branch("EL1FSR",&EL1FSR,"EL1FSR/D");
  tree->Branch("EL2FSR",&EL2FSR,"EL2FSR/D");
  tree->Branch("EL3FSR",&EL3FSR,"EL3FSR/D");
  tree->Branch("EL4FSR",&EL4FSR,"EL4FSR/D");
  tree->Branch("mvaL1",&mvaL1,"mvaL1/D");
  tree->Branch("mvaL2",&mvaL2,"mvaL2/D");
  tree->Branch("mvaL3",&mvaL3,"mvaL3/D");
  tree->Branch("mvaL4",&mvaL4,"mvaL4/D");
  tree->Branch("SipL1",&SipL1,"SipL1/D");
  tree->Branch("SipL2",&SipL2,"SipL2/D");
  tree->Branch("SipL3",&SipL3,"SipL3/D");
  tree->Branch("SipL4",&SipL4,"SipL4/D");
  tree->Branch("IPL1",&IPL1,"IPL1/D");
  tree->Branch("IPL2",&IPL2,"IPL2/D");
  tree->Branch("IPL3",&IPL3,"IPL3/D");
  tree->Branch("IPL4",&IPL4,"IPL4/D");
  tree->Branch("dIPL1",&dIPL1,"dIPL1/D");
  tree->Branch("dIPL2",&dIPL2,"dIPL2/D");
  tree->Branch("dIPL3",&dIPL3,"dIPL3/D");
  tree->Branch("dIPL4",&dIPL4,"dIPL4/D");
  tree->Branch("EZ1",&EZ1,"EZ1/D");
  tree->Branch("EZ2",&EZ2,"EZ2/D");

  tree->Branch("pTL1",&pTL1,"pTL1/D");
  tree->Branch("pTL2",&pTL2,"pTL2/D");
  tree->Branch("pTL3",&pTL3,"pTL3/D");
  tree->Branch("pTL4",&pTL4,"pTL4/D");
  tree->Branch("pXL1",&pXL1,"pXL1/D");
  tree->Branch("pXL2",&pXL2,"pXL2/D");
  tree->Branch("pXL3",&pXL3,"pXL3/D");
  tree->Branch("pXL4",&pXL4,"pXL4/D");
  tree->Branch("pYL1",&pYL1,"pYL1/D");
  tree->Branch("pYL2",&pYL2,"pYL2/D");
  tree->Branch("pYL3",&pYL3,"pYL3/D");
  tree->Branch("pYL4",&pYL4,"pYL4/D");
  tree->Branch("pZL1",&pZL1,"pZL1/D");
  tree->Branch("pZL2",&pZL2,"pZL2/D");
  tree->Branch("pZL3",&pZL3,"pZL3/D");
  tree->Branch("pZL4",&pZL4,"pZL4/D");

  tree->Branch("pTL1FSR",&pTL1FSR,"pTL1FSR/D");
  tree->Branch("pTL2FSR",&pTL2FSR,"pTL2FSR/D");
  tree->Branch("pTL3FSR",&pTL3FSR,"pTL3FSR/D");
  tree->Branch("pTL4FSR",&pTL4FSR,"pTL4FSR/D");
  tree->Branch("pXL1FSR",&pXL1FSR,"pXL1FSR/D");
  tree->Branch("pXL2FSR",&pXL2FSR,"pXL2FSR/D");
  tree->Branch("pXL3FSR",&pXL3FSR,"pXL3FSR/D");
  tree->Branch("pXL4FSR",&pXL4FSR,"pXL4FSR/D");
  tree->Branch("pYL1FSR",&pYL1FSR,"pYL1FSR/D");
  tree->Branch("pYL2FSR",&pYL2FSR,"pYL2FSR/D");
  tree->Branch("pYL3FSR",&pYL3FSR,"pYL3FSR/D");
  tree->Branch("pYL4FSR",&pYL4FSR,"pYL4FSR/D");
  tree->Branch("pZL1FSR",&pZL1FSR,"pZL1FSR/D");
  tree->Branch("pZL2FSR",&pZL2FSR,"pZL2FSR/D");
  tree->Branch("pZL3FSR",&pZL3FSR,"pZL3FSR/D");
  tree->Branch("pZL4FSR",&pZL4FSR,"pZL4FSR/D");

  tree->Branch("pTZ1",&pTZ1,"pTZ1/D");
  tree->Branch("pTZ2",&pTZ2,"pTZ2/D");
  tree->Branch("pXZ1",&pXZ1,"pXZ1/D");
  tree->Branch("pXZ2",&pXZ2,"pXZ2/D");
  tree->Branch("pYZ1",&pYZ1,"pYZ1/D");
  tree->Branch("pYZ2",&pYZ2,"pYZ2/D");
  tree->Branch("pZZ1",&pZZ1,"pZZ1/D");
  tree->Branch("pZZ2",&pZZ2,"pZZ2/D");
  tree->Branch("chargeL1",&chargeL1,"chargeL1/I");
  tree->Branch("chargeL2",&chargeL2,"chargeL2/I");
  tree->Branch("chargeL3",&chargeL3,"chargeL3/I");
  tree->Branch("chargeL4",&chargeL4,"chargeL4/I");
  tree->Branch("etaL1",&etaL1,"etaL1/D");
  tree->Branch("etaL2",&etaL2,"etaL2/D");
  tree->Branch("etaL3",&etaL3,"etaL3/D");
  tree->Branch("etaL4",&etaL4,"etaL4/D");
  tree->Branch("phiL1",&phiL1,"phiL1/D");
  tree->Branch("phiL2",&phiL2,"phiL2/D");
  tree->Branch("phiL3",&phiL3,"phiL3/D");
  tree->Branch("phiL4",&phiL4,"phiL4/D");
  tree->Branch("isoTrackL1",&isoTrackL1,"isoTrackL1/D");
  tree->Branch("isoTrackL2",&isoTrackL2,"isoTrackL2/D");
  tree->Branch("isoTrackL3",&isoTrackL3,"isoTrackL3/D");
  tree->Branch("isoTrackL4",&isoTrackL4,"isoTrackL4/D");
  tree->Branch("isoEcalL1",&isoEcalL1,"isoEcalL1/D");
  tree->Branch("isoEcalL2",&isoEcalL2,"isoEcalL2/D");
  tree->Branch("isoEcalL3",&isoEcalL3,"isoEcalL3/D");
  tree->Branch("isoEcalL4",&isoEcalL4,"isoEcalL4/D");
  tree->Branch("isoHcalL1",&isoHcalL1,"isoHcalL1/D");
  tree->Branch("isoHcalL2",&isoHcalL2,"isoHcalL2/D");
  tree->Branch("isoHcalL3",&isoHcalL3,"isoHcalL3/D");
  tree->Branch("isoHcalL4",&isoHcalL4,"isoHcalL4/D");
  tree->Branch("isoNHL1",&isoNHL1,"isoNHL1/D");
  tree->Branch("isoNHL2",&isoNHL2,"isoNHL2/D");
  tree->Branch("isoNHL3",&isoNHL3,"isoNHL3/D");
  tree->Branch("isoNHL4",&isoNHL4,"isoNHL4/D");
  tree->Branch("isoCHL1",&isoCHL1,"isoCHL1/D");
  tree->Branch("isoCHL2",&isoCHL2,"isoCHL2/D");
  tree->Branch("isoCHL3",&isoCHL3,"isoCHL3/D");
  tree->Branch("isoCHL4",&isoCHL4,"isoCHL4/D");
  tree->Branch("isoPhotL1",&isoPhotL1,"isoPhotL1/D");
  tree->Branch("isoPhotL2",&isoPhotL2,"isoPhotL2/D");
  tree->Branch("isoPhotL3",&isoPhotL3,"isoPhotL3/D");
  tree->Branch("isoPhotL4",&isoPhotL4,"isoPhotL4/D");
  tree->Branch("RelIsoL1",&RelIsoL1,"RelIsoL1/D");
  tree->Branch("RelIsoL2",&RelIsoL2,"RelIsoL2/D");
  tree->Branch("RelIsoL3",&RelIsoL3,"RelIsoL3/D");
  tree->Branch("RelIsoL4",&RelIsoL4,"RelIsoL4/D");
  tree->Branch("RelIsoUCL1",&RelIsoUCL1,"RelIsoUCL1/D");
  tree->Branch("RelIsoUCL2",&RelIsoUCL2,"RelIsoUCL2/D");
  tree->Branch("RelIsoUCL3",&RelIsoUCL3,"RelIsoUCL3/D");
  tree->Branch("RelIsoUCL4",&RelIsoUCL4,"RelIsoUCL4/D");
  tree->Branch("muRho",&muRho,"muRho/D");
  tree->Branch("elRho",&elRho,"elRho/D");
  tree->Branch("worstIso",&worstIso,"worstIso/D");
  tree->Branch("worstSIP",&worstSIP,"worstSIP/D");
  tree->Branch("helcosthetaZ1",&cosTheta1,"helcosthetaZ1/F");
  tree->Branch("helcosthetaZ2",&cosTheta2,"helcosthetaZ2/F");
  tree->Branch("helphi",&Phi,"helphi/F");
  tree->Branch("costhetastar",&cosThetaStar,"costhetastar/F");
  tree->Branch("phistarZ1",&phiStar1,"phistarZ1/F");
  //tree->Branch("phistarZ2",&phiStar2,"phistarZ2/D");
  //tree->Branch("Phi1",&Phi1,"Phi1/D");
  //tree->Branch("phiStar12",&phiStar12,"phiStar12/D");
  //tree->Branch("Phi2",&Phi2,"Phi2/D");
  tree->Branch("nJets",&nJets,"nJets/I");
  tree->Branch("nVtx",&nVtx,"nVtx/I");
  tree->Branch("nPhotons",&nPhotons,"nPhotons/I");
  tree->Branch("metVal",&metVal,"metVal/D");
  tree->Branch("extraLep_id",extraLep_id,"extraLep_id[10]/I");
  tree->Branch("extraLep_pT",extraLep_pT,"extraLep_pT[10]/D");
  tree->Branch("extraLep_pX",extraLep_pX,"extraLep_pX[10]/D");
  tree->Branch("extraLep_pY",extraLep_pY,"extraLep_pY[10]/D");
  tree->Branch("extraLep_pZ",extraLep_pZ,"extraLep_pZ[10]/D");
  tree->Branch("extraLep_e",extraLep_e,"extraLep_e[10]/D");
  tree->Branch("extraLep_iso",extraLep_iso,"extraLep_iso[10]/D");
  tree->Branch("extraLep_chIso",extraLep_chIso,"extraLep_chIso[10]/D");
  tree->Branch("extraLep_nhIso",extraLep_nhIso,"extraLep_nhIso[10]/D");
  tree->Branch("extraLep_phIso",extraLep_phIso,"extraLep_phIso[10]/D");
  tree->Branch("extraLep_eta",extraLep_eta,"extraLep_eta[10]/D");
  tree->Branch("extraLep_phi",extraLep_phi,"extraLep_phi[10]/D");
  tree->Branch("extraLep_sip",extraLep_sip,"extraLep_sip[10]/D");

  tree->Branch("massErrorUF",&m4lErrUF,"massErrorUF/D");
  tree->Branch("massErrorUFCorr",&m4lErrUFCorr,"massErrorUFCorr/D");
  tree->Branch("massErrorUCSD",&m4lErrUCSD,"massErrorUCSD/D");
  tree->Branch("massErrorUCSDCorr",&m4lErrUCSDCorr,"massErrorUCSDCorr/D");

  tree->Branch("FSR_Z1",&FSR_Z1,"FSR_Z1/O");
  tree->Branch("FSR_Z2",&FSR_Z2,"FSR_Z2/O");
  tree->Branch("FSRPhot1_Pt",&FSRPhot1_Pt,"FSRPhot1_Pt/D");
  tree->Branch("FSRPhot2_Pt",&FSRPhot2_Pt,"FSRPhot2_Pt/D");
  tree->Branch("FSRPhot1_eta",&FSRPhot1_eta,"FSRPhot1_eta/D");
  tree->Branch("FSRPhot2_eta",&FSRPhot2_eta,"FSRPhot2_eta/D");
  tree->Branch("FSRPhot1_phi",&FSRPhot1_phi,"FSRPhot1_phi/D");
  tree->Branch("FSRPhot2_phi",&FSRPhot2_phi,"FSRPhot2_phi/D");
  tree->Branch("FSRPhot1_Px",&FSRPhot1_Px,"FSRPhot1_Px/D");
  tree->Branch("FSRPhot2_Px",&FSRPhot2_Px,"FSRPhot2_Px/D");
  tree->Branch("FSRPhot1_Py",&FSRPhot1_Py,"FSRPhot1_Py/D");
  tree->Branch("FSRPhot2_Py",&FSRPhot2_Py,"FSRPhot2_Py/D");
  tree->Branch("FSRPhot1_Pz",&FSRPhot1_Pz,"FSRPhot1_Pz/D");
  tree->Branch("FSRPhot2_Pz",&FSRPhot2_Pz,"FSRPhot2_Pz/D");
  tree->Branch("FSRPhot1_E",&FSRPhot1_E,"FSRPhot1_E/D");
  tree->Branch("FSRPhot2_E",&FSRPhot2_E,"FSRPhot2_E/D");

  //PDF M4L
  tree->Branch("pdfSigM4l",&pdfSigM4l,"pdfSigM4l/D");
  tree->Branch("pdfSigM4l_noFSR",&pdfSigM4l_noFSR,"pdfSigM4l_noFSR/D");
  tree->Branch("pdfBkgM4l",&pdfBkgM4l,"pdfBkgM4l/D");
  tree->Branch("pdfBkgM4l_noFSR",&pdfBkgM4l_noFSR,"pdfBkgM4l_noFSR/D");
  tree->Branch("pdfSigM4l_ScaleUp",&pdfSigM4l_ScaleUp,"pdfSigM4l_ScaleUp/D");
  tree->Branch("pdfSigM4l_ScaleUp_noFSR",&pdfSigM4l_ScaleUp_noFSR,"pdfSigM4l_ScaleUp_noFSR/D");
  tree->Branch("pdfBkgM4l_ScaleUp",&pdfBkgM4l_ScaleUp,"pdfBkgM4l_ScaleUp/D");
  tree->Branch("pdfBkgM4l_ScaleUp_noFSR",&pdfBkgM4l_ScaleUp_noFSR,"pdfBkgM4l_ScaleUp_noFSR/D");
  tree->Branch("pdfSigM4l_ScaleDown",&pdfSigM4l_ScaleDown,"pdfSigM4l_ScaleDown/D");
  tree->Branch("pdfSigM4l_ScaleDown_noFSR",&pdfSigM4l_ScaleDown_noFSR,"pdfSigM4l_ScaleDown_noFSR/D");
  tree->Branch("pdfBkgM4l_ScaleDown",&pdfBkgM4l_ScaleDown,"pdfBkgM4l_ScaleDown/D");
  tree->Branch("pdfBkgM4l_ScaleDown_noFSR",&pdfBkgM4l_ScaleDown_noFSR,"pdfBkgM4l_ScaleDown_noFSR/D");
  tree->Branch("pdfSigM4l_ResUp",&pdfSigM4l_ResUp,"pdfSigM4l_ResUp/D");
  tree->Branch("pdfSigM4l_ResUp_noFSR",&pdfSigM4l_ResUp_noFSR,"pdfSigM4l_ResUp_noFSR/D");
  tree->Branch("pdfBkgM4l_ResUp",&pdfBkgM4l_ResUp,"pdfBkgM4l_ResUp/D");
  tree->Branch("pdfBkgM4l_ResUp_noFSR",&pdfBkgM4l_ResUp_noFSR,"pdfBkgM4l_ResUp_noFSR/D");
  tree->Branch("pdfSigM4l_ResDown",&pdfSigM4l_ResDown,"pdfSigM4l_ResDown/D");
  tree->Branch("pdfSigM4l_ResDown_noFSR",&pdfSigM4l_ResDown_noFSR,"pdfSigM4l_ResDown_noFSR/D");
  tree->Branch("pdfBkgM4l_ResDown",&pdfBkgM4l_ResDown,"pdfBkgM4l_ResDown/D");
  tree->Branch("pdfBkgM4l_ResDown_noFSR",&pdfBkgM4l_ResDown_noFSR,"pdfBkgM4l_ResDown_noFSR/D");

  //MELA
  tree->Branch("melaLD",&melaLD,"melaLD/D");
  tree->Branch("mela_Sig",&mela_Sig,"mela_Sig/D");
  tree->Branch("mela_Bkg",&mela_Bkg,"mela_Bkg/D");
  /*
  tree->Branch("melaLD_Pt",&melaLD_Pt,"melaLD_Pt/F");
  tree->Branch("mela_Sig_Pt",&mela_Sig_Pt,"mela_Sig_Pt/F");
  tree->Branch("mela_Bkg_Pt",&mela_Bkg_Pt,"mela_Bkg_Pt/F");
  tree->Branch("melaLD_Y",&melaLD_Y,"melaLD_Y/F");
  tree->Branch("mela_Sig_Y",&mela_Sig_Y,"mela_Sig_Y/F");
  tree->Branch("mela_Bkg_Y",&mela_Bkg_Y,"mela_Bkg_Y/F");
  tree->Branch("melaLD_PtY",&melaLD_PtY,"melaLD_PtY/F");
  tree->Branch("mela_Sig_PtY",&mela_Sig_PtY,"mela_Sig_PtY/F");
  tree->Branch("mela_Bkg_PtY",&mela_Bkg_PtY,"mela_Bkg_PtY/F");
  */
  tree->Branch("pseudoMelaLD",&pseudoMelaLD,"pseudoMelaLD/F");
  tree->Branch("pseudoMela_SM",&pseudoMela_SM,"pseudoMela_SM/F");
  tree->Branch("pseudoMela_PS",&pseudoMela_PS,"pseudoMela_PS/F");

  tree->Branch("spin1EvMelaLD",&spin1EvMelaLD,"spin1EvMelaLD/F");
  tree->Branch("spin1EvMela_SM",&spin1EvMela_SM,"spin1EvMela_SM/F");
  tree->Branch("spin1EvMela_S1E",&spin1EvMela_S1E,"spin1EvMela_S1E/F");

  tree->Branch("spin1OddMelaLD",&spin1OddMelaLD,"spin1OddMelaLD/F");
  tree->Branch("spin1OddMela_SM",&spin1OddMela_SM,"spin1OddMela_SM/F");
  tree->Branch("spin1OddMela_S1O",&spin1OddMela_S1O,"spin1OddMela_S1O/F");

  tree->Branch("spin2MelaLD",&spin2MelaLD,"spin2MelaLD/F");
  tree->Branch("spin2Mela_SM",&spin2Mela_SM,"spin2Mela_SM/F");
  tree->Branch("spin2Mela_S2",&spin2Mela_S2,"spin2Mela_S2/F");

  /*
  tree->Branch("p0plus_melaNorm",&p0plus_melaNorm,"p0plus_melaNorm/F");
  tree->Branch("p0plus_mela",&p0plus_mela,"p0plus_mela/F");
  tree->Branch("p0minus_mela",&p0minus_mela,"p0minus_mela/F");
  tree->Branch("p0plus_VAJHU",&p0plus_VAJHU,"p0plus_VAJHU/F");
  tree->Branch("p0minus_VAJHU",&p0minus_VAJHU,"p0minus_VAJHU/F");
  tree->Branch("p0plus_VAMCFM",&p0plus_VAMCFM,"p0plus_VAMCFM/F");
  tree->Branch("p1_mela",&p1_mela,"p1_mela/F");
  tree->Branch("p1_VAJHU",&p1_VAJHU,"p1_VAJHU/F");
  tree->Branch("p2_mela",&p2_mela,"p2_mela/F");
  tree->Branch("p2_VAJHU",&p2_VAJHU,"p2_VAJHU/F");
  tree->Branch("mela_bkg_analytic",&mela_bkg_analytic,"mela_bkg_analytic/F");
  tree->Branch("mela_bkg_VAMCFM",&mela_bkg_VAMCFM,"mela_bkg_VAMCFM/F");
  tree->Branch("mela_bkg_ggzz_VAMCFM",&mela_bkg_ggzz_VAMCFM,"mela_bkg_ggzz_VAMCFM/F");
  tree->Branch("mela_bkg_VAMCFMNorm",&mela_bkg_VAMCFMNorm,"mela_bkg_VAMCFMNorm/F");
  tree->Branch("mela_p0_pt",&mela_p0_pt,"mela_p0_pt/F");
  tree->Branch("mela_p0_y",&mela_p0_y,"mela_p0_y/F");
  tree->Branch("p0plus_m4l",&p0plus_m4l,"p0plus_m4l/F");
  tree->Branch("bkg_m4l",&bkg_m4l,"bkg_m4l/F");
  */

  tree->Branch("MEKD_wPDF",&MEKD_wPDF,"MEKD_wPDF/D");
  tree->Branch("MEKD_ME_H_wPDF",&MEKD_ME_H_wPDF,"MEKD_ME_H_wPDF/D");
  tree->Branch("MEKD_ME_ZZ_wPDF",&MEKD_ME_ZZ_wPDF,"MEKD_ME_ZZ_wPDF/D");

  tree->Branch("MEKD_noPDF",&MEKD_noPDF,"MEKD_noPDF/D");
  tree->Branch("MEKD_ME_H_noPDF",&MEKD_ME_H_noPDF,"MEKD_ME_H_noPDF/D");
  tree->Branch("MEKD_ME_ZZ_noPDF",&MEKD_ME_ZZ_noPDF,"MEKD_ME_ZZ_noPDF/D");
 
  tree->Branch("MEKD_wPDF_noFSR",&MEKD_wPDF_noFSR,"MEKD_wPDF_noFSR/D");
  tree->Branch("MEKD_ME_H_wPDF_noFSR",&MEKD_ME_H_wPDF_noFSR,"MEKD_ME_H_wPDF_noFSR/D");
  tree->Branch("MEKD_ME_ZZ_wPDF_noFSR",&MEKD_ME_ZZ_wPDF_noFSR,"MEKD_ME_ZZ_wPDF_noFSR/D");
 
  tree->Branch("MEKD_noPDF_noFSR",&MEKD_noPDF_noFSR,"MEKD_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_H_noPDF_noFSR",&MEKD_ME_H_noPDF_noFSR,"MEKD_ME_H_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_ZZ_noPDF_noFSR",&MEKD_ME_ZZ_noPDF_noFSR,"MEKD_ME_ZZ_noPDF_noFSR/D");
 
  tree->Branch("interferenceWeight",&interferenceWeight,"interferenceWeight/D");
  tree->Branch("JHUKD_H_qqZZ_noPDF",&JHUKD_H_qqZZ_noPDF,"JHUKD_H_qqZZ_noPDF/D");
  tree->Branch("JHU_ME_H",&JHU_ME_H,"JHU_ME_H/D");
  tree->Branch("MCFM_ME_qqZZ",&MCFM_ME_qqZZ,"MCFM_ME_qqZZ/D");
  tree->Branch("JHUKD_H_h0M_noPDF",&JHUKD_H_h0M_noPDF,"JHUKD_H_h0M_noPDF/D");
  tree->Branch("JHU_ME_h0M",&JHU_ME_h0M,"JHU_ME_h0M/D");
  tree->Branch("JHUKD_H_h0P_noPDF",&JHUKD_H_h0P_noPDF,"JHUKD_H_h0P_noPDF/D");
  tree->Branch("JHU_ME_h0P",&JHU_ME_h0P,"JHU_ME_h0P/D");
  tree->Branch("JHUKD_H_h1M_noPDF",&JHUKD_H_h1M_noPDF,"JHUKD_H_h1M_noPDF/D");
  tree->Branch("JHU_ME_h1M",&JHU_ME_h1M,"JHU_ME_h1M/D");
  tree->Branch("JHUKD_H_h1P_noPDF",&JHUKD_H_h1P_noPDF,"JHUKD_H_h1P_noPDF/D");
  tree->Branch("JHU_ME_h1P",&JHU_ME_h1P,"JHU_ME_h1P/D");
  tree->Branch("JHUKD_H_ggh2P_noPDF",&JHUKD_H_ggh2P_noPDF,"JHUKD_H_ggh2P_noPDF/D");
  tree->Branch("JHU_ME_ggh2P",&JHU_ME_ggh2P,"JHU_ME_ggh2P/D");
  tree->Branch("JHUKD_H_qqh2P_noPDF",&JHUKD_H_qqh2P_noPDF,"JHUKD_H_qqh2P_noPDF/D");
  tree->Branch("JHU_ME_qqh2P",&JHU_ME_qqh2P,"JHU_ME_qqh2P/D");

  tree->Branch("JHUKD_H_h2hP_noPDF",&JHUKD_H_h2hP_noPDF,"JHUKD_H_h2hP_noPDF/D");
  tree->Branch("JHU_ME_h2hP",&JHU_ME_h2hP,"JHU_ME_h2hP/D");
  tree->Branch("JHUKD_H_h2hM_noPDF",&JHUKD_H_h2hM_noPDF,"JHUKD_H_h2hM_noPDF/D");
  tree->Branch("JHU_ME_h2hM",&JHU_ME_h2hM,"JHU_ME_h2hM/D");
  tree->Branch("JHUKD_H_h2bP_noPDF",&JHUKD_H_h2bP_noPDF,"JHUKD_H_h2bP_noPDF/D");
  tree->Branch("JHU_ME_h2bP",&JHU_ME_h2bP,"JHU_ME_h2bP/D");
  tree->Branch("JHUKD_H_h2P_prodInd_noPDF",&JHUKD_H_h2P_prodInd_noPDF,"JHUKD_H_h2P_prodInd_noPDF/D");
  tree->Branch("JHU_ME_h2P_prodInd",&JHU_ME_h2P_prodInd,"JHU_ME_h2P_prodInd/D");
  tree->Branch("JHUKD_H_h1P_prodInd_noPDF",&JHUKD_H_h1P_prodInd_noPDF,"JHUKD_H_h1P_prodInd_noPDF/D");
  tree->Branch("JHU_ME_h1P_prodInd",&JHU_ME_h1P_prodInd,"JHU_ME_h1P_prodInd/D");
  tree->Branch("JHUKD_H_h1M_prodInd_noPDF",&JHUKD_H_h1M_prodInd_noPDF,"JHUKD_H_h1M_prodInd_noPDF/D");
  tree->Branch("JHU_ME_h1M_prodInd",&JHU_ME_h1M_prodInd,"JHU_ME_h1M_prodInd/D");

  tree->Branch("MEKD_noPDF_noFSR",&MEKD_noPDF_noFSR,"MEKD_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_H_noPDF_noFSR",&MEKD_ME_H_noPDF_noFSR,"MEKD_ME_H_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_ZZ_noPDF_noFSR",&MEKD_ME_ZZ_noPDF_noFSR,"MEKD_ME_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_h0M_ZZ_noPDF_noFSR",&MEKD_h0M_ZZ_noPDF_noFSR,"MEKD_h0M_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h0M_noPDF_noFSR",&MEKD_ME_h0M_noPDF_noFSR,"MEKD_ME_h0M_noPDF_noFSR/D");
  tree->Branch("MEKD_h0P_ZZ_noPDF_noFSR",&MEKD_h0P_ZZ_noPDF_noFSR,"MEKD_h0P_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h0P_noPDF_noFSR",&MEKD_ME_h0P_noPDF_noFSR,"MEKD_ME_h0P_noPDF_noFSR/D");
  tree->Branch("MEKD_h1M_ZZ_noPDF_noFSR",&MEKD_h1M_ZZ_noPDF_noFSR,"MEKD_h1M_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h1M_noPDF_noFSR",&MEKD_ME_h1M_noPDF_noFSR,"MEKD_ME_h1M_noPDF_noFSR/D");
  tree->Branch("MEKD_h1P_ZZ_noPDF_noFSR",&MEKD_h1P_ZZ_noPDF_noFSR,"MEKD_h1P_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h1P_noPDF_noFSR",&MEKD_ME_h1P_noPDF_noFSR,"MEKD_ME_h1P_noPDF_noFSR/D");
  tree->Branch("MEKD_qqh2P_ZZ_noPDF_noFSR",&MEKD_qqh2P_ZZ_noPDF_noFSR,"MEKD_qqh2P_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_qqh2P_noPDF_noFSR",&MEKD_ME_qqh2P_noPDF_noFSR,"MEKD_ME_qqh2P_noPDF_noFSR/D");
  tree->Branch("MEKD_ggh2P_ZZ_noPDF_noFSR",&MEKD_ggh2P_ZZ_noPDF_noFSR,"MEKD_ggh2P_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_ggh2P_noPDF_noFSR",&MEKD_ME_ggh2P_noPDF_noFSR,"MEKD_ME_ggh2P_noPDF_noFSR/D");

  tree->Branch("MEKD_h2hP_ZZ_noPDF_noFSR",&MEKD_h2hP_ZZ_noPDF_noFSR,"MEKD_h2hP_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h2hP_noPDF_noFSR",&MEKD_ME_h2hP_noPDF_noFSR,"MEKD_ME_h2hP_noPDF_noFSR/D");
  tree->Branch("MEKD_h2hM_ZZ_noPDF_noFSR",&MEKD_h2hM_ZZ_noPDF_noFSR,"MEKD_h2hM_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h2hM_noPDF_noFSR",&MEKD_ME_h2hM_noPDF_noFSR,"MEKD_ME_h2hM_noPDF_noFSR/D");
  tree->Branch("MEKD_h2bP_ZZ_noPDF_noFSR",&MEKD_h2bP_ZZ_noPDF_noFSR,"MEKD_h2bP_ZZ_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h2bP_noPDF_noFSR",&MEKD_ME_h2bP_noPDF_noFSR,"MEKD_ME_h2bP_noPDF_noFSR/D");
  tree->Branch("MEKD_h2P_ZZ_prodInd_noPDF_noFSR",&MEKD_h2P_ZZ_prodInd_noPDF_noFSR,"MEKD_h2P_ZZ_prodInd_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h2P_prodInd_noPDF_noFSR",&MEKD_ME_h2P_prodInd_noPDF_noFSR,"MEKD_ME_h2P_prodInd_noPDF_noFSR/D");
  tree->Branch("MEKD_h1P_ZZ_prodInd_noPDF_noFSR",&MEKD_h1P_ZZ_prodInd_noPDF_noFSR,"MEKD_h1P_ZZ_prodInd_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h1P_prodInd_noPDF_noFSR",&MEKD_ME_h1P_prodInd_noPDF_noFSR,"MEKD_ME_h1P_prodInd_noPDF_noFSR/D");
  tree->Branch("MEKD_h1M_ZZ_prodInd_noPDF_noFSR",&MEKD_h1M_ZZ_prodInd_noPDF_noFSR,"MEKD_h1M_ZZ_prodInd_noPDF_noFSR/D");
  tree->Branch("MEKD_ME_h1M_prodInd_noPDF_noFSR",&MEKD_ME_h1M_prodInd_noPDF_noFSR,"MEKD_ME_h1M_prodInd_noPDF_noFSR/D");

  tree->Branch("JHUKD_H_qqZZ_noPDF_noFSR",&JHUKD_H_qqZZ_noPDF_noFSR,"JHUKD_H_qqZZ_noPDF_noFSR/D");
  tree->Branch("JHU_ME_H_noPDF_noFSR",&JHU_ME_H_noPDF_noFSR,"JHU_ME_H_noPDF_noFSR/D");
  tree->Branch("JHU_ME_ZZ_noPDF_noFSR",&JHU_ME_ZZ_noPDF_noFSR,"JHU_ME_ZZ_noPDF_noFSR/D");
  tree->Branch("JHUKD_h0M_ZZ_noPDF_noFSR",&JHUKD_h0M_ZZ_noPDF_noFSR,"JHUKD_h0M_ZZ_noPDF_noFSR/D");
  tree->Branch("JHU_ME_h0M_noPDF_noFSR",&JHU_ME_h0M_noPDF_noFSR,"JHU_ME_h0M_noPDF_noFSR/D");
  tree->Branch("JHUKD_h0P_ZZ_noPDF_noFSR",&JHUKD_h0P_ZZ_noPDF_noFSR,"JHUKD_h0P_ZZ_noPDF_noFSR/D");
  tree->Branch("JHU_ME_h0P_noPDF_noFSR",&JHU_ME_h0P_noPDF_noFSR,"JHU_ME_h0P_noPDF_noFSR/D");
  tree->Branch("JHUKD_h1M_ZZ_noPDF_noFSR",&JHUKD_h1M_ZZ_noPDF_noFSR,"JHUKD_h1M_ZZ_noPDF_noFSR/D");
  tree->Branch("JHU_ME_h1M_noPDF_noFSR",&JHU_ME_h1M_noPDF_noFSR,"JHU_ME_h1M_noPDF_noFSR/D");
  tree->Branch("JHUKD_h1P_ZZ_noPDF_noFSR",&JHUKD_h1P_ZZ_noPDF_noFSR,"JHUKD_h1P_ZZ_noPDF_noFSR/D");
  tree->Branch("JHU_ME_h1P_noPDF_noFSR",&JHU_ME_h1P_noPDF_noFSR,"JHU_ME_h1P_noPDF_noFSR/D");
  tree->Branch("JHUKD_qqh2P_ZZ_noPDF_noFSR",&JHUKD_qqh2P_ZZ_noPDF_noFSR,"JHUKD_qqh2P_ZZ_noPDF_noFSR/D");
  tree->Branch("JHU_ME_qqh2P_noPDF_noFSR",&JHU_ME_qqh2P_noPDF_noFSR,"JHU_ME_qqh2P_noPDF_noFSR/D");
  tree->Branch("JHUKD_ggh2P_ZZ_noPDF_noFSR",&JHUKD_ggh2P_ZZ_noPDF_noFSR,"JHUKD_ggh2P_ZZ_noPDF_noFSR/D");
  tree->Branch("JHU_ME_ggh2P_noPDF_noFSR",&JHU_ME_ggh2P_noPDF_noFSR,"JHU_ME_ggh2P_noPDF_noFSR/D");
  
  tree->Branch("GENidL1_S1",&GENidL1_S1,"GENidL1_S1/I");
  tree->Branch("GENidL2_S1",&GENidL2_S1,"GENidL2_S1/I");
  tree->Branch("GENidL3_S1",&GENidL3_S1,"GENidL3_S1/I");
  tree->Branch("GENidL4_S1",&GENidL4_S1,"GENidL4_S1/I");
  tree->Branch("GENEL1_S1",&GENEL1_S1,"GENEL1_S1/D");
  tree->Branch("GENEL2_S1",&GENEL2_S1,"GENEL2_S1/D");
  tree->Branch("GENEL3_S1",&GENEL3_S1,"GENEL3_S1/D");
  tree->Branch("GENEL4_S1",&GENEL4_S1,"GENEL4_S1/D");

  tree->Branch("GENidL1_S3",&GENidL1_S3,"GENidL1_S3/I");
  tree->Branch("GENidL2_S3",&GENidL2_S3,"GENidL2_S3/I");
  tree->Branch("GENidL3_S3",&GENidL3_S3,"GENidL3_S3/I");
  tree->Branch("GENidL4_S3",&GENidL4_S3,"GENidL4_S3/I");
  tree->Branch("GENEL1_S3",&GENEL1_S3,"GENEL1_S3/D");
  tree->Branch("GENEL2_S3",&GENEL2_S3,"GENEL2_S3/D");
  tree->Branch("GENEL3_S3",&GENEL3_S3,"GENEL3_S3/D");
  tree->Branch("GENEL4_S3",&GENEL4_S3,"GENEL4_S3/D");

  tree->Branch("GENEZ1",&GENEZ1,"GENEZ1/D");
  tree->Branch("GENEZ2",&GENEZ2,"GENEZ2/D");

  tree->Branch("GENMZ1",&GENMZ1,"GENMZ1/D");
  tree->Branch("GENMZ2",&GENMZ2,"GENMZ2/D");

  tree->Branch("GENpTL1_S1",&GENpTL1_S1,"GENpTL1_S1/D");
  tree->Branch("GENpTL2_S1",&GENpTL2_S1,"GENpTL2_S1/D");
  tree->Branch("GENpTL3_S1",&GENpTL3_S1,"GENpTL3_S1/D");
  tree->Branch("GENpTL4_S1",&GENpTL4_S1,"GENpTL4_S1/D");
  tree->Branch("GENpXL1_S1",&GENpXL1_S1,"GENpXL1_S1/D");
  tree->Branch("GENpXL2_S1",&GENpXL2_S1,"GENpXL2_S1/D");
  tree->Branch("GENpXL3_S1",&GENpXL3_S1,"GENpXL3_S1/D");
  tree->Branch("GENpXL4_S1",&GENpXL4_S1,"GENpXL4_S1/D");
  tree->Branch("GENpYL1_S1",&GENpYL1_S1,"GENpYL1_S1/D");
  tree->Branch("GENpYL2_S1",&GENpYL2_S1,"GENpYL2_S1/D");
  tree->Branch("GENpYL3_S1",&GENpYL3_S1,"GENpYL3_S1/D");
  tree->Branch("GENpYL4_S1",&GENpYL4_S1,"GENpYL4_S1/D");
  tree->Branch("GENpZL1_S1",&GENpZL1_S1,"GENpZL1_S1/D");
  tree->Branch("GENpZL2_S1",&GENpZL2_S1,"GENpZL2_S1/D");
  tree->Branch("GENpZL3_S1",&GENpZL3_S1,"GENpZL3_S1/D");
  tree->Branch("GENpZL4_S1",&GENpZL4_S1,"GENpZL4_S1/D");

  tree->Branch("GENpTL1_S3",&GENpTL1_S3,"GENpTL1_S3/D");
  tree->Branch("GENpTL2_S3",&GENpTL2_S3,"GENpTL2_S3/D");
  tree->Branch("GENpTL3_S3",&GENpTL3_S3,"GENpTL3_S3/D");
  tree->Branch("GENpTL4_S3",&GENpTL4_S3,"GENpTL4_S3/D");
  tree->Branch("GENpXL1_S3",&GENpXL1_S3,"GENpXL1_S3/D");
  tree->Branch("GENpXL2_S3",&GENpXL2_S3,"GENpXL2_S3/D");
  tree->Branch("GENpXL3_S3",&GENpXL3_S3,"GENpXL3_S3/D");
  tree->Branch("GENpXL4_S3",&GENpXL4_S3,"GENpXL4_S3/D");
  tree->Branch("GENpYL1_S3",&GENpYL1_S3,"GENpYL1_S3/D");
  tree->Branch("GENpYL2_S3",&GENpYL2_S3,"GENpYL2_S3/D");
  tree->Branch("GENpYL3_S3",&GENpYL3_S3,"GENpYL3_S3/D");
  tree->Branch("GENpYL4_S3",&GENpYL4_S3,"GENpYL4_S3/D");
  tree->Branch("GENpZL1_S3",&GENpZL1_S3,"GENpZL1_S3/D");
  tree->Branch("GENpZL2_S3",&GENpZL2_S3,"GENpZL2_S3/D");
  tree->Branch("GENpZL3_S3",&GENpZL3_S3,"GENpZL3_S3/D");
  tree->Branch("GENpZL4_S3",&GENpZL4_S3,"GENpZL4_S3/D");

  tree->Branch("GENpTZ1",&GENpTZ1,"GENpTZ1/D");
  tree->Branch("GENpTZ2",&GENpTZ2,"GENpTZ2/D");
  tree->Branch("GENpXZ1",&GENpXZ1,"GENpXZ1/D");
  tree->Branch("GENpXZ2",&GENpXZ2,"GENpXZ2/D");
  tree->Branch("GENpYZ1",&GENpYZ1,"GENpYZ1/D");
  tree->Branch("GENpYZ2",&GENpYZ2,"GENpYZ2/D");
  tree->Branch("GENpZZ1",&GENpZZ1,"GENpZZ1/D");
  tree->Branch("GENpZZ2",&GENpZZ2,"GENpZZ2/D");

  tree->Branch("GENchargeL1_S1",&GENchargeL1_S1,"GENchargeL1_S1/I");
  tree->Branch("GENchargeL2_S1",&GENchargeL2_S1,"GENchargeL2_S1/I");
  tree->Branch("GENchargeL3_S1",&GENchargeL3_S1,"GENchargeL3_S1/I");
  tree->Branch("GENchargeL4_S1",&GENchargeL4_S1,"GENchargeL4_S1/I");
  tree->Branch("GENetaL1_S1",&GENetaL1_S1,"GENetaL1_S1/D");
  tree->Branch("GENetaL2_S1",&GENetaL2_S1,"GENetaL2_S1/D");
  tree->Branch("GENetaL3_S1",&GENetaL3_S1,"GENetaL3_S1/D");
  tree->Branch("GENetaL4_S1",&GENetaL4_S1,"GENetaL4_S1/D");
  tree->Branch("GENphiL1_S1",&GENphiL1_S1,"GENphiL1_S1/D");
  tree->Branch("GENphiL2_S1",&GENphiL2_S1,"GENphiL2_S1/D");
  tree->Branch("GENphiL3_S1",&GENphiL3_S1,"GENphiL3_S1/D");
  tree->Branch("GENphiL4_S1",&GENphiL4_S1,"GENphiL4_S1/D");

  tree->Branch("GENchargeL1_S3",&GENchargeL1_S3,"GENchargeL1_S3/I");
  tree->Branch("GENchargeL2_S3",&GENchargeL2_S3,"GENchargeL2_S3/I");
  tree->Branch("GENchargeL3_S3",&GENchargeL3_S3,"GENchargeL3_S3/I");
  tree->Branch("GENchargeL4_S3",&GENchargeL4_S3,"GENchargeL4_S3/I");
  tree->Branch("GENetaL1_S3",&GENetaL1_S3,"GENetaL1_S3/D");
  tree->Branch("GENetaL2_S3",&GENetaL2_S3,"GENetaL2_S3/D");
  tree->Branch("GENetaL3_S3",&GENetaL3_S3,"GENetaL3_S3/D");
  tree->Branch("GENetaL4_S3",&GENetaL4_S3,"GENetaL4_S3/D");
  tree->Branch("GENphiL1_S3",&GENphiL1_S3,"GENphiL1_S3/D");
  tree->Branch("GENphiL2_S3",&GENphiL2_S3,"GENphiL2_S3/D");
  tree->Branch("GENphiL3_S3",&GENphiL3_S3,"GENphiL3_S3/D");
  tree->Branch("GENphiL4_S3",&GENphiL4_S3,"GENphiL4_S3/D");

  tree->Branch("GENMH",&GENMH,"GENMH/D");
  tree->Branch("GENM4L",&GENM4L,"GENM4L/D");

  tree->Branch("m4l_GENMatched",&m4l_GENMatched,"m4l_GENMatched/D");

  tree->Branch("idL1_GENMatched",&idL1_GENMatched,"idL1_GENMatched/I");
  tree->Branch("idL2_GENMatched",&idL2_GENMatched,"idL2_GENMatched/I");
  tree->Branch("idL3_GENMatched",&idL3_GENMatched,"idL3_GENMatched/I");
  tree->Branch("idL4_GENMatched",&idL4_GENMatched,"idL4_GENMatched/I");
  tree->Branch("pTL1_GENMatched",&pTL1_GENMatched,"pTL1_GENMatched/D");
  tree->Branch("pTL2_GENMatched",&pTL2_GENMatched,"pTL2_GENMatched/D");
  tree->Branch("pTL3_GENMatched",&pTL3_GENMatched,"pTL3_GENMatched/D");
  tree->Branch("pTL4_GENMatched",&pTL4_GENMatched,"pTL4_GENMatched/D");
  tree->Branch("pXL1_GENMatched",&pXL1_GENMatched,"pXL1_GENMatched/D");
  tree->Branch("pXL2_GENMatched",&pXL2_GENMatched,"pXL2_GENMatched/D");
  tree->Branch("pXL3_GENMatched",&pXL3_GENMatched,"pXL3_GENMatched/D");
  tree->Branch("pXL4_GENMatched",&pXL4_GENMatched,"pXL4_GENMatched/D");
  tree->Branch("pYL1_GENMatched",&pYL1_GENMatched,"pYL1_GENMatched/D");
  tree->Branch("pYL2_GENMatched",&pYL2_GENMatched,"pYL2_GENMatched/D");
  tree->Branch("pYL3_GENMatched",&pYL3_GENMatched,"pYL3_GENMatched/D");
  tree->Branch("pYL4_GENMatched",&pYL4_GENMatched,"pYL4_GENMatched/D");
  tree->Branch("pZL1_GENMatched",&pZL1_GENMatched,"pZL1_GENMatched/D");
  tree->Branch("pZL2_GENMatched",&pZL2_GENMatched,"pZL2_GENMatched/D");
  tree->Branch("pZL3_GENMatched",&pZL3_GENMatched,"pZL3_GENMatched/D");
  tree->Branch("pZL4_GENMatched",&pZL4_GENMatched,"pZL4_GENMatched/D");
  tree->Branch("EL1_GENMatched",&EL1_GENMatched,"EL1_GENMatched/D");
  tree->Branch("EL2_GENMatched",&EL2_GENMatched,"EL2_GENMatched/D");
  tree->Branch("EL3_GENMatched",&EL3_GENMatched,"EL3_GENMatched/D");
  tree->Branch("EL4_GENMatched",&EL4_GENMatched,"EL4_GENMatched/D");
  tree->Branch("chargeL1_GENMatched",&chargeL1_GENMatched,"chargeL1_GENMatched/I");
  tree->Branch("chargeL2_GENMatched",&chargeL2_GENMatched,"chargeL2_GENMatched/I");
  tree->Branch("chargeL3_GENMatched",&chargeL3_GENMatched,"chargeL3_GENMatched/I");
  tree->Branch("chargeL4_GENMatched",&chargeL4_GENMatched,"chargeL4_GENMatched/I");
  tree->Branch("etaL1_GENMatched",&etaL1_GENMatched,"etaL1_GENMatched/D");
  tree->Branch("etaL2_GENMatched",&etaL2_GENMatched,"etaL2_GENMatched/D");
  tree->Branch("etaL3_GENMatched",&etaL3_GENMatched,"etaL3_GENMatched/D");
  tree->Branch("etaL4_GENMatched",&etaL4_GENMatched,"etaL4_GENMatched/D");
  tree->Branch("phiL1_GENMatched",&phiL1_GENMatched,"phiL1_GENMatched/D");
  tree->Branch("phiL2_GENMatched",&phiL2_GENMatched,"phiL2_GENMatched/D");
  tree->Branch("phiL3_GENMatched",&phiL3_GENMatched,"phiL3_GENMatched/D");
  tree->Branch("phiL4_GENMatched",&phiL4_GENMatched,"phiL4_GENMatched/D");

  tree->Branch("EZ1_GENMatched",&EZ1_GENMatched,"EZ1_GENMatched/D");
  tree->Branch("EZ2_GENMatched",&EZ2_GENMatched,"EZ2_GENMatched/D");
  tree->Branch("pTZ1_GENMatched",&pTZ1_GENMatched,"pTZ1_GENMatched/D");
  tree->Branch("pTZ2_GENMatched",&pTZ2_GENMatched,"pTZ2_GENMatched/D");
  tree->Branch("pXZ1_GENMatched",&pXZ1_GENMatched,"pXZ1_GENMatched/D");
  tree->Branch("pXZ2_GENMatched",&pXZ2_GENMatched,"pXZ2_GENMatched/D");
  tree->Branch("pYZ1_GENMatched",&pYZ1_GENMatched,"pYZ1_GENMatched/D");
  tree->Branch("pYZ2_GENMatched",&pYZ2_GENMatched,"pYZ2_GENMatched/D");
  tree->Branch("pZZ1_GENMatched",&pZZ1_GENMatched,"pZZ1_GENMatched/D");
  tree->Branch("pZZ2_GENMatched",&pZZ2_GENMatched,"pZZ2_GENMatched/D");

  tree->Branch("VBFJet1",&VBFJet1,"VBFJet1/O");
  tree->Branch("etaVBFJet1",&etaVBFJet1,"etaVBFJet1/D");
  tree->Branch("phiVBFJet1",&phiVBFJet1,"phiVBFJet1/D");
  tree->Branch("pTVBFJet1",&pTVBFJet1,"pTVBFJet1/D");
  tree->Branch("massVBFJet1",&massVBFJet1,"massVBFJet1/D");
  tree->Branch("VBFJet2",&VBFJet2,"VBFJet2/O");
  tree->Branch("etaVBFJet2",&etaVBFJet2,"etaVBFJet2/D");
  tree->Branch("phiVBFJet2",&phiVBFJet2,"phiVBFJet2/D");
  tree->Branch("pTVBFJet2",&pTVBFJet2,"pTVBFJet2/D");
  tree->Branch("massVBFJet2",&massVBFJet2,"massVBFJet2/D");
  tree->Branch("VBFDiJetMass",&VBFDiJetMass,"VBFDiJetMass/D");
  tree->Branch("VBFDeltaEta",&VBFDeltaEta,"VBFDeltaEta/D");
  tree->Branch("FisherDiscrim",&FisherDiscrim,"FisherDiscrim/D");
}



void UFHZZ4LAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                   std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, 
                                   std::vector<pat::Jet> selectedVBFJets, std::vector<pat::Jet> correctedVBFJets)
{
  using namespace pat;
  using namespace std;

  //Tree Variable
  if( RecoFourMuEvent ){ finalState = 1;}
  if( RecoFourEEvent  ){ finalState = 2;}
  if( RecoTwoETwoMuEvent ){ finalState = 3;}
  if( RecoTwoMuTwoEEvent ){ finalState = 4;}
  Run = iEvent.id().run();
  Event = iEvent.id().event();
  LumiSect = iEvent.id().luminosityBlock();
  mass4l = HiggsCandVec.M();
  m4lNoFSR = HiggsCandVecNoFSR.M();
  if(RecoFourMuEvent){mass4mu = HiggsCandVec.M();}
  else{mass4mu = -1;}
  if(RecoFourEEvent){ mass4e = HiggsCandVec.M();}
  else{ mass4e = -1;}
  if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){ mass2e2mu = HiggsCandVec.M(); }
  else{ mass2e2mu = -1;}
  pT4l = HiggsCandVec.Pt();
  massZ1 = Z1Vec.M();
  massZ2 = Z2Vec.M();
  eta4l = HiggsCandVec.Eta();
  phi4l = HiggsCandVec.Phi();
  muRho = muonRho;
  elRho = elecRho;
  worstSIP = highestSip;
  pTZ1 = Z1Vec.Pt();
  pXZ1 = Z1Vec.Px();
  pYZ1 = Z1Vec.Py();
  pZZ1 = Z1Vec.Pz();
  pTZ2 = Z2Vec.Pt();
  pXZ2 = Z2Vec.Px();
  pYZ2 = Z2Vec.Py();
  pZZ2 = Z2Vec.Pz();
  EZ1 = Z1Vec.E();
  EZ2 = Z2Vec.E();
  m4lErrUF = massErrorUF; m4lErrUFCorr = massErrorUFCorr;   m4lErrUFADCorr = massErrorUFADCorr;
  m4lErrUCSD = massErrorUCSD; 
  m4lErrUCSDCorr = massErrorUCSDCorr; 

  MC_weight = scaleWeight*eventWeight;
  pTL1FSR = Lep1.Pt();
  pXL1FSR = Lep1.Px();
  pYL1FSR = Lep1.Py();
  pZL1FSR = Lep1.Pz();
  EL1FSR = Lep1.E();
  pTL2FSR = Lep2.Pt();
  pXL2FSR = Lep2.Px();
  pYL2FSR = Lep2.Py();
  pZL2FSR = Lep2.Pz();
  EL2FSR = Lep2.E();
  pTL3FSR = Lep3.Pt();
  pXL3FSR = Lep3.Px();
  pYL3FSR = Lep3.Py();
  pZL3FSR = Lep3.Pz();
  EL3FSR = Lep3.E();
  pTL4FSR = Lep4.Pt();
  pXL4FSR = Lep4.Px();
  pYL4FSR = Lep4.Py();
  pZL4FSR = Lep4.Pz();
  EL4FSR = Lep4.E();

  if(RecoFourMuEvent) tmpEvent = "4mu";
  if(RecoTwoMuTwoEEvent) tmpEvent = "2mu2e";
  if(RecoTwoETwoMuEvent)  tmpEvent = "2e2mu";
  if(RecoFourEEvent)  tmpEvent = "4e";

  if( RecoFourMuEvent )
  {
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);
    missingHitsL1 = -1;

    idL2 = selectedMuons[1].pdgId();
    mvaL2 = selectedMuons[1].isPFMuon();
    EL2 = selectedMuons[1].energy();
    SipL2 = helper.getSIP3D(selectedMuons[1]);
    IPL2 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL2 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL2 = selectedMuons[1].pt();
    pXL2 = selectedMuons[1].px();
    pYL2 = selectedMuons[1].py();
    pZL2 = selectedMuons[1].pz();
    chargeL2 = selectedMuons[1].charge();
    etaL2 = selectedMuons[1].eta();
    phiL2 = selectedMuons[1].phi();
    isoNHL2 = selectedMuons[1].neutralHadronIso();
    isoCHL2 = selectedMuons[1].chargedHadronIso();
    isoPhotL2 = selectedMuons[1].photonIso();
    isoTrackL2 = selectedMuons[1].trackIso();
    isoEcalL2 = selectedMuons[1].ecalIso();
    isoHcalL2 = selectedMuons[1].hcalIso();
    RelIsoL2 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL2 = helper.pfIso(selectedMuons[1],0);
    missingHitsL2 = -1;

    idL3 = selectedMuons[2].pdgId();
    mvaL3 = selectedMuons[2].isPFMuon();
    EL3 = selectedMuons[2].energy();
    SipL3 = helper.getSIP3D(selectedMuons[2]);
    IPL3 = fabs(selectedMuons[2].dB(pat::Muon::PV3D));
    dIPL3 = selectedMuons[2].edB(pat::Muon::PV3D);
    pTL3 = selectedMuons[2].pt();
    pXL3 = selectedMuons[2].px();
    pYL3 = selectedMuons[2].py();
    pZL3 = selectedMuons[2].pz();
    chargeL3 = selectedMuons[2].charge();
    etaL3 = selectedMuons[2].eta();
    phiL3 = selectedMuons[2].phi();
    isoNHL3 = selectedMuons[2].neutralHadronIso();
    isoCHL3 = selectedMuons[2].chargedHadronIso();
    isoPhotL3 = selectedMuons[2].photonIso();
    isoTrackL3 = selectedMuons[2].trackIso();
    isoEcalL3 = selectedMuons[2].ecalIso();
    isoHcalL3 = selectedMuons[2].hcalIso();
    RelIsoL3 = helper.pfIso(selectedMuons[2],muonRho);
    RelIsoUCL3 = helper.pfIso(selectedMuons[2],0);
    missingHitsL3 = -1;

    idL4 = selectedMuons[3].pdgId();
    mvaL4 = selectedMuons[3].isPFMuon();
    IPL4 = fabs(selectedMuons[3].dB(pat::Muon::PV3D));
    dIPL4 = selectedMuons[3].edB(pat::Muon::PV3D);
    SipL4 = helper.getSIP3D(selectedMuons[3]);
    EL4 = selectedMuons[3].energy();
    pTL4 = selectedMuons[3].pt();
    pXL4 = selectedMuons[3].px();
    pYL4 = selectedMuons[3].py();
    pZL4 = selectedMuons[3].pz();
    chargeL4 = selectedMuons[3].charge();
    etaL4 = selectedMuons[3].eta();
    phiL4 = selectedMuons[3].phi();
    isoNHL4 = selectedMuons[3].neutralHadronIso();
    isoCHL4 = selectedMuons[3].chargedHadronIso();
    isoPhotL4 = selectedMuons[3].photonIso();
    isoTrackL4 = selectedMuons[3].trackIso();
    isoEcalL4 = selectedMuons[3].ecalIso();
    isoHcalL4 = selectedMuons[3].hcalIso();
    RelIsoL4 =  helper.pfIso(selectedMuons[3],muonRho);
    RelIsoUCL4 =  helper.pfIso(selectedMuons[3],0);
    missingHitsL4 = -1;
  }
  else if( RecoFourEEvent )
  {
    idL1 = selectedElectrons[0].pdgId();
    mvaL1 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID) ;
    EL1 = selectedElectrons[0].energy();
    SipL1 = helper.getSIP3D(selectedElectrons[0]);
    IPL1 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL1 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL1 = selectedElectrons[0].pt();
    pXL1 = selectedElectrons[0].px();
    pYL1 = selectedElectrons[0].py();
    pZL1 = selectedElectrons[0].pz();
    chargeL1 = selectedElectrons[0].charge();
    etaL1 = selectedElectrons[0].eta();
    phiL1 = selectedElectrons[0].phi();
    isoNHL1 = selectedElectrons[0].neutralHadronIso();
    isoCHL1 = selectedElectrons[0].chargedHadronIso();
    isoPhotL1 = selectedElectrons[0].photonIso();
    isoTrackL1 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL1 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL1 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL1 =  helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL1 = helper.pfIso(selectedElectrons[0],0);
    //missingHitsL1 = selectedElectrons[0].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
    ecalSeedL1 = selectedElectrons[0].ecalDrivenSeed();

    idL2 = selectedElectrons[1].pdgId();
    mvaL2 = elecID=="noEID" ? -100 : selectedElectrons[1].electronID(elecID);
    EL2 = selectedElectrons[1].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[1]);
    IPL2 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[1].pt();
    pXL2 = selectedElectrons[1].px();
    pYL2 = selectedElectrons[1].py();
    pZL2 = selectedElectrons[1].pz();
    chargeL2 = selectedElectrons[1].charge();
    etaL2 = selectedElectrons[1].eta();
    phiL2 = selectedElectrons[1].phi();
    isoNHL2 = selectedElectrons[1].neutralHadronIso();
    isoCHL2 = selectedElectrons[1].chargedHadronIso();
    isoPhotL2 = selectedElectrons[1].photonIso();
    isoTrackL2 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[1],0);
    //missingHitsL2 = selectedElectrons[1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
    ecalSeedL2 = selectedElectrons[1].ecalDrivenSeed();

    idL3 = selectedElectrons[2].pdgId();
    mvaL3 = elecID=="noEID" ? -100 : selectedElectrons[2].electronID(elecID);
    EL3 = selectedElectrons[2].energy();
    SipL3 = helper.getSIP3D(selectedElectrons[2]);
    IPL3 = fabs(selectedElectrons[2].dB(pat::Electron::PV3D));
    dIPL3 = selectedElectrons[2].edB(pat::Electron::PV3D);
    pTL3 = selectedElectrons[2].pt();
    pXL3 = selectedElectrons[2].px();
    pYL3 = selectedElectrons[2].py();
    pZL3 = selectedElectrons[2].pz();
    chargeL3 = selectedElectrons[2].charge();
    etaL3 = selectedElectrons[2].eta();
    phiL3 = selectedElectrons[2].phi();
    isoNHL3 = selectedElectrons[2].neutralHadronIso();
    isoCHL3 = selectedElectrons[2].chargedHadronIso();
    isoPhotL3 = selectedElectrons[2].photonIso();
    isoTrackL3 = selectedElectrons[2].dr03TkSumPt();
    isoEcalL3 = selectedElectrons[2].dr03EcalRecHitSumEt();
    isoHcalL3 = selectedElectrons[2].dr03HcalTowerSumEt();
    RelIsoL3 = helper.pfIso(selectedElectrons[2],elecRho);
    RelIsoUCL3 = helper.pfIso(selectedElectrons[2],0);
    //missingHitsL3 = selectedElectrons[2].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
    ecalSeedL3 = selectedElectrons[2].ecalDrivenSeed();

    idL4 = selectedElectrons[3].pdgId();
    mvaL4 = elecID=="noEID" ? -100 : selectedElectrons[3].electronID(elecID);
    EL4 = selectedElectrons[3].energy();
    SipL4 = helper.getSIP3D(selectedElectrons[3]);
    IPL4 = fabs(selectedElectrons[3].dB(pat::Electron::PV3D));
    dIPL4 = selectedElectrons[3].edB(pat::Electron::PV3D);
    pTL4 = selectedElectrons[3].pt();
    pXL4 = selectedElectrons[3].px();
    pYL4 = selectedElectrons[3].py();
    pZL4 = selectedElectrons[3].pz();
    chargeL4 = selectedElectrons[3].charge();
    etaL4 = selectedElectrons[3].eta();
    phiL4 = selectedElectrons[3].phi();
    isoNHL4 = selectedElectrons[3].neutralHadronIso();
    isoCHL4 = selectedElectrons[3].chargedHadronIso();
    isoPhotL4 = selectedElectrons[3].photonIso();
    isoTrackL4 = selectedElectrons[3].dr03TkSumPt();
    isoEcalL4 = selectedElectrons[3].dr03EcalRecHitSumEt();
    isoHcalL4 = selectedElectrons[3].dr03HcalTowerSumEt();
    RelIsoL4 = helper.pfIso(selectedElectrons[3],elecRho);
    RelIsoUCL4 = helper.pfIso(selectedElectrons[3],0);
    //missingHitsL4 = selectedElectrons[3].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
    ecalSeedL4 = selectedElectrons[3].ecalDrivenSeed();
  }
  else if( RecoTwoETwoMuEvent )
  {
    idL1 = selectedElectrons[0].pdgId();
    mvaL1 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL1 = selectedElectrons[0].energy();
    SipL1 = helper.getSIP3D(selectedElectrons[0]);
    IPL1 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL1 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL1 = selectedElectrons[0].pt();
    pXL1 = selectedElectrons[0].px();
    pYL1 = selectedElectrons[0].py();
    pZL1 = selectedElectrons[0].pz();
    chargeL1 = selectedElectrons[0].charge();
    etaL1 = selectedElectrons[0].eta();
    phiL1 = selectedElectrons[0].phi();
    isoNHL1 = selectedElectrons[0].neutralHadronIso();
    isoCHL1 = selectedElectrons[0].chargedHadronIso();
    isoPhotL1 = selectedElectrons[0].photonIso();
    isoTrackL1 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL1 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL1 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL1 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL1 = helper.pfIso(selectedElectrons[0],0);
    //missingHitsL1 = selectedElectrons[0].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); //for miniAOD
    ecalSeedL1 = selectedElectrons[0].ecalDrivenSeed();

    idL2 = selectedElectrons[1].pdgId();
    mvaL2 = elecID=="noEID" ? -100 : selectedElectrons[1].electronID(elecID);
    EL2 = selectedElectrons[1].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[1]);
    IPL2 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[1].pt();
    pXL2 = selectedElectrons[1].px();
    pYL2 = selectedElectrons[1].py();
    pZL2 = selectedElectrons[1].pz();
    chargeL2 = selectedElectrons[1].charge();
    etaL2 = selectedElectrons[1].eta();
    phiL2 = selectedElectrons[1].phi();
    isoNHL2 = selectedElectrons[1].neutralHadronIso();
    isoCHL2 = selectedElectrons[1].chargedHadronIso();
    isoPhotL2 = selectedElectrons[1].photonIso();
    isoTrackL2 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[1],0);
    //missingHitsL2 = selectedElectrons[1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
    ecalSeedL2 = selectedElectrons[1].ecalDrivenSeed();

    idL3 = selectedMuons[0].pdgId();
    mvaL3 = selectedMuons[0].isPFMuon();
    EL3 = selectedMuons[0].energy();
    SipL3 = helper.getSIP3D(selectedMuons[0]);
    IPL3 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL3 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL3 = selectedMuons[0].pt();
    pXL3 = selectedMuons[0].px();
    pYL3 = selectedMuons[0].py();
    pZL3 = selectedMuons[0].pz();
    chargeL3 = selectedMuons[0].charge();
    etaL3 = selectedMuons[0].eta();
    phiL3 = selectedMuons[0].phi();
    isoNHL3 = selectedMuons[0].neutralHadronIso();
    isoCHL3 = selectedMuons[0].chargedHadronIso();
    isoPhotL3 = selectedMuons[0].photonIso();
    isoTrackL3 = selectedMuons[0].trackIso();
    isoEcalL3 = selectedMuons[0].ecalIso();
    isoHcalL3 = selectedMuons[0].hcalIso();
    RelIsoL3 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL3 = helper.pfIso(selectedMuons[0],0);
    missingHitsL3 = -1;

    idL4 = selectedMuons[1].pdgId();
    mvaL4 = selectedMuons[1].isPFMuon();
    EL4 = selectedMuons[1].energy();
    SipL4 = helper.getSIP3D(selectedMuons[1]);
    IPL4 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL4 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL4 = selectedMuons[1].pt();
    pXL4 = selectedMuons[1].px();
    pYL4 = selectedMuons[1].py();
    pZL4 = selectedMuons[1].pz();
    chargeL4 = selectedMuons[1].charge();
    etaL4 = selectedMuons[1].eta();
    phiL4 = selectedMuons[1].phi();
    isoNHL4 = selectedMuons[1].neutralHadronIso();
    isoCHL4 = selectedMuons[1].chargedHadronIso();
    isoPhotL4 = selectedMuons[1].photonIso();
    isoTrackL4 = selectedMuons[1].trackIso();
    isoEcalL4 = selectedMuons[1].ecalIso();
    isoHcalL4 = selectedMuons[1].hcalIso();
    RelIsoL4 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL4 = helper.pfIso(selectedMuons[1],0);
    missingHitsL4 = -1;
  }
  else if( RecoTwoMuTwoEEvent )
  {
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);
    missingHitsL1 = -1;

    idL2 = selectedMuons[1].pdgId();
    mvaL2 = selectedMuons[1].isPFMuon();
    EL2 = selectedMuons[1].energy();
    SipL2 = helper.getSIP3D(selectedMuons[1]);
    IPL2 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL2 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL2 = selectedMuons[1].pt();
    pXL2 = selectedMuons[1].px();
    pYL2 = selectedMuons[1].py();
    pZL2 = selectedMuons[1].pz();
    chargeL2 = selectedMuons[1].charge();
    etaL2 = selectedMuons[1].eta();
    phiL2 = selectedMuons[1].phi();
    isoNHL2 = selectedMuons[1].neutralHadronIso();
    isoCHL2 = selectedMuons[1].chargedHadronIso();
    isoPhotL2 = selectedMuons[1].photonIso();
    isoTrackL2 = selectedMuons[1].trackIso();
    isoEcalL2 = selectedMuons[1].ecalIso();
    isoHcalL2 = selectedMuons[1].hcalIso();
    RelIsoL2 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL2 = helper.pfIso(selectedMuons[1],0);
    missingHitsL2 = -1;

    idL3 = selectedElectrons[0].pdgId();
    mvaL3 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL3 = selectedElectrons[0].energy();
    SipL3 = helper.getSIP3D(selectedElectrons[0]);
    IPL3 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL3 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL3 = selectedElectrons[0].pt();
    pXL3 = selectedElectrons[0].px();
    pYL3 = selectedElectrons[0].py();
    pZL3 = selectedElectrons[0].pz();
    chargeL3 = selectedElectrons[0].charge();
    etaL3 = selectedElectrons[0].eta();
    phiL3 = selectedElectrons[0].phi();
    isoNHL3 = selectedElectrons[0].neutralHadronIso();
    isoCHL3 = selectedElectrons[0].chargedHadronIso();
    isoPhotL3 = selectedElectrons[0].photonIso();
    isoTrackL3 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL3 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL3 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL3 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL3 = helper.pfIso(selectedElectrons[0],0);
    //missingHitsL3 = selectedElectrons[0].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); // for miniAOD
    ecalSeedL3 = selectedElectrons[0].ecalDrivenSeed();

    idL4 = selectedElectrons[1].pdgId();
    mvaL4 = elecID=="noEID" ? -100 : selectedElectrons[1].electronID(elecID);
    EL4 = selectedElectrons[1].energy();
    SipL4 = helper.getSIP3D(selectedElectrons[1]);
    IPL4 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL4 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL4 = selectedElectrons[1].pt();
    pXL4 = selectedElectrons[1].px();
    pYL4 = selectedElectrons[1].py();
    pZL4 = selectedElectrons[1].pz();
    chargeL4 = selectedElectrons[1].charge();
    etaL4 = selectedElectrons[1].eta();
    phiL4 = selectedElectrons[1].phi();
    isoNHL4 = selectedElectrons[1].neutralHadronIso();
    isoCHL4 = selectedElectrons[1].chargedHadronIso();
    isoPhotL4 = selectedElectrons[1].photonIso();
    isoTrackL4 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL4 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL4 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL4 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL4 = helper.pfIso(selectedElectrons[1],0);
    //missingHitsL4 = selectedElectrons[1].gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); //for miniAOD
    ecalSeedL4 = selectedElectrons[1].ecalDrivenSeed();
  } 
  else if (RecoFourMixEvent==1) // mueemu
  {      
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);

    idL2 = selectedElectrons[0].pdgId();
    mvaL2 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL2 = selectedElectrons[0].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[0]);
    IPL2 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[0].pt();
    pXL2 = selectedElectrons[0].px();
    pYL2 = selectedElectrons[0].py();
    pZL2 = selectedElectrons[0].pz();
    chargeL2 = selectedElectrons[0].charge();
    etaL2 = selectedElectrons[0].eta();
    phiL2 = selectedElectrons[0].phi();
    isoNHL2 = selectedElectrons[0].neutralHadronIso();
    isoCHL2 = selectedElectrons[0].chargedHadronIso();
    isoPhotL2 = selectedElectrons[0].photonIso();
    isoTrackL2 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[0],0);

    idL3 = selectedElectrons[1].pdgId();
    mvaL3 = elecID=="noEID" ? -100 : selectedElectrons[1].electronID(elecID);
    EL3 = selectedElectrons[1].energy();
    SipL3 = helper.getSIP3D(selectedElectrons[1]);
    IPL3 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL3 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL3 = selectedElectrons[1].pt();
    pXL3 = selectedElectrons[1].px();
    pYL3 = selectedElectrons[1].py();
    pZL3 = selectedElectrons[1].pz();
    chargeL3 = selectedElectrons[1].charge();
    etaL3 = selectedElectrons[1].eta();
    phiL3 = selectedElectrons[1].phi();
    isoNHL3 = selectedElectrons[1].neutralHadronIso();
    isoCHL3 = selectedElectrons[1].chargedHadronIso();
    isoPhotL3 = selectedElectrons[1].photonIso();
    isoTrackL3 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL3 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL3 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL3 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL3 = helper.pfIso(selectedElectrons[1],0);

    idL4 = selectedMuons[1].pdgId();
    mvaL4 = selectedMuons[1].isPFMuon();
    EL4 = selectedMuons[1].energy();
    SipL4 = helper.getSIP3D(selectedMuons[1]);
    IPL4 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL4 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL4 = selectedMuons[1].pt();
    pXL4 = selectedMuons[1].px();
    pYL4 = selectedMuons[1].py();
    pZL4 = selectedMuons[1].pz();
    chargeL4 = selectedMuons[1].charge();
    etaL4 = selectedMuons[1].eta();
    phiL4 = selectedMuons[1].phi();
    isoNHL4 = selectedMuons[1].neutralHadronIso();
    isoCHL4 = selectedMuons[1].chargedHadronIso();
    isoPhotL4 = selectedMuons[1].photonIso();
    isoTrackL4 = selectedMuons[1].trackIso();
    isoEcalL4 = selectedMuons[1].ecalIso();
    isoHcalL4 = selectedMuons[1].hcalIso();
    RelIsoL4 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL4 = helper.pfIso(selectedMuons[1],0); 
  } 
  else if (RecoFourMixEvent==2) // muemue
  {      
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);

    idL2 = selectedElectrons[0].pdgId();
    mvaL2 =  elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL2 = selectedElectrons[0].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[0]);
    IPL2 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[0].pt();
    pXL2 = selectedElectrons[0].px();
    pYL2 = selectedElectrons[0].py();
    pZL2 = selectedElectrons[0].pz();
    chargeL2 = selectedElectrons[0].charge();
    etaL2 = selectedElectrons[0].eta();
    phiL2 = selectedElectrons[0].phi();
    isoNHL2 = selectedElectrons[0].neutralHadronIso();
    isoCHL2 = selectedElectrons[0].chargedHadronIso();
    isoPhotL2 = selectedElectrons[0].photonIso();
    isoTrackL2 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[0],0);

    idL3 = selectedMuons[1].pdgId();
    mvaL3 = selectedMuons[1].isPFMuon();
    EL3 = selectedMuons[1].energy();
    SipL3 = helper.getSIP3D(selectedMuons[1]);
    IPL3 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL3 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL3 = selectedMuons[1].pt();
    pXL3 = selectedMuons[1].px();
    pYL3 = selectedMuons[1].py();
    pZL3 = selectedMuons[1].pz();
    chargeL3 = selectedMuons[1].charge();
    etaL3 = selectedMuons[1].eta();
    phiL3 = selectedMuons[1].phi();
    isoNHL3 = selectedMuons[1].neutralHadronIso();
    isoCHL3 = selectedMuons[1].chargedHadronIso();
    isoPhotL3 = selectedMuons[1].photonIso();
    isoTrackL3 = selectedMuons[1].trackIso();
    isoEcalL3 = selectedMuons[1].ecalIso();
    isoHcalL3 = selectedMuons[1].hcalIso();
    RelIsoL3 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL3 = helper.pfIso(selectedMuons[1],0);

    idL4 = selectedElectrons[1].pdgId();
    mvaL4 = elecID=="noEID" ? -100 : selectedElectrons[1].electronID(elecID);
    EL4 = selectedElectrons[1].energy();
    SipL4 = helper.getSIP3D(selectedElectrons[1]);
    IPL4 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL4 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL4 = selectedElectrons[1].pt();
    pXL4 = selectedElectrons[1].px();
    pYL4 = selectedElectrons[1].py();
    pZL4 = selectedElectrons[1].pz();
    chargeL4 = selectedElectrons[1].charge();
    etaL4 = selectedElectrons[1].eta();
    phiL4 = selectedElectrons[1].phi();
    isoNHL4 = selectedElectrons[1].neutralHadronIso();
    isoCHL4 = selectedElectrons[1].chargedHadronIso();
    isoPhotL4 = selectedElectrons[1].photonIso();
    isoTrackL4 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL4 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL4 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL4 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL4 = helper.pfIso(selectedElectrons[1],0); 
  } 
  else if (RecoFourMixEvent==3) // mueee
  {      
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);

    idL2 = selectedElectrons[0].pdgId();
    mvaL2 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL2 = selectedElectrons[0].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[0]);
    IPL2 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[0].pt();
    pXL2 = selectedElectrons[0].px();
    pYL2 = selectedElectrons[0].py();
    pZL2 = selectedElectrons[0].pz();
    chargeL2 = selectedElectrons[0].charge();
    etaL2 = selectedElectrons[0].eta();
    phiL2 = selectedElectrons[0].phi();
    isoNHL2 = selectedElectrons[0].neutralHadronIso();
    isoCHL2 = selectedElectrons[0].chargedHadronIso();
    isoPhotL2 = selectedElectrons[0].photonIso();
    isoTrackL2 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[0],0);

    idL3 = selectedElectrons[1].pdgId();
    mvaL3 = selectedElectrons[1].electronID(elecID);
    EL3 = selectedElectrons[1].energy();
    SipL3 = helper.getSIP3D(selectedElectrons[1]);
    IPL3 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL3 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL3 = selectedElectrons[1].pt();
    pXL3 = selectedElectrons[1].px();
    pYL3 = selectedElectrons[1].py();
    pZL3 = selectedElectrons[1].pz();
    chargeL3 = selectedElectrons[1].charge();
    etaL3 = selectedElectrons[1].eta();
    phiL3 = selectedElectrons[1].phi();
    isoNHL3 = selectedElectrons[1].neutralHadronIso();
    isoCHL3 = selectedElectrons[1].chargedHadronIso();
    isoPhotL3 = selectedElectrons[1].photonIso();
    isoTrackL3 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL3 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL3 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL3 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL3 = helper.pfIso(selectedElectrons[1],0);

    idL4 = selectedElectrons[2].pdgId();
    mvaL4 = elecID=="noEID" ? -100 : selectedElectrons[2].electronID(elecID);
    EL4 = selectedElectrons[2].energy();
    SipL4 = helper.getSIP3D(selectedElectrons[2]);
    IPL4 = fabs(selectedElectrons[2].dB(pat::Electron::PV3D));
    dIPL4 = selectedElectrons[2].edB(pat::Electron::PV3D);
    pTL4 = selectedElectrons[2].pt();
    pXL4 = selectedElectrons[2].px();
    pYL4 = selectedElectrons[2].py();
    pZL4 = selectedElectrons[2].pz();
    chargeL4 = selectedElectrons[2].charge();
    etaL4 = selectedElectrons[2].eta();
    phiL4 = selectedElectrons[2].phi();
    isoNHL4 = selectedElectrons[2].neutralHadronIso();
    isoCHL4 = selectedElectrons[2].chargedHadronIso();
    isoPhotL4 = selectedElectrons[2].photonIso();
    isoTrackL4 = selectedElectrons[2].dr03TkSumPt();
    isoEcalL4 = selectedElectrons[2].dr03EcalRecHitSumEt();
    isoHcalL4 = selectedElectrons[2].dr03HcalTowerSumEt();
    RelIsoL4 = helper.pfIso(selectedElectrons[2],elecRho);
    RelIsoUCL4 = helper.pfIso(selectedElectrons[2],0); 
  } 
  else if (RecoFourMixEvent==4) // muemumu
  {      
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);

    idL2 = selectedElectrons[0].pdgId();
    mvaL2 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL2 = selectedElectrons[0].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[0]);
    IPL2 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[0].pt();
    pXL2 = selectedElectrons[0].px();
    pYL2 = selectedElectrons[0].py();
    pZL2 = selectedElectrons[0].pz();
    chargeL2 = selectedElectrons[0].charge();
    etaL2 = selectedElectrons[0].eta();
    phiL2 = selectedElectrons[0].phi();
    isoNHL2 = selectedElectrons[0].neutralHadronIso();
    isoCHL2 = selectedElectrons[0].chargedHadronIso();
    isoPhotL2 = selectedElectrons[0].photonIso();
    isoTrackL2 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[0],0);

    idL3 = selectedMuons[1].pdgId();
    mvaL3 = selectedMuons[1].isPFMuon();
    EL3 = selectedMuons[1].energy();
    SipL3 = helper.getSIP3D(selectedMuons[1]);
    IPL3 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL3 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL3 = selectedMuons[1].pt();
    pXL3 = selectedMuons[1].px();
    pYL3 = selectedMuons[1].py();
    pZL3 = selectedMuons[1].pz();
    chargeL3 = selectedMuons[1].charge();
    etaL3 = selectedMuons[1].eta();
    phiL3 = selectedMuons[1].phi();
    isoNHL3 = selectedMuons[1].neutralHadronIso();
    isoCHL3 = selectedMuons[1].chargedHadronIso();
    isoPhotL3 = selectedMuons[1].photonIso();
    isoTrackL3 = selectedMuons[1].trackIso();
    isoEcalL3 = selectedMuons[1].ecalIso();
    isoHcalL3 = selectedMuons[1].hcalIso();
    RelIsoL3 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL3 = helper.pfIso(selectedMuons[1],0);

    idL4 = selectedMuons[2].pdgId();
    mvaL4 = selectedMuons[2].isPFMuon();
    EL4 = selectedMuons[2].energy();
    SipL4 = helper.getSIP3D(selectedMuons[2]);
    IPL4 = fabs(selectedMuons[2].dB(pat::Muon::PV3D));
    dIPL4 = selectedMuons[2].edB(pat::Muon::PV3D);
    pTL4 = selectedMuons[2].pt();
    pXL4 = selectedMuons[2].px();
    pYL4 = selectedMuons[2].py();
    pZL4 = selectedMuons[2].pz();
    chargeL4 = selectedMuons[2].charge();
    etaL4 = selectedMuons[2].eta();
    phiL4 = selectedMuons[2].phi();
    isoNHL4 = selectedMuons[2].neutralHadronIso();
    isoCHL4 = selectedMuons[2].chargedHadronIso();
    isoPhotL4 = selectedMuons[2].photonIso();
    isoTrackL4 = selectedMuons[2].trackIso();
    isoEcalL4 = selectedMuons[2].ecalIso();
    isoHcalL4 = selectedMuons[2].hcalIso();
    RelIsoL4 = helper.pfIso(selectedMuons[2],muonRho);
    RelIsoUCL4 = helper.pfIso(selectedMuons[2],0);
  } 
  else if (RecoFourMixEvent==5) // eeemu
  {      
    idL1 = selectedElectrons[0].pdgId();
    mvaL1 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL1 = selectedElectrons[0].energy();
    SipL1 = helper.getSIP3D(selectedElectrons[0]);
    IPL1 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL1 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL1 = selectedElectrons[0].pt();
    pXL1 = selectedElectrons[0].px();
    pYL1 = selectedElectrons[0].py();
    pZL1 = selectedElectrons[0].pz();
    chargeL1 = selectedElectrons[0].charge();
    etaL1 = selectedElectrons[0].eta();
    phiL1 = selectedElectrons[0].phi();
    isoNHL1 = selectedElectrons[0].neutralHadronIso();
    isoCHL1 = selectedElectrons[0].chargedHadronIso();
    isoPhotL1 = selectedElectrons[0].photonIso();
    isoTrackL1 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL1 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL1 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL1 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL1 = helper.pfIso(selectedElectrons[0],0);

    idL2 = selectedElectrons[1].pdgId();
    mvaL2 = elecID=="noEID" ? -100 : selectedElectrons[1].electronID(elecID);
    EL2 = selectedElectrons[1].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[1]);
    IPL2 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[1].pt();
    pXL2 = selectedElectrons[1].px();
    pYL2 = selectedElectrons[1].py();
    pZL2 = selectedElectrons[1].pz();
    chargeL2 = selectedElectrons[1].charge();
    etaL2 = selectedElectrons[1].eta();
    phiL2 = selectedElectrons[1].phi();
    isoNHL2 = selectedElectrons[1].neutralHadronIso();
    isoCHL2 = selectedElectrons[1].chargedHadronIso();
    isoPhotL2 = selectedElectrons[1].photonIso();
    isoTrackL2 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[1],0);

    idL3 = selectedElectrons[2].pdgId();
    mvaL3 = elecID=="noEID" ? -100 : selectedElectrons[2].electronID(elecID);
    EL3 = selectedElectrons[2].energy();
    SipL3 = helper.getSIP3D(selectedElectrons[2]);
    IPL3 = fabs(selectedElectrons[2].dB(pat::Electron::PV3D));
    dIPL3 = selectedElectrons[2].edB(pat::Electron::PV3D);
    pTL3 = selectedElectrons[2].pt();
    pXL3 = selectedElectrons[2].px();
    pYL3 = selectedElectrons[2].py();
    pZL3 = selectedElectrons[2].pz();
    chargeL3 = selectedElectrons[2].charge();
    etaL3 = selectedElectrons[2].eta();
    phiL3 = selectedElectrons[2].phi();
    isoNHL3 = selectedElectrons[2].neutralHadronIso();
    isoCHL3 = selectedElectrons[2].chargedHadronIso();
    isoPhotL3 = selectedElectrons[2].photonIso();
    isoTrackL3 = selectedElectrons[2].dr03TkSumPt();
    isoEcalL3 = selectedElectrons[2].dr03EcalRecHitSumEt();
    isoHcalL3 = selectedElectrons[2].dr03HcalTowerSumEt();
    RelIsoL3 = helper.pfIso(selectedElectrons[2],elecRho);
    RelIsoUCL3 = helper.pfIso(selectedElectrons[2],0);

    idL4 = selectedMuons[0].pdgId();
    mvaL4 = selectedMuons[0].isPFMuon();
    EL4 = selectedMuons[0].energy();
    SipL4 = helper.getSIP3D(selectedMuons[0]);
    IPL4 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL4 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL4 = selectedMuons[0].pt();
    pXL4 = selectedMuons[0].px();
    pYL4 = selectedMuons[0].py();
    pZL4 = selectedMuons[0].pz();
    chargeL4 = selectedMuons[0].charge();
    etaL4 = selectedMuons[0].eta();
    phiL4 = selectedMuons[0].phi();
    isoNHL4 = selectedMuons[0].neutralHadronIso();
    isoCHL4 = selectedMuons[0].chargedHadronIso();
    isoPhotL4 = selectedMuons[0].photonIso();
    isoTrackL4 = selectedMuons[0].trackIso();
    isoEcalL4 = selectedMuons[0].ecalIso();
    isoHcalL4 = selectedMuons[0].hcalIso();
    RelIsoL4 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL4 = helper.pfIso(selectedMuons[0],0); 
  } 
  else if (RecoFourMixEvent==6) // eemue
  {      
    idL1 = selectedElectrons[0].pdgId();
    mvaL1 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL1 = selectedElectrons[0].energy();
    SipL1 = helper.getSIP3D(selectedElectrons[0]);
    IPL1 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL1 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL1 = selectedElectrons[0].pt();
    pXL1 = selectedElectrons[0].px();
    pYL1 = selectedElectrons[0].py();
    pZL1 = selectedElectrons[0].pz();
    chargeL1 = selectedElectrons[0].charge();
    etaL1 = selectedElectrons[0].eta();
    phiL1 = selectedElectrons[0].phi();
    isoNHL1 = selectedElectrons[0].neutralHadronIso();
    isoCHL1 = selectedElectrons[0].chargedHadronIso();
    isoPhotL1 = selectedElectrons[0].photonIso();
    isoTrackL1 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL1 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL1 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL1 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL1 = helper.pfIso(selectedElectrons[0],0);

    idL2 = selectedElectrons[1].pdgId();
    mvaL2 = elecID=="noEID" ? -100 : selectedElectrons[1].electronID(elecID);
    EL2 = selectedElectrons[1].energy();
    SipL2 = helper.getSIP3D(selectedElectrons[1]);
    IPL2 = fabs(selectedElectrons[1].dB(pat::Electron::PV3D));
    dIPL2 = selectedElectrons[1].edB(pat::Electron::PV3D);
    pTL2 = selectedElectrons[1].pt();
    pXL2 = selectedElectrons[1].px();
    pYL2 = selectedElectrons[1].py();
    pZL2 = selectedElectrons[1].pz();
    chargeL2 = selectedElectrons[1].charge();
    etaL2 = selectedElectrons[1].eta();
    phiL2 = selectedElectrons[1].phi();
    isoNHL2 = selectedElectrons[1].neutralHadronIso();
    isoCHL2 = selectedElectrons[1].chargedHadronIso();
    isoPhotL2 = selectedElectrons[1].photonIso();
    isoTrackL2 = selectedElectrons[1].dr03TkSumPt();
    isoEcalL2 = selectedElectrons[1].dr03EcalRecHitSumEt();
    isoHcalL2 = selectedElectrons[1].dr03HcalTowerSumEt();
    RelIsoL2 = helper.pfIso(selectedElectrons[1],elecRho);
    RelIsoUCL2 = helper.pfIso(selectedElectrons[1],0);

    idL3 = selectedMuons[0].pdgId();
    mvaL3 = selectedMuons[0].isPFMuon();
    EL3 = selectedMuons[0].energy();
    SipL3 = helper.getSIP3D(selectedMuons[0]);
    IPL3 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL3 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL3 = selectedMuons[0].pt();
    pXL3 = selectedMuons[0].px();
    pYL3 = selectedMuons[0].py();
    pZL3 = selectedMuons[0].pz();
    chargeL3 = selectedMuons[0].charge();
    etaL3 = selectedMuons[0].eta();
    phiL3 = selectedMuons[0].phi();
    isoNHL3 = selectedMuons[0].neutralHadronIso();
    isoCHL3 = selectedMuons[0].chargedHadronIso();
    isoPhotL3 = selectedMuons[0].photonIso();
    isoTrackL3 = selectedMuons[0].trackIso();
    isoEcalL3 = selectedMuons[0].ecalIso();
    isoHcalL3 = selectedMuons[0].hcalIso();
    RelIsoL3 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL3 = helper.pfIso(selectedMuons[0],0);

    idL4 = selectedElectrons[2].pdgId();
    mvaL4 = elecID=="noEID" ? -100 : selectedElectrons[2].electronID(elecID);
    EL4 = selectedElectrons[2].energy();
    SipL4 = helper.getSIP3D(selectedElectrons[2]);
    IPL4 = fabs(selectedElectrons[2].dB(pat::Electron::PV3D));
    dIPL4 = selectedElectrons[2].edB(pat::Electron::PV3D);
    pTL4 = selectedElectrons[2].pt();
    pXL4 = selectedElectrons[2].px();
    pYL4 = selectedElectrons[2].py();
    pZL4 = selectedElectrons[2].pz();
    chargeL4 = selectedElectrons[2].charge();
    etaL4 = selectedElectrons[2].eta();
    phiL4 = selectedElectrons[2].phi();
    isoNHL4 = selectedElectrons[2].neutralHadronIso();
    isoCHL4 = selectedElectrons[2].chargedHadronIso();
    isoPhotL4 = selectedElectrons[2].photonIso();
    isoTrackL4 = selectedElectrons[2].dr03TkSumPt();
    isoEcalL4 = selectedElectrons[2].dr03EcalRecHitSumEt();
    isoHcalL4 = selectedElectrons[2].dr03HcalTowerSumEt();
    RelIsoL4 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL4 = helper.pfIso(selectedElectrons[0],0);
  } 
  else if (RecoFourMixEvent==7) // mumuemu
  {      
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);

    idL2 = selectedMuons[1].pdgId();
    mvaL2 = selectedMuons[1].isPFMuon();
    EL2 = selectedMuons[1].energy();
    SipL2 = helper.getSIP3D(selectedMuons[1]);
    IPL2 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL2 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL2 = selectedMuons[1].pt();
    pXL2 = selectedMuons[1].px();
    pYL2 = selectedMuons[1].py();
    pZL2 = selectedMuons[1].pz();
    chargeL2 = selectedMuons[1].charge();
    etaL2 = selectedMuons[1].eta();
    phiL2 = selectedMuons[1].phi();
    isoNHL2 = selectedMuons[1].neutralHadronIso();
    isoCHL2 = selectedMuons[1].chargedHadronIso();
    isoPhotL2 = selectedMuons[1].photonIso();
    isoTrackL2 = selectedMuons[1].trackIso();
    isoEcalL2 = selectedMuons[1].ecalIso();
    isoHcalL2 = selectedMuons[1].hcalIso();
    RelIsoL2 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL2 = helper.pfIso(selectedMuons[1],0);

    idL3 = selectedElectrons[0].pdgId();
    mvaL3 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL3 = selectedElectrons[0].energy();
    SipL3 = helper.getSIP3D(selectedElectrons[0]);
    IPL3 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL3 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL3 = selectedElectrons[0].pt();
    pXL3 = selectedElectrons[0].px();
    pYL3 = selectedElectrons[0].py();
    pZL3 = selectedElectrons[0].pz();
    chargeL3 = selectedElectrons[0].charge();
    etaL3 = selectedElectrons[0].eta();
    phiL3 = selectedElectrons[0].phi();
    isoNHL3 = selectedElectrons[0].neutralHadronIso();
    isoCHL3 = selectedElectrons[0].chargedHadronIso();
    isoPhotL3 = selectedElectrons[0].photonIso();
    isoTrackL3 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL3 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL3 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL3 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL3 = helper.pfIso(selectedElectrons[0],0);

    idL4 = selectedMuons[2].pdgId();
    mvaL4 = selectedMuons[2].isPFMuon();
    EL4 = selectedMuons[2].energy();
    SipL4 = helper.getSIP3D(selectedMuons[2]);
    IPL4 = fabs(selectedMuons[2].dB(pat::Muon::PV3D));
    dIPL4 = selectedMuons[2].edB(pat::Muon::PV3D);
    pTL4 = selectedMuons[2].pt();
    pXL4 = selectedMuons[2].px();
    pYL4 = selectedMuons[2].py();
    pZL4 = selectedMuons[2].pz();
    chargeL4 = selectedMuons[2].charge();
    etaL4 = selectedMuons[2].eta();
    phiL4 = selectedMuons[2].phi();
    isoNHL4 = selectedMuons[2].neutralHadronIso();
    isoCHL4 = selectedMuons[2].chargedHadronIso();
    isoPhotL4 = selectedMuons[2].photonIso();
    isoTrackL4 = selectedMuons[2].trackIso();
    isoEcalL4 = selectedMuons[2].ecalIso();
    isoHcalL4 = selectedMuons[2].hcalIso();
    RelIsoL4 = helper.pfIso(selectedMuons[2],muonRho);
    RelIsoUCL4 = helper.pfIso(selectedMuons[2],0);
  } 
  else if (RecoFourMixEvent==8) // mumumue
  {      
    idL1 = selectedMuons[0].pdgId();
    mvaL1 = selectedMuons[0].isPFMuon();
    EL1 = selectedMuons[0].energy();
    SipL1 = helper.getSIP3D(selectedMuons[0]);
    IPL1 = fabs(selectedMuons[0].dB(pat::Muon::PV3D));
    dIPL1 = selectedMuons[0].edB(pat::Muon::PV3D);
    pTL1 = selectedMuons[0].pt();
    pXL1 = selectedMuons[0].px();
    pYL1 = selectedMuons[0].py();
    pZL1 = selectedMuons[0].pz();
    chargeL1 = selectedMuons[0].charge();
    etaL1 = selectedMuons[0].eta();
    phiL1 = selectedMuons[0].phi();
    isoNHL1 = selectedMuons[0].neutralHadronIso();
    isoCHL1 = selectedMuons[0].chargedHadronIso();
    isoPhotL1 = selectedMuons[0].photonIso();
    isoTrackL1 = selectedMuons[0].trackIso();
    isoEcalL1 = selectedMuons[0].ecalIso();
    isoHcalL1 = selectedMuons[0].hcalIso();
    RelIsoL1 = helper.pfIso(selectedMuons[0],muonRho);
    RelIsoUCL1 = helper.pfIso(selectedMuons[0],0);

    idL2 = selectedMuons[1].pdgId();
    mvaL2 = selectedMuons[1].isPFMuon();
    EL2 = selectedMuons[1].energy();
    SipL2 = helper.getSIP3D(selectedMuons[1]);
    IPL2 = fabs(selectedMuons[1].dB(pat::Muon::PV3D));
    dIPL2 = selectedMuons[1].edB(pat::Muon::PV3D);
    pTL2 = selectedMuons[1].pt();
    pXL2 = selectedMuons[1].px();
    pYL2 = selectedMuons[1].py();
    pZL2 = selectedMuons[1].pz();
    chargeL2 = selectedMuons[1].charge();
    etaL2 = selectedMuons[1].eta();
    phiL2 = selectedMuons[1].phi();
    isoNHL2 = selectedMuons[1].neutralHadronIso();
    isoCHL2 = selectedMuons[1].chargedHadronIso();
    isoPhotL2 = selectedMuons[1].photonIso();
    isoTrackL2 = selectedMuons[1].trackIso();
    isoEcalL2 = selectedMuons[1].ecalIso();
    isoHcalL2 = selectedMuons[1].hcalIso();
    RelIsoL2 = helper.pfIso(selectedMuons[1],muonRho);
    RelIsoUCL2 = helper.pfIso(selectedMuons[1],0);

    idL3 = selectedMuons[2].pdgId();
    mvaL3 = selectedMuons[2].isPFMuon();
    EL3 = selectedMuons[2].energy();
    SipL3 = helper.getSIP3D(selectedMuons[2]);
    IPL3 = fabs(selectedMuons[2].dB(pat::Muon::PV3D));
    dIPL3 = selectedMuons[2].edB(pat::Muon::PV3D);
    pTL3 = selectedMuons[2].pt();
    pXL3 = selectedMuons[2].px();
    pYL3 = selectedMuons[2].py();
    pZL3 = selectedMuons[2].pz();
    chargeL3 = selectedMuons[2].charge();
    etaL3 = selectedMuons[2].eta();
    phiL3 = selectedMuons[2].phi();
    isoNHL3 = selectedMuons[2].neutralHadronIso();
    isoCHL3 = selectedMuons[2].chargedHadronIso();
    isoPhotL3 = selectedMuons[2].photonIso();
    isoTrackL3 = selectedMuons[2].trackIso();
    isoEcalL3 = selectedMuons[2].ecalIso();
    isoHcalL3 = selectedMuons[2].hcalIso();
    RelIsoL3 = helper.pfIso(selectedMuons[2],muonRho);
    RelIsoUCL3 = helper.pfIso(selectedMuons[2],0);

    idL4 = selectedElectrons[0].pdgId();
    mvaL4 = elecID=="noEID" ? -100 : selectedElectrons[0].electronID(elecID);
    EL4 = selectedElectrons[0].energy();
    SipL4 = helper.getSIP3D(selectedElectrons[0]);
    IPL4 = fabs(selectedElectrons[0].dB(pat::Electron::PV3D));
    dIPL4 = selectedElectrons[0].edB(pat::Electron::PV3D);
    pTL4 = selectedElectrons[0].pt();
    pXL4 = selectedElectrons[0].px();
    pYL4 = selectedElectrons[0].py();
    pZL4 = selectedElectrons[0].pz();
    chargeL4 = selectedElectrons[0].charge();
    etaL4 = selectedElectrons[0].eta();
    phiL4 = selectedElectrons[0].phi();
    isoNHL4 = selectedElectrons[0].neutralHadronIso();
    isoCHL4 = selectedElectrons[0].chargedHadronIso();
    isoPhotL4 = selectedElectrons[0].photonIso();
    isoTrackL4 = selectedElectrons[0].dr03TkSumPt();
    isoEcalL4 = selectedElectrons[0].dr03EcalRecHitSumEt();
    isoHcalL4 = selectedElectrons[0].dr03HcalTowerSumEt();
    RelIsoL4 = helper.pfIso(selectedElectrons[0],elecRho);
    RelIsoUCL4 = helper.pfIso(selectedElectrons[0],0);
  } 
  else
  {      
    idL1 = -999;
    EL1 = -999.99;
    SipL1 = -999.99;
    pTL1 = -999.99;
    pXL1 = -999.99;
    pYL1 = -999.99;
    pZL1 = -999.99;
    chargeL1 = -999;
    etaL1 = -999.99;
    phiL1 = -999.99;
    isoTrackL1 = -999.99;
    isoEcalL1 = -999.99;
    isoHcalL1 = -999.99;
    RelIsoL1 = -999.99;
    RelIsoUCL1 = -999.99;

    idL2 = -999;
    EL2 = -999.99;
    SipL2 = -999.99;
    pTL2 = -999.99;
    pXL2 = -999.99;
    pYL2 = -999.99;
    pZL2 = -999.99;
    chargeL2 = -999;
    etaL2 = -999.99;
    phiL2 = -999.99;
    isoTrackL2 = -999.99;
    isoEcalL2 = -999.99;
    isoHcalL2 = -999.99;
    RelIsoL2 = -999.99;
    RelIsoUCL2 = -999.99;

    idL3 = -999;
    EL3 = -999.99;
    SipL3 = -999.99;
    pTL3 = -999.99;
    pXL3 = -999.99;
    pYL3 = -999.99;
    pZL3 = -999.99;
    chargeL3 = -999;
    etaL3 = -999.99;
    phiL3 = -999.99;
    isoTrackL3 = -999.99;
    isoEcalL3 = -999.99;
    isoHcalL3 = -999.99;
    RelIsoL3 = -999.99;
    RelIsoUCL3 = -999.99;
          
    idL4 = -999;
    EL4 = -999.99;
    SipL4 = -999.99;
    pTL4 = -999.99;
    pXL4 = -999.99;
    pYL4 = -999.99;
    pZL4 = -999.99;
    chargeL4 = -999;
    etaL4 = -999.99;
    phiL4 = -999.99;
    isoTrackL4 = -999.99;
    isoEcalL4 = -999.99;
    isoHcalL4 = -999.99;
    RelIsoL4 = -999.99;
    RelIsoUCL4 = -999.99; 
  } 

  // Hengne: There are still 4 combinations leading with electron that are not considered above:
  // e m m m 
  // e m m e
  // e m e m
  // e m e e
  // Why not consider them?


  double tempDeltaR = 999;
  vector<pat::Jet> finalVBFJets;
  for( unsigned int k = 0; k < selectedVBFJets.size(); k++)
  {
    bool isDeltaR = true;
    //check for overlap with leptons and gammas
    tempDeltaR = deltaR(selectedVBFJets[k].eta(),selectedVBFJets[k].phi(),etaL1,phiL1);
    if (tempDeltaR < 0.5) isDeltaR = false;
    tempDeltaR = deltaR(selectedVBFJets[k].eta(),selectedVBFJets[k].phi(),etaL2,phiL2);
    if (tempDeltaR < 0.5) isDeltaR = false;
    tempDeltaR = deltaR(selectedVBFJets[k].eta(),selectedVBFJets[k].phi(),etaL3,phiL3);
    if (tempDeltaR < 0.5) isDeltaR = false;
    tempDeltaR = deltaR(selectedVBFJets[k].eta(),selectedVBFJets[k].phi(),etaL4,phiL4);
    if (tempDeltaR < 0.5) isDeltaR = false;
    
    if (isDeltaR)
    {
      if (correctedVBFJets[k].pt() > 30 && fabs(selectedVBFJets[k].eta()) < 4.7)
      {
        finalVBFJets.push_back(correctedVBFJets[k]);
      }	  
    }
  }


  //VBF                                                                                                                                            
  if(finalVBFJets.size() >= 1){VBFJet1 = true;}
  if(finalVBFJets.size() == 2){VBFJet2 = true;}

  if(VBFJet1)
  {
    etaVBFJet1 = finalVBFJets[0].eta();
    pTVBFJet1 = finalVBFJets[0].pt();
    phiVBFJet1 = finalVBFJets[0].phi();
    massVBFJet1 = finalVBFJets[0].mass();
  }
  if(VBFJet2)
  {
    etaVBFJet2 = finalVBFJets[1].eta();
    pTVBFJet2 = finalVBFJets[1].pt();
    phiVBFJet2 = finalVBFJets[1].phi();
    massVBFJet2 = finalVBFJets[1].mass();
    VBFDiJetMass = (finalVBFJets[0].p4()+finalVBFJets[1].p4()).M();
    VBFDeltaEta = fabs(finalVBFJets[0].eta()-finalVBFJets[1].eta());
    // OLD MORIOND --- FisherDiscrim = 0.09407*fabs(VBFDeltaEta) + 4.1581e-4*VBFDiJetMass;
    FisherDiscrim = 0.18*fabs(VBFDeltaEta) + 1.92e-4*VBFDiJetMass;
  }
  
}




void UFHZZ4LAna::setGENVariables(std::vector<reco::GenParticle> Higgs, 
		     std::vector<reco::GenParticle> Zs, 
		     std::vector<reco::GenParticle> leptonsS1, std::vector<reco::GenParticle> leptonsS3)
{

  if( Higgs.size() == 1) GENMH = Higgs[0].mass();
  if( leptonsS1.size() == 4)
  {
    GENidL1_S1 = leptonsS1[0].pdgId(); GENidL2_S1 = leptonsS1[1].pdgId(); 
    GENidL3_S1 = leptonsS1[2].pdgId(); GENidL4_S1 = leptonsS1[3].pdgId();
    GENpTL1_S1 = leptonsS1[0].pt(); GENpTL2_S1 = leptonsS1[1].pt(); 
    GENpTL3_S1 = leptonsS1[2].pt(); GENpTL4_S1 = leptonsS1[3].pt();
    GENpXL1_S1 = leptonsS1[0].px(); GENpXL2_S1 = leptonsS1[1].px(); 
    GENpXL3_S1 = leptonsS1[2].px(); GENpXL4_S1 = leptonsS1[3].px();
    GENpYL1_S1 = leptonsS1[0].py(); GENpYL2_S1 = leptonsS1[1].py(); 
    GENpYL3_S1 = leptonsS1[2].py(); GENpYL4_S1 = leptonsS1[3].py();
    GENpZL1_S1 = leptonsS1[0].pz(); GENpZL2_S1 = leptonsS1[1].pz(); 
    GENpZL3_S1 = leptonsS1[2].pz(); GENpZL4_S1 = leptonsS1[3].pz();
    GENEL1_S1 = leptonsS1[0].energy(); GENEL2_S1 = leptonsS1[1].energy(); 
    GENEL3_S1 = leptonsS1[2].energy(); GENEL4_S1 = leptonsS1[3].energy();
    GENchargeL1_S1 = leptonsS1[0].charge(); GENchargeL2_S1 = leptonsS1[1].charge();
    GENchargeL3_S1 = leptonsS1[2].charge(); GENchargeL4_S1 = leptonsS1[3].charge();
    GENetaL1_S1 = leptonsS1[0].eta(); GENetaL2_S1 = leptonsS1[1].eta();
    GENetaL3_S1 = leptonsS1[2].eta(); GENetaL4_S1 = leptonsS1[3].eta();
    GENphiL1_S1 = leptonsS1[0].phi(); GENphiL2_S1 = leptonsS1[1].phi();
    GENphiL3_S1 = leptonsS1[2].phi(); GENphiL4_S1 = leptonsS1[3].phi();
    GENM4L = (leptonsS1[0].p4()+leptonsS1[1].p4()+leptonsS1[2].p4()+leptonsS1[3].p4()).M();
  }

  if( Zs.size() == 2)
  {
    GENEZ1 = Zs[0].energy(); GENEZ2 = Zs[1].energy();
    GENpTZ1 = Zs[0].pt();    GENpTZ2 = Zs[1].pt();
    GENpXZ1 = Zs[0].px();    GENpXZ2 = Zs[1].px();
    GENpYZ1 = Zs[0].py();    GENpYZ2 = Zs[1].py();
    GENpZZ1 = Zs[0].pz();    GENpZZ2 = Zs[1].pz();
    GENMZ1  = Zs[0].mass();  GENMZ2  = Zs[1].mass();
  }
  if( leptonsS3.size() == 4)
  {
    GENidL1_S3 = leptonsS3[0].pdgId(); GENidL2_S3 = leptonsS3[1].pdgId(); 
    GENidL3_S3 = leptonsS3[2].pdgId(); GENidL4_S3 = leptonsS3[3].pdgId();
    GENpTL1_S3 = leptonsS3[0].pt(); GENpTL2_S3 = leptonsS3[1].pt(); 
    GENpTL3_S3 = leptonsS3[2].pt(); GENpTL4_S3 = leptonsS3[3].pt();
    GENpXL1_S3 = leptonsS3[0].px(); GENpXL2_S3 = leptonsS3[1].px(); 
    GENpXL3_S3 = leptonsS3[2].px(); GENpXL4_S3 = leptonsS3[3].px();
    GENpYL1_S3 = leptonsS3[0].py(); GENpYL2_S3 = leptonsS3[1].py(); 
    GENpYL3_S3 = leptonsS3[2].py(); GENpYL4_S3 = leptonsS3[3].py();
    GENpZL1_S3 = leptonsS3[0].pz(); GENpZL2_S3 = leptonsS3[1].pz(); 
    GENpZL3_S3 = leptonsS3[2].pz(); GENpZL4_S3 = leptonsS3[3].pz();
    GENEL1_S3 = leptonsS3[0].energy(); GENEL2_S3 = leptonsS3[1].energy(); 
    GENEL3_S3 = leptonsS3[2].energy(); GENEL4_S3 = leptonsS3[3].energy();
    GENchargeL1_S3 = leptonsS3[0].charge(); GENchargeL2_S3 = leptonsS3[1].charge();
    GENchargeL3_S3 = leptonsS3[2].charge(); GENchargeL4_S3 = leptonsS3[3].charge();
    GENetaL1_S3 = leptonsS3[0].eta(); GENetaL2_S3 = leptonsS3[1].eta();
    GENetaL3_S3 = leptonsS3[2].eta(); GENetaL4_S3 = leptonsS3[3].eta();
    GENphiL1_S3 = leptonsS3[0].phi(); GENphiL2_S3 = leptonsS3[1].phi();
    GENphiL3_S3 = leptonsS3[2].phi(); GENphiL4_S3 = leptonsS3[3].phi();
  }
}


void UFHZZ4LAna::setGENMatchedVariables(std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons)
{
  TLorentzVector m1,m2,m3,m4;
  //Muons
  for(unsigned int i = 0; i < selectedMuons.size(); i++)
  {
    for(unsigned int j = 0; j < selectedMuons[i].genParticleRefs().size(); j++)
    {
      if( selectedMuons[i].genParticle(j)->status() != 1 ) continue;
      if( abs(selectedMuons[i].genParticle(j)->pdgId()) != 13 ) continue;
      if( !genAna.IsMotherZ(selectedMuons[i].genParticle(j)) ) continue;
      double genZpx, genZpy, genZpz, genZE;
      TLorentzVector m; int mPdgId=0;
      TLorentzVector *genZ1 = new TLorentzVector();
      genAna.getStatusThree(selectedMuons[i].genParticle(j),m,13,mPdgId);
      genAna.getMotherZ(selectedMuons[i].genParticle(j),genZpx,genZpy,genZpz,genZE);
      genZ1->SetPxPyPzE(genZpx,genZpy,genZpz,genZE);
      int charge;
      if(mPdgId > 0) charge = -1;
      else charge = 1;
      if(m.Pt()<=0.0) continue;
      if(i == 0)
      {
        if(!RecoTwoETwoMuEvent)
        {
          m1 = m;
          idL1_GENMatched = mPdgId; pTL1_GENMatched = m.Pt(); 
          pXL1_GENMatched = m.Px(); pYL1_GENMatched = m.Py(); pZL1_GENMatched = m.Pz(); EL1_GENMatched = m.Energy(); 
          chargeL1_GENMatched = charge;
          etaL1_GENMatched = m.Eta(); phiL1_GENMatched = m.Phi();
          if(genZ1 != NULL)
          {
            EZ1_GENMatched = genZ1->Energy();  pTZ1_GENMatched = genZ1->Pt(); pXZ1_GENMatched = genZ1->Px(); 
            pYZ1_GENMatched = genZ1->Py();   pZZ1_GENMatched = genZ1->Pz();
          } 
        }
        else
        {
          m3 = m;
          idL3_GENMatched = mPdgId; pTL3_GENMatched = m.Pt();
          pXL3_GENMatched = m.Px(); pYL3_GENMatched = m.Py(); pZL3_GENMatched = m.Pz(); EL3_GENMatched = m.Energy();
          chargeL3_GENMatched = charge; 
          etaL3_GENMatched = m.Eta(); phiL3_GENMatched = m.Phi();
          if(genZ1 != NULL)
          {
            EZ2_GENMatched = genZ1->Energy();  pTZ2_GENMatched = genZ1->Pt(); pXZ2_GENMatched = genZ1->Px();
            pYZ2_GENMatched = genZ1->Py();   pZZ2_GENMatched = genZ1->Pz();
          }
        }
      }
      if(i == 1)
      {
        if(!RecoTwoETwoMuEvent)
        {
          m2 = m;
          idL2_GENMatched = mPdgId; pTL2_GENMatched = m.Pt(); 
          pXL2_GENMatched = m.Px(); pYL2_GENMatched = m.Py(); pZL2_GENMatched = m.Pz(); EL2_GENMatched = m.Energy(); 
          chargeL2_GENMatched = charge; 
          etaL2_GENMatched = m.Eta(); phiL2_GENMatched = m.Phi();
        }
        else
        {
          m4 = m;
          idL4_GENMatched = mPdgId; pTL4_GENMatched = m.Pt();
          pXL4_GENMatched = m.Px(); pYL4_GENMatched = m.Py(); pZL4_GENMatched = m.Pz(); EL4_GENMatched = m.Energy();
          chargeL4_GENMatched = charge;
          etaL4_GENMatched = m.Eta(); phiL4_GENMatched = m.Phi();
        }
      }
      if(i == 2)
      {
        m3 = m;
        idL3_GENMatched = mPdgId; pTL3_GENMatched = m.Pt(); 
        pXL3_GENMatched = m.Px(); pYL3_GENMatched = m.Py(); pZL3_GENMatched = m.Pz(); EL3_GENMatched = m.Energy(); 
        chargeL3_GENMatched = charge;
        etaL3_GENMatched = m.Eta(); phiL3_GENMatched = m.Phi();
        if(genZ1 != NULL)
        {
          EZ2_GENMatched = genZ1->Energy();  pTZ2_GENMatched = genZ1->Pt(); pXZ2_GENMatched = genZ1->Px();
          pYZ2_GENMatched = genZ1->Py(); pZZ2_GENMatched = genZ1->Pz();
        }
      }
      if(i == 3)
      {
        m4 = m;
        idL4_GENMatched = mPdgId; pTL4_GENMatched = m.Pt(); 
        pXL4_GENMatched = m.Px(); pYL4_GENMatched = m.Py(); pZL4_GENMatched = m.Pz(); EL4_GENMatched = m.Energy(); 
        chargeL4_GENMatched = charge; 
        etaL4_GENMatched = m.Eta(); phiL4_GENMatched = m.Phi();
      }
    }
  }

  //Electrons
  for(unsigned int i = 0; i < selectedElectrons.size(); i++)
  {
    for(unsigned int j = 0; j < selectedElectrons[i].genParticleRefs().size(); j++)
    {
      if( selectedElectrons[i].genParticle(j)->status() != 1 ) continue;
      if( abs(selectedElectrons[i].genParticle(j)->pdgId()) != 11 ) continue;
      if( !genAna.IsMotherZ(selectedElectrons[i].genParticle(j)) ) continue;
      double genZpx, genZpy, genZpz, genZE;
      TLorentzVector m; int mPdgId=0; 
      TLorentzVector *genZ1 = new TLorentzVector();
      genAna.getStatusThree(selectedElectrons[i].genParticle(j),m,11,mPdgId);
      genAna.getMotherZ(selectedElectrons[i].genParticle(j),genZpx,genZpy,genZpz,genZE);
      genZ1->SetPxPyPzE(genZpx,genZpy,genZpz,genZE);
      int charge;
      if(mPdgId > 0) charge = -1;
      else charge = 1;
      if(m.Pt()<=0.0) continue;
      if(i == 0)
      {
        if(!RecoTwoMuTwoEEvent)
        {
          m1 = m;
          idL1_GENMatched = mPdgId; pTL1_GENMatched = m.Pt(); 
          pXL1_GENMatched = m.Px(); pYL1_GENMatched = m.Py(); pZL1_GENMatched = m.Pz(); EL1_GENMatched = m.Energy(); 
          chargeL1_GENMatched = charge; 
          etaL1_GENMatched = m.Eta(); phiL1_GENMatched = m.Phi();
          if(genZ1 != NULL)
          {
            EZ1_GENMatched = genZ1->Energy();  pTZ1_GENMatched = genZ1->Pt(); pXZ1_GENMatched = genZ1->Px(); 
            pYZ1_GENMatched = genZ1->Py(); pZZ1_GENMatched = genZ1->Pz();
          } 
        }
        else 
        {
          m3 = m;
          idL3_GENMatched = mPdgId; pTL3_GENMatched = m.Pt();
          pXL3_GENMatched = m.Px(); pYL3_GENMatched = m.Py(); pZL3_GENMatched = m.Pz(); EL3_GENMatched = m.Energy();
          chargeL3_GENMatched = charge; 
          etaL3_GENMatched = m.Eta(); phiL3_GENMatched = m.Phi();
          if(genZ1 != NULL)
          {
            EZ2_GENMatched = genZ1->Energy();  pTZ2_GENMatched = genZ1->Pt(); pXZ2_GENMatched = genZ1->Px();
            pYZ2_GENMatched = genZ1->Py(); pZZ2_GENMatched = genZ1->Pz();
          }
        }
      }
      if(i == 1)
      {
        if(!RecoTwoMuTwoEEvent)
        {
          m2 = m;
          idL2_GENMatched = mPdgId; pTL2_GENMatched = m.Pt(); 
          pXL2_GENMatched = m.Px(); pYL2_GENMatched = m.Py(); pZL2_GENMatched = m.Pz(); EL2_GENMatched = m.Energy(); 
          chargeL2_GENMatched = charge; 
          etaL2_GENMatched = m.Eta(); phiL2_GENMatched = m.Phi();
        }
        else
        {
          m4 = m;
          idL4_GENMatched = mPdgId; pTL4_GENMatched = m.Pt();
          pXL4_GENMatched = m.Px(); pYL4_GENMatched = m.Py(); pZL4_GENMatched = m.Pz(); EL4_GENMatched = m.Energy();
          chargeL4_GENMatched = charge; 
          etaL4_GENMatched = m.Eta(); phiL4_GENMatched = m.Phi();
        }
      }
      if(i == 2)
      {
        m3 = m;
        idL3_GENMatched = mPdgId; pTL3_GENMatched = m.Pt(); 
        pXL3_GENMatched = m.Px(); pYL3_GENMatched = m.Py(); pZL3_GENMatched = m.Pz(); EL3_GENMatched = m.Energy(); 
        chargeL3_GENMatched = charge; 
        etaL3_GENMatched = m.Eta(); phiL3_GENMatched = m.Phi();
        if(genZ1 != NULL)
        {
          EZ2_GENMatched = genZ1->Energy();  pTZ2_GENMatched = genZ1->Pt(); pXZ2_GENMatched = genZ1->Px();
          pYZ2_GENMatched = genZ1->Py(); pZZ2_GENMatched = genZ1->Pz();
        }
      }
      if(i == 3)
      {
        m4 = m;
        idL4_GENMatched = mPdgId; pTL4_GENMatched = m.Pt(); 
        pXL4_GENMatched = m.Px(); pYL4_GENMatched = m.Py(); pZL4_GENMatched = m.Pz(); EL4_GENMatched = m.Energy(); 
        chargeL4_GENMatched = charge;
        etaL4_GENMatched = m.Eta(); phiL4_GENMatched = m.Phi();
      }	  
    }
  }

  m4l_GENMatched = (m1+m2+m3+m4).M();

}	      





void UFHZZ4LAna::bookResolutionHistograms()
{
  //Resolution
  histContainer_["ptrelerrorForZ_mu_barrel"]=fs->make<TH1F>("ptrelerrorForZ_mu_barrel","muon pt relative error; #sigma_{p_{T}}^{RECO}/p_{T}; N Events", 
							    1000, -1, 1);
  histContainer_["ptrelerrorForZ_mu_endcap"]=fs->make<TH1F>("ptrelerrorForZ_mu_endcap","muon pt relative error; #sigma_{p_{T}}^{RECO}/p_{T}; N Events", 
							    1000, -1, 1);
  histContainer_["ptrelerrorForZ_mu"]=fs->make<TH1F>("ptrelerrorForZ_mu","muon pt relative error; #sigma_{p_{T}}^{RECO}/p_{T}; N Events", 1000, -1, 1);
  histContainer2D_["ptrelerrorForZ_mu_vs_eta"]=fs->make<TH2F>("ptrelerrorForZ_mu_vs_eta","muon pt relative error; #sigma_{p_{T}}^{RECO}/p_{T}; N Events", 
							      1000, -1, 1, 24, 0, 2.4);
}

void UFHZZ4LAna::fillResolutionHistograms(edm::Handle<edm::View<pat::Muon> > muons)
{
  //  For Muon Pt Resolution  start
  /********** M U O N  C U T S **********/
  double muPtCut = 5;
  double muEtaCut = 2.4;
  double muChi2Cut = 10.0;
  int muMuonHits = 0;
  int muNumMatches = 0;
  for(edm::View<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu)
  {
    if( mu->isGlobalMuon() != 1 || mu->isTrackerMuon() != 1) continue;
    if(mu->normChi2() >= muChi2Cut ||  mu->globalTrack()->hitPattern().numberOfValidMuonHits() <= muMuonHits)  continue;
    if( mu->numberOfMatches() <= muNumMatches ) continue;
    double iso =  (mu->trackIso() + mu->ecalIso() + mu->hcalIso() )/mu->pt() ;
    if ( iso > 0.15 ) continue;
    if( mu->pt() < muPtCut || abs(mu->eta()) > muEtaCut ) continue;
    for(edm::View<pat::Muon>::const_iterator mu2=mu+1; mu2!=muons->end(); ++mu2)
    {
      if( mu2->isGlobalMuon() != 1 || mu2->isTrackerMuon() != 1) continue;
      if(mu2->normChi2() >= muChi2Cut ||  mu2->globalTrack()->hitPattern().numberOfValidMuonHits() <= muMuonHits)  continue;
      if( mu2->numberOfMatches() <= muNumMatches ) continue;
      double iso2 =  (mu2->trackIso() + mu2->ecalIso() + mu2->hcalIso() )/mu2->pt() ;
      if ( iso2 > 0.15 ) continue;
      if( mu2->pt() < muPtCut || abs(mu2->eta()) > muEtaCut ) continue;
      double mass = (mu->p4()+mu2->p4()).M();
      if(mass<95 and mass>85) 
      {
        double pterr = mu->muonBestTrack()->ptError(); // miniAOD
        double recpt = mu->pt();
        double eta = mu->eta();
        histContainer_["ptrelerrorForZ_mu"]->Fill(pterr/recpt,eventWeight);
        histContainer2D_["ptrelerrorForZ_mu_vs_eta"]->Fill(pterr/recpt, fabs(eta));
        if(abs(eta)>1.2) histContainer_["ptrelerrorForZ_mu_endcap"]->Fill(pterr/recpt,eventWeight);
        else histContainer_["ptrelerrorForZ_mu_barrel"]->Fill(pterr/recpt,eventWeight);
      }
    }
  }
  //  For Muon Pt Resolution  end

}



////////////Muons/////////////////
bool UFHZZ4LAna::findZ(std::vector<pat::PFParticle> photons, std::vector<double> deltaRVec, 
                       pat::Muon &muon1, pat::Muon &muon2,int taken1, int &taken2, int &assocMuon, 
                       math::XYZTLorentzVector &ZVec, math::XYZTLorentzVector &photVec, bool &foundPhoton)
{

  using namespace std;
  using namespace pat;
  using namespace reco;

  double deltaR1, deltaR2, smallestDeltaR, totalSmallestDeltaR = 999;
  double photHighestPt = 0;
  math::XYZTLorentzVector mll, mllgam;
  bool foundPhot = false;
  double massDiffPhot = 999, massDiffNoPhot = 999;
  double coneSize = 0.4;
  double muon1Iso = 999, muon2Iso = 999;
  double assocMuonTmp = 999;
  double photIsoCut = 1.0;
  bool foundZ = false;
  double isoVetoMuons = 0.00;

  if( !photons.empty() && photons.size() > 0 && doFsrRecovery)
  {
    for(int i=0; i<(int)photons.size(); i++)
    {
      if( taken1 == i ) continue;
      //pt, eta checks
      if( photons[i].pt() < 4 ) continue;
      if( photons[i].eta() > 2.4 ) continue;
      if( (photons[i].userFloat("fsrPhotonPFIsoChHad03pt02")+photons[i].userFloat("fsrPhotonPFIsoNHad03")
          +photons[i].userFloat("fsrPhotonPFIsoPhoton03")
          +photons[i].userFloat("fsrPhotonPFIsoChHadPU03pt02"))/photons[i].pt() > photIsoCut) continue;
      //calc both deltaRs
      deltaR1 = deltaR(muon1.eta(),muon1.phi(),photons[i].eta(),photons[i].phi());
      deltaR2 = deltaR(muon2.eta(),muon2.phi(),photons[i].eta(),photons[i].phi());
      //associate with closest lepton
      if( deltaR1 < deltaR2 ){ assocMuonTmp = 1; smallestDeltaR = deltaR1;}
      else{ assocMuonTmp = 2; smallestDeltaR = deltaR2;}
      if( smallestDeltaR > 0.5 || smallestDeltaR > deltaRVec[i] ) continue;
      if( photons[i].pt() < photHighestPt ) continue;
      //calc P vectors
      mllgam = muon1.p4()+muon2.p4()+photons[i].p4();
      mll = muon1.p4()+muon2.p4();
      //check inv mass
      if( mllgam.M() < 4 || mllgam.M() > 100) continue;
      massDiffPhot = fabs(mllgam.M() - Zmass);
      massDiffNoPhot = fabs(mll.M() - Zmass);

      //if its smaller with phot, keep phot
      if( massDiffPhot < massDiffNoPhot )
      {
        //check iso cone
        if( deltaR1 < coneSize && deltaR1 > isoVetoMuons && deltaR1 > 1e-06){muon1Iso = helper.pfIsoFSR(muon1,muonRho,photons[i].pt());}
        else{muon1Iso = helper.pfIso(muon1,muonRho);}
        if( deltaR2 < coneSize && deltaR2 > isoVetoMuons && deltaR2 > 1e-06){muon2Iso = helper.pfIsoFSR(muon2,muonRho,photons[i].pt());}
        else{muon2Iso = helper.pfIso(muon2,muonRho);}

        if(muon1Iso < isoCut && muon2Iso < isoCut)
        {
          foundPhot = true;
          photHighestPt = photons[i].pt();
          photVec = photons[i].p4();
          ZVec = mllgam;
          assocMuon = assocMuonTmp;
          taken2 = i;
          foundZ = true;
        }
      }
    }
    if(!foundPhot)
    {
      bool useDR = false, usePT = false;
      for(unsigned int i = 0; i < photons.size(); i++)
      {
        if( taken1 == (int)i ) continue;
        //pt, eta checks
        if( photons[i].pt() < 2 ) continue;
        if( photons[i].eta() > 2.4 ) continue;
        if( photons[i].pt() < 4 ) useDR = true;
        if( photons[i].pt() > 4 ) usePT = true;
        if( usePT && photons[i].pt() < photHighestPt) continue;
        //calc both deltaRs
        deltaR1 = deltaR(muon1.eta(),muon1.phi(),photons[i].eta(),photons[i].phi());
        deltaR2 = deltaR(muon2.eta(),muon2.phi(),photons[i].eta(),photons[i].phi());
        //associate with closest lepton
        if( deltaR1 < deltaR2 ){ assocMuonTmp = 1; smallestDeltaR = deltaR1;}
        else{ assocMuonTmp = 2; smallestDeltaR = deltaR2;}
        if( smallestDeltaR > 0.07  || smallestDeltaR > deltaRVec[i] ) continue;
        if( smallestDeltaR > totalSmallestDeltaR && useDR ) continue;
        //calc P vectors
        mllgam = muon1.p4()+muon2.p4()+photons[i].p4();
        mll = muon1.p4()+muon2.p4();
        //check inv mass
        if( mllgam.M() < 4 || mllgam.M() > 100) continue;
        massDiffPhot = fabs(mllgam.M() - Zmass);
        massDiffNoPhot = fabs(mll.M() - Zmass);
   
        //if its smaller with phot, keep phot
        if( massDiffPhot < massDiffNoPhot )
        {
          if( deltaR1 < coneSize && deltaR1 > isoVetoMuons && deltaR1 > 1e-06 ){muon1Iso = helper.pfIsoFSR(muon1,muonRho,photons[i].pt());}
          else{muon1Iso = helper.pfIso(muon1,muonRho);}
          if( deltaR2 < coneSize && deltaR2 > isoVetoMuons && deltaR2 > 1e-06 ){muon2Iso = helper.pfIsoFSR(muon2,muonRho,photons[i].pt());}
          else{muon2Iso = helper.pfIso(muon2,muonRho);}
		
          if(muon1Iso < isoCut && muon2Iso < isoCut)
          {
            foundPhot = true;
            if(useDR) totalSmallestDeltaR = smallestDeltaR;
            if(usePT) photHighestPt = photons[i].pt();
            photVec = photons[i].p4();
            ZVec = mllgam;
            assocMuon = assocMuonTmp;
            taken2 = (int)i;
            foundZ = true;
          }
        }
        useDR = false;
        usePT = false;
      }
    }
    if(!foundPhot)
    {
      muon1Iso = helper.pfIso(muon1,muonRho);
      muon2Iso = helper.pfIso(muon2,muonRho);
      mll = muon1.p4()+muon2.p4();

      if( muon1Iso < isoCut && muon2Iso < isoCut )
      {
        ZVec = muon1.p4()+muon2.p4();
        taken2 = 999;
        assocMuon = 999;
        foundPhot = false;
        foundZ = true;
      }
    }
  }
  else{
    muon1Iso = helper.pfIso(muon1,muonRho);
    muon2Iso = helper.pfIso(muon2,muonRho);
    mll = muon1.p4()+muon2.p4();

    if( muon1Iso < isoCut && muon2Iso < isoCut )
    {
      ZVec = muon1.p4()+muon2.p4();
      taken2 = 999;
      assocMuon = 999;
      foundPhot = false;
      foundZ = true;
    }
  }	  
  
  foundPhoton = foundPhot;
  return foundZ;
} 


////////////Electrons/////////////////
bool UFHZZ4LAna::findZ(std::vector<pat::PFParticle> photons, std::vector<double> deltaRVec, 
                       pat::Electron &electron1, pat::Electron &electron2,int taken1, int &taken2, int &assocElec, 
                       math::XYZTLorentzVector &ZVec, math::XYZTLorentzVector &photVec, bool &foundPhoton) 
{

  using namespace std;
  using namespace pat;
  using namespace reco;

  double deltaR1, deltaR2, smallestDeltaR, totalSmallestDeltaR = 999;
  double photHighestPt = 0;
  math::XYZTLorentzVector mll, mllgam;
  bool foundPhot = false;
  double massDiffPhot = 999, massDiffNoPhot = 999;
  double coneSize = 0.4;
  double elec1Iso = 999, elec2Iso = 999;
  double assocElecTmp = 999;
  double photIsoCut = 1.0; 
  bool foundZ = false;
  foundPhoton = false;
  double isoVetoEE = 0.08;
  double isoVetoE1 = 0, isoVetoE2 = 0;
  if(electron1.superCluster()->eta() > 1.479) isoVetoE1 = isoVetoEE;
  if(electron2.superCluster()->eta() > 1.479) isoVetoE2 = isoVetoEE;

  if( !photons.empty() && photons.size() > 0 && doFsrRecovery)
  {
    for(int i=0; i<(int)photons.size(); i++)
    {
      if( taken1==i ) continue;
      //pt, eta checks
      if( photons[i].pt() < 4 ) continue;
      if( photons[i].eta() > 2.4 ) continue;
      if( (photons[i].userFloat("fsrPhotonPFIsoChHad03pt02")+photons[i].userFloat("fsrPhotonPFIsoNHad03")
          +photons[i].userFloat("fsrPhotonPFIsoPhoton03")
          +photons[i].userFloat("fsrPhotonPFIsoChHadPU03pt02"))/photons[i].pt() > photIsoCut) continue;
      //calc both deltaRs
      deltaR1 = deltaR(electron1.eta(),electron1.phi(),photons[i].eta(),photons[i].phi());
      deltaR2 = deltaR(electron2.eta(),electron2.phi(),photons[i].eta(),photons[i].phi());
      //associate with closest lepton
      if( deltaR1 < deltaR2 ){ assocElecTmp = 1; smallestDeltaR = deltaR1;}
      else{ assocElecTmp = 2; smallestDeltaR = deltaR2;}
      if( smallestDeltaR > 0.5  || smallestDeltaR > deltaRVec[i]) continue;
      if( photons[i].pt() < photHighestPt ) continue;
      //calc P vectors
      mllgam = electron1.p4()+electron2.p4()+photons[i].p4();
      mll = electron1.p4()+electron2.p4();
      //check inv mass
      if( mllgam.M() < 4 || mllgam.M() > 100) continue;
      massDiffPhot = fabs(mllgam.M() - Zmass);
      massDiffNoPhot = fabs(mll.M() - Zmass);

      //if its smaller with phot, keep phot
      if( massDiffPhot < massDiffNoPhot )
      {
        //check  iso cone
        if( deltaR1 < coneSize && deltaR1 > isoVetoE1){elec1Iso = helper.pfIsoFSR(electron1,elecRho,photons[i].pt());}
        else{elec1Iso = helper.pfIso(electron1,elecRho);}
        //check electron2's iso cone
        if( deltaR2 < coneSize  && deltaR2 > isoVetoE2){elec2Iso = helper.pfIsoFSR(electron2,elecRho,photons[i].pt());}
        else{elec2Iso = helper.pfIso(electron2,elecRho);}
	      
        if(elec1Iso < isoCut && elec2Iso < isoCut)
        {
          foundPhot = true;
          photHighestPt = photons[i].pt();
          photVec = photons[i].p4();
          ZVec = mllgam;
          assocElec = assocElecTmp;
          taken2 = (int)i;
          foundZ = true;
        }
      }
    }
    if(!foundPhot)
    {
      bool useDR = false, usePT = false;
      for(int i=0; i<(int)photons.size(); i++)
      {
        if( taken1==i ) continue;
        //pt, eta checks
        if( photons[i].pt() < 2 ) continue;
        if( photons[i].eta() > 2.4 ) continue;
        if( photons[i].pt() < 4 ) useDR = true;
        if( photons[i].pt() > 4 ) usePT = true;
        if( usePT && photons[i].pt() < photHighestPt) continue;
        //calc both deltaRs
        deltaR1 = deltaR(electron1.eta(),electron1.phi(),photons[i].eta(),photons[i].phi());
        deltaR2 = deltaR(electron2.eta(),electron2.phi(),photons[i].eta(),photons[i].phi());
        //associate with closest lepton
        if( deltaR1 < deltaR2 ){ assocElecTmp = 1; smallestDeltaR = deltaR1;}
        else{ assocElecTmp = 2; smallestDeltaR = deltaR2;}
        if( smallestDeltaR > 0.07  || smallestDeltaR > deltaRVec[i]) continue;
        if( smallestDeltaR > totalSmallestDeltaR ) continue;
        //calc P vectors
        mllgam = electron1.p4()+electron2.p4()+photons[i].p4();
        mll = electron1.p4()+electron2.p4();
        //check inv mass
        if( mllgam.M() < 4 || mllgam.M() > 100) continue;
        massDiffPhot = fabs(mllgam.M() - Zmass);
        massDiffNoPhot = fabs(mll.M() - Zmass);
        //if its smaller with phot, keep phot
        if( massDiffPhot < massDiffNoPhot )
        {
          //check electron1's iso cone
          if( deltaR1 < coneSize && deltaR1 > isoVetoE1 ){elec1Iso = helper.pfIsoFSR(electron1,elecRho,photons[i].pt());}
          else{elec1Iso = helper.pfIso(electron1,elecRho);}
          if( deltaR2 < coneSize && deltaR2 > isoVetoE2 ){elec2Iso = helper.pfIsoFSR(electron2,elecRho,photons[i].pt());}
          else{elec2Iso = helper.pfIso(electron2,elecRho);}
          if(elec1Iso < isoCut && elec2Iso < isoCut)
          {
            foundPhot = true;
            if(useDR) totalSmallestDeltaR = smallestDeltaR;
            if(usePT) photHighestPt = photons[i].pt();
            photHighestPt = photons[i].pt();
            photVec = photons[i].p4();
            ZVec = mllgam;
            assocElec = assocElecTmp;
            taken2 = (int)i;
            foundZ = true;
          }
        }
        useDR = false;
        usePT = false;
      }
    }
    if(!foundPhot)
    {
      elec1Iso = helper.pfIso(electron1,elecRho);
      elec2Iso = helper.pfIso(electron2,elecRho);
      mll = electron1.p4()+electron2.p4();
      if( elec1Iso < isoCut && elec2Iso < isoCut )
      {
        ZVec = electron1.p4()+electron2.p4();
        taken2 = 999;
        assocElec = 999;
        foundPhot = false;
        foundZ = true;
      }
    }
  }
  else
  {
    elec1Iso = helper.pfIso(electron1,elecRho);
    elec2Iso = helper.pfIso(electron2,elecRho);
    mll = electron1.p4()+electron2.p4();
    if( elec1Iso < isoCut && elec2Iso < isoCut )
    {
      ZVec = electron1.p4()+electron2.p4();
      taken2 = 999;
      assocElec = 999;
      foundPhot = false;
      foundZ = true;
    }
  }	  
  
  foundPhoton = foundPhot;
  return foundZ;
} 


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void UFHZZ4LAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(UFHZZ4LAna);
