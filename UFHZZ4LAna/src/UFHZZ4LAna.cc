// -*- C++ -*-
//
// Package:    UFHZZ4LAna
// Class:      UFHZZ4LAna
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
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <set>

#define PI 3.14159

// user include files 
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

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
//#include "ZZMatrixElement/MELA/interface/Mela.h"  // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/PseudoMELA.h" // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/SpinOneEvenMELA.h" // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/SpinOneOddMELA.h" // removed for miniAOD
//#include "ZZMatrixElement/MELA/interface/SpinTwoMinimalMELA.h" // removed for miniAOD
#include "ZZMatrixElement/MELA/src/computeAngles.h" // removed for miniAOD

//MEKD
//#include "ZZMatrixElement/MEKD/interface/MEKD.h" // removed for miniAOD
//#include "ZZMatrixElement/MEKD/interface/MEKD_MG.h" // removed for miniAOD
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h" // removed for miniAOD

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

//Mass Error
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMassErr.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPerLepResolution.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LDiLepResolution.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LResolution.h"
//GEN
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LGENAna.h"
//VBF Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJets.h"

// Jet energy correction
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <vector>

// Kinematic Fit
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

// EWK corrections
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/EwkCorrections.h"

//
// class declaration
//
using namespace MEMNames;
using namespace EwkCorrections;

class UFHZZ4LAna : public edm::EDAnalyzer {
public:
    explicit UFHZZ4LAna(const edm::ParameterSet&);
    ~UFHZZ4LAna();
  
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    static bool sortByPt( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.pt() > p2.pt()); };
  
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
  
    virtual void beginRun(edm::Run const&, const edm::EventSetup& iSetup);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup);
  
    void findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons, const edm::Event& iEvent);

    void bookResolutionHistograms();
    void fillResolutionHistograms(edm::Handle<edm::View<pat::Muon> > muons);
  
    //MELA
    HZZ4LAngles angles;
    //Helper Class
    HZZ4LHelper helper;
    //Sip Class
    HZZ4LSipAna sipAna;
    //Mass Err
    HZZ4LMassErr massErr;
    //GEN
    HZZ4LGENAna genAna;
    //VBF
    HZZ4LJets jetHelper;
    //PU Reweighting
    edm::LumiReWeighting *lumiWeight;
    HZZ4LPileUp pileUp;
    //JES Uncertainties
    JetCorrectionUncertainty *jecunc;

    // kfactors
    TSpline3 *kFactor_ggzz;
    std::vector<std::vector<float> > tableEwk;

    //Saved Events Trees
    TTree *passedEventsTree_All;

    void bookPassedEventTree(TString treeName, TTree *tree);
    void setTreeVariables( const edm::Event&, const edm::EventSetup&, 
                           std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, std::vector<pat::Jet> goodJets);
    void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                         edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                         edm::Handle<edm::View<reco::GenJet> > genJets);
    bool mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts);

    // -------------------------
    // RECO level information
    // -------------------------

    // Event Variables
    ULong64_t Run, Event, LumiSect;
    int nVtx;
    int finalState;
    std::string triggersPassed;
    bool passedTrig, passedFullSelection, passedZ4lSelection, passedQCDcut;

    // Event Weights
    float genWeight, pileupWeight, dataMCWeight, eventWeight;
    float k_ggZZ, k_qqZZ_qcd, k_qqZZ_ewk;
    // pdf weights                                                                   
    vector<float> pdfWeights;
    vector<int> pdfWeightIDs;
    float pdfRMSup, pdfRMSdown, pdfENVup, pdfENVdown;
    // lepton variables
    vector<double> lep_pt;
    vector<double> lep_pterr;
    vector<double> lep_eta;
    vector<double> lep_phi;
    vector<double> lep_mass;
    vector<double> lepFSR_pt;
    vector<double> lepFSR_eta;
    vector<double> lepFSR_phi;
    vector<double> lepFSR_mass;
    int lep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z2 sub
    float pTL1, pTL2, pTL3, pTL4;
    float etaL1, etaL2, etaL3, etaL4;
    int idL1, idL2, idL3, idL4;
    float pTL1FSR, pTL2FSR, pTL3FSR, pTL4FSR;
    vector<int> lep_genindex; //position of lepton in GENlep_p4 (if gen matched, -1 if not gen matched)
    vector<int> lep_id;
    vector<float> lep_mva;
    vector<int> lep_ecalDriven;
    vector<int> lep_tightId;
    vector<int> lep_tightIdSUS;
    vector<float> lep_Sip;
    vector<float> lep_IP;
    vector<float> lep_isoNH;
    vector<float> lep_isoCH;
    vector<float> lep_isoPhot;
    vector<float> lep_isoPU;
    vector<float> lep_isoPUcorr;
    vector<float> lep_RelIso;
    vector<float> lep_RelIsoNoFSR;
    vector<float> lep_MiniIso;
    vector<float> lep_ptRatio;
    vector<float> lep_ptRel;
    vector<int> lep_missingHits;
    vector<string> lep_filtersMatched; // for each lepton, all filters it is matched to
    int nisoleptons;
    double muRho, elRho, rhoSUS;

    // Higgs candidate variables
    vector<double> H_pt;
    vector<double> H_eta;
    vector<double> H_phi;
    vector<double> H_mass;
    vector<double> H_noFSR_pt;
    vector<double> H_noFSR_eta;
    vector<double> H_noFSR_phi;
    vector<double> H_noFSR_mass;
    float mass4l, mass4l_noFSR, mass4e, mass4mu, mass2e2mu, pT4l, eta4l, phi4l, rapidity4l;
    float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;

    // kin fit
    float mass4lREFIT, massZ1REFIT, massZ1Err, mass4lErr, mass4lErrREFIT;

    // Z candidate variables
    vector<double> Z_pt;
    vector<double> Z_eta;
    vector<double> Z_phi;
    vector<double> Z_mass;
    vector<double> Z_noFSR_pt;
    vector<double> Z_noFSR_eta;
    vector<double> Z_noFSR_phi;
    vector<double> Z_noFSR_mass;
    int Z_Hindex[2]; // position of Z1 and Z2 in Z_p4
    float massZ1, massZ2, pTZ1, pTZ2;

    // MET
    float met;
    float met_phi;

    // Jets
    vector<double> jet_pt;
    vector<double> jet_eta;
    vector<double> jet_phi;
    vector<double> jet_mass;
    vector<float> jet_pumva, jet_csvv2;
    vector<int>   jet_isbtag;
    vector<double> jet_jesup_pt;
    vector<double> jet_jesup_eta;
    vector<double> jet_jesup_phi;
    vector<double> jet_jesup_mass;
    vector<double> jet_jesdn_pt;
    vector<double> jet_jesdn_eta;
    vector<double> jet_jesdn_phi;
    vector<double> jet_jesdn_mass;
    vector<double> jet_jerup_pt;
    vector<double> jet_jerup_eta;
    vector<double> jet_jerup_phi;
    vector<double> jet_jerup_mass;
    vector<double> jet_jerdn_pt;
    vector<double> jet_jerdn_eta;
    vector<double> jet_jerdn_phi;
    vector<double> jet_jerdn_mass;    
    int njets_pt30_eta4p7;
    int njets_pt30_eta4p7_jesup;
    int njets_pt30_eta4p7_jesdn;
    int njets_pt30_eta4p7_jerup;
    int njets_pt30_eta4p7_jerdn;
    int nbjets_pt30_eta4p7;
    int nvjets_pt40_eta2p4;
    float pt_leadingjet_pt30_eta4p7;
    float pt_leadingjet_pt30_eta4p7_jesup;
    float pt_leadingjet_pt30_eta4p7_jesdn;
    float pt_leadingjet_pt30_eta4p7_jerup;
    float pt_leadingjet_pt30_eta4p7_jerdn;
    float absrapidity_leadingjet_pt30_eta4p7;
    float absrapidity_leadingjet_pt30_eta4p7_jesup;
    float absrapidity_leadingjet_pt30_eta4p7_jesdn;
    float absrapidity_leadingjet_pt30_eta4p7_jerup;
    float absrapidity_leadingjet_pt30_eta4p7_jerdn;
    float absdeltarapidity_hleadingjet_pt30_eta4p7;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jesup;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jerup;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn;
    float DijetMass, DijetDEta;
    float DijetFisher;

    // FSR Photons
    int nFSRPhotons;
    vector<int> fsrPhotons_lepindex; 
    vector<double> fsrPhotons_pt;
    vector<double> fsrPhotons_pterr;
    vector<double> fsrPhotons_eta;
    vector<double> fsrPhotons_phi;
    vector<double> fsrPhotons_mass;
    vector<float> fsrPhotons_dR;
    vector<float> fsrPhotons_iso;
    vector<float> allfsrPhotons_dR;
    vector<float> allfsrPhotons_pt;
    vector<float> allfsrPhotons_iso;

    // Z4l? FIXME
    float theta12, theta13, theta14;  
    float minM3l, Z4lmaxP, minDeltR, m3l_soft;
    float minMass2Lep, maxMass2Lep;
    float thetaPhoton, thetaPhotonZ;

    //Resolution
    float massErrorUCSD, massErrorUCSDCorr, massErrorUF, massErrorUFCorr, massErrorUFADCorr;
    HZZ4LPerLepResolution * PerLepReso;
    HZZ4LDiLepResolution * DiLepReso;
    HZZ4LResolution * FourLepReso;
    
    // Event Category
    int EventCat;

    // -------------------------
    // GEN level information
    // -------------------------

    //Event variables
    int GENfinalState;
    bool passedFiducialSelection;

    // lepton variables
    vector<double> GENlep_pt; 
    vector<double> GENlep_eta; 
    vector<double> GENlep_phi; 
    vector<double> GENlep_mass; 
    vector<int> GENlep_id; 
    vector<int> GENlep_status; 
    vector<int> GENlep_MomId; 
    vector<int> GENlep_MomMomId;
    int GENlep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
    vector<float> GENlep_isoCH; 
    vector<float> GENlep_isoNH; 
    vector<float> GENlep_isoPhot;
    vector<float> GENlep_RelIso; 


    // Higgs candidate variables (calculated using selected gen leptons)
    vector<double> GENH_pt; 
    vector<double> GENH_eta; 
    vector<double> GENH_phi; 
    vector<double> GENH_mass; 
    float GENmass4l, GENmass4e, GENmass4mu, GENmass2e2mu, GENpT4l, GENeta4l, GENrapidity4l;
    float GENMH; //mass directly from gen particle with id==25
    float GENcosTheta1, GENcosTheta2, GENcosThetaStar, GENPhi, GENPhi1;

    // Z candidate variables
    vector<double> GENZ_pt; 
    vector<double> GENZ_eta; 
    vector<double> GENZ_phi; 
    vector<double> GENZ_mass; 
    vector<int> GENZ_DaughtersId; vector<int> GENZ_MomId;
    float  GENmassZ1, GENmassZ2, GENpTZ1, GENpTZ2, GENdPhiZZ, GENmassZZ, GENpTZZ;

    // Higgs variables directly from GEN particle
    float GENHmass;

    // Jets
    vector<double> GENjet_pt; 
    vector<double> GENjet_eta; 
    vector<double> GENjet_phi; 
    vector<double> GENjet_mass; 
    int GENnjets_pt30_eta4p7;
    float GENpt_leadingjet_pt30_eta4p7;
    float GENabsrapidity_leadingjet_pt30_eta4p7;
    float GENabsdeltarapidity_hleadingjet_pt30_eta4p7;

    //KinZfitter
    KinZfitter *kinZfitter;

    // MEM
    MEMs*  combinedMEM;

    double me_0plus_JHU, me_qqZZ_MCFM, p0plus_m4l, bkg_m4l;
    double D_bkg_kin, D_bkg;   

    double bkg_VAMCFM;
    double p0minus_VAJHU, Dgg10_VAMCFM;
    double phjj_VAJHU, pvbf_VAJHU;

    double D_g4;
    double Djet_VAJHU;

    // mass err
    std::map<TString,TH1F*> hContainer_;
    std::map<TString,TH2F*> hContainer2D_;
    std::map<TString,TH3F*> hContainer3D_;

    // a vector<float> for each vector<double>
    vector<float>  lep_pt_float;
    vector<float>  lep_pterr_float;
    vector<float>  lep_eta_float;
    vector<float>  lep_phi_float;
    vector<float>  lep_mass_float;
    vector<float>  lepFSR_pt_float;
    vector<float>  lepFSR_eta_float;
    vector<float>  lepFSR_phi_float;
    vector<float>  lepFSR_mass_float;
    vector<float>  H_pt_float;
    vector<float>  H_eta_float;
    vector<float>  H_phi_float;
    vector<float>  H_mass_float;
    vector<float>  H_noFSR_pt_float;
    vector<float>  H_noFSR_eta_float;
    vector<float>  H_noFSR_phi_float;
    vector<float>  H_noFSR_mass_float;
    vector<float>  Z_pt_float;
    vector<float>  Z_eta_float;
    vector<float>  Z_phi_float;
    vector<float>  Z_mass_float;
    vector<float>  Z_noFSR_pt_float;
    vector<float>  Z_noFSR_eta_float;
    vector<float>  Z_noFSR_phi_float;
    vector<float>  Z_noFSR_mass_float;
    vector<float>  jet_pt_float;
    vector<float>  jet_eta_float;
    vector<float>  jet_phi_float;
    vector<float>  jet_mass_float;
    vector<float>  jet_jesup_pt_float;
    vector<float>  jet_jesup_eta_float;
    vector<float>  jet_jesup_phi_float;
    vector<float>  jet_jesup_mass_float;
    vector<float>  jet_jesdn_pt_float;
    vector<float>  jet_jesdn_eta_float;
    vector<float>  jet_jesdn_phi_float;
    vector<float>  jet_jesdn_mass_float;
    vector<float>  jet_jerup_pt_float;
    vector<float>  jet_jerup_eta_float;
    vector<float>  jet_jerup_phi_float;
    vector<float>  jet_jerup_mass_float;
    vector<float>  jet_jerdn_pt_float;
    vector<float>  jet_jerdn_eta_float;
    vector<float>  jet_jerdn_phi_float;
    vector<float>  jet_jerdn_mass_float;
    vector<float>  fsrPhotons_pt_float;
    vector<float>  fsrPhotons_pterr_float;
    vector<float>  fsrPhotons_eta_float;
    vector<float>  fsrPhotons_phi_float;
    vector<float>  fsrPhotons_mass_float;
    vector<float>  GENlep_pt_float;
    vector<float>  GENlep_eta_float;
    vector<float>  GENlep_phi_float;
    vector<float>  GENlep_mass_float;
    vector<float>  GENH_pt_float;
    vector<float>  GENH_eta_float;
    vector<float>  GENH_phi_float;
    vector<float>  GENH_mass_float;
    vector<float>  GENZ_pt_float;
    vector<float>  GENZ_eta_float;
    vector<float>  GENZ_phi_float;
    vector<float>  GENZ_mass_float;
    vector<float>  GENjet_pt_float;
    vector<float>  GENjet_eta_float;
    vector<float>  GENjet_phi_float;
    vector<float>  GENjet_mass_float;

    // Global Variables but not stored in the tree
    vector<double> lep_ptreco;
    vector<int> lep_ptid;
    vector<int> lep_ptindex;
    vector<pat::Muon> recoMuons;
    vector<pat::Electron> recoElectrons;
    vector<pat::PFParticle> fsrPhotons; 
    TLorentzVector HVec, HVecNoFSR, Z1Vec, Z2Vec;
    TLorentzVector GENZ1Vec, GENZ2Vec;
    bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;
    bool foundHiggsCandidate;
    bool firstEntry;

    // hist container
    std::map<std::string,TH1F*> histContainer_;
    std::map<std::string,TH2F*> histContainer2D_;
 
    //Input edm
    edm::EDGetTokenT<edm::View<pat::Photon> > photonSrc_;
    edm::EDGetTokenT<edm::View<pat::Electron> > elecSrc_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muonSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
    edm::EDGetTokenT<edm::View<pat::MET> > metSrc_;
    //edm::InputTag triggerSrc_;
    edm::EDGetTokenT<edm::TriggerResults> triggerSrc_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<std::vector<reco::Conversion> > conversionSrc_;
    edm::EDGetTokenT<double> muRhoSrc_;
    edm::EDGetTokenT<double> elRhoSrc_;
    edm::EDGetTokenT<double> rhoSrcSUS_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSrc_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsSrc_;
    edm::EDGetTokenT<edm::View<pat::PFParticle> > fsrPhotonsSrc_;
    edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsSrc_;
    edm::EDGetTokenT<GenEventInfoProduct> generatorSrc_;
    edm::EDGetTokenT<LHEEventProduct> lheInfoSrc_;
    edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;

    // Configuration
    const float Zmass;
    float mZ1Low, mZ2Low;
    float mZ1High, mZ2High;
    float m4lLowCut;
    float jetpt_cut, jeteta_cut;
    std::string elecID;
    bool isMC, isSignal;
    float mH;
    float crossSection, filterEff, Lumi;
    bool weightEvents;
    float isoCutEl, isoCutMu; 
    double isoConeSizeEl, isoConeSizeMu;
    float sip3dCut;
    float leadingPtCut, subleadingPtCut;
    float genIsoCutEl, genIsoCutMu;
    double genIsoConeSizeEl, genIsoConeSizeMu;
    float _elecPtCut, _muPtCut;
    float BTagCut;
    bool reweightForPU;
    bool interactiveRun;
    std::string PUVersion;
    bool doFsrRecovery, doPUJetID;
    int jetIDLevel;
    bool doSUSYSelection;
    bool bStudyResolution;
    bool bStudyDiLeptonResolution;
    bool bStudyFourLeptonResolution;
    std::vector<std::string> triggerList;
    int skimLooseLeptons;
    bool verbose;
    // register to the TFileService
    edm::Service<TFileService> fs;

    // Counters
    float nEventsTotal;
    float sumWeightsTotal;


};


UFHZZ4LAna::UFHZZ4LAna(const edm::ParameterSet& iConfig) :
    histContainer_(),
    histContainer2D_(),
    photonSrc_(consumes<edm::View<pat::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc"))),
    elecSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
    muonSrc_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
    jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc"))),
    metSrc_(consumes<edm::View<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>("metSrc"))),
    triggerSrc_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    vertexSrc_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
    beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
    conversionSrc_(consumes<std::vector<reco::Conversion> >(iConfig.getUntrackedParameter<edm::InputTag>("conversionSrc"))),
    muRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc"))),
    elRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc"))),
    rhoSrcSUS_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("rhoSrcSUS"))),
    pileupSrc_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc"))),
    pfCandsSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandsSrc"))),
    fsrPhotonsSrc_(consumes<edm::View<pat::PFParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("fsrPhotonsSrc"))),
    prunedgenParticlesSrc_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc"))),
    packedgenParticlesSrc_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("packedgenParticlesSrc"))),
    genJetsSrc_(consumes<edm::View<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genJetsSrc"))),
    generatorSrc_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorSrc"))),
    lheInfoSrc_(consumes<LHEEventProduct>(iConfig.getUntrackedParameter<edm::InputTag>("lheInfoSrc"))),
    lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""))),
    Zmass(91.1876),
    mZ1Low(iConfig.getUntrackedParameter<double>("mZ1Low",40.0)),
    mZ2Low(iConfig.getUntrackedParameter<double>("mZ2Low",12.0)), // was 12
    mZ1High(iConfig.getUntrackedParameter<double>("mZ1High",120.0)),
    mZ2High(iConfig.getUntrackedParameter<double>("mZ2High",120.0)),
    m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",70.0)),
    jetpt_cut(iConfig.getUntrackedParameter<double>("jetpt_cut",10.0)),
    jeteta_cut(iConfig.getUntrackedParameter<double>("eta_cut",4.7)),
    elecID(iConfig.getUntrackedParameter<std::string>("elecID","NonTrig")),
    isMC(iConfig.getUntrackedParameter<bool>("isMC",true)),
    isSignal(iConfig.getUntrackedParameter<bool>("isSignal",false)),
    mH(iConfig.getUntrackedParameter<double>("mH",0.0)),
    crossSection(iConfig.getUntrackedParameter<double>("CrossSection",1.0)),
    filterEff(iConfig.getUntrackedParameter<double>("FilterEff",1.0)),
    Lumi(iConfig.getUntrackedParameter<double>("Lumi",1.0)),
    weightEvents(iConfig.getUntrackedParameter<bool>("weightEvents",false)),
    isoCutEl(iConfig.getUntrackedParameter<double>("isoCutEl",0.35)),
    isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",0.35)),
    isoConeSizeEl(iConfig.getUntrackedParameter<double>("isoConeSizeEl",0.3)),
    isoConeSizeMu(iConfig.getUntrackedParameter<double>("isoConeSizeMu",0.3)),
    sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
    leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
    subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
    genIsoCutEl(iConfig.getUntrackedParameter<double>("genIsoCutEl",0.35)), 
    genIsoCutMu(iConfig.getUntrackedParameter<double>("genIsoCutMu",0.35)), 
    genIsoConeSizeEl(iConfig.getUntrackedParameter<double>("genIsoConeSizeEl",0.3)), 
    genIsoConeSizeMu(iConfig.getUntrackedParameter<double>("genIsoConeSizeMu",0.3)), 
    _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",7)),
    _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",5)),
    BTagCut(iConfig.getUntrackedParameter<double>("BTagCut",0.89)),
    reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
    interactiveRun(iConfig.getUntrackedParameter<bool>("interactiveRun",false)),
    PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Fall15_76X")),
    doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",true)),
    doPUJetID(iConfig.getUntrackedParameter<bool>("doPUJetID",false)),
    jetIDLevel(iConfig.getUntrackedParameter<int>("jetIDLevel",1)),
    doSUSYSelection(iConfig.getUntrackedParameter<bool>("doSUSYSelection",false)),
    bStudyResolution(iConfig.getUntrackedParameter<bool>("bStudyResolution",false)),
    bStudyDiLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyDiLeptonResolution",false)),
    bStudyFourLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyFourLeptonResolution",false)),
    triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList")),
    skimLooseLeptons(iConfig.getUntrackedParameter<int>("skimLooseLeptons",2)),    
    verbose(iConfig.getUntrackedParameter<bool>("verbose",false))

{
  
    if(!isMC){reweightForPU = false;}

    nEventsTotal=0.0;
    sumWeightsTotal=0.0;
    histContainer_["NEVENTS"]=fs->make<TH1F>("nEvents","nEvents in Sample",2,0,2);
    histContainer_["SUMWEIGHTS"]=fs->make<TH1F>("sumWeights","sum Weights of Sample",2,0,2);
    histContainer_["NVTX"]=fs->make<TH1F>("nVtx","Number of Vertices",36,-0.5,35.5);
    histContainer_["NVTX_RW"]=fs->make<TH1F>("nVtx_ReWeighted","Number of Vertices",36,-0.5,35.5);
    histContainer_["NINTERACT"]=fs->make<TH1F>("nInteractions","Number of True Interactions",61,-0.5,60.5);
    histContainer_["NINTERACT_RW"]=fs->make<TH1F>("nInteraction_ReWeighted","Number of True Interactions",61,-0.5,60.5);

    passedEventsTree_All = new TTree("passedEvents","passedEvents");

    if(bStudyResolution) {
        PerLepReso = new HZZ4LPerLepResolution();
        if(bStudyDiLeptonResolution) {  DiLepReso = new HZZ4LDiLepResolution(); }
        if(bStudyFourLeptonResolution)  { FourLepReso = new HZZ4LResolution(); }
    }

    edm::FileInPath jecfileInPath("UFHZZAnalysisRun2/UFHZZ4LAna/data/Summer15_25nsV3M3_DATA_UncertaintySources_AK4PFchs.txt");
    jecunc = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecfileInPath.fullPath().c_str(),"Total")));

    edm::FileInPath kfacfileInPath("UFHZZAnalysisRun2/UFHZZ4LAna/data/Kfactor_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
    TFile *fKFactor = TFile::Open(kfacfileInPath.fullPath().c_str());
    kFactor_ggzz = (TSpline3*) fKFactor->Get("sp_Kfactor");
    fKFactor->Close();
    delete fKFactor;

    tableEwk = readFile_and_loadEwkTable("ZZBG");   
   
    kinZfitter = new KinZfitter(!isMC);
    combinedMEM = new MEMs(13.0,125,"CTEQ6L",false);
 
    
}



UFHZZ4LAna::~UFHZZ4LAna()
{
    //destructor --- don't do anything here  
}



// ------------ method called for each event  ------------
void
UFHZZ4LAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;
    using namespace std;
    using namespace pat;
    using namespace trigger;
    using namespace MEMNames;
    using namespace EwkCorrections;

    nEventsTotal += 1.0;

    Run = iEvent.id().run();
    Event = iEvent.id().event();
    LumiSect = iEvent.id().luminosityBlock();

    if (verbose) {
       cout<<"Run: " << Run << ",Event: " << Event << ",LumiSect: "<<LumiSect<<endl;
    }

    // ======= Get Collections ======= //
    if (verbose) {cout<<"getting collections"<<endl;}

    // trigger collection
    edm::Handle<edm::TriggerResults> trigger;
    iEvent.getByToken(triggerSrc_,trigger);
    const edm::TriggerNames trigNames = iEvent.triggerNames(*trigger);

    // trigger Objects
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);

    // vertex collection
    edm::Handle<reco::VertexCollection> vertex;
    iEvent.getByToken(vertexSrc_,vertex);

    // photon collection 
    edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByToken(photonSrc_,photons);
  
    // electron collection
    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByToken(elecSrc_,electrons);
    if (verbose) cout<<electrons->size()<<" total electrons in the collection"<<endl;

    // muon collection
    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByToken(muonSrc_,muons);
    if (verbose) cout<<muons->size()<<" total muons in the collection"<<endl;

    // met collection 
    edm::Handle<edm::View<pat::MET> > mets;
    iEvent.getByToken(metSrc_,mets);
    
    // Rho Correction
    edm::Handle<double> eventRhoMu;
    iEvent.getByToken(muRhoSrc_,eventRhoMu);
    muRho = *eventRhoMu;

    edm::Handle<double> eventRhoE;
    iEvent.getByToken(elRhoSrc_,eventRhoE);
    elRho = *eventRhoE;

    edm::Handle<double> eventRhoSUS;
    iEvent.getByToken(rhoSrcSUS_,eventRhoSUS);
    rhoSUS = *eventRhoSUS;

    // Conversions
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    iEvent.getByToken(conversionSrc_, theConversions);
 
    // Beam Spot
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpotSrc_,beamSpot);
    const reco::BeamSpot BS = *beamSpot;

    // Particle Flow Cands
    edm::Handle<pat::PackedCandidateCollection> pfCands;
    iEvent.getByToken(pfCandsSrc_,pfCands);

    // Photons
    edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;
    iEvent.getByToken(fsrPhotonsSrc_,photonsForFsr);
  
    // Jets
    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByToken(jetSrc_,jets);

    // GEN collections
    edm::Handle<reco::GenParticleCollection> prunedgenParticles;
    iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);

    edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles;
    iEvent.getByToken(packedgenParticlesSrc_, packedgenParticles);
    
    edm::Handle<edm::View<reco::GenJet> > genJets;
    iEvent.getByToken(genJetsSrc_, genJets);
    
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(generatorSrc_,genEventInfo);
    
    //vector<edm::Handle<LHEEventProduct> > lheInfos;
    //iEvent.getManyByType(lheInfos); // using this method because the label is not always the same (e.g. "source" in the ttH sample)

    edm::Handle<LHEEventProduct> lheInfo;
    iEvent.getByToken(lheInfoSrc_, lheInfo);

    // ============ Initialize Variables ============= //

    // Event Variables
    if (verbose) {cout<<"clear variables"<<endl;}
    nVtx = -1.0;
    finalState = -1;
    triggersPassed="";
    passedTrig=false; passedFullSelection=false; passedZ4lSelection=false; passedQCDcut=false;

    // Event Weights
    genWeight=1.0; pileupWeight=1.0; dataMCWeight=1.0; eventWeight=1.0;
    k_ggZZ=1.0; k_qqZZ_qcd = 1.1; /*initialized to inclusive k-factor*/ k_qqZZ_ewk = 1.0;

    pdfWeights.clear();
    pdfRMSup=1.0; pdfRMSdown=1.0; pdfENVup=1.0; pdfENVdown=1.0;

    //lepton variables
    lep_pt.clear(); lep_pterr.clear(); lep_eta.clear(); lep_phi.clear(); lep_mass.clear(); 
    lepFSR_pt.clear(); lepFSR_eta.clear(); lepFSR_phi.clear(); lepFSR_mass.clear(); 
    for (int i=0; i<4; ++i) {lep_Hindex[i]=-1;}
    pTL1=-1.0; pTL2=-1.0; pTL3=-1.0; pTL4=-1.0;
    etaL1=9999.0; etaL2=9999.0; etaL3=9999.0; etaL4=9999.0;
    idL1=9999; idL2=9999; idL3=9999; idL4=9999;
    pTL1FSR=-1.0; pTL2FSR=-1.0; pTL3FSR=-1.0; pTL4FSR=-1.0;
    lep_genindex.clear(); lep_id.clear(); 
    lep_mva.clear(); lep_ecalDriven.clear(); lep_tightId.clear(); lep_tightIdSUS.clear();
    lep_Sip.clear(); lep_IP.clear(); 
    lep_isoNH.clear(); lep_isoCH.clear(); lep_isoPhot.clear(); lep_isoPU.clear(); lep_isoPUcorr.clear(); 
    lep_RelIso.clear(); lep_RelIsoNoFSR.clear(); lep_MiniIso.clear();
    lep_ptRatio.clear(); lep_ptRel.clear();
    lep_missingHits.clear();
    lep_filtersMatched.clear();    
    nisoleptons=0;

    // Higgs candidate variables
    H_pt.clear(); H_eta.clear(); H_phi.clear(); H_mass.clear(); 
    H_noFSR_pt.clear(); H_noFSR_eta.clear(); H_noFSR_phi.clear(); H_noFSR_mass.clear(); 
    mass4l=-1.0; mass4l_noFSR=-1.0; mass4e=-1.0; mass4mu=-1.0; mass2e2mu=-1.0; pT4l=-1.0; eta4l=9999.0; phi4l=9999.0; rapidity4l=9999.0;
    cosTheta1=9999.0; cosTheta2=9999.0; cosThetaStar=9999.0; Phi=9999.0; Phi1=9999.0;

    // kin fitter
    mass4lREFIT = -999.0; massZ1REFIT = -999.0; massZ1Err = -999.0; mass4lErr = -999.0; mass4lErrREFIT = -999.0;

    // Z candidate variables
    Z_pt.clear(); Z_eta.clear(); Z_phi.clear(); Z_mass.clear(); 
    Z_noFSR_pt.clear(); Z_noFSR_eta.clear(); Z_noFSR_phi.clear(); Z_noFSR_mass.clear(); 
    for (int i=0; i<2; ++i) {Z_Hindex[i]=-1;}
    massZ1=-1.0; massZ2=-1.0; pTZ1=-1.0; pTZ2=-1.0;

    // MET
    met=-1.0;
    met_phi=9999.0;

    // Jets
    jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_mass.clear(); 
    jet_jesup_pt.clear(); jet_jesup_eta.clear(); jet_jesup_phi.clear(); jet_jesup_mass.clear(); 
    jet_jesdn_pt.clear(); jet_jesdn_eta.clear(); jet_jesdn_phi.clear(); jet_jesdn_mass.clear(); 
    jet_jerup_pt.clear(); jet_jerup_eta.clear(); jet_jerup_phi.clear(); jet_jerup_mass.clear(); 
    jet_jerdn_pt.clear(); jet_jerdn_eta.clear(); jet_jerdn_phi.clear(); jet_jerdn_mass.clear(); 
    jet_pumva.clear(); jet_csvv2.clear(); jet_isbtag.clear();

    njets_pt30_eta4p7=0;
    njets_pt30_eta4p7_jesup=0;
    njets_pt30_eta4p7_jesdn=0;
    njets_pt30_eta4p7_jerup=0;
    njets_pt30_eta4p7_jerdn=0;

    nbjets_pt30_eta4p7=0; nvjets_pt40_eta2p4=0;

    pt_leadingjet_pt30_eta4p7=-1.0;
    pt_leadingjet_pt30_eta4p7_jesup=-1.0;
    pt_leadingjet_pt30_eta4p7_jesdn=-1.0;
    pt_leadingjet_pt30_eta4p7_jerup=-1.0;
    pt_leadingjet_pt30_eta4p7_jerdn=-1.0;

    absrapidity_leadingjet_pt30_eta4p7=-1.0;
    absrapidity_leadingjet_pt30_eta4p7_jesup=-1.0;
    absrapidity_leadingjet_pt30_eta4p7_jesdn=-1.0;
    absrapidity_leadingjet_pt30_eta4p7_jerup=-1.0;
    absrapidity_leadingjet_pt30_eta4p7_jerdn=-1.0;

    absdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
    absdeltarapidity_hleadingjet_pt30_eta4p7_jesup=-1.0;
    absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn=-1.0;
    absdeltarapidity_hleadingjet_pt30_eta4p7_jerup=-1.0;
    absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn=-1.0;

    DijetMass=-1.0; DijetDEta=9999.0;
    DijetFisher=9999.0;

    // FSR Photons
    nFSRPhotons=0;
    fsrPhotons_lepindex.clear(); fsrPhotons_pt.clear(); fsrPhotons_pterr.clear(); fsrPhotons_eta.clear(); fsrPhotons_phi.clear();
    fsrPhotons_dR.clear(); fsrPhotons_iso.clear();
    allfsrPhotons_dR.clear(); allfsrPhotons_pt.clear(); allfsrPhotons_iso.clear();

    // Z4l? FIXME
    theta12=9999.0; theta13=9999.0; theta14=9999.0;
    minM3l=-1.0; Z4lmaxP=-1.0; minDeltR=9999.0; m3l_soft=-1.0;
    minMass2Lep=-1.0; maxMass2Lep=-1.0;
    thetaPhoton=9999.0; thetaPhotonZ=9999.0;

    // -------------------------
    // GEN level information
    // ------------------------- 

    //Event variables
    GENfinalState=-1;
    passedFiducialSelection=false;

    // lepton variables
    GENlep_pt.clear(); GENlep_eta.clear(); GENlep_phi.clear(); GENlep_mass.clear();
    GENlep_id.clear();
    GENlep_status.clear();
    GENlep_MomId.clear();
    GENlep_MomMomId.clear();
    for (int i=0; i<4; ++i) {GENlep_Hindex[i]=-1;};//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
    GENlep_isoCH.clear();
    GENlep_isoNH.clear();
    GENlep_isoPhot.clear();
    GENlep_RelIso.clear();

    // Higgs candidate variables (calculated using selected gen leptons)
    GENH_pt.clear(); GENH_eta.clear(); GENH_phi.clear(); GENH_mass.clear();
    GENmass4l=-1.0; GENmassZ1=-1.0; GENmassZ2=-1.0; GENpT4l=-1.0; GENeta4l=9999.0; GENrapidity4l=9999.0; GENMH=-1.0;
    GENcosTheta1=9999.0; GENcosTheta2=9999.0; GENcosThetaStar=9999.0; GENPhi=9999.0; GENPhi1=9999.0;

    // Z candidate variables
    GENZ_DaughtersId.clear(); GENZ_MomId.clear();
    GENZ_pt.clear(); GENZ_eta.clear(); GENZ_phi.clear(); GENZ_mass.clear();
    GENmassZ1=-1.0; GENmassZ2=-1.0; GENpTZ1=-1.0; GENpTZ2=-1.0, GENdPhiZZ=9999.0, GENmassZZ=-1.0, GENpTZZ=-1.0;

    // Higgs variables directly from GEN particle
    GENHmass=-1.0;

    // Jets
    GENjet_pt.clear(); GENjet_eta.clear(); GENjet_phi.clear(); GENjet_mass.clear(); 
    GENnjets_pt30_eta4p7=0;
    GENpt_leadingjet_pt30_eta4p7=-1.0;
    GENabsrapidity_leadingjet_pt30_eta4p7=-1.0;
    GENabsdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;

    if (verbose) {cout<<"clear other variables"<<endl; }
    // Resolution
    massErrorUCSD=-1.0; massErrorUCSDCorr=-1.0; massErrorUF=-1.0; massErrorUFCorr=-1.0; massErrorUFADCorr=-1.0;

    // Event Category
    EventCat=-1;

    // Global variables not stored in tree
    lep_ptreco.clear(); lep_ptid.clear(); lep_ptindex.clear();
    recoMuons.clear(); recoElectrons.clear(); fsrPhotons.clear();
    HVec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    HVecNoFSR.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    Z1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    Z2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    GENZ1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    GENZ2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    RecoFourMuEvent = false; RecoFourEEvent = false;
    RecoTwoETwoMuEvent = false; RecoTwoMuTwoEEvent = false;
    foundHiggsCandidate = false;

    // Float vectors
    lep_pt_float.clear(); lep_pterr_float.clear(); lep_eta_float.clear(); lep_phi_float.clear(); lep_mass_float.clear();
    lepFSR_pt_float.clear(); lepFSR_eta_float.clear(); lepFSR_phi_float.clear(); lepFSR_mass_float.clear();
    H_pt_float.clear(); H_eta_float.clear(); H_phi_float.clear(); H_mass_float.clear();
    H_noFSR_pt_float.clear(); H_noFSR_eta_float.clear(); H_noFSR_phi_float.clear(); H_noFSR_mass_float.clear();
    Z_pt_float.clear(); Z_eta_float.clear(); Z_phi_float.clear(); Z_mass_float.clear();
    Z_noFSR_pt_float.clear(); Z_noFSR_eta_float.clear(); Z_noFSR_phi_float.clear(); Z_noFSR_mass_float.clear();
    jet_pt_float.clear(); jet_eta_float.clear(); jet_phi_float.clear(); jet_mass_float.clear();
    jet_jesup_pt_float.clear(); jet_jesup_eta_float.clear(); jet_jesup_phi_float.clear(); jet_jesup_mass_float.clear();
    jet_jesdn_pt_float.clear(); jet_jesdn_eta_float.clear(); jet_jesdn_phi_float.clear(); jet_jesdn_mass_float.clear();
    jet_jerup_pt_float.clear(); jet_jerup_eta_float.clear(); jet_jerup_phi_float.clear(); jet_jerup_mass_float.clear();
    jet_jerdn_pt_float.clear(); jet_jerdn_eta_float.clear(); jet_jerdn_phi_float.clear();  jet_jerdn_mass_float.clear();
    fsrPhotons_pt_float.clear(); fsrPhotons_pterr_float.clear(); fsrPhotons_eta_float.clear(); fsrPhotons_phi_float.clear(); fsrPhotons_mass_float.clear();

    // ====================== Do Analysis ======================== //

    std::map<int, TLorentzVector> fsrmap;
    vector<reco::Candidate*> selectedLeptons;
    std::map<unsigned int, TLorentzVector> selectedFsrMap;

    fsrmap.clear();
    selectedFsrMap.clear();
    selectedLeptons.clear();

    if (verbose) cout<<"start pileup reweighting"<<endl;
    // PU information
    if(isMC && reweightForPU) {        
        edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByToken(pileupSrc_, PupInfo);      

        if (verbose) cout<<"got pileup info"<<endl;

        std::vector<PileupSummaryInfo>::const_iterator PVI;      
        int npv = -1;
        for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            int BX = PVI->getBunchCrossing();
            if(BX == 0) { npv = PVI->getTrueNumInteractions(); continue;}
        }        
        if (verbose) cout<<"N true interations = "<<npv<<endl;
        pileupWeight = pileUp.getPUWeight(npv,PUVersion);
        if (verbose) cout<<"pileup weight = "<<pileupWeight<<", filling histograms"<<endl;
        histContainer_["NINTERACT"]->Fill(npv);
        histContainer_["NINTERACT_RW"]->Fill(npv,pileupWeight);
    } else { pileupWeight = 1.0;}   

    if (verbose) {cout<<"finished pileup reweighting"<<endl; }
    
    if(isMC) {
        float tmpWeight = genEventInfo->weight();
        genWeight = (tmpWeight > 0 ? 1.0 : -1.0);
        double rms = 0.0;

        //if(lheInfos.size()>0){            
        //edm::Handle<LHEEventProduct> lheInfo = lheInfos.front();
        if(lheInfo.isValid()){
            //for(unsigned int i = 0; i < pdfWeightIDs.size(); i++) {
            for(unsigned int i = 0; i < lheInfo->weights().size(); i++) {
                tmpWeight = genEventInfo->weight();
                //tmpWeight *= lheInfo->weights()[pdfWeightIDs[i]].wgt/lheInfo->originalXWGTUP();
                tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
                pdfWeights.push_back(tmpWeight);
                // NNPDF30 variations
                if (i>=9 && i<=109) {
                    rms += tmpWeight*tmpWeight;
                    if (tmpWeight>pdfENVup) pdfENVup=tmpWeight;
                    if (tmpWeight<pdfENVdown) pdfENVdown=tmpWeight;
                }
            }
            pdfRMSup=sqrt(rms/100.0); pdfRMSdown=1.0/pdfRMSup;
            if (verbose) cout<<"pdfRMSup "<<pdfRMSup<<" pdfRMSdown "<<pdfRMSdown<<endl;
        }
        
        if (verbose) cout<<"setting gen variables"<<endl;       
        setGENVariables(prunedgenParticles,packedgenParticles,genJets); 
        if (verbose) { cout<<"finshed setting gen variables"<<endl;  }

        if (int(GENZ_pt.size()) == 2) {        
            GENZ1Vec.SetPtEtaPhiM(GENZ_pt[0], GENZ_eta[0], GENZ_phi[0], GENZ_mass[0]);
            GENZ2Vec.SetPtEtaPhiM(GENZ_pt[1], GENZ_eta[1], GENZ_phi[1], GENZ_mass[1]);
            vector<vector<float> > tableEwk;
            tableEwk = readFile_and_loadEwkTable("ZZBG");   
            k_qqZZ_ewk = getEwkCorrections(prunedgenParticles, tableEwk, genEventInfo, GENZ1Vec, GENZ2Vec);
            k_ggZZ = kFactor_ggzz->Eval((GENZ1Vec+GENZ2Vec).M());
            if (verbose) cout<<"ZZmass: "<< (GENZ1Vec+GENZ2Vec).M() <<"k factor qqZZ ewk: "<<k_qqZZ_ewk<<" ggZZ qcd"<<k_ggZZ << endl;

        }

    }
    sumWeightsTotal += genWeight;

    eventWeight = pileupWeight*genWeight;

    unsigned int _tSize = trigger->size();
    // create a string with all passing trigger names
    for (unsigned int i=0; i<_tSize; ++i) {
        std::string triggerName = trigNames.triggerName(i);
        if (trigger->accept(i)) triggersPassed += triggerName; 
    }
    if (firstEntry) cout<<"triggersPassed: "<<triggersPassed<<endl;
    firstEntry = false;
    // check if any of the triggers in the user list have passed
    for (unsigned int i=0; i<triggerList.size(); ++i) {
        if (strstr(triggersPassed.c_str(),triggerList.at(i).c_str())) passedTrig=true;
    }

    if (verbose) cout<<"checking PV"<<endl;       
    const reco::Vertex *PV = 0;
    int theVertex = -1;
    for (unsigned int i=0; i<vertex->size(); i++) {
        PV = &(vertex->at(i));        
        if (PV->chi2()==0 && PV->ndof()==0) continue;
        if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
        theVertex=(int)i; break;
    }        
 
    if(theVertex >= 0)  {

        if (verbose) cout<<"good PV "<<theVertex<<endl;       

        //N Vertex 
        if (verbose) {cout<<"fill nvtx histogram"<<endl;}
        nVtx = vertex->size();
        histContainer_["NVTX"]->Fill(nVtx);
        histContainer_["NVTX_RW"]->Fill(nVtx,pileupWeight);

        //MET
        if (verbose) {cout<<"get met value"<<endl;}
        if (!mets->empty()) {
            met = (*mets)[0].et();
            met_phi = (*mets)[0].phi();
        }

        if (verbose) cout<<"start lepton analysis"<<endl;           
        vector<pat::Muon> AllMuons;
        vector<pat::Electron> AllElectrons;  
        AllMuons = helper.goodLooseMuons2012(muons,_muPtCut);
        AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut);

        helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);
        if (verbose) cout<<AllMuons.size()<<" loose muons "<<AllElectrons.size()<<" loose elctrons"<<endl;
        recoMuons = helper.goodMuons2015_noIso_noPf(AllMuons,_muPtCut,PV,sip3dCut);
        recoElectrons = helper.goodElectrons2015_noIso_noBdt(AllElectrons,_elecPtCut,elecID,PV,iEvent,sip3dCut);

        //sort electrons and muons by pt
        if (verbose) cout<<recoMuons.size()<<" good muons and "<<recoElectrons.size()<<" good electrons to be sorted"<<endl;
        if (verbose) cout<<"start pt-sorting leptons"<<endl;
        if (verbose) cout<<"adding muons to sorted list"<<endl;           
        for(unsigned int i = 0; i < recoMuons.size(); i++) {
            if (lep_ptreco.size()==0 || recoMuons[i].pt()<lep_ptreco[lep_ptreco.size()-1]) {
                lep_ptreco.push_back(recoMuons[i].pt());
                lep_ptid.push_back(recoMuons[i].pdgId());
                lep_ptindex.push_back(i);
                continue;
            }
            for (unsigned int j=0; j<lep_ptreco.size(); j++) {
                if (recoMuons[i].pt()>lep_ptreco[j]) {
                    lep_ptreco.insert(lep_ptreco.begin()+j,recoMuons[i].pt());
                    lep_ptid.insert(lep_ptid.begin()+j,recoMuons[i].pdgId());
                    lep_ptindex.insert(lep_ptindex.begin()+j,i);
                    break;
                }
            }
        }
        if (verbose) cout<<"adding electrons to sorted list"<<endl;           
        for(unsigned int i = 0; i < recoElectrons.size(); i++) {
            if (lep_ptreco.size()==0 || recoElectrons[i].pt()<lep_ptreco[lep_ptreco.size()-1]) {
                lep_ptreco.push_back(recoElectrons[i].pt());
                lep_ptid.push_back(recoElectrons[i].pdgId());
                lep_ptindex.push_back(i);
                continue;
            }
            for (unsigned int j=0; j<lep_ptreco.size(); j++) {
                if (recoElectrons[i].pt()>lep_ptreco[j]) {
                    lep_ptreco.insert(lep_ptreco.begin()+j,recoElectrons[i].pt());
                    lep_ptid.insert(lep_ptid.begin()+j,recoElectrons[i].pdgId());
                    lep_ptindex.insert(lep_ptindex.begin()+j,i);
                    break;
                }
            }
        }

        for(unsigned int i = 0; i < lep_ptreco.size(); i++) {
            
            if (verbose) cout<<"sorted lepton "<<i<<" pt "<<lep_ptreco[i]<<" id "<<lep_ptid[i]<<" index "<<lep_ptindex[i]<<endl;

            if (abs(lep_ptid[i])==11) {
                lep_id.push_back(recoElectrons[lep_ptindex[i]].pdgId());
                lep_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
                double perr = 0.0;
                if (recoElectrons[lep_ptindex[i]].ecalDriven()) {
                    perr = recoElectrons[lep_ptindex[i]].p4Error(reco::GsfElectron::P4_COMBINATION);
                }
                else {
                    double ecalEnergy = recoElectrons[lep_ptindex[i]].correctedEcalEnergy();
                    double err2 = 0.0;
                    if (recoElectrons[lep_ptindex[i]].isEB()) {
                        err2 += (5.24e-02*5.24e-02)/ecalEnergy;
                        err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
                        err2 += 1.00e-02*1.00e-02;
                    } else if (recoElectrons[lep_ptindex[i]].isEE()) {
                        err2 += (1.46e-01*1.46e-01)/ecalEnergy;
                        err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
                        err2 += 1.94e-03*1.94e-03;
                    }
                    perr = ecalEnergy * sqrt(err2);
                }
                double pterr = perr*recoElectrons[lep_ptindex[i]].pt()/recoElectrons[lep_ptindex[i]].p();
                lep_pterr.push_back(pterr);
                lep_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
                lep_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
                lep_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
                lepFSR_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
                lepFSR_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
                lepFSR_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
                lepFSR_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
                if (isoConeSizeEl==0.4) {
                    lep_RelIso.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
                    lep_RelIsoNoFSR.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
                    lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].chargedHadronIso());
                    lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].neutralHadronIso());
                    lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].photonIso());
                    lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].puChargedHadronIso());
                    lep_isoPUcorr.push_back(helper.getPUIso(recoElectrons[lep_ptindex[i]],elRho));
                } else if (isoConeSizeEl==0.3) {
                    lep_RelIso.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
                    lep_RelIsoNoFSR.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
                    lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumChargedHadronPt);
                    lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumNeutralHadronEt);
                    lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPhotonEt);
                    lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPUPt);
                    lep_isoPUcorr.push_back(helper.getPUIso03(recoElectrons[lep_ptindex[i]],elRho));
                }
                lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoElectrons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
                lep_Sip.push_back(helper.getSIP3D(recoElectrons[lep_ptindex[i]]));           
                lep_mva.push_back(recoElectrons[lep_ptindex[i]].electronID(elecID)); 
                lep_ecalDriven.push_back(recoElectrons[lep_ptindex[i]].ecalDriven()); 
                lep_tightId.push_back(helper.passTight_BDT_Id(recoElectrons[lep_ptindex[i]],elecID));           
                lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoElectrons[lep_ptindex[i]],elecID,PV,BS,theConversions));           
                lep_ptRatio.push_back(helper.ptRatio(recoElectrons[lep_ptindex[i]],jets,isMC));           
                lep_ptRel.push_back(helper.ptRel(recoElectrons[lep_ptindex[i]],jets,isMC));           
                lep_genindex.push_back(-1.0);
            }
            if (abs(lep_ptid[i])==13) {            
                lep_id.push_back(recoMuons[lep_ptindex[i]].pdgId());
                lep_pt.push_back(recoMuons[lep_ptindex[i]].pt());
                lep_pterr.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->ptError());
                lep_eta.push_back(recoMuons[lep_ptindex[i]].eta());
                lep_phi.push_back(recoMuons[lep_ptindex[i]].phi());
                lep_mass.push_back(recoMuons[lep_ptindex[i]].mass());
                lepFSR_pt.push_back(recoMuons[lep_ptindex[i]].pt());
                lepFSR_eta.push_back(recoMuons[lep_ptindex[i]].eta());
                lepFSR_phi.push_back(recoMuons[lep_ptindex[i]].phi());
                lepFSR_mass.push_back(recoMuons[lep_ptindex[i]].mass());
                if (isoConeSizeMu==0.4) {
                    lep_RelIso.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
                    lep_RelIsoNoFSR.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
                    lep_isoCH.push_back(recoMuons[lep_ptindex[i]].chargedHadronIso());
                    lep_isoNH.push_back(recoMuons[lep_ptindex[i]].neutralHadronIso());
                    lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].photonIso());
                    lep_isoPU.push_back(recoMuons[lep_ptindex[i]].puChargedHadronIso());
                    lep_isoPUcorr.push_back(helper.getPUIso(recoMuons[lep_ptindex[i]],muRho));
                } else if (isoConeSizeMu==0.3) {
                    lep_RelIso.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
                    lep_RelIsoNoFSR.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
                    lep_isoCH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumChargedHadronPt);
                    lep_isoNH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumNeutralHadronEt);
                    lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPhotonEt);
                    lep_isoPU.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPUPt);
                    lep_isoPUcorr.push_back(helper.getPUIso03(recoMuons[lep_ptindex[i]],muRho));
                }
                lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoMuons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
                lep_Sip.push_back(helper.getSIP3D(recoMuons[lep_ptindex[i]]));            
                lep_mva.push_back(recoMuons[lep_ptindex[i]].isPFMuon());  
                lep_ecalDriven.push_back(0);  
                lep_tightId.push_back(recoMuons[lep_ptindex[i]].isPFMuon());         
                lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoMuons[lep_ptindex[i]],PV));
                lep_ptRatio.push_back(helper.ptRatio(recoMuons[lep_ptindex[i]],jets,isMC));           
                lep_ptRel.push_back(helper.ptRel(recoMuons[lep_ptindex[i]],jets,isMC));                      
                lep_genindex.push_back(-1.0);
            }
            if (verbose) {cout<<" RelIso: "<<lep_RelIso[i]<<" isoCH: "<<lep_isoCH[i]<<" isoNH: "<<lep_isoNH[i]
                              <<" isoPhot: "<<lep_isoPhot[i]<<" isoPUcorr: "<<lep_isoPUcorr[i]<<" Sip: "<<lep_Sip[i]
                              <<" MiniIso: "<<lep_MiniIso[i]<<" ptRatio: "<<lep_ptRatio[i]<<" ptRel: "<<lep_ptRel[i]<<endl;}
        }


        if (verbose) cout<<"start trigger matching"<<endl;
        // trigger Matching
        for(unsigned int i = 0; i < lep_pt.size(); i++) {

            TLorentzVector reco;
            reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

            double reco_eta = reco.Eta();
            double reco_phi = reco.Phi();

            std::string filtersMatched = "";

            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    double hlt_eta = obj.eta();
                    double hlt_phi = obj.phi();
                    double dR =  deltaR(reco_eta,reco_phi,hlt_eta,hlt_phi); 
                    if (dR<0.5) {
                        for (unsigned h = 0; h < obj.filterLabels().size(); ++h) filtersMatched += obj.filterLabels()[h];
                    }
            }

            //if (verbose) cout<<"Trigger matching lep id: "<<lep_id[i]<<" pt: "<<reco.Pt()<<" filters: "<<filtersMatched<<endl;
            lep_filtersMatched.push_back(filtersMatched);

        }
            
        // GEN matching
        if(isMC) {
            if (verbose) cout<<"begin gen matching"<<endl;
            // for each reco lepton find the nearest gen lepton with same ID
            for(unsigned int i = 0; i < lep_pt.size(); i++) {

                double minDr=9999.0;
                
                TLorentzVector reco, gen;
                reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

                for (unsigned int j = 0; j < GENlep_id.size(); j++) {

                    if (GENlep_id[j]!=lep_id[i]) continue;

                    gen.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);
                    double thisDr = deltaR(reco.Eta(),reco.Phi(),gen.Eta(),gen.Phi());

                    if (thisDr<minDr && thisDr<0.5) {
                        lep_genindex[i]=j;
                        minDr=thisDr;
                    }

                } // all gen leptons
 
            } // all reco leptons
            if (verbose) {cout<<"finished gen matching"<<endl;}
        } //isMC

        if( (recoMuons.size() + recoElectrons.size()) >= (uint)skimLooseLeptons ) {
            if (verbose) cout<<"found two leptons"<<endl;

            // Mass Resolution Study
            /*
            if((recoMuons.size() + recoElectrons.size()) >=2 && bStudyResolution) {
                if (verbose) cout<<"begin 2lep mass resolution study"<<endl;
                vector<pat::Muon> recoIsoMuons;
                vector<pat::Electron> recoIsoElectrons;	    
                recoIsoMuons = helper.goodMuons2012_Iso(AllMuons,_muPtCut, muRho, isoCutMu, PV);
                recoIsoElectrons = helper.goodElectrons2012_Iso(AllElectrons,_elecPtCut, elRho, isoCutEl, elecID,PV);	   
                PerLepReso->fillHistograms(hContainer_, hContainer2D_, hContainer3D_, recoIsoElectrons, recoIsoMuons, eventWeight, !isMC);
                if(bStudyDiLeptonResolution) { DiLepReso->fillHistograms(hContainer_, hContainer2D_, recoIsoElectrons, recoIsoMuons, eventWeight, !isMC); }
                if (verbose) {cout<<"finished 2lep mass resolution study"<<endl;}// cout << "" << endl;}
            }
            */

            // FSR Photons
            if(doFsrRecovery) {
                
                if (verbose) cout<<"checking "<<photonsForFsr->size()<<" fsr photon candidates"<<endl;
                
                // try to find an fsr photon for each lepton
                unsigned int Nleptons = lep_pt.size();
                for (unsigned int i=0; i<Nleptons; i++) {
                        
                    TLorentzVector thisLep;
                    thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

                    double minDrOEt2 = 999.0; double selectedPhotonIso=9999.0; double selectedPhotonDr=9999.0;
                    bool selected=false; pat::PFParticle selectedPhoton; 
                    for(edm::View<pat::PFParticle>::const_iterator phot=photonsForFsr->begin(); phot!=photonsForFsr->end(); ++phot) {


                        // preselection
                        if (fabs(phot->eta()) > 2.4) continue;
                        if (phot->pt()<2.0) continue;
                        double fsrDr = deltaR(thisLep.Eta(), thisLep.Phi(), phot->eta(), phot->phi());
                        if (fsrDr>0.5) continue;

                        // check super cluster veto against all electrons for each photon
                        // at the same time check that this is the closest lepton for this photon
                        bool matched=false;
                        bool closest=true;                    
                        for (unsigned int j=0; j<Nleptons; j++) {                                            

                            TLorentzVector otherLep;
                            otherLep.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

                            double fsrDrOther = deltaR(otherLep.Eta(), otherLep.Phi(), phot->eta(), phot->phi());
                            if (j!=i && fsrDrOther<fsrDr) {closest=false;}

                            if ( abs(lep_id[(int)j])==11) {                            
                                double fsrDr_eSC = deltaR(recoElectrons[lep_ptindex[j]].superCluster()->eta(), recoElectrons[lep_ptindex[j]].superCluster()->phi(), phot->eta(), phot->phi());
                                double fsrDeltaPhi_eSC = fabs(deltaPhi(phot->phi(),recoElectrons[lep_ptindex[j]].superCluster()->phi()));
                                double fsrDeltaEta_eSC = fabs(phot->eta()-recoElectrons[lep_ptindex[j]].superCluster()->eta());                            
                                if ( fsrDr_eSC<0.15 || (fsrDeltaPhi_eSC<2.0 && fsrDeltaEta_eSC<0.05) ) matched=true;
                            }

                        }
                        if (matched) continue;
                        if (!closest) continue;

                        // comput iso, dR/ET^2
                        double photoniso = (phot->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02")+phot->userFloat("fsrPhotonPFIsoNHadPhoton03"))/phot->pt();
                        double fsrDrOEt2 = fsrDr/(phot->pt()*phot->pt());

                        // fill all fsr photons before iso, dR/Et^2 cuts
                        allfsrPhotons_dR.push_back(fsrDr);
                        allfsrPhotons_iso.push_back(photoniso);
                        allfsrPhotons_pt.push_back(phot->pt());

                        // require photon iso, dR/Et^2
                        if (photoniso>1.8) continue;
                        if (fsrDrOEt2>0.012) continue;

                        // this photon is now a good one, check if it is the best one
                        if ( verbose) cout<<"fsr photon cand, pt: "<<phot->pt()<<" eta: "<<phot->eta()<<" phi: "<<phot->phi()
                                          <<" isoCHPUNoPU: "<<phot->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02")
                                          <<" isoNHPhoton: "<<phot->userFloat("fsrPhotonPFIsoNHadPhoton03")
                                          <<" photoniso: "<<photoniso<<" DrOEt2: "<< fsrDrOEt2 <<"minDrOEt2: "<<minDrOEt2<<endl;                        
                        if( fsrDrOEt2 < minDrOEt2 ) {
                            selected = true;
                            selectedPhoton=(*phot);
                            selectedPhotonIso=photoniso;
                            selectedPhotonDr=fsrDr;
                            minDrOEt2 = fsrDrOEt2;
                        }
                        
                    } // all photons

                    if (selected) {
                        nFSRPhotons++;
                        fsrPhotons.push_back(selectedPhoton); 
                        fsrPhotons_pt.push_back(selectedPhoton.pt());
                        double perr = PFEnergyResolution().getEnergyResolutionEm(selectedPhoton.energy(), selectedPhoton.eta());
                        double pterr = perr*selectedPhoton.pt()/selectedPhoton.p();
                        fsrPhotons_pterr.push_back(pterr);
                        fsrPhotons_eta.push_back(selectedPhoton.eta());
                        fsrPhotons_phi.push_back(selectedPhoton.phi());
                        fsrPhotons_lepindex.push_back((int)i);
                        fsrPhotons_dR.push_back(selectedPhotonDr);
                        fsrPhotons_iso.push_back(selectedPhotonIso);
                        TLorentzVector phofsr;
                        phofsr.SetPtEtaPhiM(selectedPhoton.pt(),selectedPhoton.eta(),selectedPhoton.phi(),0.0);
                        TLorentzVector lepfsr;
                        lepfsr = thisLep+phofsr;
                        lepFSR_pt[i] = lepfsr.Pt();
                        lepFSR_eta[i] = lepfsr.Eta();
                        lepFSR_phi[i] = lepfsr.Phi();
                        lepFSR_mass[i] = lepfsr.M();

                        fsrmap[i] = phofsr;
 
                    }

                } // all leptons

                // subtract selected photons from all leptons isolations
                for (unsigned int i=0; i<Nleptons; i++) {
                        
                    TLorentzVector lep_nofsr;
                    lep_nofsr.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

                    double isoFSR=0.0;
                    for (unsigned int j=0; j<fsrPhotons.size(); j++) {
                        double fsrDr = deltaR(lep_nofsr.Eta(),lep_nofsr.Phi(),fsrPhotons[j].eta(),fsrPhotons[j].phi());
                        bool isoVeto=true;
                        if (abs(lep_id[i])==13 && fsrDr>0.01) isoVeto=false;
                        if (abs(lep_id[i])==11 && (abs(recoElectrons[lep_ptindex[i]].superCluster()->eta())<1.479 || fsrDr>0.08)) isoVeto=false;
                        if (fsrDr<((abs(lep_id[i])==11)?isoConeSizeEl:isoConeSizeMu) && !isoVeto) isoFSR += fsrPhotons[j].pt();
                    }

                    double RelIsoNoFSR = (lep_isoCH[i]+std::max(lep_isoNH[i]+lep_isoPhot[i]-lep_isoPUcorr[i]-isoFSR,0.0))/lep_nofsr.Pt();
                    lep_RelIsoNoFSR[i] = RelIsoNoFSR;               
                    if (verbose) cout<<"lep pt: "<<lep_nofsr.Pt()<<" eta: "<<lep_nofsr.Eta()<<" phi: "<<lep_nofsr.Phi()<<" RelIsoNoFSR: "<<RelIsoNoFSR<<" tightId? "<<lep_tightId[i]<<endl;
                }
                if (verbose) {cout<<"finished filling fsr photon candidates"<<endl;}
            } // doFsrRecovery

            // creat vectors for selected objects
            vector<pat::Muon> selectedMuons;
            vector<pat::Electron> selectedElectrons;
                        
            if (verbose) cout<<"begin looking for higgs candidate"<<endl;                    
            findHiggsCandidate(selectedMuons,selectedElectrons,iEvent);
            if (verbose) {cout<<"found higgs candidate? "<<foundHiggsCandidate<<endl; }
            
            // Jets
            vector<pat::Jet> goodJets;
            //double tempDeltaR = -999;
            
            if (verbose) cout<<"begin filling jet candidates"<<endl;                                        
            for(unsigned int i = 0; i < jets->size(); ++i) {
                
                const pat::Jet & jet = jets->at(i);
                
                //JetID ID
                if (verbose) cout<<"checking jetid..."<<endl;
                float jpumva=0.;
                bool passPU=true;
                jpumva=jet.userFloat("pileupJetId:fullDiscriminant");
                if (verbose) cout<< " jet pu mva  "<<jpumva <<endl;
                if(jet.pt()>20){
                    if(abs(jet.eta())>3.){
                        if(jpumva<=-0.45)passPU=false;
                    }else if(abs(jet.eta())>2.75){
                        if(jpumva<=-0.55)passPU=false;
                    }else if(abs(jet.eta())>2.5){
                        if(jpumva<=-0.6)passPU=false;
                    }else if(jpumva<=-0.63)passPU=false;
                }else{
                    if(abs(jet.eta())>3.){
                        if(jpumva<=-0.95)passPU=false;
                    }else if(abs(jet.eta())>2.75){
                        if(jpumva<=-0.94)passPU=false;
                    }else if(abs(jet.eta())>2.5){
                        if(jpumva<=-0.96)passPU=false;
                    }else if(jpumva<=-0.95)passPU=false;
                }
                
                if (verbose) cout<<"pt: "<<jet.pt()<<" eta: "<<jet.eta()<<" passPU: "<<passPU
                                 <<" jetid: "<<jetHelper.patjetID(jet)<<endl;
                
                if( jetHelper.patjetID(jet)>=jetIDLevel && (passPU || !doPUJetID) ) {
                    
                    if (verbose) cout<<"passed pf jet id and pu jet id"<<endl;
                    
                    // apply loose pt cut here (10 GeV cut is already applied in MINIAOD) since we are before JES/JER corrections
                    if(jet.pt() > 10.0 && fabs(jet.eta()) < jeteta_cut) {
                        
                        // apply scale factor for PU Jets by demoting 1-data/MC % of jets jets in certain pt/eta range 
                        // Configured now that SF is 1.0
                        if (verbose) cout<<"adding pu jet scale factors..."<<endl;       
                        bool dropit=false;
                        if (abs(jet.eta())>3.0 && isMC) {
                            TRandom3 rand;
                            rand.SetSeed(abs(static_cast<int>(sin(jet.phi())*100000)));
                            float coin = rand.Uniform(1.); 
                            if (jet.pt()>=20.0 && jet.pt()<36.0 && coin>1.0) dropit=true;
                            if (jet.pt()>=36.0 && jet.pt()<50.0 && coin>1.0) dropit=true;
                            if (jet.pt()>=50.0 && coin>1.0) dropit=true;
                        }                                        
                        
                        if (!dropit) {
                            if (verbose) cout<<"adding jet candidate, pt: "<<jet.pt()<<" eta: "<<jet.eta()<<endl;
                            goodJets.push_back(jet);
                        } // pu jet scale factor
                        
                    } // pass loose pt cut 
                } // pass loose pf jet id and pu jet id
            } // all jets
                        
            if( foundHiggsCandidate ){
                
                for(unsigned int i = 0; i<4;i++){
                    
                    int index = lep_Hindex[i];
                    if(fsrmap[index].Pt()!=0){
                        if (verbose) cout<<"find a fsr photon for "<<i<<" th Higgs lepton"<<endl;
                        selectedFsrMap[i] = fsrmap[index];
                    }
                }
                
                
                if (verbose) cout<<"storing H_p4_noFSR"<<endl; 
                math::XYZTLorentzVector tmpHVec;
                if (verbose) cout<<"selectedMuons "<<selectedMuons.size()<<" selectedElectrons "<<selectedElectrons.size()<<endl;
                if(RecoFourMuEvent) {
                    tmpHVec = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedMuons[2].p4() + selectedMuons[3].p4();
                    HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                    
                    reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedMuons[0]);
                    selectedLeptons.push_back(c1);
                    reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedMuons[1]);
                    selectedLeptons.push_back(c2);
                    reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedMuons[2]);
                    selectedLeptons.push_back(c3);
                    reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedMuons[3]);
                    selectedLeptons.push_back(c4);
                    
                }
                else if(RecoFourEEvent) {
                    tmpHVec = selectedElectrons[0].p4() + selectedElectrons[1].p4() + selectedElectrons[2].p4() + selectedElectrons[3].p4();
                    HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                    
                    reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedElectrons[0]);
                    selectedLeptons.push_back(c1);
                    reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedElectrons[1]);
                    selectedLeptons.push_back(c2);
                    reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedElectrons[2]);
                    selectedLeptons.push_back(c3);
                    reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedElectrons[3]);
                    selectedLeptons.push_back(c4);
                    
                }
                else if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){
                    tmpHVec = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedElectrons[0].p4() + selectedElectrons[1].p4();
                    HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                    
                    if(RecoTwoETwoMuEvent){
                        
                        reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedElectrons[0]);
                        selectedLeptons.push_back(c1);
                        reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedElectrons[1]);
                        selectedLeptons.push_back(c2);
                        reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedMuons[0]);
                        selectedLeptons.push_back(c3);
                        reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedMuons[1]);
                        selectedLeptons.push_back(c4);
                    }
                    else{
                        
                        reco::Candidate *c1 = dynamic_cast<reco::Candidate* >(&selectedMuons[0]);
                        selectedLeptons.push_back(c1);
                        reco::Candidate *c2 = dynamic_cast<reco::Candidate* >(&selectedMuons[1]);
                        selectedLeptons.push_back(c2);
                        reco::Candidate *c3 = dynamic_cast<reco::Candidate* >(&selectedElectrons[0]);
                        selectedLeptons.push_back(c3);
                        reco::Candidate *c4 = dynamic_cast<reco::Candidate* >(&selectedElectrons[1]);
                        selectedLeptons.push_back(c4);
                    }
                    
                    
                }
                H_noFSR_pt.push_back(HVecNoFSR.Pt());
                H_noFSR_eta.push_back(HVecNoFSR.Eta());
                H_noFSR_phi.push_back(HVecNoFSR.Phi());
                H_noFSR_mass.push_back(HVecNoFSR.M());
                
                if(Z2Vec.M() > 0) {
                    passedZ4lSelection = true;
                    if(Z2Vec.M() > mZ2Low && passedTrig) passedFullSelection = true;
                }
                
                
            } // found higgs candidate 
            else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed higgs candidate"<<endl;}


            //Set All the Variables for Saved Trees (after finding higgs candidate)
            if (verbose) cout<<"begin setting tree variables"<<endl;
            setTreeVariables(iEvent, iSetup, selectedMuons, selectedElectrons, recoMuons, recoElectrons, goodJets);
            if (verbose) cout<<"finshed setting tree variables"<<endl;


            // Comput Matrix Elelements (after filling jets)
            if (foundHiggsCandidate) {

                if (verbose) cout<<"Kin fitter begin with lep size "<<selectedLeptons.size()<<" fsr size "<<selectedFsrMap.size()<<endl;

                kinZfitter->Setup(selectedLeptons, selectedFsrMap);
                kinZfitter->KinRefitZ1();
                
                mass4lREFIT = (float)kinZfitter->GetRefitM4l();
                mass4lErrREFIT = (float)kinZfitter->GetRefitM4lErrFullCov();
                mass4lErr = (float)kinZfitter->GetM4lErr();
                massZ1REFIT = (float)kinZfitter->GetRefitMZ1();
                massZ1Err = (float)kinZfitter->GetMZ1Err();
                
                if (verbose) cout<<"mass4l "<<mass4l<<" mass4lREFIT "<<mass4lREFIT<<" massErr "<<mass4lErr<<" massErrREFIT "<<mass4lErrREFIT<<" massZ1REFIT "<<massZ1REFIT<<endl;

                int tmpIdL1,tmpIdL2,tmpIdL3,tmpIdL4;

                TLorentzVector  L11P4, L12P4, L21P4, L22P4, J1P4,  J2P4;
                TLorentzVector Lep1, Lep2, Lep3, Lep4,  Jet1, Jet2;

                TLorentzVector nullFourVector(0, 0, 0, 0);//hualin                                                                                                                                                                                   
                Lep1.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[0]],lepFSR_eta[lep_Hindex[0]],lepFSR_phi[lep_Hindex[0]],lepFSR_mass[lep_Hindex[0]]);
                Lep2.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[1]],lepFSR_eta[lep_Hindex[1]],lepFSR_phi[lep_Hindex[1]],lepFSR_mass[lep_Hindex[1]]);
                Lep3.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[2]],lepFSR_eta[lep_Hindex[2]],lepFSR_phi[lep_Hindex[2]],lepFSR_mass[lep_Hindex[2]]);
                Lep4.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[3]],lepFSR_eta[lep_Hindex[3]],lepFSR_phi[lep_Hindex[3]],lepFSR_mass[lep_Hindex[3]]);
                
                if (njets_pt30_eta4p7 > 0) {
                    J1P4.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_mass[0]);
                }
                if (njets_pt30_eta4p7 > 1) {
                    J2P4.SetPtEtaPhiM(jet_pt[1],jet_eta[1],jet_phi[1],jet_mass[1]);
                }
                
                L11P4.SetPxPyPzE(Lep1.Px(),Lep1.Py(),Lep1.Pz(),Lep1.E()); tmpIdL1 = idL1;
                L12P4.SetPxPyPzE(Lep2.Px(),Lep2.Py(),Lep2.Pz(),Lep2.E()); tmpIdL2 = idL2;
                L21P4.SetPxPyPzE(Lep3.Px(),Lep3.Py(),Lep3.Pz(),Lep3.E()); tmpIdL3 = idL3;
                L22P4.SetPxPyPzE(Lep4.Px(),Lep4.Py(),Lep4.Pz(),Lep4.E()); tmpIdL4 = idL4;
                
                vector<TLorentzVector> P4s;
                vector<int> tmpIDs;
                vector<TLorentzVector> partPprod;
                vector<int> partIdprod;

                P4s.push_back(L11P4); P4s.push_back(L12P4);
                P4s.push_back(L21P4); P4s.push_back(L22P4);

                tmpIDs.push_back(tmpIdL1); tmpIDs.push_back(tmpIdL2);
                tmpIDs.push_back(tmpIdL3); tmpIDs.push_back(tmpIdL4);

                partPprod.push_back(L11P4); partPprod.push_back(L12P4);
                partPprod.push_back(L21P4); partPprod.push_back(L22P4);
                partPprod.push_back(njets_pt30_eta4p7 > 0 ? J1P4 : nullFourVector);
                partPprod.push_back(njets_pt30_eta4p7 > 1 ? J2P4 : nullFourVector);

                partIdprod.push_back(tmpIdL1); partIdprod.push_back(tmpIdL2);
                partIdprod.push_back(tmpIdL3); partIdprod.push_back(tmpIdL4);
                partIdprod.push_back(0); partIdprod.push_back(0);

                combinedMEM->computePm4l(P4s,tmpIDs,MEMNames::kNone,p0plus_m4l,bkg_m4l);
                combinedMEM->computeME(MEMNames::kSMHiggs     ,MEMNames::kJHUGen    ,P4s, tmpIDs,me_0plus_JHU); // higgs, vector algebra, JHUgen
                combinedMEM->computeME(MEMNames::kqqZZ        ,MEMNames::kMCFM      ,P4s, tmpIDs,me_qqZZ_MCFM); // background, vector algebra, MCFM
                combinedMEM->computeME (MEMNames::k0minus, MEMNames::kJHUGen, P4s, tmpIDs, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME
                combinedMEM->computeME (MEMNames::kggHZZ_10, MEMNames::kMCFM, P4s, tmpIDs, Dgg10_VAMCFM); // Direct calculation of Dgg (D^kin for off-shell) from MCFM MEs
                
                if (njets_pt30_eta4p7>=2){
                    combinedMEM->computeME (MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partPprod, partIdprod, phjj_VAJHU); // SM gg->H+2j
                    combinedMEM->computeME (MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_VAJHU); // SM VBF->H
                    Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
                } else {
                    Djet_VAJHU = -1;
                }
                
                D_bkg_kin = me_0plus_JHU / (me_0plus_JHU + me_qqZZ_MCFM);
                D_bkg = me_0plus_JHU * p0plus_m4l / (me_0plus_JHU * p0plus_m4l + me_qqZZ_MCFM * bkg_m4l); // superMELA 
                D_g4 = me_0plus_JHU / ( me_0plus_JHU + p0minus_VAJHU ); // D_0- 
                
                mela::computeAngles(P4s[0], tmpIDs[0], P4s[1], tmpIDs[1], P4s[2], tmpIDs[2], P4s[3], tmpIDs[3], cosThetaStar,cosTheta1,cosTheta2,Phi,Phi1);

                if (verbose) cout<<"D_bkg_kin: "<<D_bkg_kin<< ", D_bkg: " << D_bkg << ", Dgg: " << Dgg10_VAMCFM << ", HJJ_VBF: " << Djet_VAJHU << " ,D0-: " << D_g4 << endl;
                if (verbose) cout<<"cosThetaStar: "<<cosThetaStar<< ", cosTheta1: " << cosTheta1 << ", cosTheta2: " << cosTheta2 << ", Phi: " << Phi << " , Phi1: " << Phi1 << endl;
                
            }

            // fill the vector<float>
            lep_pt_float.assign(lep_pt.begin(),lep_pt.end());
            lep_pterr_float.assign(lep_pterr.begin(),lep_pterr.end());
            lep_eta_float.assign(lep_eta.begin(),lep_eta.end());
            lep_phi_float.assign(lep_phi.begin(),lep_phi.end());
            lep_mass_float.assign(lep_mass.begin(),lep_mass.end());
            lepFSR_pt_float.assign(lepFSR_pt.begin(),lepFSR_pt.end());
            lepFSR_eta_float.assign(lepFSR_eta.begin(),lepFSR_eta.end());
            lepFSR_phi_float.assign(lepFSR_phi.begin(),lepFSR_phi.end());
            lepFSR_mass_float.assign(lepFSR_mass.begin(),lepFSR_mass.end());
            H_pt_float.assign(H_pt.begin(),H_pt.end());
            H_eta_float.assign(H_eta.begin(),H_eta.end());
            H_phi_float.assign(H_phi.begin(),H_phi.end());
            H_mass_float.assign(H_mass.begin(),H_mass.end());
            H_noFSR_pt_float.assign(H_noFSR_pt.begin(),H_noFSR_pt.end());
            H_noFSR_eta_float.assign(H_noFSR_eta.begin(),H_noFSR_eta.end());
            H_noFSR_phi_float.assign(H_noFSR_phi.begin(),H_noFSR_phi.end());
            H_noFSR_mass_float.assign(H_noFSR_mass.begin(),H_noFSR_mass.end());
            Z_pt_float.assign(Z_pt.begin(),Z_pt.end());
            Z_eta_float.assign(Z_eta.begin(),Z_eta.end());
            Z_phi_float.assign(Z_phi.begin(),Z_phi.end());
            Z_mass_float.assign(Z_mass.begin(),Z_mass.end());
            Z_noFSR_pt_float.assign(Z_noFSR_pt.begin(),Z_noFSR_pt.end());
            Z_noFSR_eta_float.assign(Z_noFSR_eta.begin(),Z_noFSR_eta.end());
            Z_noFSR_phi_float.assign(Z_noFSR_phi.begin(),Z_noFSR_phi.end());
            Z_noFSR_mass_float.assign(Z_noFSR_mass.begin(),Z_noFSR_mass.end());
            jet_pt_float.assign(jet_pt.begin(),jet_pt.end());
            jet_eta_float.assign(jet_eta.begin(),jet_eta.end());
            jet_phi_float.assign(jet_phi.begin(),jet_phi.end());
            jet_mass_float.assign(jet_mass.begin(),jet_mass.end());
            jet_jesup_pt_float.assign(jet_jesup_pt.begin(),jet_jesup_pt.end());
            jet_jesup_eta_float.assign(jet_jesup_eta.begin(),jet_jesup_eta.end());
            jet_jesup_phi_float.assign(jet_jesup_phi.begin(),jet_jesup_phi.end());
            jet_jesup_mass_float.assign(jet_jesup_mass.begin(),jet_jesup_mass.end());
            jet_jesdn_pt_float.assign(jet_jesdn_pt.begin(),jet_jesdn_pt.end());
            jet_jesdn_eta_float.assign(jet_jesdn_eta.begin(),jet_jesdn_eta.end());
            jet_jesdn_phi_float.assign(jet_jesdn_phi.begin(),jet_jesdn_phi.end());
            jet_jesdn_mass_float.assign(jet_jesdn_mass.begin(),jet_jesdn_mass.end());
            jet_jerup_pt_float.assign(jet_jerup_pt.begin(),jet_jerup_pt.end());
            jet_jerup_eta_float.assign(jet_jerup_eta.begin(),jet_jerup_eta.end());
            jet_jerup_phi_float.assign(jet_jerup_phi.begin(),jet_jerup_phi.end());
            jet_jerup_mass_float.assign(jet_jerup_mass.begin(),jet_jerup_mass.end());
            jet_jerdn_pt_float.assign(jet_jerdn_pt.begin(),jet_jerdn_pt.end());
            jet_jerdn_eta_float.assign(jet_jerdn_eta.begin(),jet_jerdn_eta.end());
            jet_jerdn_phi_float.assign(jet_jerdn_phi.begin(),jet_jerdn_phi.end());
            jet_jerdn_mass_float.assign(jet_jerdn_mass.begin(),jet_jerdn_mass.end());
            fsrPhotons_pt_float.assign(fsrPhotons_pt.begin(),fsrPhotons_pt.end());
            fsrPhotons_pterr_float.assign(fsrPhotons_pterr.begin(),fsrPhotons_pterr.end());
            fsrPhotons_eta_float.assign(fsrPhotons_eta.begin(),fsrPhotons_eta.end());
            fsrPhotons_phi_float.assign(fsrPhotons_phi.begin(),fsrPhotons_phi.end());
            fsrPhotons_mass_float.assign(fsrPhotons_mass.begin(),fsrPhotons_mass.end());

            if (!isMC) passedEventsTree_All->Fill();        

        } //if 2 lepID
        else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed  2 ID"<<endl;}
    } //primary vertex,notDuplicate
    else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed primary vertex"<<endl;}
    

    GENlep_pt_float.clear(); GENlep_pt_float.assign(GENlep_pt.begin(),GENlep_pt.end());
    GENlep_eta_float.clear(); GENlep_eta_float.assign(GENlep_eta.begin(),GENlep_eta.end());
    GENlep_phi_float.clear(); GENlep_phi_float.assign(GENlep_phi.begin(),GENlep_phi.end());
    GENlep_mass_float.clear(); GENlep_mass_float.assign(GENlep_mass.begin(),GENlep_mass.end());
    GENH_pt_float.clear(); GENH_pt_float.assign(GENH_pt.begin(),GENH_pt.end());
    GENH_eta_float.clear(); GENH_eta_float.assign(GENH_eta.begin(),GENH_eta.end());
    GENH_phi_float.clear(); GENH_phi_float.assign(GENH_phi.begin(),GENH_phi.end());
    GENH_mass_float.clear(); GENH_mass_float.assign(GENH_mass.begin(),GENH_mass.end());
    GENZ_pt_float.clear(); GENZ_pt_float.assign(GENZ_pt.begin(),GENZ_pt.end());
    GENZ_eta_float.clear(); GENZ_eta_float.assign(GENZ_eta.begin(),GENZ_eta.end());
    GENZ_phi_float.clear(); GENZ_phi_float.assign(GENZ_phi.begin(),GENZ_phi.end());
    GENZ_mass_float.clear(); GENZ_mass_float.assign(GENZ_mass.begin(),GENZ_mass.end());
    GENjet_pt_float.clear(); GENjet_pt_float.assign(GENjet_pt.begin(),GENjet_pt.end());
    GENjet_eta_float.clear(); GENjet_eta_float.assign(GENjet_eta.begin(),GENjet_eta.end());
    GENjet_phi_float.clear(); GENjet_phi_float.assign(GENjet_phi.begin(),GENjet_phi.end());
    GENjet_mass_float.clear(); GENjet_mass_float.assign(GENjet_mass.begin(),GENjet_mass.end());

    if (isMC) passedEventsTree_All->Fill();

}



// ------------ method called once each job just before starting event loop  ------------
void 
UFHZZ4LAna::beginJob()
{
    using namespace edm;
    using namespace std;
    using namespace pat;

    bookPassedEventTree("passedEvents", passedEventsTree_All);

    bookResolutionHistograms();
    if(bStudyResolution){
        PerLepReso->bookHistograms(fs, hContainer_, hContainer2D_, hContainer3D_);
        if(bStudyDiLeptonResolution) { DiLepReso->bookHistograms(fs, hContainer_, hContainer2D_); }
        if(bStudyFourLeptonResolution) { FourLepReso->bookHistograms(fs, hContainer_, hContainer2D_);}
    }

    firstEntry = true;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
UFHZZ4LAna::endJob() 
{
    histContainer_["NEVENTS"]->SetBinContent(1,nEventsTotal);
    histContainer_["NEVENTS"]->GetXaxis()->SetBinLabel(1,"N Events in Sample");
    histContainer_["SUMWEIGHTS"]->SetBinContent(1,sumWeightsTotal);
    histContainer_["SUMWEIGHTS"]->GetXaxis()->SetBinLabel(1,"sum Weights in Sample");
}

void
UFHZZ4LAna::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    massErr.init(iSetup);
    if (isMC) {
        pdfWeightIDs.clear();
        edm::Handle<LHERunInfoProduct> run;
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        try {
            iRun.getByToken( lheRunInfoToken_, run );
            LHERunInfoProduct myLHERunInfoProduct = *(run.product());
            for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
                std::cout << iter->tag() << std::endl;
                std::vector<std::string> lines = iter->lines();
                for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                    std::string pdfid=lines.at(iLine);
                    if (pdfid.substr(1,6)=="weight" && pdfid.substr(8,2)=="id") {
                        std::cout<<pdfid<<std::endl;
                        std::string pdf_weight_id = pdfid.substr(12,4);
                        int pdf_weightid=atoi(pdf_weight_id.c_str());
                        std::cout<<"parsed id: "<<pdf_weightid<<std::endl;
                        pdfWeightIDs.push_back(pdf_weightid);
                    }
                }
            }
        }
        catch(...) {
            std::cout<<"No LHERunInfoProduct"<<std::endl;
        }
    }

}


// ------------ method called when ending the processing of a run  ------------
void 
UFHZZ4LAna::endRun(const edm::Run& iRun, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
UFHZZ4LAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
UFHZZ4LAna::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup)
{
    using namespace edm;
    using namespace std;
    // Keep track of all the events run over
    edm::Handle<MergeableCounter> numEventsCounter;
    lumiSeg.getByLabel("nEventsTotal", numEventsCounter);    
    if(numEventsCounter.isValid()) {
        std::cout<<"numEventsCounter->value "<<numEventsCounter->value<<endl;
        nEventsTotal += numEventsCounter->value;        
    }
}

// ============================ UF Functions ============================= //


//Find Z1,Z2, and Higgs candidate
//Pass good leptons for analysis as candMuons and candElectrons
//Pass empty vectors of pat leptons as selectedMuons and selectedElectrons
// these will be filled in the function and then useable for more analysis.
void
UFHZZ4LAna::findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons,const edm::Event& iEvent )
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    const double Zmass = 91.1876;

    unsigned int Nlep = lepFSR_pt.size();
    if (verbose) cout<<Nlep<<" leptons in total"<<endl;

    // First, make all Z candidates including any FSR photons
    int n_Zs=0;
    vector<int> Z_lepindex1;
    vector<int> Z_lepindex2;

    for(unsigned int i=0; i<Nlep; i++){
        if (!(lep_tightId[i])) continue; // checking tight lepton ID
        if (lep_RelIsoNoFSR[i]>((abs(lep_id[i])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
        for(unsigned int j=i+1; j<Nlep; j++){
            if (!(lep_tightId[j])) continue; // checking tight lepton ID
            if (lep_RelIsoNoFSR[j]>((abs(lep_id[j])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed

            // same flavor opposite charge
            if((lep_id[i]+lep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
            lj.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

            TLorentzVector lifsr, ljfsr;
            lifsr.SetPtEtaPhiM(lepFSR_pt[i],lepFSR_eta[i],lepFSR_phi[i],lepFSR_mass[i]);
            ljfsr.SetPtEtaPhiM(lepFSR_pt[j],lepFSR_eta[j],lepFSR_phi[j],lepFSR_mass[j]);

            TLorentzVector liljfsr = lifsr+ljfsr;

            if (verbose) {
                cout<<"OSSF pair: i="<<i<<" id1="<<lep_id[i]<<" j="<<j<<" id2="<<lep_id[j]<<" pt1: "<<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;    
            }
            
            TLorentzVector Z, Z_noFSR;
            Z = lifsr+ljfsr;
            Z_noFSR = li+lj;
            
            if (verbose) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;

            if (Z.M()>0.0) {
                n_Zs++;
                Z_pt.push_back(Z.Pt());
                Z_eta.push_back(Z.Eta());
                Z_phi.push_back(Z.Phi());
                Z_mass.push_back(Z.M());
                Z_noFSR_pt.push_back(Z_noFSR.Pt());
                Z_noFSR_eta.push_back(Z_noFSR.Eta());
                Z_noFSR_phi.push_back(Z_noFSR.Phi());
                Z_noFSR_mass.push_back(Z_noFSR.M());
                Z_lepindex1.push_back(i);
                Z_lepindex2.push_back(j);
                if (verbose) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
            }

        } // lep i
    } // lep j
    
    if( (recoMuons.size() + recoElectrons.size()) < 4 ) return;
        
    if (verbose) cout<<"found four leptons"<<endl;     

    bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
    for(unsigned int i =0; i<recoMuons.size(); i++) {
        if(recoMuons[i].pdgId()<0) Nmm = Nmm+1;
        if(recoMuons[i].pdgId()>0) Nmp = Nmp+1;
    }
    for(unsigned int i =0; i<recoElectrons.size(); i++) {
        if(recoElectrons[i].pdgId()<0) Nem = Nem+1;
        if(recoElectrons[i].pdgId()>0) Nep = Nep+1;
    }
    
    if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
    if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
    if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu
    
    // four proper charge flavor combination
    if(!properLep_ID) return;

    // Consider all ZZ candidates
    double minZ1DeltaM=9999.9;
    double maxZ2SumPt=0.0;
    for (int i=0; i<n_Zs; i++) {
        for (int j=i+1; j<n_Zs; j++) {
 
            int i1 = Z_lepindex1[i]; int i2 = Z_lepindex2[i];                            
            int j1 = Z_lepindex1[j]; int j2 = Z_lepindex2[j];                            

            if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;

            TLorentzVector lep_i1, lep_i2, lep_j1, lep_j2;
            lep_i1.SetPtEtaPhiM(lepFSR_pt[i1],lepFSR_eta[i1],lepFSR_phi[i1],lepFSR_mass[i1]);
            lep_i2.SetPtEtaPhiM(lepFSR_pt[i2],lepFSR_eta[i2],lepFSR_phi[i2],lepFSR_mass[i2]);
            lep_j1.SetPtEtaPhiM(lepFSR_pt[j1],lepFSR_eta[j1],lepFSR_phi[j1],lepFSR_mass[j1]);
            lep_j2.SetPtEtaPhiM(lepFSR_pt[j2],lepFSR_eta[j2],lepFSR_phi[j2],lepFSR_mass[j2]);

            TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr, lep_j2_nofsr;
            lep_i1_nofsr.SetPtEtaPhiM(lep_pt[i1],lep_eta[i1],lep_phi[i1],lep_mass[i1]);
            lep_i2_nofsr.SetPtEtaPhiM(lep_pt[i2],lep_eta[i2],lep_phi[i2],lep_mass[i2]);
            lep_j1_nofsr.SetPtEtaPhiM(lep_pt[j1],lep_eta[j1],lep_phi[j1],lep_mass[j1]);
            lep_j2_nofsr.SetPtEtaPhiM(lep_pt[j2],lep_eta[j2],lep_phi[j2],lep_mass[j2]);

            TLorentzVector Zi, Zj;
            Zi.SetPtEtaPhiM(Z_pt[i],Z_eta[i],Z_phi[i],Z_mass[i]);
            Zj.SetPtEtaPhiM(Z_pt[j],Z_eta[j],Z_phi[j],Z_mass[j]);
            
            if (verbose) {cout<<"ZZ candidate Zi->M() "<<Zi.M()<<" Zj->M() "<<Zj.M()<<endl;}

            TLorentzVector Z1, Z2;
            int Z1index, Z2index;
            int Z1_lepindex[2] = {0,0};
            int Z2_lepindex[2] = {0,0};
            double Z1DeltaM, Z2SumPt;

            if (abs(Zi.M()-Zmass)<abs(Zj.M()-Zmass)) { 
                Z1index = i; Z2index = j;
                Z1 = Zi; Z2 = Zj;                 
                if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
                else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }                
                if (lep_j1.Pt()>lep_j2.Pt()) { Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2; } 
                else { Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1; }                
                Z1DeltaM = abs(Zi.M()-Zmass);                
                Z2SumPt = lep_j1_nofsr.Pt()+lep_j2_nofsr.Pt();
            }
            else { 
                Z1index = j; Z2index = i;
                Z1 = Zj; Z2 = Zi; 
                if (lep_j1.Pt()>lep_j2.Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
                else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
                if (lep_i1.Pt()>lep_i2.Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
                else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
                Z1DeltaM = abs(Zj.M()-Zmass);
                Z2SumPt = lep_i1_nofsr.Pt()+lep_i2_nofsr.Pt();
            }
           
            // Check Leading and Subleading pt Cut
            vector<double> allPt;
            allPt.push_back(lep_i1_nofsr.Pt()); allPt.push_back(lep_i2_nofsr.Pt());
            allPt.push_back(lep_j1_nofsr.Pt()); allPt.push_back(lep_j2_nofsr.Pt());
            std::sort(allPt.begin(), allPt.end());
            if (verbose) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
            if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;
            
            // Check dR(li,lj)>0.02 for any i,j
            vector<double> alldR;
            alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi(),lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi()));
            alldR.push_back(deltaR(lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));
            alldR.push_back(deltaR(lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));            
            if (verbose) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
            if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

            // Check M(l+,l-)>4.0 GeV for any OS pair
            // Do not include FSR photons
            vector<double> allM;
            TLorentzVector i1i2;
            i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
            TLorentzVector j1j2;
            j1j2 = (lep_j1_nofsr)+(lep_j2_nofsr); allM.push_back(j1j2.M());            

            if (lep_id[i1]*lep_id[j1]<0) {
                TLorentzVector i1j1;
                i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
                TLorentzVector i2j2;
                i2j2 = (lep_i2_nofsr)+(lep_j2_nofsr); allM.push_back(i2j2.M());
            } else {
                TLorentzVector i1j2;
                i1j2 = (lep_i1_nofsr)+(lep_j2_nofsr); allM.push_back(i1j2.M());
                TLorentzVector i2j1;
                i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
            }
            if (verbose) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
            if (*min_element(allM.begin(),allM.end())<4.0) {passedQCDcut=false; continue;}

            // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
            // only for 4mu or 4e ZZ candidates
            bool passSmartCut=true;
            if ( abs(lep_id[i1])==abs(lep_id[j1])) {
                TLorentzVector Za, Zb;
                if (lep_id[i1]==lep_id[j1]) {                  
                    Za = (lep_i1)+(lep_j2);
                    Zb = (lep_i2)+(lep_j1);                    
                } else {
                    Za = (lep_i1)+(lep_j1);
                    Zb = (lep_i2)+(lep_j2);
                }                
                if ( abs(Za.M()-Zmass)<abs(Zb.M()-Zmass) ) {
                    if (verbose) cout<<"abs(Za.M()-Zmass)-abs(Z1.M()-Zmass): "<<abs(Za.M()-Zmass)-abs(Z1.M()-Zmass)<<" Zb.M(): "<<Zb.M()<<endl;
                    if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) && Zb.M()<mZ2Low ) passSmartCut=false;
                }
                else {
                    if (verbose) cout<<"abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass): "<<abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass)<<" Za.M(): "<<Za.M()<<endl;
                    if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) && Za.M()<mZ2Low ) passSmartCut=false;
                }

            }
            if (!passSmartCut) continue;
            
            if (verbose) cout<<" massZ1: "<<Z1.M()<<" massZ2: "<<Z2.M()<<endl;
            if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) || (Z2.M() < mZ2Low) || (Z2.M() > mZ2High) ) continue;

            if (verbose) cout<<"good ZZ candidate, Z1DeltaM: "<<Z1DeltaM<<" minZ1DeltaM: "<<minZ1DeltaM<<" Z2SumPt: "<<Z2SumPt<<" maxZ2SumPt: "<<maxZ2SumPt<<endl;

            // Check if this candidate has the best Z1 and highest scalar sum of Z2 lepton pt            

            if ( Z1DeltaM<=minZ1DeltaM ) {

                minZ1DeltaM = Z1DeltaM;

                if (Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt) continue;

                Z_Hindex[0] = Z1index;
                lep_Hindex[0] = Z1_lepindex[0];
                lep_Hindex[1] = Z1_lepindex[1];

                maxZ2SumPt = Z2SumPt;
                Z_Hindex[1] = Z2index;
                lep_Hindex[2] = Z2_lepindex[0];
                lep_Hindex[3] = Z2_lepindex[1];

                Z1Vec = Z1;
                Z2Vec = Z2;
                HVec = Z1+Z2;

                massZ1 = Z1Vec.M();
                massZ2 = Z2Vec.M();
                mass4l = HVec.M();

                if (verbose) cout<<" new best candidate: mass4l: "<<HVec.M()<<endl;
                if (HVec.M()>m4lLowCut) foundHiggsCandidate=true;

            }

            if (verbose) cout<<"Z_Hindex[0]: "<<Z_Hindex[0]<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                             <<"Z_Hindex[1]: "<<Z_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;

        } // Zj
    } // Zi


    if(foundHiggsCandidate) {
        
        if (verbose) cout<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;

        if (verbose) cout<<" lep_id[lep_Hindex[0]]: "<<lep_id[lep_Hindex[0]]<<" lep_id[lep_Hindex[1]]: "<<lep_id[lep_Hindex[1]]
                         <<" lep_id[lep_Hindex[2]]: "<<lep_id[lep_Hindex[2]]<<" lep_id[lep_Hindex[3]]: "<<lep_id[lep_Hindex[3]]<<endl;

        if ( abs(lep_id[lep_Hindex[0]])==13 && abs(lep_id[lep_Hindex[2]])==13 ) {
            RecoFourMuEvent = true;
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[0]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[1]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[2]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[3]]]);
        }
        else if ( abs(lep_id[lep_Hindex[0]])==13 && abs(lep_id[lep_Hindex[2]])==11 ) {
            RecoTwoMuTwoEEvent = true;
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[0]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[1]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[2]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[3]]]);
        }
        else if ( abs(lep_id[lep_Hindex[0]])==11 && abs(lep_id[lep_Hindex[2]])==13 ) {
            RecoTwoETwoMuEvent = true;
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[0]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[1]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[2]]]);
            selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[3]]]);
        }
        else if ( abs(lep_id[lep_Hindex[0]])==11 && abs(lep_id[lep_Hindex[2]])==11 ) {
            RecoFourEEvent = true;
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[0]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[1]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[2]]]);
            selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[3]]]);
        }


    }

}

void UFHZZ4LAna::bookPassedEventTree(TString treeName, TTree *tree)
{     

    using namespace edm;
    using namespace pat;
    using namespace std;


    // -------------------------                                                                                                                                                                        
    // RECO level information                                                                                                                                                                           
    // -------------------------                                                                                                                                                                        
    // Event variables
    tree->Branch("Run",&Run,"Run/l");
    tree->Branch("Event",&Event,"Event/l");
    tree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    tree->Branch("nVtx",&nVtx,"nVtx/I");
    tree->Branch("finalState",&finalState,"finalState/I");
    tree->Branch("triggersPassed",&triggersPassed);
    tree->Branch("passedTrig",&passedTrig,"passedTrig/O");
    tree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    tree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    tree->Branch("passedQCDcut",&passedQCDcut,"passedQCDcut/O");
    tree->Branch("genWeight",&genWeight,"genWeight/F");
    tree->Branch("k_ggZZ",&k_ggZZ,"k_ggZZ/F");
    tree->Branch("k_qqZZ_qcd",&k_qqZZ_qcd,"k_qqZZ_qcd/F");
    tree->Branch("k_qqZZ_ewk",&k_qqZZ_ewk,"k_qqZZ_ewk/F");
    tree->Branch("pdfWeights",&pdfWeights);
    tree->Branch("pdfWeightIDs",&pdfWeightIDs);
    tree->Branch("pdfRMSup",&pdfRMSup,"pdfRMSup/F");
    tree->Branch("pdfRMSdown",&pdfRMSdown,"pdfRMSdown/F");
    tree->Branch("pdfENVup",&pdfENVup,"pdfENVup/F");
    tree->Branch("pdfENVdown",&pdfENVdown,"pdfENVdown/F");
    tree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
    tree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F");
    tree->Branch("eventWeight",&eventWeight,"eventWeight/F");
    tree->Branch("crossSection",&crossSection,"crossSection/F");
    tree->Branch("filterEff",&filterEff,"filterEff/F");

    // Lepton variables
    tree->Branch("lep_id",&lep_id);
    tree->Branch("lep_pt",&lep_pt_float);
    tree->Branch("lep_pterr",&lep_pterr_float);
    tree->Branch("lep_eta",&lep_eta_float);
    tree->Branch("lep_phi",&lep_phi_float);
    tree->Branch("lep_mass",&lep_mass_float);
    tree->Branch("lepFSR_pt",&lepFSR_pt_float);
    tree->Branch("lepFSR_eta",&lepFSR_eta_float);
    tree->Branch("lepFSR_phi",&lepFSR_phi_float);
    tree->Branch("lepFSR_mass",&lepFSR_mass_float);
    tree->Branch("lep_Hindex",&lep_Hindex,"lep_Hindex[4]/I");
    tree->Branch("lep_genindex",&lep_genindex);
    tree->Branch("lep_missingHits",&lep_missingHits);
    tree->Branch("lep_mva",&lep_mva);
    tree->Branch("lep_ecalDriven",&lep_ecalDriven);
    tree->Branch("lep_tightId",&lep_tightId);
    tree->Branch("lep_tightIdSUS",&lep_tightIdSUS);
    tree->Branch("lep_Sip",&lep_Sip);
    tree->Branch("lep_IP",&lep_IP);
    tree->Branch("lep_isoNH",&lep_isoNH);
    tree->Branch("lep_isoCH",&lep_isoCH);
    tree->Branch("lep_isoPhot",&lep_isoPhot);
    tree->Branch("lep_isoPU",&lep_isoPU);
    tree->Branch("lep_isoPUcorr",&lep_isoPUcorr);
    tree->Branch("lep_RelIso",&lep_RelIso);
    tree->Branch("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
    tree->Branch("lep_MiniIso",&lep_MiniIso);
    tree->Branch("lep_ptRatio",&lep_ptRatio);
    tree->Branch("lep_ptRel",&lep_ptRel);
    tree->Branch("lep_filtersMatched",&lep_filtersMatched);
    tree->Branch("nisoleptons",&nisoleptons,"nisoleptons/I");
    tree->Branch("muRho",&muRho,"muRho/F");
    tree->Branch("elRho",&elRho,"elRho/F");
    tree->Branch("pTL1",&pTL1,"pTL1/F");
    tree->Branch("pTL2",&pTL2,"pTL2/F");
    tree->Branch("pTL3",&pTL3,"pTL3/F");
    tree->Branch("pTL4",&pTL4,"pTL4/F");
    tree->Branch("idL1",&idL1,"idL1/I");
    tree->Branch("idL2",&idL2,"idL2/I");
    tree->Branch("idL3",&idL3,"idL3/I");
    tree->Branch("idL4",&idL4,"idL4/I");
    tree->Branch("etaL1",&etaL1,"etaL1/F");
    tree->Branch("etaL2",&etaL2,"etaL2/F");
    tree->Branch("etaL3",&etaL3,"etaL3/F");
    tree->Branch("etaL4",&etaL4,"etaL4/F");
    tree->Branch("pTL1FSR",&pTL1FSR,"pTL1FSR/F");
    tree->Branch("pTL2FSR",&pTL2FSR,"pTL2FSR/F");
    tree->Branch("pTL3FSR",&pTL3FSR,"pTL3FSR/F");
    tree->Branch("pTL4FSR",&pTL4FSR,"pTL4FSR/F");

    //Higgs Candidate Variables
    tree->Branch("H_pt",&H_pt_float);
    tree->Branch("H_eta",&H_eta_float);
    tree->Branch("H_phi",&H_phi_float);
    tree->Branch("H_mass",&H_mass_float);
    tree->Branch("H_noFSR_pt",&H_noFSR_pt_float);
    tree->Branch("H_noFSR_eta",&H_noFSR_eta_float);
    tree->Branch("H_noFSR_phi",&H_noFSR_phi_float);
    tree->Branch("H_noFSR_mass",&H_noFSR_mass_float);
    tree->Branch("mass4l",&mass4l,"mass4l/F");
    tree->Branch("mass4l_noFSR",&mass4l_noFSR,"mass4l_noFSR/F");
    tree->Branch("mass4lErr",&mass4lErr,"mass4lErr/F");
    tree->Branch("mass4lREFIT",&mass4lREFIT,"mass4lREFIT/F");
    tree->Branch("mass4lErrREFIT",&mass4lErrREFIT,"mass4lErrREFIT/F");
    tree->Branch("massZ1REFIT",&massZ1REFIT,"massZ1REFIT/F");
    tree->Branch("massZ1Err",&massZ1Err,"massZ1Err/F");
    tree->Branch("mass4mu",&mass4mu,"mass4mu/F");
    tree->Branch("mass4e",&mass4e,"mass4e/F");
    tree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/F");
    tree->Branch("pT4l",&pT4l,"pT4l/F");
    tree->Branch("eta4l",&eta4l,"eta4l/F");
    tree->Branch("phi4l",&phi4l,"phi4l/F");
    tree->Branch("rapidity4l",&rapidity4l,"rapidity4l/F");
    tree->Branch("cosTheta1",&cosTheta1,"cosTheta1/F");
    tree->Branch("cosTheta2",&cosTheta2,"cosTheta2/F");
    tree->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
    tree->Branch("Phi",&Phi,"Phi/F");
    tree->Branch("Phi1",&Phi1,"Phi1/F");

    // Z candidate variables
    tree->Branch("Z_pt",&Z_pt_float);
    tree->Branch("Z_eta",&Z_eta_float);
    tree->Branch("Z_phi",&Z_phi_float);
    tree->Branch("Z_mass",&Z_mass_float);
    tree->Branch("Z_noFSR_pt",&Z_noFSR_pt_float);
    tree->Branch("Z_noFSR_eta",&Z_noFSR_eta_float);
    tree->Branch("Z_noFSR_phi",&Z_noFSR_phi_float);
    tree->Branch("Z_noFSR_mass",&Z_noFSR_mass_float);
    tree->Branch("Z_Hindex",&Z_Hindex,"Z_Hindex[2]/I");
    tree->Branch("massZ1",&massZ1,"massZ1/F");
    tree->Branch("massZ2",&massZ2,"massZ2/F");  
    tree->Branch("pTZ1",&pTZ1,"pTZ1/F");
    tree->Branch("pTZ2",&pTZ2,"pTZ2/F");

    // MET
    tree->Branch("met",&met,"met/F");
    tree->Branch("met_phi",&met_phi,"met_phi/F");

    // Jets
    tree->Branch("jet_pt",&jet_pt_float);
    tree->Branch("jet_eta",&jet_eta_float);
    tree->Branch("jet_phi",&jet_phi_float);
    tree->Branch("jet_mass",&jet_mass_float);
    tree->Branch("jet_jesup_pt",&jet_jesup_pt_float);
    tree->Branch("jet_jesup_eta",&jet_jesup_eta_float);
    tree->Branch("jet_jesup_phi",&jet_jesup_phi_float);
    tree->Branch("jet_jesup_mass",&jet_jesup_mass_float);
    tree->Branch("jet_jesdn_pt",&jet_jesdn_pt_float);
    tree->Branch("jet_jesdn_eta",&jet_jesdn_eta_float);
    tree->Branch("jet_jesdn_phi",&jet_jesdn_phi_float);
    tree->Branch("jet_jesdn_mass",&jet_jesdn_mass_float);
    tree->Branch("jet_jerup_pt",&jet_jerup_pt_float);
    tree->Branch("jet_jerup_eta",&jet_jerup_eta_float);
    tree->Branch("jet_jerup_phi",&jet_jerup_phi_float);
    tree->Branch("jet_jerup_mass",&jet_jerup_mass_float);
    tree->Branch("jet_jerdn_pt",&jet_jerdn_pt_float);
    tree->Branch("jet_jerdn_eta",&jet_jerdn_eta_float);
    tree->Branch("jet_jerdn_phi",&jet_jerdn_phi_float);
    tree->Branch("jet_jerdn_mass",&jet_jerdn_mass_float);
    tree->Branch("jet_pumva",&jet_pumva);
    tree->Branch("jet_csvv2",&jet_csvv2);
    tree->Branch("jet_isbtag",&jet_isbtag);
    tree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
    tree->Branch("njets_pt30_eta4p7_jesup",&njets_pt30_eta4p7_jesup,"njets_pt30_eta4p7_jesup/I");
    tree->Branch("njets_pt30_eta4p7_jesdn",&njets_pt30_eta4p7_jesdn,"njets_pt30_eta4p7_jesdn/I");
    tree->Branch("njets_pt30_eta4p7_jerup",&njets_pt30_eta4p7_jerup,"njets_pt30_eta4p7_jerup/I");
    tree->Branch("njets_pt30_eta4p7_jerdn",&njets_pt30_eta4p7_jerdn,"njets_pt30_eta4p7_jerdn/I");
    tree->Branch("nbjets_pt30_eta4p7",&nbjets_pt30_eta4p7,"nbjets_pt30_eta4p7/I");
    tree->Branch("nvjets_pt40_eta2p4",&nvjets_pt40_eta2p4,"nvjets_pt40_eta2p4/I");
    tree->Branch("pt_leadingjet_pt30_eta4p7",&pt_leadingjet_pt30_eta4p7,"pt_leadingjet_pt30_eta4p7/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jesup",&pt_leadingjet_pt30_eta4p7_jesup,"pt_leadingjet_pt30_eta4p7_jesup/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jesdn",&pt_leadingjet_pt30_eta4p7_jesdn,"pt_leadingjet_pt30_eta4p7_jesdn/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jerup",&pt_leadingjet_pt30_eta4p7_jerup,"pt_leadingjet_pt30_eta4p7_jerup/F");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jerdn",&pt_leadingjet_pt30_eta4p7_jerdn,"pt_leadingjet_pt30_eta4p7_jerdn/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7",&absrapidity_leadingjet_pt30_eta4p7,"absrapidity_leadingjet_pt30_eta4p7/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesup",&absrapidity_leadingjet_pt30_eta4p7_jesup,"absrapidity_leadingjet_pt30_eta4p7_jesup/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesdn",&absrapidity_leadingjet_pt30_eta4p7_jesdn,"absrapidity_leadingjet_pt30_eta4p7_jesdn/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerup",&absrapidity_leadingjet_pt30_eta4p7_jerup,"absrapidity_leadingjet_pt30_eta4p7_jerup/F");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerdn",&absrapidity_leadingjet_pt30_eta4p7_jerdn,"absrapidity_leadingjet_pt30_eta4p7_jerdn/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7",&absdeltarapidity_hleadingjet_pt30_eta4p7,"absdeltarapidity_hleadingjet_pt30_eta4p7/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesup/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerup/F");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn/F");
    tree->Branch("DijetMass",&DijetMass,"DijetMass/F");
    tree->Branch("DijetDEta",&DijetDEta,"DijetDEta/F");
    tree->Branch("DijetFisher",&DijetFisher,"DijetFisher/F");

    // FSR Photons
    tree->Branch("nFSRPhotons",&nFSRPhotons,"nFSRPhotons/I");
    tree->Branch("allfsrPhotons_dR",&allfsrPhotons_dR);
    tree->Branch("allfsrPhotons_iso",&allfsrPhotons_iso);
    tree->Branch("allfsrPhotons_pt",&allfsrPhotons_pt);
    tree->Branch("fsrPhotons_lepindex",&fsrPhotons_lepindex);
    tree->Branch("fsrPhotons_pt",&fsrPhotons_pt_float);
    tree->Branch("fsrPhotons_pterr",&fsrPhotons_pterr_float);
    tree->Branch("fsrPhotons_eta",&fsrPhotons_eta_float);
    tree->Branch("fsrPhotons_phi",&fsrPhotons_phi_float);
    tree->Branch("fsrPhotons_dR",&fsrPhotons_dR);
    tree->Branch("fsrPhotons_iso",&fsrPhotons_iso);

    // Z4l? FIXME
    tree->Branch("theta12",&theta12,"theta12/F"); 
    tree->Branch("theta13",&theta13,"theta13/F"); 
    tree->Branch("theta14",&theta14,"theta14/F");
    tree->Branch("minM3l",&minM3l,"minM3l/F"); 
    tree->Branch("Z4lmaxP",&Z4lmaxP,"Z4lmaxP/F"); 
    tree->Branch("minDeltR",&minDeltR,"minDeltR/F"); 
    tree->Branch("m3l_soft",&m3l_soft,"m3l_soft/F");
    tree->Branch("minMass2Lep",&minMass2Lep,"minMass2Lep/F"); 
    tree->Branch("maxMass2Lep",&maxMass2Lep,"maxMass2Lep/F");
    tree->Branch("thetaPhoton",&thetaPhoton,"thetaPhoton/F"); 
    tree->Branch("thetaPhotonZ",&thetaPhotonZ,"thetaPhotonZ/F");

    // Resolution
    tree->Branch("massErrorUF",&massErrorUF,"massErrorUF/F");
    tree->Branch("massErrorUFCorr",&massErrorUFCorr,"massErrorUFCorr/F");
    tree->Branch("massErrorUCSD",&massErrorUCSD,"massErrorUCSD/F");
    tree->Branch("massErrorUCSDCorr",&massErrorUCSDCorr,"massErrorUCSDCorr/F");

    // Event Category
    tree->Branch("EventCat",&EventCat,"EventCat/I");

    // -------------------------                                                                                                                                                                        
    // GEN level information                                                                                                                                                                            
    // -------------------------                                                                                                                                                                        
    //Event variables
    tree->Branch("GENfinalState",&GENfinalState,"GENfinalState/I");
    tree->Branch("passedFiducialSelection",&passedFiducialSelection,"passedFiducialSelection/O");

    // lepton variables
    tree->Branch("GENlep_pt",&GENlep_pt_float);
    tree->Branch("GENlep_eta",&GENlep_eta_float);
    tree->Branch("GENlep_phi",&GENlep_phi_float);
    tree->Branch("GENlep_mass",&GENlep_mass_float);
    tree->Branch("GENlep_id",&GENlep_id);
    tree->Branch("GENlep_status",&GENlep_status);
    tree->Branch("GENlep_MomId",&GENlep_MomId);
    tree->Branch("GENlep_MomMomId",&GENlep_MomMomId);
    tree->Branch("GENlep_Hindex",&GENlep_Hindex,"GENlep_Hindex[4]/I");
    tree->Branch("GENlep_isoCH",&GENlep_isoCH);
    tree->Branch("GENlep_isoNH",&GENlep_isoNH);
    tree->Branch("GENlep_isoPhot",&GENlep_isoPhot);
    tree->Branch("GENlep_RelIso",&GENlep_RelIso);

    // Higgs candidate variables (calculated using selected gen leptons)
    tree->Branch("GENH_pt",&GENH_pt_float);
    tree->Branch("GENH_eta",&GENH_eta_float);
    tree->Branch("GENH_phi",&GENH_phi_float);
    tree->Branch("GENH_mass",&GENH_mass_float);
    tree->Branch("GENmass4l",&GENmass4l,"GENmass4l/F");
    tree->Branch("GENmass4mu",&GENmass4mu,"GENmass4mu/F");
    tree->Branch("GENmass4e",&GENmass4e,"GENmass4e/F");
    tree->Branch("GENmass2e2mu",&GENmass2e2mu,"GENmass2e2mu/F");
    tree->Branch("GENpT4l",&GENpT4l,"GENpT4l/F");
    tree->Branch("GENeta4l",&GENeta4l,"GENeta4l/F");
    tree->Branch("GENrapidity4l",&GENrapidity4l,"GENrapidity4l/F");
    tree->Branch("GENcosTheta1",&GENcosTheta1,"GENcosTheta1/F");
    tree->Branch("GENcosTheta2",&GENcosTheta2,"GENcosTheta2/F");
    tree->Branch("GENcosThetaStar",&GENcosThetaStar,"GENcosThetaStar/F");
    tree->Branch("GENPhi",&GENPhi,"GENPhi/F");
    tree->Branch("GENPhi1",&GENPhi1,"GENPhi1/F");
    tree->Branch("GENMH",&GENMH,"GENMH/F");

    // Z candidate variables
    tree->Branch("GENZ_pt",&GENZ_pt_float);
    tree->Branch("GENZ_eta",&GENZ_eta_float);
    tree->Branch("GENZ_phi",&GENZ_phi_float);
    tree->Branch("GENZ_mass",&GENZ_mass_float);
    tree->Branch("GENZ_DaughtersId",&GENZ_DaughtersId); 
    tree->Branch("GENZ_MomId",&GENZ_MomId);
    tree->Branch("GENmassZ1",&GENmassZ1,"GENmassZ1/F");
    tree->Branch("GENmassZ2",&GENmassZ2,"GENmassZ2/F");  
    tree->Branch("GENpTZ1",&GENpTZ1,"GENpTZ1/F");
    tree->Branch("GENpTZ2",&GENpTZ2,"GENpTZ2/F");
    tree->Branch("GENdPhiZZ",&GENdPhiZZ,"GENdPhiZZ/F");
    tree->Branch("GENmassZZ",&GENmassZZ,"GENmassZZ/F");
    tree->Branch("GENpTZZ",&GENpTZZ,"GENpTZZ/F");

    // Higgs variables directly from GEN particle
    tree->Branch("GENHmass",&GENHmass,"GENHmass/F");

    // Jets
    tree->Branch("GENjet_pt",&GENjet_pt_float);
    tree->Branch("GENjet_eta",&GENjet_eta_float);
    tree->Branch("GENjet_phi",&GENjet_phi_float);
    tree->Branch("GENjet_mass",&GENjet_mass_float);
    tree->Branch("GENnjets_pt30_eta4p7",&GENnjets_pt30_eta4p7,"GENnjets_pt30_eta4p7/I");
    tree->Branch("GENpt_leadingjet_pt30_eta4p7",&GENpt_leadingjet_pt30_eta4p7,"GENpt_leadingjet_pt30_eta4p7/F");
    tree->Branch("GENabsrapidity_leadingjet_pt30_eta4p7",&GENabsrapidity_leadingjet_pt30_eta4p7,"GENabsrapidity_leadingjet_pt30_eta4p7/F");
    tree->Branch("GENabsdeltarapidity_hleadingjet_pt30_eta4p7",&GENabsdeltarapidity_hleadingjet_pt30_eta4p7,"GENabsdeltarapidity_hleadingjet_pt30_eta4p7/F");

    //ME
    tree->Branch("me_0plus_JHU", &me_0plus_JHU, "me_0plus_JHU/D");
    tree->Branch("me_qqZZ_MCFM", &me_qqZZ_MCFM, "me_qqZZ_MCFM/D");
    tree->Branch("p0plus_m4l", &p0plus_m4l, "p0plus_m4l/D");
    tree->Branch("bkg_m4l", &bkg_m4l, "bkg_m4l/D");
    tree->Branch("D_bkg_kin", &D_bkg_kin, "D_bkg_kin/D");
    tree->Branch("D_bkg", &D_bkg, "D_bkg/D");
    tree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM, "Dgg10_VAMCFM/D");
    tree->Branch("Djet_VAJHU", &Djet_VAJHU, "Djet_VAJHU/D");
    tree->Branch("D_g4", &D_g4, "D_g4/D");

}

void UFHZZ4LAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                   std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, 
                                   std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, std::vector<pat::Jet> goodJets)
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    // Jet Info
    double tempDeltaR = 999.0;
    for( unsigned int k = 0; k < goodJets.size(); k++) {

        if (verbose) cout<<"jet pt: "<<goodJets[k].pt()<<" eta: "<<goodJets[k].eta()<<" phi: "<<goodJets[k].phi()<<endl;

        bool isDeltaR_eta4p7 = true;

        // check overlap with isolated leptons
        unsigned int Nleptons = lep_pt.size();
        for (unsigned int i=0; i<Nleptons; i++) {
            if (abs(lep_id[i])==13 && lep_RelIsoNoFSR[i]>isoCutMu) continue;
            if (abs(lep_id[i])==11 && lep_RelIsoNoFSR[i]>isoCutEl) continue;
            if (!(lep_tightId[i])) continue;
            TLorentzVector thisLep;
            thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),thisLep.Eta(),thisLep.Phi());
            if (tempDeltaR<0.4) {
                isDeltaR_eta4p7 = false;
            }
        }

        // check overlap with fsr photons
        unsigned int N = fsrPhotons_pt.size();
        for(unsigned int i=0; i<N; i++) {
            
            // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton
            if (!lep_tightId[fsrPhotons_lepindex[i]]) continue;
            double RelIsoNoFSR=lep_RelIsoNoFSR[fsrPhotons_lepindex[i]];
            if (RelIsoNoFSR>((abs(lep_id[fsrPhotons_lepindex[i]])==11) ? isoCutEl : isoCutMu)) continue;

            TLorentzVector pho;
            pho.SetPtEtaPhiM(fsrPhotons_pt[i],fsrPhotons_eta[i],fsrPhotons_phi[i],0.0);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),pho.Eta(),pho.Phi());
            if (tempDeltaR<0.4) {
                isDeltaR_eta4p7 = false;
            }
        }

        
        /*
        // JER Smearing
        double factor = 1.0;
        double factorup = 1.0;
        double factordn = 1.0;
        
        if ( abs(goodJets[k].eta()) < 0.5 ) {
            factor = 1.079;
            factorup = 1.105;
            factordn = 1.053;
        }
        else if ( abs(goodJets[k].eta()) < 1.1 && abs(goodJets[k].eta()) >= 0.5 ) {
            factor = 1.099;
            factorup = 1.127;
            factordn = 1.071;
        }
        else if ( abs(goodJets[k].eta()) < 1.7 && abs(goodJets[k].eta()) >= 1.1 ) {
            factor = 1.121;
            factorup = 1.150;
            factordn = 1.092;
        }
        else if ( abs(goodJets[k].eta()) < 2.3 && abs(goodJets[k].eta()) >= 1.7 ) {
            factor = 1.208;
            factorup = 1.254;
            factordn = 1.162;
        }
        else if ( abs(goodJets[k].eta()) < 2.8 && abs(goodJets[k].eta()) >= 2.3 ) {
            factor = 1.254;
            factorup = 1.316;
            factordn = 1.192;
        }
        else if ( abs(goodJets[k].eta()) < 3.2 && abs(goodJets[k].eta()) >= 2.8 ) {
            factor = 1.395;
            factorup = 1.458;
            factordn = 1.332;
        }
        else if ( abs(goodJets[k].eta()) < 5.0 && abs(goodJets[k].eta()) >= 3.2 ) {
            factor = 1.056;
            factorup = 1.247;
            factordn = 0.865;
        }
            
        double pt_jer = goodJets[k].pt();
        double pt_jerup = goodJets[k].pt();
        double pt_jerdn = goodJets[k].pt();
        const reco::GenJet * genJet = goodJets[k].genJet();
        if (genJet && genJet->pt()>15. && (abs(genJet->pt()/pt_jer-1)<0.5)) {
            double gen_pt = genJet->pt();
            pt_jer = max(0.0,gen_pt+factor*(goodJets[k].pt()-gen_pt));
            pt_jerup = max(0.0,gen_pt+factorup*(goodJets[k].pt()-gen_pt));
            pt_jerdn = max(0.0,gen_pt+factordn*(goodJets[k].pt()-gen_pt));
        }

        double jercorrection = pt_jer/goodJets[k].pt();
        double jercorrectionup = pt_jerup/goodJets[k].pt();
        double jercorrectiondn = pt_jerdn/goodJets[k].pt();
        */

        // FIXME: For now, no smearing
        double jercorrection = 1.0;
        double jercorrectionup = 1.0;
        double jercorrectiondn = 1.0;

        double jetPx_jer = jercorrection * goodJets[k].px();
        double jetPy_jer = jercorrection * goodJets[k].py();
        double jetPz_jer = jercorrection * goodJets[k].pz();
        double jetE_jer = sqrt(jetPx_jer*jetPx_jer + jetPy_jer*jetPy_jer + jetPz_jer*jetPz_jer + goodJets[k].mass()*goodJets[k].mass());          
        TLorentzVector *jet_jer = new TLorentzVector(jetPx_jer,jetPy_jer,jetPz_jer,jetE_jer);

        double jetPx_jerup = jercorrectionup * goodJets[k].px();
        double jetPy_jerup = jercorrectionup * goodJets[k].py();
        double jetPz_jerup = jercorrectionup * goodJets[k].pz();
        double jetE_jerup = sqrt(jetPx_jerup*jetPx_jerup + jetPy_jerup*jetPy_jerup + jetPz_jerup*jetPz_jerup + goodJets[k].mass()*goodJets[k].mass());          
        TLorentzVector *jet_jerup = new TLorentzVector(jetPx_jerup,jetPy_jerup,jetPz_jerup,jetE_jerup);

        double jetPx_jerdn = jercorrectiondn * goodJets[k].px();
        double jetPy_jerdn = jercorrectiondn * goodJets[k].py();
        double jetPz_jerdn = jercorrectiondn * goodJets[k].pz();
        double jetE_jerdn = sqrt(jetPx_jerdn*jetPx_jerdn + jetPy_jerdn*jetPy_jerdn + jetPz_jerdn*jetPz_jerdn + goodJets[k].mass()*goodJets[k].mass());          
        TLorentzVector *jet_jerdn = new TLorentzVector(jetPx_jerdn,jetPy_jerdn,jetPz_jerdn,jetE_jerdn);

        //std::cout<<"Jet nominal: "<<goodJets[k].pt()<<" JER corrected: "<<jet_jer->Pt()<<" JER up: "<<jet_jerup->Pt()<<" JER dn: "<<jet_jerdn->Pt()<<std::endl;

        jecunc->setJetPt(jet_jer->Pt());
        jecunc->setJetEta(goodJets[k].eta());
        double jecunc_up = 1.0+jecunc->getUncertainty(true);
        jecunc->setJetPt(jet_jer->Pt());
        jecunc->setJetEta(goodJets[k].eta());
        double jecunc_dn = 1.0-jecunc->getUncertainty(false);

        if (jet_jer->Pt() > 30.0 && fabs(goodJets[k].eta())<4.7) {
            if (isDeltaR_eta4p7) { 
                njets_pt30_eta4p7++;
                if (jet_jer->Pt() > pt_leadingjet_pt30_eta4p7) {
                    pt_leadingjet_pt30_eta4p7 = jet_jer->Pt();
                    absrapidity_leadingjet_pt30_eta4p7 = jet_jer->Rapidity(); //take abs later
                }
                jet_pt.push_back(jet_jer->Pt());
                jet_eta.push_back(jet_jer->Eta());
                jet_phi.push_back(jet_jer->Phi());
                jet_mass.push_back(jet_jer->M());
                jet_pumva.push_back(goodJets[k].userFloat("pileupJetId:fullDiscriminant"));
                jet_csvv2.push_back(goodJets[k].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                jet_isbtag.push_back(1 ? goodJets[k].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>BTagCut : 0);
                if (goodJets[k].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>BTagCut) nbjets_pt30_eta4p7++;
            }
        }
        
        // JER up
        if (jet_jerup->Pt() > 30.0 && fabs(jet_jerup->Eta())<4.7) {
            if (isDeltaR_eta4p7) {
                njets_pt30_eta4p7_jerup++;
                if (jet_jerup->Pt() > pt_leadingjet_pt30_eta4p7_jerup) {
                    pt_leadingjet_pt30_eta4p7_jerup = jet_jerup->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jerup = jet_jerup->Rapidity(); //take abs later
                }
                TLorentzVector jet_jerup(jetPx_jerup,jetPy_jerup,jetPz_jerup,jetE_jerup);
                jet_jerup_pt.push_back(jet_jerup.Pt());
                jet_jerup_eta.push_back(jet_jerup.Eta());
                jet_jerup_phi.push_back(jet_jerup.Phi());
                jet_jerup_mass.push_back(jet_jerup.M());

            }
        }

        // JER dn
        if (jet_jerdn->Pt() > 30.0 && fabs(jet_jerdn->Eta())<4.7) {
            if (isDeltaR_eta4p7) {
                njets_pt30_eta4p7_jerdn++;
                if (jet_jerdn->Pt() > pt_leadingjet_pt30_eta4p7_jerdn) {
                    pt_leadingjet_pt30_eta4p7_jerdn = jet_jerdn->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jerdn = jet_jerdn->Rapidity(); //take abs later
                }
                TLorentzVector jet_jerdn(jetPx_jerdn,jetPy_jerdn,jetPz_jerdn,jetE_jerdn);
                jet_jerdn_pt.push_back(jet_jerdn.Pt());
                jet_jerdn_eta.push_back(jet_jerdn.Eta());
                jet_jerdn_phi.push_back(jet_jerdn.Phi());
                jet_jerdn_mass.push_back(jet_jerdn.M());
            }
        }

        double jetPx_jesup = jecunc_up * jet_jer->Px();
        double jetPy_jesup = jecunc_up * jet_jer->Py();
        double jetPz_jesup = jecunc_up * jet_jer->Pz();
        double jetE_jesup = sqrt(jetPx_jesup*jetPx_jesup + jetPy_jesup*jetPy_jesup + jetPz_jesup*jetPz_jesup + jet_jer->M()*jet_jer->M());
        TLorentzVector *jet_jesup = new TLorentzVector(jetPx_jesup,jetPy_jesup,jetPz_jesup,jetE_jesup);
            
        if (jet_jesup->Pt() > 30.0 && fabs(jet_jesup->Eta())<4.7) {
            if (isDeltaR_eta4p7) {
                njets_pt30_eta4p7_jesup++;
                if (jet_jesup->Pt() > pt_leadingjet_pt30_eta4p7_jesup) {
                    pt_leadingjet_pt30_eta4p7_jesup = jet_jesup->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jesup = jet_jesup->Rapidity(); //take abs later
                }                             
                TLorentzVector jet_jesup(jetPx_jesup,jetPy_jesup,jetPz_jesup,jetE_jesup);
                jet_jesup_pt.push_back(jet_jesup.Pt());
                jet_jesup_eta.push_back(jet_jesup.Eta());
                jet_jesup_phi.push_back(jet_jesup.Phi());
                jet_jesup_mass.push_back(jet_jesup.M());
            }
        }

        double jetPx_jesdn = jecunc_dn * jet_jer->Px();
        double jetPy_jesdn = jecunc_dn * jet_jer->Py();
        double jetPz_jesdn = jecunc_dn * jet_jer->Pz();
        double jetE_jesdn = sqrt(jetPx_jesdn*jetPx_jesdn + jetPy_jesdn*jetPy_jesdn + jetPz_jesdn*jetPz_jesdn + jet_jer->M()*jet_jer->M());
        TLorentzVector *jet_jesdn = new TLorentzVector(jetPx_jesdn,jetPy_jesdn,jetPz_jesdn,jetE_jesdn);

        if (jet_jesdn->Pt() > 30.0 && fabs(jet_jesdn->Eta())<4.7) {
            if (isDeltaR_eta4p7) {
                njets_pt30_eta4p7_jesdn++;
                if (jet_jesdn->Pt() > pt_leadingjet_pt30_eta4p7_jesdn) {
                    pt_leadingjet_pt30_eta4p7_jesdn = jet_jesdn->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jesdn = jet_jesdn->Rapidity(); //take abs later
                }                
                TLorentzVector jet_jesdn(jetPx_jesdn,jetPy_jesdn,jetPz_jesdn,jetE_jesdn);
                jet_jesdn_pt.push_back(jet_jesdn.Pt());
                jet_jesdn_eta.push_back(jet_jesdn.Eta());
                jet_jesdn_phi.push_back(jet_jesdn.Phi());
                jet_jesdn_mass.push_back(jet_jesdn.M());
            }
        }
        
    } // loop over jets

    if(njets_pt30_eta4p7>1){
        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_mass[0]);
        jet2.SetPtEtaPhiM(jet_pt[1],jet_eta[1],jet_phi[1],jet_mass[1]);
        TLorentzVector Dijet;
        Dijet = jet1+jet2; 
        DijetMass = Dijet.M();
        DijetDEta = fabs(jet1.Eta()-jet2.Eta());
        // OLD MORIOND --- FisherDiscrim = 0.09407*fabs(VBFDeltaEta) + 4.1581e-4*VBFDiJetMass;
        DijetFisher = 0.18*fabs(DijetDEta) + 1.92e-4*DijetMass;
    }

    // Double loop over jets, for V-jet tagging
    for (int i=0; i<njets_pt30_eta4p7; i++) {
        for (int j=i+1; j<njets_pt30_eta4p7; j++) {
            if (i==j) continue;
            TLorentzVector ijet, jjet;
            ijet.SetPtEtaPhiM(jet_pt[i],jet_eta[i],jet_phi[i],jet_mass[i]);
            jjet.SetPtEtaPhiM(jet_pt[j],jet_eta[j],jet_phi[j],jet_mass[j]);
            if (ijet.Pt()<40.0 || abs(ijet.Eta())>2.4) continue;
            if (jjet.Pt()<40.0 || abs(jjet.Eta())>2.4) continue;
            TLorentzVector Dijet;
            Dijet = ijet+jjet;
            double mass = Dijet.M();
            if (mass > 60 && mass < 120) nvjets_pt40_eta2p4++;            
        }
    }

    // Higgs Variables
    if( RecoFourMuEvent ){ finalState = 1;}
    if( RecoFourEEvent  ){ finalState = 2;}
    if( RecoTwoETwoMuEvent ){ finalState = 3;}
    if( RecoTwoMuTwoEEvent ){ finalState = 4;}

    H_pt.push_back(HVec.Pt());
    H_eta.push_back(HVec.Eta());
    H_phi.push_back(HVec.Phi());
    H_mass.push_back(HVec.M());
    mass4l = HVec.M();
    mass4l_noFSR = HVecNoFSR.M();

    if(RecoFourMuEvent){mass4mu = HVec.M();}
    else{mass4mu = -1;}
    if(RecoFourEEvent){ mass4e = HVec.M();}
    else{ mass4e = -1;}
    if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){ mass2e2mu = HVec.M(); }
    else{ mass2e2mu = -1;}

    pT4l = HVec.Pt(); eta4l = HVec.Eta();
    rapidity4l = HVec.Rapidity(); phi4l = HVec.Phi();

    pTZ1 = Z1Vec.Pt(); pTZ2 = Z2Vec.Pt();
    massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M();

    if (njets_pt30_eta4p7>0) absdeltarapidity_hleadingjet_pt30_eta4p7 = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7);
    if (njets_pt30_eta4p7_jesup>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jesup = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jesup);
    if (njets_pt30_eta4p7_jesdn>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jesdn);
    if (njets_pt30_eta4p7_jerup>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jerup = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jerup);
    if (njets_pt30_eta4p7_jerdn>0) absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn = fabs(rapidity4l-absrapidity_leadingjet_pt30_eta4p7_jerdn);
    
    if (njets_pt30_eta4p7>0) absrapidity_leadingjet_pt30_eta4p7 = fabs(absrapidity_leadingjet_pt30_eta4p7);
    if (njets_pt30_eta4p7_jesup>0) absrapidity_leadingjet_pt30_eta4p7_jesup = fabs(absrapidity_leadingjet_pt30_eta4p7_jesup);
    if (njets_pt30_eta4p7_jesdn>0) absrapidity_leadingjet_pt30_eta4p7_jesdn = fabs(absrapidity_leadingjet_pt30_eta4p7_jesdn);
    if (njets_pt30_eta4p7_jerup>0) absrapidity_leadingjet_pt30_eta4p7_jerup = fabs(absrapidity_leadingjet_pt30_eta4p7_jerup);
    if (njets_pt30_eta4p7_jerdn>0) absrapidity_leadingjet_pt30_eta4p7_jerdn = fabs(absrapidity_leadingjet_pt30_eta4p7_jerdn);

    if (foundHiggsCandidate) {

        TLorentzVector Lep1FSR, Lep2FSR, Lep3FSR, Lep4FSR;
        Lep1FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[0]],lepFSR_eta[lep_Hindex[0]],lepFSR_phi[lep_Hindex[0]],lepFSR_mass[lep_Hindex[0]]);
        Lep2FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[1]],lepFSR_eta[lep_Hindex[1]],lepFSR_phi[lep_Hindex[1]],lepFSR_mass[lep_Hindex[1]]);
        Lep3FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[2]],lepFSR_eta[lep_Hindex[2]],lepFSR_phi[lep_Hindex[2]],lepFSR_mass[lep_Hindex[2]]);
        Lep4FSR.SetPtEtaPhiM(lepFSR_pt[lep_Hindex[3]],lepFSR_eta[lep_Hindex[3]],lepFSR_phi[lep_Hindex[3]],lepFSR_mass[lep_Hindex[3]]);
        pTL1FSR = Lep1FSR.Pt();
        pTL2FSR = Lep2FSR.Pt();
        pTL3FSR = Lep3FSR.Pt();
        pTL4FSR = Lep4FSR.Pt();

        TLorentzVector Lep1, Lep2, Lep3, Lep4;
        Lep1.SetPtEtaPhiM(lep_pt[lep_Hindex[0]],lep_eta[lep_Hindex[0]],lep_phi[lep_Hindex[0]],lep_mass[lep_Hindex[0]]);
        Lep2.SetPtEtaPhiM(lep_pt[lep_Hindex[1]],lep_eta[lep_Hindex[1]],lep_phi[lep_Hindex[1]],lep_mass[lep_Hindex[1]]);
        Lep3.SetPtEtaPhiM(lep_pt[lep_Hindex[2]],lep_eta[lep_Hindex[2]],lep_phi[lep_Hindex[2]],lep_mass[lep_Hindex[2]]);
        Lep4.SetPtEtaPhiM(lep_pt[lep_Hindex[3]],lep_eta[lep_Hindex[3]],lep_phi[lep_Hindex[3]],lep_mass[lep_Hindex[3]]);

        idL1 = lep_id[lep_Hindex[0]]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta();
        idL2 = lep_id[lep_Hindex[1]]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta();
        idL3 = lep_id[lep_Hindex[2]]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta();       
        idL4 = lep_id[lep_Hindex[3]]; pTL4 = Lep4.Pt(); etaL4 = Lep4.Eta();

    }

    for(unsigned int i = 0; i < lep_pt.size(); i++) {
        if ((int)i==lep_Hindex[0] || (int)i==lep_Hindex[1] || (int)i==lep_Hindex[2] || (int)i==lep_Hindex[3]) { nisoleptons++; }
        else {
            if (abs(lep_id[i])==11 && lep_tightId[i]==1 && lep_RelIsoNoFSR[i]<isoCutEl) { nisoleptons++; }
            else if (abs(lep_id[i])==13 && lep_tightId[i]==1 && lep_RelIsoNoFSR[i]<isoCutMu) { nisoleptons++; }
        }
    }

    // Event Categories
    if (nisoleptons==4 && njets_pt30_eta4p7>1 && nbjets_pt30_eta4p7<2 && DijetFisher>0.5) {EventCat=2;}
    else if (nisoleptons==4 && ( (nvjets_pt40_eta2p4>0 && pT4l>mass4l) || (njets_pt30_eta4p7==2 && nbjets_pt30_eta4p7==2) )) {EventCat=4;}
    else if (nisoleptons>4 && njets_pt30_eta4p7<3 && nbjets_pt30_eta4p7==0) { EventCat=3;}
    else if ((nisoleptons>4) || (njets_pt30_eta4p7>2 && nbjets_pt30_eta4p7>0)) {EventCat=5;}
    else if (njets_pt30_eta4p7>0) {EventCat=1;}
    else {EventCat=0;}


}




void UFHZZ4LAna::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                                 edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                                 edm::Handle<edm::View<reco::GenJet> > genJets)
{

    reco::GenParticleCollection::const_iterator genPart;
    int j = -1;
    int nGENLeptons=0;

    if (verbose) cout<<"begin looping on gen particles"<<endl;
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
        j++;

        if (abs(genPart->pdgId())==11  || abs(genPart->pdgId())==13 || abs(genPart->pdgId())==15) {

            if (!(genPart->status()==1 || abs(genPart->pdgId())==15)) continue;
            if (!(genAna.MotherID(&prunedgenParticles->at(j))==23 || abs(genAna.MotherID(&prunedgenParticles->at(j)))==24) ) continue;

            nGENLeptons++;
            if (verbose) cout<<"found a gen lepton: id "<<genPart->pdgId()<<" pt: "<<genPart->pt()<<" eta: "<<genPart->eta()<<endl;

            // Collect FSR photons
            TLorentzVector lep_dressed;            
            lep_dressed.SetPtEtaPhiE(genPart->pt(),genPart->eta(),genPart->phi(),genPart->energy());
            set<int> gen_fsrset;
            for(size_t k=0; k<packedgenParticles->size();k++){
                if( (*packedgenParticles)[k].status() != 1) continue; // stable particles only
                if( (*packedgenParticles)[k].pdgId() != 22) continue; // only photons
                double this_dR_lgamma = deltaR(genPart->eta(), genPart->phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
                bool idmatch=false;
                if ((*packedgenParticles)[k].mother(0)->pdgId()==genPart->pdgId() ) idmatch=true;
                const reco::Candidate * mother = (*packedgenParticles)[k].mother(0);
                for(size_t m=0;m<mother->numberOfMothers();m++) {
                    if ( (*packedgenParticles)[k].mother(m)->pdgId() == genPart->pdgId() ) idmatch=true;
                }
                if (!idmatch) continue;
                if(this_dR_lgamma<((abs(genPart->pdgId())==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {
                    gen_fsrset.insert(k);
                    TLorentzVector gamma;
                    gamma.SetPtEtaPhiE((*packedgenParticles)[k].pt(),(*packedgenParticles)[k].eta(),(*packedgenParticles)[k].phi(),(*packedgenParticles)[k].energy());
                    lep_dressed = lep_dressed+gamma;
                }
            } // Dressed leptons loop
            if (verbose) cout<<"gen lep pt "<<genPart->pt()<< " dressed pt: " << lep_dressed.Pt()<<endl;
  
            GENlep_id.push_back( genPart->pdgId() );
            GENlep_status.push_back(genPart->status());
            GENlep_pt.push_back( lep_dressed.Pt() );
            GENlep_eta.push_back( lep_dressed.Eta() );
            GENlep_phi.push_back( lep_dressed.Phi() );
            GENlep_mass.push_back( lep_dressed.M() );
            GENlep_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
            GENlep_MomMomId.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j)));

            TLorentzVector thisLep;
            thisLep.SetPtEtaPhiM(lep_dressed.Pt(),lep_dressed.Eta(),lep_dressed.Phi(),lep_dressed.M());
            // GEN iso calculation
            double this_GENiso=0.0;
            double this_GENneutraliso=0.0;
            double this_GENchargediso=0.0;
            if (verbose) cout<<"gen iso calculation"<<endl;
            for(size_t k=0; k<packedgenParticles->size();k++){
                if( (*packedgenParticles)[k].status() != 1 ) continue; // stable particles only         
                if (abs((*packedgenParticles)[k].pdgId())==12 || abs((*packedgenParticles)[k].pdgId())==14 || abs((*packedgenParticles)[k].pdgId())==16) continue; // exclude neutrinos
                if ((abs((*packedgenParticles)[k].pdgId())==11 || abs((*packedgenParticles)[k].pdgId())==13)) continue; // exclude leptons
                if (gen_fsrset.find(k)!=gen_fsrset.end()) continue; // exclude particles which were selected as fsr photons
                double this_dRvL = deltaR(thisLep.Eta(), thisLep.Phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
                if(this_dRvL<((abs(genPart->pdgId())==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {
                    //if (verbose) cout<<"adding to geniso id: "<<(*packedgenParticles)[k].pdgId()<<" status: "<<(*packedgenParticles)[k].status()<<" pt: "<<(*packedgenParticles)[k].pt()<<" dR: "<<this_dRvL<<endl;
                    this_GENiso = this_GENiso + (*packedgenParticles)[k].pt();
                    if ((*packedgenParticles)[k].charge()==0) this_GENneutraliso = this_GENneutraliso + (*packedgenParticles)[k].pt();
                    if ((*packedgenParticles)[k].charge()!=0) this_GENchargediso = this_GENchargediso + (*packedgenParticles)[k].pt();
                }
            } // GEN iso loop
            this_GENiso = this_GENiso/thisLep.Pt();
            if (verbose) cout<<"gen lep pt: "<<thisLep.Pt()<<" rel iso: "<<this_GENiso<<endl;
            GENlep_RelIso.push_back(this_GENiso);
            // END GEN iso calculation

        } // leptons
        
        if (genPart->pdgId()==25) {
            GENMH=genPart->mass();
            GENH_pt.push_back(genPart->pt());
            GENH_eta.push_back(genPart->eta());
            GENH_phi.push_back(genPart->phi());
            GENH_mass.push_back(genPart->mass());
        }

        
        if (genPart->pdgId()==23 && (genPart->status()>=20 && genPart->status()<30) ) {
            const reco::Candidate *Zdau0=genPart->daughter(0);
            int ZdauId = fabs(Zdau0->pdgId());
            if (fabs(Zdau0->pdgId())==23) {
                int ndau = genPart->numberOfDaughters();
                for (int d=0; d<ndau; d++) {
                    const reco::Candidate *Zdau=genPart->daughter(d);
                    if (verbose) cout<<"ZDau "<<d<<" id "<<fabs(Zdau->pdgId())<<endl;
                    if (fabs(Zdau->pdgId())<17) {
                        ZdauId = fabs(Zdau->pdgId());
                        break;
                    }
                }
            }
            if (verbose) cout<<"GENZ status "<<genPart->status()<<" MomId: "<<genAna.MotherID(&prunedgenParticles->at(j))<< "DauId: "<<ZdauId<<endl;

            if (Zdau0) GENZ_DaughtersId.push_back(ZdauId);           
            GENZ_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));                
            GENZ_pt.push_back(genPart->pt());
            GENZ_eta.push_back(genPart->eta());
            GENZ_phi.push_back(genPart->phi());
            GENZ_mass.push_back(genPart->mass());
        }

    }

    if (GENlep_pt.size()>=4) {
        
        unsigned int L1_nocuts=99; unsigned int L2_nocuts=99; unsigned int L3_nocuts=99; unsigned int L4_nocuts=99;
        bool passedFiducialSelectionNoCuts = mZ1_mZ2(L1_nocuts, L2_nocuts, L3_nocuts, L4_nocuts, false);
        if (passedFiducialSelectionNoCuts) {
            TLorentzVector Z1_1, Z1_2, Z2_1, Z2_2;
            Z1_1.SetPtEtaPhiM(GENlep_pt[L1_nocuts],GENlep_eta[L1_nocuts],GENlep_phi[L1_nocuts],GENlep_mass[L1_nocuts]);
            Z1_2.SetPtEtaPhiM(GENlep_pt[L2_nocuts],GENlep_eta[L2_nocuts],GENlep_phi[L2_nocuts],GENlep_mass[L2_nocuts]);
            Z2_1.SetPtEtaPhiM(GENlep_pt[L3_nocuts],GENlep_eta[L3_nocuts],GENlep_phi[L3_nocuts],GENlep_mass[L3_nocuts]);
            Z2_2.SetPtEtaPhiM(GENlep_pt[L4_nocuts],GENlep_eta[L4_nocuts],GENlep_phi[L4_nocuts],GENlep_mass[L4_nocuts]);         
            GENdPhiZZ = deltaPhi((Z1_1+Z1_2).Phi(),(Z2_1+Z2_2).Phi());            
            GENmassZZ = (Z1_1+Z1_2+Z2_1+Z2_2).M();            
            GENpTZZ = (Z1_1+Z1_2+Z2_1+Z2_2).Pt();            
            int genfs;
            if (abs(GENlep_id[L1_nocuts])==abs(GENlep_id[L3_nocuts])) genfs=1;
            else genfs=2;
            k_qqZZ_qcd = helper.kfactor_qqZZ_qcd_dPhi(GENdPhiZZ,genfs);
            if (verbose) cout<<"qcd kfactor qqZZ: "<<k_qqZZ_qcd<<endl;
        }
    }    
    
   
    /////// DO THE FIDUCIAL VOLUME CALCULATION //////////////
    if (verbose) cout<<"begin fiducial volume calculation"<<endl;
    passedFiducialSelection=false;
    int nFiducialLeptons = 0;
    int nFiducialPtLead=0;
    int nFiducialPtSublead=0;
    
    for (unsigned int i=0; i<GENlep_id.size(); ++i) {    

        TLorentzVector thisLep;
        thisLep.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
        
        if ( ( (abs(GENlep_id[i]) == 13 && thisLep.Pt() > 5.0 && abs(thisLep.Eta()) < 2.4)
               || (abs(GENlep_id[i]) == 11 && thisLep.Pt() > 7.0 && abs(thisLep.Eta()) < 2.5) )
             && GENlep_RelIso[i]<((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu) ) {
            nFiducialLeptons++;
            if (verbose) cout<<nFiducialLeptons<<" fiducial leptons, id;"<<GENlep_id[i]<<" pt: "<<thisLep.Pt()<<" eta: "<<thisLep.Eta()<<endl; 
            if (thisLep.Pt()>leadingPtCut) nFiducialPtLead++;
            if (thisLep.Pt()>subleadingPtCut) nFiducialPtSublead++;
        }                
    }
    
    if (nFiducialLeptons>=4 && nFiducialPtLead>=1 && nFiducialPtSublead>=2) {
                                
        // START FIDUCIAL EVENT TOPOLOGY CUTS
        unsigned int L1=99; unsigned int L2=99; unsigned int L3=99; unsigned int L4=99;
        GENmass4l = -1.0; GENmass4e = -1.0; GENmass4mu = -1.0; GENmass2e2mu = -1.0;
        GENmassZ1 = -1.0; GENmassZ2 = -1.0; GENpT4l = -1.0; GENeta4l = 999.; GENrapidity4l = 999.;
        passedFiducialSelection = mZ1_mZ2(L1, L2, L3, L4, true);      
                
        GENlep_Hindex[0] = L1; GENlep_Hindex[1] = L2; GENlep_Hindex[2] = L3; GENlep_Hindex[3] = L4;

        if (passedFiducialSelection) {
            
            TLorentzVector LS3_Z1_1, LS3_Z1_2, LS3_Z2_1, LS3_Z2_2;
            LS3_Z1_1.SetPtEtaPhiM(GENlep_pt[L1],GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1]);
            LS3_Z1_2.SetPtEtaPhiM(GENlep_pt[L2],GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2]);
            LS3_Z2_1.SetPtEtaPhiM(GENlep_pt[L3],GENlep_eta[L3],GENlep_phi[L3],GENlep_mass[L3]);
            LS3_Z2_2.SetPtEtaPhiM(GENlep_pt[L4],GENlep_eta[L4],GENlep_phi[L4],GENlep_mass[L4]);
            
            GENmass4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).M();
          
            if (abs(GENlep_id[L1])==11 && abs(GENlep_id[L3])==11) {GENmass4e = GENmass4l;};
            if (abs(GENlep_id[L1])==13 && abs(GENlep_id[L3])==13) {GENmass4mu = GENmass4l;};
            if ( (abs(GENlep_id[L1])==11 || abs(GENlep_id[L1])==13) &&
                 (abs(GENlep_id[L3])==11 || abs(GENlep_id[L3])==13) &&
                 (abs(GENlep_id[L1])!=abs(GENlep_id[L3]) ) ) {GENmass2e2mu = GENmass4l;};
            GENpT4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Pt();
            GENeta4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Eta();
            GENrapidity4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Rapidity();
            GENmassZ1 = (LS3_Z1_1+LS3_Z1_2).M();
            GENmassZ2 = (LS3_Z2_1+LS3_Z2_2).M();
                    
            int tmpIdL1,tmpIdL2,tmpIdL3,tmpIdL4;
            TLorentzVector GENL11P4, GENL12P4, GENL21P4, GENL22P4;
            if(GENlep_id[L1] < 0){ GENL11P4.SetPxPyPzE(LS3_Z1_1.Px(),LS3_Z1_1.Py(),LS3_Z1_1.Pz(),LS3_Z1_1.E()); tmpIdL1 = GENlep_id[L1];}
            else{ GENL11P4.SetPxPyPzE(LS3_Z1_2.Px(),LS3_Z1_2.Py(),LS3_Z1_2.Pz(),LS3_Z1_2.E()); tmpIdL1 = GENlep_id[L2];}
            if(GENlep_id[L2] > 0){ GENL12P4.SetPxPyPzE(LS3_Z1_2.Px(),LS3_Z1_2.Py(),LS3_Z1_2.Pz(),LS3_Z1_2.E()); tmpIdL2 = GENlep_id[L2];}
            else{ GENL12P4.SetPxPyPzE(LS3_Z1_1.Px(),LS3_Z1_1.Py(),LS3_Z1_1.Pz(),LS3_Z1_1.E()); tmpIdL2 = GENlep_id[L1];}
            if(GENlep_id[L3] < 0){ GENL21P4.SetPxPyPzE(LS3_Z2_1.Px(),LS3_Z2_1.Py(),LS3_Z2_1.Pz(),LS3_Z2_1.E()); tmpIdL3 = GENlep_id[L3];}
            else{ GENL21P4.SetPxPyPzE(LS3_Z2_2.Px(),LS3_Z2_2.Py(),LS3_Z2_2.Pz(),LS3_Z2_2.E()); tmpIdL3 = GENlep_id[L4];}
            if(GENlep_id[L4] > 0) { GENL22P4.SetPxPyPzE(LS3_Z2_2.Px(),LS3_Z2_2.Py(),LS3_Z2_2.Pz(),LS3_Z2_2.E()); tmpIdL4 = GENlep_id[L4];}
            else{ GENL22P4.SetPxPyPzE(LS3_Z2_1.Px(),LS3_Z2_1.Py(),LS3_Z2_1.Pz(),LS3_Z2_1.E()); tmpIdL4 = GENlep_id[L3];}
                    
            vector<TLorentzVector> GENP4s;
            vector<int> GENtmpIDs;
            GENP4s.push_back(GENL11P4); 
            GENP4s.push_back(GENL12P4); 
            GENP4s.push_back(GENL21P4); 
            GENP4s.push_back(GENL22P4); 
            GENtmpIDs.push_back(tmpIdL1);
            GENtmpIDs.push_back(tmpIdL2);
            GENtmpIDs.push_back(tmpIdL3);
            GENtmpIDs.push_back(tmpIdL4);
                    
            //mela::computeAngles(GENP4s[0], GENtmpIDs[0], GENP4s[1], GENtmpIDs[1], GENP4s[2], GENtmpIDs[2], GENP4s[3], GENtmpIDs[3], GENcosThetaStar,GENcosTheta1,GENcosTheta2,GENPhi,GENPhi1);
            
        }
        
        //cout<<"passedFiducialTopology? "<<passedFiducialTopology<<" GENmZ1Z2 = "<<GENmZ1Z2<<endl;
                
        bool passedMassOS = true; bool passedElMuDeltaR = true; bool passedDeltaR = true;            
        unsigned int N=GENlep_pt.size();
        for(unsigned int i = 0; i<N; i++) {
            for(unsigned int j = i+1; j<N; j++) {
                        
                // only consider the leptons from Z1 and Z2
                if (!(i==L1 || i==L2 || i==L3 || i==L4)) continue; 
                if (!(j==L1 || j==L2 || j==L3 || j==L4)) continue;
                        
                TLorentzVector li, lj;
                li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
                lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);

                TLorentzVector mll = li+lj;
                        
                if(GENlep_id[i]*GENlep_id[j]<0) {
                    if(mll.M()<=4) { passedMassOS = false; break; }
                }
                        
                if(abs(GENlep_id[i]) != abs(GENlep_id[j])) {
                    double deltaR = li.DeltaR(lj);
                    if(deltaR<=0.02) { passedElMuDeltaR = false; break; }
                }
                double deltaRll = li.DeltaR(lj);
                if(deltaRll<=0.02) { passedDeltaR = false; break; }
            }
        }
        
        if(passedMassOS==false || passedElMuDeltaR==false || passedDeltaR==false) passedFiducialSelection=false;
                
        if (passedFiducialSelection) {

            // DO GEN JETS
            if (verbose) cout<<"begin filling gen jets"<<endl;    
            edm::View<reco::GenJet>::const_iterator genjet;            
            for(genjet = genJets->begin(); genjet != genJets->end(); genjet++) {
                                        
                double pt = genjet->pt();  double eta = genjet->eta();
                if (pt<30.0 || abs(eta)>4.7) continue;

                bool inDR_pt30_eta4p7 = false;
                unsigned int N=GENlep_pt.size();
                for(unsigned int i = 0; i<N; i++) {
                    //if (GENlep_status[i]!=1) continue;
                    if (!(abs(GENlep_id[i])==11 || abs(GENlep_id[i])==13)) continue;
                    TLorentzVector genlep;
                    genlep.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
                    double dR = deltaR(genlep.Eta(), genlep.Phi(), genjet->eta(),genjet->phi());                        
                    if(dR<0.5) {
                        inDR_pt30_eta4p7=true;
                    }                                
                }
                
                if (verbose) cout<<"check overlap of gen jet with gen leptons"<<endl;
                // count number of gen jets which no gen leptons are inside its cone
                if (!inDR_pt30_eta4p7) { 
                    GENnjets_pt30_eta4p7++;
                    GENjet_pt.push_back(genjet->pt());
                    GENjet_eta.push_back(genjet->eta());
                    GENjet_phi.push_back(genjet->phi());
                    GENjet_mass.push_back(genjet->mass());
                    if (pt>GENpt_leadingjet_pt30_eta4p7) {
                        GENpt_leadingjet_pt30_eta4p7=pt;
                        GENabsrapidity_leadingjet_pt30_eta4p7=genjet->rapidity(); //take abs later
                    }
                }

            }// loop over gen jets

            if (GENnjets_pt30_eta4p7>0) GENabsdeltarapidity_hleadingjet_pt30_eta4p7 = fabs(GENrapidity4l-GENabsrapidity_leadingjet_pt30_eta4p7);
            if (GENnjets_pt30_eta4p7>0) GENabsrapidity_leadingjet_pt30_eta4p7 = fabs(GENabsrapidity_leadingjet_pt30_eta4p7);
            
        } //passedFiducialSelection

    } // 4 fiducial leptons

}

bool UFHZZ4LAna::mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts)
{

    double offshell = 999.0; bool findZ1 = false; bool passZ1 = false;

    L1 = 0; L2 = 0;

    unsigned int N = GENlep_pt.size();

    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<N; j++){


            if((GENlep_id[i]+GENlep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);

            if (verbose) cout<<"gen lep i id: "<<GENlep_id[i]<<" pt: "<<li.Pt()<<" lep j id: "<<GENlep_id[j]<<" pt: "<<lj.Pt()<<endl;

            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
                
                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;                
            }

            TLorentzVector mll = li+lj;
            if (verbose) cout<<"gen mass ij: "<<mll.M()<<endl;

            if(abs(mll.M()-Zmass)<offshell){
                double mZ1 = mll.M();
                L1 = i; L2 = j; findZ1 = true; offshell = abs(mZ1-Zmass);          
            }
        }    
    }

    TLorentzVector l1, l2;
    l1.SetPtEtaPhiM(GENlep_pt[L1],GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1]);
    l2.SetPtEtaPhiM(GENlep_pt[L2],GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2]);
    TLorentzVector ml1l2 = l1+l2;

    if(ml1l2.M()>40 && ml1l2.M()<120 && findZ1) passZ1 = true;
    if (!makeCuts) passZ1 = true;

    double pTL34 = 0.0; bool findZ2 = false;

    for(unsigned int i=0; i<N; i++){
        if(i==L1 || i==L2) continue; // can not be the lep from Z1
        for(unsigned int j=i+1; j<N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;            

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);
            TLorentzVector Z2 = li+lj;

            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
                
                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
            }

            if ( ( li.Pt()+lj.Pt() ) >=pTL34 ) { // choose high sum pT pair satisfy the following selection
                double mZ2 = Z2.M();
                if( (mZ2>12 && mZ2<120) || (!makeCuts) ) {
                    L3 = i; L4 = j; findZ2 = true; pTL34 = li.Pt()+lj.Pt();
                } else {
                    // still assign L3 and L4 to this pair if we don't have a passing Z2 yet
                    if (findZ2 == false) {L3 = i; L4 = j;}
                }
            }
            
        } // lj
    } // li

    if(passZ1 && findZ2) return true;
    else return false;
    
}



void UFHZZ4LAna::bookResolutionHistograms()
{

    using namespace edm;
    using namespace pat;
    using namespace std;

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

    using namespace edm;
    using namespace pat;
    using namespace std;

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
            if(mass<95 and mass>85) {
	    
                double pterr = mu->innerTrack()->ptError();
                double recpt = mu->pt();
                double eta = mu->eta();
                histContainer_["ptrelerrorForZ_mu"]->Fill(pterr/recpt,eventWeight);
                histContainer2D_["ptrelerrorForZ_mu_vs_eta"]->Fill(pterr/recpt, fabs(eta));
                if(abs(eta)>1.2){
                    histContainer_["ptrelerrorForZ_mu_endcap"]->Fill(pterr/recpt,eventWeight);
                }else{
                    histContainer_["ptrelerrorForZ_mu_barrel"]->Fill(pterr/recpt,eventWeight);
                }
            }
        }
    }
    //  For Muon Pt Resolution  end
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
UFHZZ4LAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(UFHZZ4LAna);
