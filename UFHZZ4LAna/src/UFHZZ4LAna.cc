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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
//#include "CMGTools/External/interface/PileupJetIdAlgo.h" // removed for miniAOD
//#include "CMGTools/External/interface/PileupJetIdentifier.h" // removed for miniAOD

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
//Iso
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LIsoEff.h"
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

//
// class declaration
//
using namespace MEMNames;

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
  
    void findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons,
                            std::vector< pat::PFParticle > &selectedFsrPhotons,const edm::Event& iEvent);

    void bookResolutionHistograms();
    void fillResolutionHistograms(edm::Handle<edm::View<pat::Muon> > muons);
  
    //MELA
    HZZ4LAngles angles;
    //Helper Class
    HZZ4LHelper helper;
    //Sip Class
    HZZ4LSipAna sipAna;
    //Iso Class
    HZZ4LIsoEff *isoEff;
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

    //Saved Events Trees
    TTree *passedEventsTree_All;
    void bookPassedEventTree(TString treeName, TTree *tree);
    void setTreeVariables( const edm::Event&, const edm::EventSetup&, 
                           std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, std::vector<pat::Jet> goodJets);
    void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                         edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                         edm::Handle<edm::View<reco::GenJet> > genJets);
    bool mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4);

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
    double genWeight, pileupWeight, dataMCWeight, eventWeight;

    // lepton variables
    TClonesArray *lep_p4;
    TClonesArray *lep_p4_FSR; //lep p4 also including any recovered FSR photons
    int lep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z2 sub
    double pTL1, pTL2, pTL3, pTL4;
    double etaL1, etaL2, etaL3, etaL4;
    int idL1, idL2, idL3, idL4;
    double pTL1FSR, pTL2FSR, pTL3FSR, pTL4FSR;
    vector<int> lep_genindex; //position of lepton in GENlep_p4 (if gen matched, -1 if not gen matched)
    vector<int> lep_id;
    vector<double> lep_mva;
    vector<int>   lep_tightId;
    vector<double> lep_Sip;
    vector<double> lep_IP;
    vector<double> lep_isoNH;
    vector<double> lep_isoCH;
    vector<double> lep_isoPhot;
    vector<double> lep_isoPU;
    vector<double> lep_isoPUcorr;
    vector<double> lep_RelIso;
    vector<int> lep_missingHits;
    int nisoleptons;
    double muRho, elRho;

    // Higgs candidate variables
    TClonesArray *H_p4;
    TClonesArray *H_p4_noFSR;
    double mass4l, mass4l_noFSR, mass4e, mass4mu, mass2e2mu, pT4l, eta4l, phi4l, rapidity4l;
    float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;

    // Z candidate variables
    TClonesArray *Z_p4;
    TClonesArray *Z_p4_noFSR;
    int Z_Hindex[2]; // position of Z1 and Z2 in Z_p4
    double massZ1, massZ2, pTZ1, pTZ2;

    // MET
    TClonesArray *met_p4;
    double met;

    // Jets
    TClonesArray *jet_p4;
    vector<double> jet_pumva, jet_csvv2;
    vector<int> jet_isbtag;
    TClonesArray *jet_p4_jesup;
    TClonesArray *jet_p4_jesdn;
    TClonesArray *jet_p4_jerup;
    TClonesArray *jet_p4_jerdn;
    
    int njets_pt30_eta4p7;
    int njets_pt30_eta4p7_jesup;
    int njets_pt30_eta4p7_jesdn;
    int njets_pt30_eta4p7_jerup;
    int njets_pt30_eta4p7_jerdn;

    int nbjets_pt30_eta4p7, nvjets_pt40_eta2p4;

    double pt_leadingjet_pt30_eta4p7;
    double pt_leadingjet_pt30_eta4p7_jesup;
    double pt_leadingjet_pt30_eta4p7_jesdn;
    double pt_leadingjet_pt30_eta4p7_jerup;
    double pt_leadingjet_pt30_eta4p7_jerdn;

    double absrapidity_leadingjet_pt30_eta4p7;
    double absrapidity_leadingjet_pt30_eta4p7_jesup;
    double absrapidity_leadingjet_pt30_eta4p7_jesdn;
    double absrapidity_leadingjet_pt30_eta4p7_jerup;
    double absrapidity_leadingjet_pt30_eta4p7_jerdn;

    double absdeltarapidity_hleadingjet_pt30_eta4p7;
    double absdeltarapidity_hleadingjet_pt30_eta4p7_jesup;
    double absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn;
    double absdeltarapidity_hleadingjet_pt30_eta4p7_jerup;
    double absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn;

    double DijetMass, DijetDEta;
    double DijetFisher;

    // FSR Photons
    TClonesArray *phofsr_p4;
    int nFSRPhotons;
    bool FSR_Z1, FSR_Z2;
    vector<int> phofsr_lepindex;

    // Z4l? FIXME
    double theta12, theta13, theta14;  
    double minM3l, Z4lmaxP, minDeltR, m3l_soft;
    double minMass2Lep, maxMass2Lep;
    double thetaPhoton, thetaPhotonZ;

    //Resolution
    double massErrorUCSD, massErrorUCSDCorr, massErrorUF, massErrorUFCorr, massErrorUFADCorr;
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
    TClonesArray *GENlep_p4;
    vector<int> GENlep_id; 
    vector<int> GENlep_status; 
    vector<int> GENlep_MomId; 
    vector<int> GENlep_MomMomId;
    int GENlep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
    vector<double> GENlep_isoCH; 
    vector<double> GENlep_isoNH; 
    vector<double> GENlep_isoPhot;
    vector<double> GENlep_RelIso; 


    // Higgs candidate variables (calculated using selected gen leptons)
    TClonesArray *GENH_p4;
    double GENmass4l, GENmass4e, GENmass4mu, GENmass2e2mu, GENpT4l, GENeta4l, GENrapidity4l;
    double GENMH; //mass directly from gen particle with id==25
    float GENcosTheta1, GENcosTheta2, GENcosThetaStar, GENPhi, GENPhi1;

    // Z candidate variables
    TClonesArray *GENZ_p4;
    vector<int> GENZ_DaughtersId; vector<int> GENZ_MomId;
    double  GENmassZ1, GENmassZ2, GENpTZ1, GENpTZ2;

    // Higgs variables directly from GEN particle
    double GENHmass;

    // Jets
    TClonesArray *GENjet_p4;
    int GENnjets_pt30_eta4p7;
    double GENpt_leadingjet_pt30_eta4p7;
    double GENabsrapidity_leadingjet_pt30_eta4p7;
    double GENabsdeltarapidity_hleadingjet_pt30_eta4p7;

    // MEM
    MEMs*  combinedMEM;

    double me_0plus_JHU, me_qqZZ_MCFM, p0plus_m4l, bkg_m4l;
    double D_bkg_kin, D_bkg;   

    double bkg_VAMCFM;
    double p0minus_VAJHU, Dgg10_VAMCFM;
    double phjj_VAJHU, pvbf_VAJHU;

    double D_g4;
    double Djet_VAJHU;

    // ----------member data -------------------

 
    // mass err
    std::map<TString,TH1F*> hContainer_;
    std::map<TString,TH2F*> hContainer2D_;
    std::map<TString,TH3F*> hContainer3D_;

    // Global Variables but not stored in the tree
    vector<double> lep_pt;
    vector<int> lep_ptid;
    vector<int> lep_ptindex;
    vector<pat::Muon> recoMuons;
    vector<pat::Electron> recoElectrons;
    vector<pat::PFParticle> fsrPhotons; 
    vector<int> fsrPhotons_lepindex; 
    vector<double> fsrPhotons_deltaR;
    vector<double> newfsrPhotons_dR;
    vector<double> newfsrPhotons_pt;
    vector<double> newfsrPhotons_iso;
    TLorentzVector HVec, HVecNoFSR, Z1Vec, Z2Vec;
    bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;
    bool foundHiggsCandidate;
    bool firstEntry;

    // hist container
    std::map<std::string,TH1F*> histContainer_;
    std::map<std::string,TH2F*> histContainer2D_;
 
    //Input tags
    edm::InputTag photonSrc_;
    edm::InputTag elecSrc_;
    edm::InputTag muonSrc_;
    edm::InputTag jetSrc_;
    edm::InputTag metSrc_;
    edm::InputTag triggerSrc_;
    edm::InputTag vertexSrc_;
    edm::InputTag muRhoSrc_;
    edm::InputTag elRhoSrc_;
    edm::InputTag pileupSrc_;

    // Configuration
    const double Zmass;
    double mZ1Low, mZ2Low;
    double mZ1High, mZ2High;
    double m4lLowCut;
    double jetpt_cut, jeteta_cut;
    std::string elecID;
    bool isMC, isSignal;
    double mH;
    double crossSection, filterEff, Lumi;
    bool weightEvents;
    double isoCutEl, isoCutMu, sip3dCut;
    double leadingPtCut, subleadingPtCut;
    double genIsoCut;
    double genFSRDrCut;
    double _elecPtCut, _muPtCut;
    double BTagCut;
    bool reweightForPU;
    bool interactiveRun;
    std::string PUVersion;
    bool doFsrRecovery, doPUJetID;
    bool doSUSYSelection;
    bool bStudyResolution;
    bool bStudyDiLeptonResolution;
    bool bStudyFourLeptonResolution;
    std::vector<std::string> triggerList;
    bool verbose;
    // register to the TFileService
    edm::Service<TFileService> fs;

    // Counters
    double nEventsTotal;
    double sumWeightsTotal;
};


UFHZZ4LAna::UFHZZ4LAna(const edm::ParameterSet& iConfig) :
    histContainer_(),
    histContainer2D_(),
    photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
    elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
    muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
    jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" )),
    metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc" )),
    triggerSrc_(iConfig.getUntrackedParameter<edm::InputTag>("triggerSrc")),
    vertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")),
    muRhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc")),
    elRhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc")),
    pileupSrc_(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc")),
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
    isoCutEl(iConfig.getUntrackedParameter<double>("isoCutEl",0.5)),
    isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",0.4)),
    sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
    leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
    subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
    genIsoCut(iConfig.getUntrackedParameter<double>("genIsoCut",0.5)), 
    genFSRDrCut(iConfig.getUntrackedParameter<double>("genFSRDrCut",0.4)), 
    _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",7)),
    _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",5)),
    BTagCut(iConfig.getUntrackedParameter<double>("BTagCut",0.814)),
    reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
    interactiveRun(iConfig.getUntrackedParameter<bool>("interactiveRun",false)),
    PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Fall15_74X")),
    doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",true)),
    doPUJetID(iConfig.getUntrackedParameter<bool>("doPUJetID",true)),
    doSUSYSelection(iConfig.getUntrackedParameter<bool>("doSUSYSelection",false)),
    bStudyResolution(iConfig.getUntrackedParameter<bool>("bStudyResolution",false)),
    bStudyDiLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyDiLeptonResolution",false)),
    bStudyFourLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyFourLeptonResolution",false)),
    triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList")),
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
   
    isoEff = new HZZ4LIsoEff("");
    
    if(bStudyResolution) {
        PerLepReso = new HZZ4LPerLepResolution();
        if(bStudyDiLeptonResolution) {  DiLepReso = new HZZ4LDiLepResolution(); }
        if(bStudyFourLeptonResolution)  { FourLepReso = new HZZ4LResolution(); }
    }

    jecunc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("UFHZZAnalysisRun2/UFHZZ4LAna/hists/Summer13_V5_DATA_UncertaintySources_AK5PF.txt","Total")));

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
    using namespace MEMNames;

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
    iEvent.getByLabel(triggerSrc_,trigger);
    const edm::TriggerNames trigNames = iEvent.triggerNames(*trigger);

    // vertex collection
    edm::Handle<reco::VertexCollection> vertex;
    iEvent.getByLabel(vertexSrc_,vertex);
    const reco::Vertex *PV = 0;
    if(!vertex->empty() && vertex->size() > 0) PV = &(vertex->at(0));

    // photon collection 
    edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByLabel(photonSrc_,photons);
  
    // electron collection
    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByLabel(elecSrc_,electrons);
    if (verbose) cout<<electrons->size()<<" total electrons in the collection"<<endl;

    // muon collection
    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByLabel(muonSrc_,muons);
    if (verbose) cout<<muons->size()<<" total muons in the collection"<<endl;

    // met collection 
    edm::Handle<edm::View<pat::MET> > mets;
    iEvent.getByLabel(metSrc_,mets);
    
    // beamspot collection
    edm::Handle<reco::BeamSpot> beamspot;
    iEvent.getByLabel("offlineBeamSpot",beamspot);
  
    // Rho Correction
    edm::Handle<double> eventRhoMu;
    iEvent.getByLabel(muRhoSrc_,eventRhoMu);
    muRho = *eventRhoMu;

    edm::Handle<double> eventRhoE;
    iEvent.getByLabel(elRhoSrc_,eventRhoE);
    elRho = *eventRhoE;

    // Conversions
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    iEvent.getByLabel("reducedEgamma","reducedConversions", theConversions);
 
    // Beam Spot
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByLabel("offlineBeamSpot",beamSpot);
    const reco::BeamSpot BS = *beamSpot;

    // Particle Flow Cands
    edm::Handle<pat::PackedCandidateCollection> pfCands;
    iEvent.getByLabel("packedPFCandidates",pfCands);

    // Photons
    edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;
    iEvent.getByLabel("boostedFsrPhotons",photonsForFsr);
  
    // Jets
    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel(jetSrc_,jets);
  
    // GEN collections
    edm::Handle<reco::GenParticleCollection> prunedgenParticles;
    iEvent.getByLabel("prunedGenParticles", prunedgenParticles);

    edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles;
    iEvent.getByLabel("packedGenParticles", packedgenParticles);
    
    edm::Handle<edm::View<reco::GenJet> > genJets;
    iEvent.getByLabel("slimmedGenJets", genJets);
    
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByLabel("generator",genEventInfo);
    
    // ============ Initialize Variables ============= //

    // Event Variables
    if (verbose) {cout<<"clear variables"<<endl;}
    nVtx = -1.0;
    finalState = -1;
    triggersPassed="";
    passedTrig=false; passedFullSelection=false; passedZ4lSelection=false; passedQCDcut=false;

    // Event Weights
    genWeight=1.0; pileupWeight=1.0; dataMCWeight=1.0; eventWeight=1.0;

    //lepton variables
    if (lep_p4->GetLast()!=-1) lep_p4->Clear(); 
    if (lep_p4_FSR->GetLast()!=-1) lep_p4_FSR->Clear();
    for (int i=0; i<4; ++i) {lep_Hindex[i]=-1;}
    pTL1=-1.0; pTL2=-1.0; pTL3=-1.0; pTL4=-1.0;
    etaL1=9999.0; etaL2=9999.0; etaL3=9999.0; etaL4=9999.0;
    idL1=9999; idL2=9999; idL3=9999; idL4=9999;
    pTL1FSR=-1.0; pTL2FSR=-1.0; pTL3FSR=-1.0; pTL4FSR=-1.0;
    lep_genindex.clear(); lep_id.clear(); 
    lep_mva.clear(); lep_tightId.clear();
    lep_Sip.clear(); lep_IP.clear(); 
    lep_isoNH.clear(); lep_isoCH.clear(); lep_isoPhot.clear(); lep_isoPU.clear(); lep_isoPUcorr.clear(); lep_RelIso.clear();
    lep_missingHits.clear();
    nisoleptons=0;

    // Higgs candidate variables
    if (H_p4->GetLast()!=-1) H_p4->Clear();
    if (H_p4_noFSR->GetLast()!=-1) H_p4_noFSR->Clear();
    mass4l=-1.0; mass4l_noFSR=-1.0; mass4e=-1.0; mass4mu=-1.0; mass2e2mu=-1.0; pT4l=-1.0; eta4l=9999.0; phi4l=9999.0; rapidity4l=9999.0;
    cosTheta1=9999.0; cosTheta2=9999.0; cosThetaStar=9999.0; Phi=9999.0; Phi1=9999.0;

    // Z candidate variables
    if (Z_p4->GetLast()!=-1) Z_p4->Clear();
    for (int i=0; i<2; ++i) {Z_Hindex[i]=-1;}
    massZ1=-1.0; massZ2=-1.0; pTZ1=-1.0; pTZ2=-1.0;

    // MET
    if (met_p4->GetLast()!=-1) met_p4->Clear();
    met=-1.0;

    // Jets
    if (jet_p4->GetLast()!=-1) jet_p4->Clear();
    if (jet_p4_jesup->GetLast()!=-1) jet_p4_jesup->Clear();
    if (jet_p4_jesdn->GetLast()!=-1) jet_p4_jesdn->Clear();
    if (jet_p4_jerup->GetLast()!=-1) jet_p4_jerup->Clear();
    if (jet_p4_jerdn->GetLast()!=-1) jet_p4_jerdn->Clear();
    
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
    if (phofsr_p4->GetLast()!=-1) phofsr_p4->Clear();
    nFSRPhotons=0;
    FSR_Z1=false; FSR_Z2=false;
    phofsr_lepindex.clear();

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
    if (GENlep_p4->GetLast()!=-1) GENlep_p4->Clear();
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
    if (GENH_p4->GetLast()!=-1) GENH_p4->Clear();
    GENmass4l=-1.0; GENmassZ1=-1.0; GENmassZ2=-1.0; GENpT4l=-1.0; GENeta4l=9999.0; GENrapidity4l=9999.0; GENMH=-1.0;
    GENcosTheta1=9999.0; GENcosTheta2=9999.0; GENcosThetaStar=9999.0; GENPhi=9999.0; GENPhi1=9999.0;

    // Z candidate variables
    if (GENZ_p4->GetLast()!=-1) GENZ_p4->Clear();
    GENZ_DaughtersId.clear(); GENZ_MomId.clear();
    GENmassZ1=-1.0; GENmassZ2=-1.0; GENpTZ1=-1.0; GENpTZ2=-1.0;

    // Higgs variables directly from GEN particle
    GENHmass=-1.0;

    // Jets
    if (GENjet_p4->GetLast()!=-1) GENjet_p4->Clear();
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
    lep_pt.clear(); lep_ptid.clear(); lep_ptindex.clear();
    recoMuons.clear(); recoElectrons.clear(); 
    fsrPhotons.clear(); fsrPhotons_lepindex.clear(); fsrPhotons_deltaR.clear();
    newfsrPhotons_dR.clear(); newfsrPhotons_pt.clear(); newfsrPhotons_iso.clear();
    HVec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    HVecNoFSR.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    Z1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    Z2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
    RecoFourMuEvent = false; RecoFourEEvent = false;
    RecoTwoETwoMuEvent = false; RecoTwoMuTwoEEvent = false;
    foundHiggsCandidate = false;

    // ====================== Do Analysis ======================== //

    if (verbose) cout<<"start pileup reweighting"<<endl;
    // PU information
    if(isMC && reweightForPU) {        
        edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByLabel(pileupSrc_, PupInfo);      

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
        double tmpWeight = genEventInfo->weight();
        genWeight = (tmpWeight > 0 ? 1.0 : -1.0);
        if (verbose) cout<<"setting gen variables"<<endl;       
        setGENVariables(prunedgenParticles,packedgenParticles,genJets); 
        if (verbose) { cout<<"finshed setting gen variables"<<endl;  }
    }
    sumWeightsTotal += genWeight;

    eventWeight = pileupWeight*genWeight;

    //Check for duplicate events in data
    if (verbose) cout<<"checking duplicates"<<endl;       
    bool notDuplicateEvent = true;
    /*
    std::vector<ULong64_t> runVec, lumiVec, eventVec;
    if(!isMC) { 
        for (unsigned int n = 0; n < runVec.size(); n++)  {
            if(Run == runVec[n] && LumiSect == lumiVec[n] && Event == eventVec[n]){notDuplicateEvent = false;}
        }      
        if (notDuplicateEvent) {
            runVec.push_back(Run);
            lumiVec.push_back(LumiSect);
            eventVec.push_back(Event);
        }      
    }
    if (verbose) { cout<<"finished checking duplicates"<<endl;}
    */

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

    bool isFake = (PV->chi2()==0 && PV->ndof()==0);

    if( notDuplicateEvent && !vertex->empty() && !isFake
        && PV->ndof()>4 && PV->position().Rho()<=2.0
        && fabs(PV->position().Z())<=24.0 ) {
    
        //N Vertex 
        if (verbose) {cout<<"fill nvtx histogram"<<endl;}
        nVtx = vertex->size();
        histContainer_["NVTX"]->Fill(nVtx);
        histContainer_["NVTX_RW"]->Fill(nVtx,pileupWeight);

        //MET
        if (verbose) {cout<<"get met value"<<endl;}
        if (!mets->empty()) {
            met = (*mets)[0].et();
            new ( (*met_p4)[0] ) TLorentzVector((*mets)[0].px(),(*mets)[0].py(),0.0,(*mets)[0].et());
        }

        if (verbose) cout<<"start lepton analysis"<<endl;           
        vector<pat::Muon> AllMuons;
        vector<pat::Electron> AllElectrons;  
        AllMuons = helper.goodLooseMuons2012(muons,_muPtCut);
        AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut);

        if (doSUSYSelection) {
            helper.cleanOverlappingLeptons_SUSY(AllMuons,AllElectrons,PV);
            if (verbose) cout<<AllMuons.size()<<" loose muons "<<AllElectrons.size()<<" loose elctrons"<<endl;
            recoMuons = helper.goodMuons2015_SUSY(AllMuons,_muPtCut,PV,pfCands,jets,sip3dCut,muRho,0.22,true);
            recoElectrons = helper.goodElectrons2015_SUSY(AllElectrons,_elecPtCut,elecID,PV,BS,pfCands,theConversions,jets,sip3dCut,elRho,0.14,true);
        } else { 
            helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);
            if (verbose) cout<<AllMuons.size()<<" loose muons "<<AllElectrons.size()<<" loose elctrons"<<endl;
            recoMuons = helper.goodMuons2015_noIso_noPf(AllMuons,_muPtCut,PV,sip3dCut);
            recoElectrons = helper.goodElectrons2015_noIso_noBdt(AllElectrons,_elecPtCut,elecID,PV,iEvent,sip3dCut);
        }


        //sort electrons and muons by pt
        if (verbose) cout<<recoMuons.size()<<" good muons and "<<recoElectrons.size()<<" good electrons to be sorted"<<endl;
        if (verbose) cout<<"start pt-sorting leptons"<<endl;
        if (verbose) cout<<"adding muons to sorted list"<<endl;           
        for(unsigned int i = 0; i < recoMuons.size(); i++) {
            if (lep_pt.size()==0 || recoMuons[i].pt()<lep_pt[lep_pt.size()-1]) {
                lep_pt.push_back(recoMuons[i].pt());
                lep_ptid.push_back(recoMuons[i].pdgId());
                lep_ptindex.push_back(i);
                continue;
            }
            for (unsigned int j=0; j<lep_pt.size(); j++) {
                if (recoMuons[i].pt()>lep_pt[j]) {
                    lep_pt.insert(lep_pt.begin()+j,recoMuons[i].pt());
                    lep_ptid.insert(lep_ptid.begin()+j,recoMuons[i].pdgId());
                    lep_ptindex.insert(lep_ptindex.begin()+j,i);
                    break;
                }
            }
        }
        if (verbose) cout<<"adding electrons to sorted list"<<endl;           
        for(unsigned int i = 0; i < recoElectrons.size(); i++) {
            if (lep_pt.size()==0 || recoElectrons[i].pt()<lep_pt[lep_pt.size()-1]) {
                lep_pt.push_back(recoElectrons[i].pt());
                lep_ptid.push_back(recoElectrons[i].pdgId());
                lep_ptindex.push_back(i);
                continue;
            }
            for (unsigned int j=0; j<lep_pt.size(); j++) {
                if (recoElectrons[i].pt()>lep_pt[j]) {
                    lep_pt.insert(lep_pt.begin()+j,recoElectrons[i].pt());
                    lep_ptid.insert(lep_ptid.begin()+j,recoElectrons[i].pdgId());
                    lep_ptindex.insert(lep_ptindex.begin()+j,i);
                    break;
                }
            }
        }

        for(unsigned int i = 0; i < lep_pt.size(); i++) {
            
            if (verbose) cout<<"sorted lepton "<<i<<" pt "<<lep_pt[i]<<" id "<<lep_ptid[i]<<" index "<<lep_ptindex[i]<<endl;

            if (abs(lep_ptid[i])==11) {
                new ( (*lep_p4)[i] ) TLorentzVector(recoElectrons[lep_ptindex[i]].px(),recoElectrons[lep_ptindex[i]].py(),recoElectrons[lep_ptindex[i]].pz(),recoElectrons[lep_ptindex[i]].energy());
                new ( (*lep_p4_FSR)[i] ) TLorentzVector(recoElectrons[lep_ptindex[i]].px(),recoElectrons[lep_ptindex[i]].py(),recoElectrons[lep_ptindex[i]].pz(),recoElectrons[lep_ptindex[i]].energy());
                lep_id.push_back(recoElectrons[lep_ptindex[i]].pdgId());
                lep_RelIso.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
                lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].chargedHadronIso());
                lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].neutralHadronIso());
                lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].photonIso());
                lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].puChargedHadronIso());
                lep_isoPUcorr.push_back(helper.getPUIso(recoElectrons[lep_ptindex[i]],elRho));
                lep_Sip.push_back(helper.getSIP3D(recoElectrons[lep_ptindex[i]]));           
                lep_mva.push_back(recoElectrons[lep_ptindex[i]].electronID(elecID)); 
                lep_tightId.push_back(helper.passTight_BDT_Id(recoElectrons[lep_ptindex[i]],elecID));           
                lep_genindex.push_back(-1.0);
            }
            if (abs(lep_ptid[i])==13) {            
                new ( (*lep_p4)[i] ) TLorentzVector(recoMuons[lep_ptindex[i]].px(),recoMuons[lep_ptindex[i]].py(),recoMuons[lep_ptindex[i]].pz(),recoMuons[lep_ptindex[i]].energy());
                new ( (*lep_p4_FSR)[i] ) TLorentzVector(recoMuons[lep_ptindex[i]].px(),recoMuons[lep_ptindex[i]].py(),recoMuons[lep_ptindex[i]].pz(),recoMuons[lep_ptindex[i]].energy());
                lep_id.push_back(recoMuons[lep_ptindex[i]].pdgId());
                lep_RelIso.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
                lep_isoCH.push_back(recoMuons[lep_ptindex[i]].chargedHadronIso());
                lep_isoNH.push_back(recoMuons[lep_ptindex[i]].neutralHadronIso());
                lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].photonIso());
                lep_isoPU.push_back(recoMuons[lep_ptindex[i]].puChargedHadronIso());
                lep_isoPUcorr.push_back(helper.getPUIso(recoMuons[lep_ptindex[i]],muRho));
                lep_Sip.push_back(helper.getSIP3D(recoMuons[lep_ptindex[i]]));            
                lep_mva.push_back(recoMuons[lep_ptindex[i]].isPFMuon());  
                lep_tightId.push_back(recoMuons[lep_ptindex[i]].isPFMuon());         
                lep_genindex.push_back(-1.0);
            }
            if (verbose) {cout<<" RelIso: "<<lep_RelIso[i]<<" isoCH: "<<lep_isoCH[i]<<" isoNH: "<<lep_isoNH[i]
                              <<" isoPhot: "<<lep_isoPhot[i]<<" isoPUcorr: "<<lep_isoPUcorr[i]<<" Sip: "<<lep_Sip[i]<<endl;}
        }

        // GEN matching
        if(isMC) {
            if (verbose) cout<<"begin gen matching"<<endl;
            // for each reco lepton find the nearest gen lepton with same ID
            for(unsigned int i = 0; i < lep_pt.size(); i++) {

                double minDr=9999.0;

                TLorentzVector *reco, *gen;
                reco = (TLorentzVector*) lep_p4->At(i);

                for (unsigned int j = 0; j < GENlep_id.size(); j++) {

                    if (GENlep_id[j]!=lep_id[i]) continue;
                    gen = (TLorentzVector*) GENlep_p4->At(j);
                    double thisDr = deltaR(reco->Eta(),reco->Phi(),gen->Eta(),gen->Phi());

                    if (thisDr<minDr && thisDr<0.5) {
                        lep_genindex[i]=j;
                        minDr=thisDr;
                    }

                } // all gen leptons
 
            } // all reco leptons
            if (verbose) {cout<<"finished gen matching"<<endl;}
        } //isMC

        if( (recoMuons.size() + recoElectrons.size()) >= 2 ) {
            if (verbose) cout<<"found two leptons"<<endl;

            // Mass Resolution Study
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
            
            // FSR Photons
            if(doFsrRecovery) {
                
                if (verbose) cout<<"checking "<<photonsForFsr->size()<<" fsr photon candidates"<<endl;
                
                for(edm::View<pat::PFParticle>::const_iterator phot=photonsForFsr->begin(); phot!=photonsForFsr->end(); ++phot) {
                    
                    bool matched = false; double minDeltaRFSR = 999.0; int lepindex=-1;
                    
                    if (fabs(phot->eta()) > 2.4) continue;
                    if (phot->pt()<2.0) continue;
                    
                    unsigned int Nleptons = lep_p4->GetLast()+1;
                    for (unsigned int i=0; i<Nleptons; i++) {
                        
                        if (matched) continue;
                        
                        TLorentzVector *thisLep;
                        thisLep = (TLorentzVector*) lep_p4->At(i);
                        
                        double fsrDr = deltaR(thisLep->Eta(), thisLep->Phi(), phot->eta(), phot->phi());
                        
                        if ( abs(lep_id[(int)i])==11) {
                            
                            double fsrDr_eSC = deltaR(recoElectrons[lep_ptindex[i]].superCluster()->eta(), recoElectrons[lep_ptindex[i]].superCluster()->phi(), phot->eta(), phot->phi());
                            double fsrDeltaPhi_eSC = fabs(deltaPhi(phot->phi(),recoElectrons[lep_ptindex[i]].superCluster()->phi()));
                            double fsrDeltaEta_eSC = fabs(phot->eta()-recoElectrons[lep_ptindex[i]].superCluster()->eta());
                            
                            if ( fsrDr_eSC<0.15 || (fsrDeltaPhi_eSC<2.0 && fsrDeltaEta_eSC<0.05) ) { 
                                matched=true;
                                continue;
                            }
                            
                        }
                        
                        if( fsrDr < minDeltaRFSR ) {
                            minDeltaRFSR = fsrDr;
                            lepindex = (int)i;
                        }
                        
                    }
                    
                    if (matched) continue;
                    
                    double photoniso = (phot->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02")+phot->userFloat("fsrPhotonPFIsoNHadPhoton03"))/phot->pt();
                    
                    if ( verbose) cout<<"fsr photon cand, pt: "<<phot->pt()<<" eta: "<<phot->eta()<<" phi: "<<phot->phi()
                                      <<" isoCHPUNoPU: "<<phot->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02")
                                      <<" isoNHPhoton: "<<phot->userFloat("fsrPhotonPFIsoNHadPhoton03")
                                      <<" photoniso: "<<photoniso<<" mindDeltaRFSR: "<<minDeltaRFSR<<endl;
                    
                    if (minDeltaRFSR < 0.5) {
                        newfsrPhotons_dR.push_back(minDeltaRFSR);
                        newfsrPhotons_iso.push_back(photoniso);
                        newfsrPhotons_pt.push_back(phot->pt());
                    }

                    if ( minDeltaRFSR < 0.07 || (phot->pt() > 4.0 && minDeltaRFSR < 0.5 && photoniso<1.0) ) {
                        if (verbose) cout<<" keeping this photon "<<endl;
                        fsrPhotons.push_back(*phot); 
                        fsrPhotons_lepindex.push_back(lepindex);
                        fsrPhotons_deltaR.push_back(minDeltaRFSR);
                    }
                    
                } // loop over fsr photon candidates
                if (verbose) {cout<<"finished filling fsr photon candidates"<<endl;}
            } // doFsrRecovery
            

            // creat vectors for selected objects
            vector<pat::Muon> selectedMuons;
            vector<pat::Electron> selectedElectrons;
            vector<pat::PFParticle> selectedFsrPhotons;
            
            if (verbose) cout<<"begin looking for higgs candidate"<<endl;                    
            findHiggsCandidate(selectedMuons,selectedElectrons,selectedFsrPhotons,iEvent);
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
                
                if( jetHelper.patjetID(jet)==1 && (passPU || !doPUJetID) ) {
                    
                    if (verbose) cout<<"passed pf jet id and pu jet id"<<endl;
                    if (verbose) cout<<"checking overlap with fsr photons..."<<endl;                                        
                    
                    bool isDeltaR_FSR = true;
                    /*
                      for(unsigned int phIndex = 0; phIndex < selectedFsrPhotons.size(); phIndex++) {
                      tempDeltaR = deltaR(jet.eta(),jet.phi(),selectedFsrPhotons[phIndex].eta(),selectedFsrPhotons[phIndex].phi());
                      if (tempDeltaR < 0.4) isDeltaR_FSR = false;
                      }
                    */
                    
                    // apply loose pt cut here (10 GeV cut is already applied in MINIAOD) since we are before JES/JER corrections
                    if(jet.pt() > 10.0 && fabs(jet.eta()) < jeteta_cut && isDeltaR_FSR) {
                        
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
                        
                    } // pass deltaR jet/fsr photons		    
                } // pass loose pf jet id and pu jet id
            } // all jets
            
            
            if( foundHiggsCandidate ){
                
                //M4L Error
                /*
                  if (verbose) cout<<"getting mass errors"<<endl;  
                  massErrorUCSD = massErr.getMassResolution(selectedElectrons, selectedMuons, selectedFsrPhotons);
                  massErrorUCSDCorr = massErr.getMassResolutionCorr(selectedElectrons, selectedMuons, selectedFsrPhotons, true, !isMC);
                  if(RecoFourMuEvent) {
                  massErrorUF = massErr.calc4muErr(selectedMuons,selectedFsrPhotons,false,!isMC);
                  massErrorUFCorr = massErr.calc4muErr(selectedMuons,selectedFsrPhotons,true,!isMC);
                  }
                  else if (RecoFourEEvent) {
                  massErrorUF = massErr.calc4eErr(selectedElectrons,selectedFsrPhotons,false,!isMC);
                  massErrorUFCorr = massErr.calc4eErr(selectedElectrons,selectedFsrPhotons,true,!isMC);                      
                  }
                  else if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent) { 
                  massErrorUF = massErr.calc2e2muErr(selectedElectrons,selectedMuons,selectedFsrPhotons,false,!isMC);
                  massErrorUFCorr = massErr.calc2e2muErr(selectedElectrons,selectedMuons,selectedFsrPhotons,true,!isMC);               
                  }
                */
                
                if (verbose) cout<<"storing H_p4_noFSR"<<endl; 
                if(RecoFourMuEvent) {
                    math::XYZTLorentzVector tmpHVec;
                    tmpHVec = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedMuons[2].p4() + selectedMuons[3].p4();
                    HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                    new ( (*H_p4_noFSR)[0] ) TLorentzVector(HVecNoFSR.Px(),HVecNoFSR.Py(),HVecNoFSR.Pz(),HVecNoFSR.E());
                }
                else if(RecoFourEEvent) {
                    math::XYZTLorentzVector tmpHVec;
                    tmpHVec = selectedElectrons[0].p4() + selectedElectrons[1].p4() + selectedElectrons[2].p4() + selectedElectrons[3].p4();
                    HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                    new ( (*H_p4_noFSR)[0] ) TLorentzVector(HVecNoFSR.Px(),HVecNoFSR.Py(),HVecNoFSR.Pz(),HVecNoFSR.E());
                }
                else if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){
                    math::XYZTLorentzVector tmpHVec;
                    tmpHVec = selectedMuons[0].p4() + selectedMuons[1].p4() + selectedElectrons[0].p4() + selectedElectrons[1].p4();
                    HVecNoFSR.SetPtEtaPhiM(tmpHVec.Pt(),tmpHVec.Eta(),tmpHVec.Phi(),tmpHVec.M());
                    new ( (*H_p4_noFSR)[0] ) TLorentzVector(HVecNoFSR.Px(),HVecNoFSR.Py(),HVecNoFSR.Pz(),HVecNoFSR.E());
                }
                
                if(Z2Vec.M() > 0) {
                    passedZ4lSelection = true;
                    if(Z2Vec.M() > mZ2Low && passedTrig) passedFullSelection = true;
                }
                
                if(bStudyResolution && bStudyFourLeptonResolution){
                    if (passedFullSelection) {
                        FourLepReso->fillHistograms(hContainer_, hContainer2D_, selectedElectrons, selectedMuons, selectedFsrPhotons, eventWeight,!isMC);
                    }
                } // bStudyFourLepResolution
                
            } // found higgs candidate 
            else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed higgs candidate"<<endl;}

            //Set All the Variables for Saved Trees (after finding higgs candidate) 
            if (verbose) cout<<"begin setting tree variables"<<endl; 
            setTreeVariables(iEvent, iSetup, selectedMuons, selectedElectrons, recoMuons, recoElectrons, goodJets);
            if (verbose) cout<<"finshed setting tree variables"<<endl;                            

            // Comput Matrix Elelements (after filling jets)
            if (foundHiggsCandidate) {

                int tmpIdL1,tmpIdL2,tmpIdL3,tmpIdL4;

                TLorentzVector  L11P4, L12P4, L21P4, L22P4, J1P4,  J2P4;
                TLorentzVector *Lep1, *Lep2, *Lep3, *Lep4,  *Jet1, *Jet2;

                TLorentzVector nullFourVector(0, 0, 0, 0);//hualin                                                                                                                                                                                   

                Lep1 = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[0]);
                Lep2 = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[1]);
                Lep3 = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[2]);
                Lep4 = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[3]);

                if (njets_pt30_eta4p7 > 0) {
                    Jet1 = (TLorentzVector*) jet_p4->At(0);
                    J1P4.SetPxPyPzE(Jet1->Px(),Jet1->Py(),Jet1->Pz(),Jet1->E());
                }
                if (njets_pt30_eta4p7 > 1) {
                    Jet2 = (TLorentzVector*) jet_p4->At(1);
                    J2P4.SetPxPyPzE(Jet2->Px(),Jet2->Py(),Jet2->Pz(),Jet2->E());
                }

                L11P4.SetPxPyPzE(Lep1->Px(),Lep1->Py(),Lep1->Pz(),Lep1->E()); tmpIdL1 = idL1;
                L12P4.SetPxPyPzE(Lep2->Px(),Lep2->Py(),Lep2->Pz(),Lep2->E()); tmpIdL2 = idL2;
                L21P4.SetPxPyPzE(Lep3->Px(),Lep3->Py(),Lep3->Pz(),Lep3->E()); tmpIdL3 = idL3;
                L22P4.SetPxPyPzE(Lep4->Px(),Lep4->Py(),Lep4->Pz(),Lep4->E()); tmpIdL4 = idL4;

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

            if (!isMC) passedEventsTree_All->Fill();        

        } //if 2 lepID
        else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed  2 ID"<<endl;}
    } //primary vertex,notDuplicate
    else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed notDuplicate"<<endl;}
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
UFHZZ4LAna::beginRun(edm::Run const&, const edm::EventSetup& iSetup)
{
    massErr.init(iSetup);
}


// ------------ method called when ending the processing of a run  ------------
void 
UFHZZ4LAna::endRun(edm::Run const&, edm::EventSetup const&)
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
UFHZZ4LAna::findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons,
                               std::vector< pat::PFParticle > &selectedFsrPhotons,const edm::Event& iEvent )
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    const double Zmass = 91.1876;

    unsigned int Nlep = lep_p4->GetLast()+1;
    if (verbose) cout<<Nlep<<" leptons in total"<<endl;

    // First, make all Z candidates including any FSR photons
    int n_Zs=0;
    vector<int> Z_lepindex1;
    vector<int> Z_lepindex2;
    vector<int> Z_fsrindex;

    for(unsigned int i=0; i<Nlep; i++){
        if (!(lep_tightId[i])) continue; // checking tight lepton ID
        for(unsigned int j=i+1; j<Nlep; j++){
            if (!(lep_tightId[j])) continue; // checking tight lepton ID

            // same flavor opposite charge
            if((lep_id[i]+lep_id[j])!=0) continue;

            TLorentzVector *li, *lj;
            li = (TLorentzVector*) lep_p4->At(i);
            lj = (TLorentzVector*) lep_p4->At(j);            

            TLorentzVector lilj = (*li)+(*lj);

            if (verbose) {
                cout<<"OSSF pair: i="<<i<<" id1="<<lep_id[i]<<" j="<<j<<" id2="<<lep_id[j]<<" pt1: "<<li->Pt()<<" pt2: "<<lj->Pt()<<" M: "<<lilj.M()<<endl;    
            }
            
            int phoindex=-1;
            int goodfsr_Npt4=0;
            double goodfsr_maxpt=0.0;
            double goodfsr_mindr=9999.9;

            for (unsigned int p=0; p<fsrPhotons.size(); p++) {

                if (!((fsrPhotons_lepindex[(int)p]==(int)i) || (fsrPhotons_lepindex[(int)p]==(int)j))) continue;

                TLorentzVector pho;
                pho.SetPtEtaPhiE(fsrPhotons[p].pt(),fsrPhotons[p].eta(),fsrPhotons[p].phi(),fsrPhotons[p].energy());

                TLorentzVector liljpho;
                liljpho = (*li)+(*lj)+pho;
                
                if (verbose) cout<<" fsr photon "<<p<<" mllgam "<<liljpho.M()<<" abs(mllgam-mz)-abs(mll-mz): "<<abs(liljpho.M()-Zmass) - abs(lilj.M()-Zmass)<<endl;
                if ( liljpho.M()<4.0 || liljpho.M()>100.0 ) continue;
                if ( abs(liljpho.M()-Zmass) > abs(lilj.M()-Zmass) ) continue;

                if (pho.Pt()>4.0) goodfsr_Npt4++;

                // if there's at least one photon with pT > 4 GeV, pick the one with highest pT 
                if (goodfsr_Npt4>0) {
                    if (pho.Pt()>goodfsr_maxpt) {
                        goodfsr_maxpt = pho.Pt();
                        phoindex=p;
                    }
                }
                // if all photons have pT < 4 GeV, pick the one that has the else  smallest dR to its closest lepton
                if (goodfsr_Npt4==0) {
                    if (fsrPhotons_deltaR[p]<goodfsr_mindr) {
                        goodfsr_mindr=fsrPhotons_deltaR[p];
                        phoindex=p;
                    }
                }
                if (verbose) cout<<"phoindex is: "<<phoindex<<endl;
            }
            
            TLorentzVector Z, Z_noFSR;
            if (phoindex>=0) {
                TLorentzVector pho;
                pho.SetPtEtaPhiE(fsrPhotons[phoindex].pt(),fsrPhotons[phoindex].eta(),fsrPhotons[phoindex].phi(),fsrPhotons[phoindex].energy());
                Z = (*li)+(*lj)+pho;
                Z_noFSR = (*li)+(*lj);
            } else {
                Z = (*li)+(*lj);
                Z_noFSR = (*li)+(*lj);
            }
            
            if (verbose) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;

            if (Z.M()>0.0) {
                n_Zs++;
                new ( (*Z_p4)[n_Zs-1] ) TLorentzVector(Z.Px(),Z.Py(),Z.Pz(),Z.E());
                new ( (*Z_p4_noFSR)[n_Zs-1] ) TLorentzVector(Z_noFSR.Px(),Z_noFSR.Py(),Z_noFSR.Pz(),Z_noFSR.E());
                Z_fsrindex.push_back(phoindex);
                Z_lepindex1.push_back(i);
                Z_lepindex2.push_back(j);
                if (verbose) cout<<" add Z_fsrindex: "<<phoindex<<" Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
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

            TLorentzVector *lep_i1, *lep_i2, *lep_j1, *lep_j2;
            lep_i1 = (TLorentzVector*) lep_p4->At(i1);
            lep_i2 = (TLorentzVector*) lep_p4->At(i2);
            lep_j1 = (TLorentzVector*) lep_p4->At(j1);
            lep_j2 = (TLorentzVector*) lep_p4->At(j2);

            TLorentzVector *Zi, *Zj;
            Zi = (TLorentzVector*) Z_p4->At(i);
            Zj = (TLorentzVector*) Z_p4->At(j);
            
            if (verbose) {cout<<"ZZ candidate Zi->M() "<<Zi->M()<<" Zj->M() "<<Zj->M()<<endl;}

            TLorentzVector Z1, Z2;
            int Z1index, Z2index;
            int Z1_lepindex[2] = {0,0};
            int Z2_lepindex[2] = {0,0};
            double Z1DeltaM, Z2SumPt;

            if (abs(Zi->M()-Zmass)<abs(Zj->M()-Zmass)) { 
                Z1index = i; Z2index = j;
                Z1 = (*Zi); Z2 = (*Zj);                 
                if (lep_i1->Pt()>lep_i2->Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
                else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }                
                if (lep_j1->Pt()>lep_j2->Pt()) { Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2; } 
                else { Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1; }                
                Z1DeltaM = abs(Zi->M()-Zmass);                
                Z2SumPt = lep_j1->Pt()+lep_j2->Pt();
            }
            else { 
                Z1index = j; Z2index = i;
                Z1 = (*Zj); Z2 = (*Zi); 
                if (lep_j1->Pt()>lep_j2->Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
                else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
                if (lep_i1->Pt()>lep_i2->Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
                else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
                Z1DeltaM = abs(Zj->M()-Zmass);
                Z2SumPt = lep_i1->Pt()+lep_i2->Pt();
            }
           
            // Check Leading and Subleading pt Cut
            vector<double> allPt;
            allPt.push_back(lep_i1->Pt()); allPt.push_back(lep_i2->Pt());
            allPt.push_back(lep_j1->Pt()); allPt.push_back(lep_j2->Pt());
            std::sort(allPt.begin(), allPt.end());
            if (verbose) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
            if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;
            
            // Check dR(li,lj)>0.02 for any i,j
            vector<double> alldR;
            alldR.push_back(deltaR(lep_i1->Eta(),lep_i1->Phi(),lep_i2->Eta(),lep_i2->Phi()));
            alldR.push_back(deltaR(lep_i1->Eta(),lep_i1->Phi(),lep_j1->Eta(),lep_j1->Phi()));
            alldR.push_back(deltaR(lep_i1->Eta(),lep_i1->Phi(),lep_j2->Eta(),lep_j2->Phi()));
            alldR.push_back(deltaR(lep_i2->Eta(),lep_i2->Phi(),lep_j1->Eta(),lep_j1->Phi()));
            alldR.push_back(deltaR(lep_i2->Eta(),lep_i2->Phi(),lep_j2->Eta(),lep_j2->Phi()));
            alldR.push_back(deltaR(lep_j1->Eta(),lep_j1->Phi(),lep_j2->Eta(),lep_j2->Phi()));            
            if (verbose) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
            if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

            // Check M(l+,l-)>4.0 GeV for any OS pair
            // Do not include FSR photons
            vector<double> allM;
            TLorentzVector i1i2;
            i1i2 = (*lep_i1)+(*lep_i2); allM.push_back(i1i2.M());
            TLorentzVector j1j2;
            j1j2 = (*lep_j1)+(*lep_j2); allM.push_back(j1j2.M());            

            if (lep_id[i1]*lep_id[j1]<0) {
                TLorentzVector i1j1;
                i1j1 = (*lep_i1)+(*lep_j1); allM.push_back(i1j1.M());
                TLorentzVector i2j2;
                i2j2 = (*lep_i2)+(*lep_j2); allM.push_back(i2j2.M());
            } else {
                TLorentzVector i1j2;
                i1j2 = (*lep_i1)+(*lep_j2); allM.push_back(i1j2.M());
                TLorentzVector i2j1;
                i2j1 = (*lep_i2)+(*lep_j1); allM.push_back(i2j1.M());
            }
            if (verbose) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
            if (*min_element(allM.begin(),allM.end())<4.0) {passedQCDcut=false; continue;}

            // Check which leptons include any FSR photons in their isolation cones
            double coneSize=0.4;
            double isoFSRi1=0.0, isoFSRi2=0.0, isoFSRj1=0.0, isoFSRj2=0.0;

            if (Z_fsrindex[i]>=0) {
                int ipho = Z_fsrindex[i];

                bool isoVeto=true;
                double dR_pho_i1 = deltaR(fsrPhotons[ipho].eta(),fsrPhotons[ipho].phi(),lep_i1->Eta(),lep_i1->Phi());
                if (abs(lep_id[i1])==13 && dR_pho_i1>0.01) isoVeto=false;
                if (abs(lep_id[i1])==11 && (abs(recoElectrons[lep_ptindex[i1]].superCluster()->eta())<1.479 || dR_pho_i1>0.08)) isoVeto=false;
                if (dR_pho_i1<coneSize && !isoVeto) isoFSRi1 += fsrPhotons[ipho].pt();

                isoVeto=true;
                double dR_pho_i2 = deltaR(fsrPhotons[ipho].eta(),fsrPhotons[ipho].phi(),lep_i2->Eta(),lep_i2->Phi());
                if (abs(lep_id[i2])==13 && dR_pho_i2>0.01) isoVeto=false;
                if (abs(lep_id[i2])==11 && (abs(recoElectrons[lep_ptindex[i2]].superCluster()->eta())<1.479 || dR_pho_i2>0.08)) isoVeto=false;
                if (dR_pho_i2<coneSize && !isoVeto) isoFSRi2 += fsrPhotons[ipho].pt();

                isoVeto=true;
                double dR_pho_j1 = deltaR(fsrPhotons[ipho].eta(),fsrPhotons[ipho].phi(),lep_j1->Eta(),lep_j1->Phi());
                if (abs(lep_id[j1])==13 && dR_pho_j1>0.01) isoVeto=false;
                if (abs(lep_id[j1])==11 && (abs(recoElectrons[lep_ptindex[j1]].superCluster()->eta())<1.479 || dR_pho_j1>0.08)) isoVeto=false;
                if (dR_pho_j1<coneSize && !isoVeto) isoFSRj1 += fsrPhotons[ipho].pt();

                isoVeto=true;
                double dR_pho_j2 = deltaR(fsrPhotons[ipho].eta(),fsrPhotons[ipho].phi(),lep_j2->Eta(),lep_j2->Phi());
                if (abs(lep_id[j2])==13 && dR_pho_j2>0.01) isoVeto=false;
                if (abs(lep_id[j2])==11 && (abs(recoElectrons[lep_ptindex[j2]].superCluster()->eta())<1.479 || dR_pho_j2>0.08)) isoVeto=false;
                if (dR_pho_j2<coneSize && !isoVeto) isoFSRj2 += fsrPhotons[ipho].pt();

            }
            if (Z_fsrindex[j]>=0) {
                int jpho = Z_fsrindex[j];

                bool isoVeto=true;
                double dR_pho_i1 = deltaR(fsrPhotons[jpho].eta(),fsrPhotons[jpho].phi(),lep_i1->Eta(),lep_i1->Phi());
                if (abs(lep_id[i1])==13 && dR_pho_i1>0.01) isoVeto=false;
                if (abs(lep_id[i1])==11 && (abs(recoElectrons[lep_ptindex[i1]].superCluster()->eta())<1.479 || dR_pho_i1>0.08)) isoVeto=false;
                if (dR_pho_i1<coneSize && !isoVeto) isoFSRi1 += fsrPhotons[jpho].pt();

                isoVeto=true;
                double dR_pho_i2 = deltaR(fsrPhotons[jpho].eta(),fsrPhotons[jpho].phi(),lep_i2->Eta(),lep_i2->Phi());
                if (abs(lep_id[i2])==13 && dR_pho_i2>0.01) isoVeto=false;
                if (abs(lep_id[i2])==11 && (abs(recoElectrons[lep_ptindex[i2]].superCluster()->eta())<1.479 || dR_pho_i2>0.08)) isoVeto=false;
                if (dR_pho_i2<coneSize && !isoVeto) isoFSRi2 += fsrPhotons[jpho].pt();

                isoVeto=true;
                double dR_pho_j1 = deltaR(fsrPhotons[jpho].eta(),fsrPhotons[jpho].phi(),lep_j1->Eta(),lep_j1->Phi());
                if (abs(lep_id[j1])==13 && dR_pho_j1>0.01) isoVeto=false;
                if (abs(lep_id[j1])==11 && (abs(recoElectrons[lep_ptindex[j1]].superCluster()->eta())<1.479 || dR_pho_j1>0.08)) isoVeto=false;
                if (dR_pho_j1<coneSize && !isoVeto) isoFSRj1 += fsrPhotons[jpho].pt();

                isoVeto=true;
                double dR_pho_j2 = deltaR(fsrPhotons[jpho].eta(),fsrPhotons[jpho].phi(),lep_j2->Eta(),lep_j2->Phi());
                if (abs(lep_id[j2])==13 && dR_pho_j2>0.01) isoVeto=false;
                if (abs(lep_id[j2])==11 && (abs(recoElectrons[lep_ptindex[j2]].superCluster()->eta())<1.479 || dR_pho_j2>0.08)) isoVeto=false;
                if (dR_pho_j2<coneSize && !isoVeto) isoFSRj2 += fsrPhotons[jpho].pt();

            }

            double isoLi1 = (lep_isoCH[i1]+std::max(lep_isoNH[i1]+lep_isoPhot[i1]-lep_isoPUcorr[i1]-isoFSRi1,0.0))/lep_i1->Pt();
            double isoLi2 = (lep_isoCH[i2]+std::max(lep_isoNH[i2]+lep_isoPhot[i2]-lep_isoPUcorr[i2]-isoFSRi2,0.0))/lep_i2->Pt();
            double isoLj1 = (lep_isoCH[j1]+std::max(lep_isoNH[j1]+lep_isoPhot[j1]-lep_isoPUcorr[j1]-isoFSRj1,0.0))/lep_j1->Pt();
            double isoLj2 = (lep_isoCH[j2]+std::max(lep_isoNH[j2]+lep_isoPhot[j2]-lep_isoPUcorr[j2]-isoFSRj2,0.0))/lep_j2->Pt();

            if (verbose) cout<<"isoLi1: "<<isoLi1<<" isoLi2: "<<isoLi2<<" isoLj1: "<<isoLj1<<" isoLj2: "<<isoLj2<<endl;
            // Check isolation cut, subtracting any FSR photons
            if (isoLi1>((abs(lep_id[i1])==11) ? isoCutEl : isoCutMu)) continue;
            if (isoLi2>((abs(lep_id[i2])==11) ? isoCutEl : isoCutMu)) continue;
            if (isoLj1>((abs(lep_id[j1])==11) ? isoCutEl : isoCutMu)) continue;
            if (isoLj2>((abs(lep_id[j2])==11) ? isoCutEl : isoCutMu)) continue;


            // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
            // only for 4mu or 4e ZZ candidates
            // consider FSR photons already available in this ZZ candidate
            bool passSmartCut=true;
            if ( abs(lep_id[i1])==abs(lep_id[j1])) {

                TLorentzVector Za, Zb;

                TLorentzVector tmpZa, tmpZb;
                int i_tmpZa, j_tmpZa, i_tmpZb, j_tmpZb;
                if (lep_id[i1]==lep_id[j1]) {                  
                    tmpZa = (*lep_i1)+(*lep_j2);
                    i_tmpZa = i1; j_tmpZa = j2;
                    tmpZb = (*lep_i2)+(*lep_j1);                    
                    i_tmpZb = i2; j_tmpZb = j1;
                } else {
                    tmpZa = (*lep_i1)+(*lep_j1);
                    i_tmpZa = i1; j_tmpZa = j1;
                    tmpZb = (*lep_i2)+(*lep_j2);
                    i_tmpZb = i2; j_tmpZb = j2;
                }

                double phopt_i=0.0;
                double phodr_i=9999.9;

                Za = tmpZa;
                Zb = tmpZb;
                
                if (Z_fsrindex[i]>=0 && fsrPhotons_lepindex[Z_fsrindex[i]]==i_tmpZa) {
                    TLorentzVector pho;
                    int p = Z_fsrindex[i];
                    pho.SetPtEtaPhiE(fsrPhotons[p].pt(),fsrPhotons[p].eta(),fsrPhotons[p].phi(),fsrPhotons[p].energy());
                    TLorentzVector tmpZa_fsr;
                    tmpZa_fsr = tmpZa+pho;                
                    if ( tmpZa_fsr.M()>4.0 && tmpZa_fsr.M()<100.0 &&  abs(tmpZa_fsr.M()-Zmass) < abs(tmpZa.M()-Zmass)) {
                        phopt_i = fsrPhotons[p].pt();
                        phodr_i = fsrPhotons_deltaR[p];
                        Za = tmpZa_fsr;
                    }
                }
                if (Z_fsrindex[j]>=0 && fsrPhotons_lepindex[Z_fsrindex[j]]==j_tmpZa) {
                    TLorentzVector pho;
                    int p = Z_fsrindex[j];
                    pho.SetPtEtaPhiE(fsrPhotons[p].pt(),fsrPhotons[p].eta(),fsrPhotons[p].phi(),fsrPhotons[p].energy());
                    TLorentzVector tmpZa_fsr;
                    tmpZa_fsr = tmpZa+pho;                
                    if ( tmpZa_fsr.M()>4.0 && tmpZa_fsr.M()<100.0 &&  abs(tmpZa_fsr.M()-Zmass) < abs(Za.M()-Zmass)) {
                        if ( (max(phopt_i,fsrPhotons[p].pt())>4.0 && fsrPhotons[p].pt()>phopt_i) 
                             || (max(phopt_i,fsrPhotons[p].pt())<4.0 && fsrPhotons_deltaR[p]<phodr_i) ) { 
                            Za = tmpZa_fsr;
                        }
                    }
                }
                phopt_i=0.0;
                phodr_i=9999.9;
                if (Z_fsrindex[i]>=0 && fsrPhotons_lepindex[Z_fsrindex[i]]==i_tmpZb) {
                    TLorentzVector pho;
                    int p = Z_fsrindex[i];
                    pho.SetPtEtaPhiE(fsrPhotons[p].pt(),fsrPhotons[p].eta(),fsrPhotons[p].phi(),fsrPhotons[p].energy());
                    TLorentzVector tmpZb_fsr;
                    tmpZb_fsr = tmpZb+pho;                
                    if ( tmpZb_fsr.M()>4.0 && tmpZb_fsr.M()<100.0 &&  abs(tmpZb_fsr.M()-Zmass) < abs(tmpZb.M()-Zmass)) {
                        phopt_i = fsrPhotons[p].pt();
                        phodr_i = fsrPhotons_deltaR[p];
                        Zb = tmpZb_fsr;
                    }
                }
                if (Z_fsrindex[j]>=0 && fsrPhotons_lepindex[Z_fsrindex[j]]==j_tmpZb) {
                    TLorentzVector pho;
                    int p = Z_fsrindex[j];
                    pho.SetPtEtaPhiE(fsrPhotons[p].pt(),fsrPhotons[p].eta(),fsrPhotons[p].phi(),fsrPhotons[p].energy());
                    TLorentzVector tmpZb_fsr;
                    tmpZb_fsr = tmpZb+pho;                
                    if ( tmpZb_fsr.M()>4.0 && tmpZb_fsr.M()<100.0 &&  abs(tmpZb_fsr.M()-Zmass) < abs(Zb.M()-Zmass)) {
                        if ( (max(phopt_i,fsrPhotons[p].pt())>4.0 && fsrPhotons[p].pt()>phopt_i) 
                             || (max(phopt_i,fsrPhotons[p].pt())<4.0 && fsrPhotons_deltaR[p]<phodr_i) ) { 
                            Zb = tmpZb_fsr;
                        }
                    }                                           
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

        TLorentzVector *lep_1, *lep_2, *lep_3, *lep_4;
        lep_1 = (TLorentzVector*) lep_p4->At(lep_Hindex[0]);
        lep_2 = (TLorentzVector*) lep_p4->At(lep_Hindex[1]);
        lep_3 = (TLorentzVector*) lep_p4->At(lep_Hindex[2]);
        lep_4 = (TLorentzVector*) lep_p4->At(lep_Hindex[3]);

        if (Z_Hindex[0]>=0 && Z_fsrindex[Z_Hindex[0]]>=0) {
            nFSRPhotons++;
            FSR_Z1 = true;
            int p = Z_fsrindex[Z_Hindex[0]];
            if (verbose) cout<<"found fsr Z1, Z_Hindex[0] "<<Z_Hindex[0]<<" Z_fsrindex[Z_Hindex[0]] "<<Z_fsrindex[Z_Hindex[0]]<<endl;
            selectedFsrPhotons.push_back(fsrPhotons[p]);
            new ( (*phofsr_p4)[nFSRPhotons-1] ) TLorentzVector(fsrPhotons[p].px(),fsrPhotons[p].py(),fsrPhotons[p].pz(),fsrPhotons[p].energy());
            phofsr_lepindex.push_back(fsrPhotons_lepindex[p]);
            //phofsr_lepindex.push_back(associatedPh1);
            TLorentzVector *pho;
            pho = (TLorentzVector*) phofsr_p4->At(nFSRPhotons-1);
            TLorentzVector lepfsr;
            if (fsrPhotons_lepindex[p]==lep_Hindex[0]) lepfsr = (*lep_1)+(*pho);
            else if (fsrPhotons_lepindex[p]==lep_Hindex[1]) lepfsr = (*lep_2)+(*pho);
            new ( (*lep_p4_FSR)[fsrPhotons_lepindex[p]] ) TLorentzVector(lepfsr.Px(),lepfsr.Py(),lepfsr.Pz(),lepfsr.E());
        }
        if (Z_Hindex[1]>=0 && Z_fsrindex[Z_Hindex[1]]>=0) {
            nFSRPhotons++;
            FSR_Z2 = true;
            int p = Z_fsrindex[Z_Hindex[1]];
            if (verbose) cout<<"found fsr Z2, Z_Hindex[1] "<<Z_Hindex[1]<<" Z_fsrindex[Z_Hindex[1]] "<<Z_fsrindex[Z_Hindex[1]]<<endl;
            selectedFsrPhotons.push_back(fsrPhotons[p]);
            new ( (*phofsr_p4)[nFSRPhotons-1] ) TLorentzVector(fsrPhotons[p].px(),fsrPhotons[p].py(),fsrPhotons[p].pz(),fsrPhotons[p].energy()); 
            phofsr_lepindex.push_back(fsrPhotons_lepindex[p]);
            TLorentzVector *pho;
            pho = (TLorentzVector*) phofsr_p4->At(nFSRPhotons-1);
            TLorentzVector lepfsr;
            if (fsrPhotons_lepindex[p]==lep_Hindex[2]) lepfsr = (*lep_3)+(*pho);
            else if (fsrPhotons_lepindex[p]==lep_Hindex[3]) lepfsr = (*lep_4)+(*pho);
            new ( (*lep_p4_FSR)[fsrPhotons_lepindex[p]] ) TLorentzVector(lepfsr.Px(),lepfsr.Py(),lepfsr.Pz(),lepfsr.E());

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
    tree->Branch("genWeight",&genWeight,"genWeight/D");
    tree->Branch("pileupWeight",&pileupWeight,"pileupWeight/D");
    tree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/D");
    tree->Branch("eventWeight",&eventWeight,"eventWeight/D");
    tree->Branch("crossSection",&crossSection,"crossSection/D");
    tree->Branch("filterEff",&filterEff,"filterEff/D");

    // Lepton variables
    lep_p4 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("lep_p4","TClonesArray", &lep_p4, 128000, 0);
    lep_p4_FSR = new TClonesArray("TLorentzVector", 10);
    tree->Branch("lep_p4_FSR","TClonesArray", &lep_p4_FSR, 128000, 0);
    tree->Branch("lep_Hindex",&lep_Hindex,"lep_Hindex[4]/I");
    tree->Branch("lep_genindex",&lep_genindex);
    tree->Branch("lep_id",&lep_id);
    tree->Branch("lep_missingHits",&lep_missingHits);
    tree->Branch("lep_mva",&lep_mva);
    tree->Branch("lep_tightId",&lep_tightId);
    tree->Branch("lep_Sip",&lep_Sip);
    tree->Branch("lep_IP",&lep_IP);
    tree->Branch("lep_isoNH",&lep_isoNH);
    tree->Branch("lep_isoCH",&lep_isoCH);
    tree->Branch("lep_isoPhot",&lep_isoPhot);
    tree->Branch("lep_isoPU",&lep_isoPU);
    tree->Branch("lep_isoPUcorr",&lep_isoPUcorr);
    tree->Branch("lep_RelIso",&lep_RelIso);
    tree->Branch("nisoleptons",&nisoleptons,"nisoleptons/I");
    tree->Branch("muRho",&muRho,"muRho/D");
    tree->Branch("elRho",&elRho,"elRho/D");
    tree->Branch("pTL1",&pTL1,"pTL1/D");
    tree->Branch("pTL2",&pTL2,"pTL2/D");
    tree->Branch("pTL3",&pTL3,"pTL3/D");
    tree->Branch("pTL4",&pTL4,"pTL4/D");
    tree->Branch("idL1",&idL1,"idL1/I");
    tree->Branch("idL2",&idL2,"idL2/I");
    tree->Branch("idL3",&idL3,"idL3/I");
    tree->Branch("idL4",&idL4,"idL4/I");
    tree->Branch("etaL1",&etaL1,"etaL1/D");
    tree->Branch("etaL2",&etaL2,"etaL2/D");
    tree->Branch("etaL3",&etaL3,"etaL3/D");
    tree->Branch("etaL4",&etaL4,"etaL4/D");
    tree->Branch("pTL1FSR",&pTL1FSR,"pTL1FSR/D");
    tree->Branch("pTL2FSR",&pTL2FSR,"pTL2FSR/D");
    tree->Branch("pTL3FSR",&pTL3FSR,"pTL3FSR/D");
    tree->Branch("pTL4FSR",&pTL4FSR,"pTL4FSR/D");

    //Higgs Candidate Variables
    H_p4 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("H_p4","TClonesArray", &H_p4, 128000, 0);
    H_p4_noFSR = new TClonesArray("TLorentzVector", 10);
    tree->Branch("H_p4_noFSR","TClonesArray", &H_p4_noFSR, 128000, 0);

    tree->Branch("mass4l",&mass4l,"mass4l/D");
    tree->Branch("mass4l_noFSR",&mass4l_noFSR,"mass4l_noFSR/D");
    tree->Branch("mass4mu",&mass4mu,"mass4mu/D");
    tree->Branch("mass4e",&mass4e,"mass4e/D");
    tree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/D");
    tree->Branch("pT4l",&pT4l,"pT4l/D");
    tree->Branch("eta4l",&eta4l,"eta4l/D");
    tree->Branch("phi4l",&phi4l,"phi4l/D");
    tree->Branch("rapidity4l",&rapidity4l,"rapidity4l/D");
    tree->Branch("cosTheta1",&cosTheta1,"cosTheta1/F");
    tree->Branch("cosTheta2",&cosTheta2,"cosTheta2/F");
    tree->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
    tree->Branch("Phi",&Phi,"Phi/F");
    tree->Branch("Phi1",&Phi1,"Phi1/F");

    // Z candidate variables
    Z_p4 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("Z_p4","TClonesArray", &Z_p4, 128000, 0);
    Z_p4_noFSR = new TClonesArray("TLorentzVector", 10);
    tree->Branch("Z_p4_noFSR","TClonesArray", &Z_p4_noFSR, 128000, 0);
    tree->Branch("Z_Hindex",&Z_Hindex,"Z_Hindex[2]/I");
    tree->Branch("massZ1",&massZ1,"massZ1/D");
    tree->Branch("massZ2",&massZ2,"massZ2/D");  
    tree->Branch("pTZ1",&pTZ1,"pTZ1/D");
    tree->Branch("pTZ2",&pTZ2,"pTZ2/D");

    // MET
    met_p4 = new TClonesArray("TLorentzVector", 1);
    tree->Branch("met_p4","TClonesArray", &met_p4, 128000, 0);
    tree->Branch("met",&met,"met/D");

    // Jets
    jet_p4 = new TClonesArray("TLorentzVector", 20);
    tree->Branch("jet_p4","TClonesArray", &jet_p4, 128000, 0);
    jet_p4_jesup = new TClonesArray("TLorentzVector", 20);
    tree->Branch("jet_p4_jesup","TClonesArray", &jet_p4_jesup, 128000, 0);
    jet_p4_jesdn = new TClonesArray("TLorentzVector", 20);
    tree->Branch("jet_p4_jesdn","TClonesArray", &jet_p4_jesdn, 128000, 0);
    jet_p4_jerup = new TClonesArray("TLorentzVector", 20);
    tree->Branch("jet_p4_jerup","TClonesArray", &jet_p4_jesup, 128000, 0);
    jet_p4_jerdn = new TClonesArray("TLorentzVector", 20);
    tree->Branch("jet_p4_jerdn","TClonesArray", &jet_p4_jesdn, 128000, 0);

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

    tree->Branch("pt_leadingjet_pt30_eta4p7",&pt_leadingjet_pt30_eta4p7,"pt_leadingjet_pt30_eta4p7/D");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jesup",&pt_leadingjet_pt30_eta4p7_jesup,"pt_leadingjet_pt30_eta4p7_jesup/D");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jesdn",&pt_leadingjet_pt30_eta4p7_jesdn,"pt_leadingjet_pt30_eta4p7_jesdn/D");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jerup",&pt_leadingjet_pt30_eta4p7_jerup,"pt_leadingjet_pt30_eta4p7_jerup/D");
    tree->Branch("pt_leadingjet_pt30_eta4p7_jerdn",&pt_leadingjet_pt30_eta4p7_jerdn,"pt_leadingjet_pt30_eta4p7_jerdn/D");

    tree->Branch("absrapidity_leadingjet_pt30_eta4p7",&absrapidity_leadingjet_pt30_eta4p7,"absrapidity_leadingjet_pt30_eta4p7/D");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesup",&absrapidity_leadingjet_pt30_eta4p7_jesup,"absrapidity_leadingjet_pt30_eta4p7_jesup/D");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesdn",&absrapidity_leadingjet_pt30_eta4p7_jesdn,"absrapidity_leadingjet_pt30_eta4p7_jesdn/D");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerup",&absrapidity_leadingjet_pt30_eta4p7_jerup,"absrapidity_leadingjet_pt30_eta4p7_jerup/D");
    tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerdn",&absrapidity_leadingjet_pt30_eta4p7_jerdn,"absrapidity_leadingjet_pt30_eta4p7_jerdn/D");

    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7",&absdeltarapidity_hleadingjet_pt30_eta4p7,"absdeltarapidity_hleadingjet_pt30_eta4p7/D");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesup/D");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn/D");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerup/D");
    tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn/D");

    tree->Branch("DijetMass",&DijetMass,"DijetMass/D");
    tree->Branch("DijetDEta",&DijetDEta,"DijetDEta/D");
    tree->Branch("DijetFisher",&DijetFisher,"DijetFisher/D");

    // FSR Photons
    phofsr_p4 = new TClonesArray("TLorentzVector", 20);
    tree->Branch("phofsr_p4","TClonesArray", &phofsr_p4, 128000, 0);    
    tree->Branch("nFSRPhotons",&nFSRPhotons,"nFSRPhotons/I");
    tree->Branch("FSR_Z1",&FSR_Z1,"FSR_Z1/O");
    tree->Branch("FSR_Z2",&FSR_Z1,"FSR_Z2/O");
    tree->Branch("phofsr_lepindex",&phofsr_lepindex);
    tree->Branch("newfsrPhotons_dR",&newfsrPhotons_dR);
    tree->Branch("newfsrPhotons_iso",&newfsrPhotons_iso);
    tree->Branch("newfsrPhotons_pt",&newfsrPhotons_pt);

    // Z4l? FIXME
    tree->Branch("theta12",&theta12,"theta12/D"); 
    tree->Branch("theta13",&theta13,"theta13/D"); 
    tree->Branch("theta14",&theta14,"theta14/D");
    tree->Branch("minM3l",&minM3l,"minM3l/D"); 
    tree->Branch("Z4lmaxP",&Z4lmaxP,"Z4lmaxP/D"); 
    tree->Branch("minDeltR",&minDeltR,"minDeltR/D"); 
    tree->Branch("m3l_soft",&m3l_soft,"m3l_soft/D");
    tree->Branch("minMass2Lep",&minMass2Lep,"minMass2Lep/D"); 
    tree->Branch("maxMass2Lep",&maxMass2Lep,"maxMass2Lep/D");
    tree->Branch("thetaPhoton",&thetaPhoton,"thetaPhoton/D"); 
    tree->Branch("thetaPhotonZ",&thetaPhotonZ,"thetaPhotonZ/D");

    // Resolution
    tree->Branch("massErrorUF",&massErrorUF,"massErrorUF/D");
    tree->Branch("massErrorUFCorr",&massErrorUFCorr,"massErrorUFCorr/D");
    tree->Branch("massErrorUCSD",&massErrorUCSD,"massErrorUCSD/D");
    tree->Branch("massErrorUCSDCorr",&massErrorUCSDCorr,"massErrorUCSDCorr/D");

    // Event Category
    tree->Branch("EventCat",&EventCat,"EventCat/I");

    // -------------------------                                                                                                                                                                        
    // GEN level information                                                                                                                                                                            
    // -------------------------                                                                                                                                                                        
    //Event variables
    tree->Branch("GENfinalState",&GENfinalState,"GENfinalState/I");
    tree->Branch("passedFiducialSelection",&passedFiducialSelection,"passedFiducialSelection/O");

    // lepton variables
    GENlep_p4 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("GENlep_p4","TClonesArray", &GENlep_p4, 128000, 0);    
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
    GENH_p4 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("GENH_p4","TClonesArray", &GENH_p4, 128000, 0);    
    tree->Branch("GENmass4l",&GENmass4l,"GENmass4l/D");
    tree->Branch("GENmass4mu",&GENmass4mu,"GENmass4mu/D");
    tree->Branch("GENmass4e",&GENmass4e,"GENmass4e/D");
    tree->Branch("GENmass2e2mu",&GENmass2e2mu,"GENmass2e2mu/D");
    tree->Branch("GENpT4l",&GENpT4l,"GENpT4l/D");
    tree->Branch("GENeta4l",&GENeta4l,"GENeta4l/D");
    tree->Branch("GENrapidity4l",&GENrapidity4l,"GENrapidity4l/D");
    tree->Branch("GENcosTheta1",&GENcosTheta1,"GENcosTheta1/F");
    tree->Branch("GENcosTheta2",&GENcosTheta2,"GENcosTheta2/F");
    tree->Branch("GENcosThetaStar",&GENcosThetaStar,"GENcosThetaStar/F");
    tree->Branch("GENPhi",&GENPhi,"GENPhi/F");
    tree->Branch("GENPhi1",&GENPhi1,"GENPhi1/F");
    tree->Branch("GENMH",&GENMH,"GENMH/D");

    // Z candidate variables
    GENZ_p4 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("GENZ_p4","TClonesArray", &GENZ_p4, 128000, 0);    
    tree->Branch("GENZ_DaughtersId",&GENZ_DaughtersId); 
    tree->Branch("GENZ_MomId",&GENZ_MomId);
    tree->Branch("GENmassZ1",&GENmassZ1,"GENmassZ1/D");
    tree->Branch("GENmassZ2",&GENmassZ2,"GENmassZ2/D");  
    tree->Branch("GENpTZ1",&GENpTZ1,"GENpTZ1/D");
    tree->Branch("GENpTZ2",&GENpTZ2,"GENpTZ2/D");

    // Higgs variables directly from GEN particle
    tree->Branch("GENHmass",&GENHmass,"GENHmass/D");

    // Jets
    GENjet_p4 = new TClonesArray("TLorentzVector", 20);
    tree->Branch("GENjet_p4","TClonesArray", &GENjet_p4, 128000, 0);    
    tree->Branch("GENnjets_pt30_eta4p7",&GENnjets_pt30_eta4p7,"GENnjets_pt30_eta4p7/I");
    tree->Branch("GENpt_leadingjet_pt30_eta4p7",&GENpt_leadingjet_pt30_eta4p7,"GENpt_leadingjet_pt30_eta4p7/D");
    tree->Branch("GENabsrapidity_leadingjet_pt30_eta4p7",&GENabsrapidity_leadingjet_pt30_eta4p7,"GENabsrapidity_leadingjet_pt30_eta4p7/D");
    tree->Branch("GENabsdeltarapidity_hleadingjet_pt30_eta4p7",&GENabsdeltarapidity_hleadingjet_pt30_eta4p7,"GENabsdeltarapidity_hleadingjet_pt30_eta4p7/D");

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
                                   std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, std::vector<pat::Jet> goodJets)
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
        unsigned int Nleptons = lep_p4->GetLast()+1;
        for (unsigned int i=0; i<Nleptons; i++) {
            if (abs(lep_id[i])==13 && lep_RelIso[i]>isoCutMu) continue;
            if (abs(lep_id[i])==11 && lep_RelIso[i]>isoCutEl) continue;
            if (!(lep_tightId[i])) continue;
            TLorentzVector *thisLep;
            thisLep = (TLorentzVector*) lep_p4->At(i);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),thisLep->Eta(),thisLep->Phi());
            if (tempDeltaR<0.4) {
                isDeltaR_eta4p7 = false;
            }
        }

        // check overlap with fsr photons
        /*
        unsigned int N = phofsr_p4->GetLast()+1;
        for(unsigned int i=0; i<N; i++) {
            TLorentzVector *pho;
            pho = (TLorentzVector*) phofsr_p4->At(i);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),pho->Eta(),pho->Phi());
            if (tempDeltaR<0.4) {
                isDeltaR_eta4p7 = false;
            }
        }
        */

        /*
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
                new ( (*jet_p4)[njets_pt30_eta4p7-1] ) TLorentzVector(jet_jer->Px(), jet_jer->Py(), jet_jer->Pz(), jet_jer->Energy());
                jet_pumva.push_back(goodJets[k].userFloat("pileupJetId:fullDiscriminant"));
                jet_csvv2.push_back(goodJets[k].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
                jet_isbtag.push_back(1 ? goodJets[k].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")>BTagCut : 0);
                if (goodJets[k].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")>BTagCut) nbjets_pt30_eta4p7++;
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
                new ( (*jet_p4_jerup)[njets_pt30_eta4p7_jerup-1] ) TLorentzVector(jetPx_jerup,jetPy_jerup,jetPz_jerup,jetE_jerup);
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
                new ( (*jet_p4_jerdn)[njets_pt30_eta4p7_jerdn-1] ) TLorentzVector(jetPx_jerdn,jetPy_jerdn,jetPz_jerdn,jetE_jerdn);
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
                new ( (*jet_p4_jesup)[njets_pt30_eta4p7_jesup-1] ) TLorentzVector(jetPx_jesup,jetPy_jesup,jetPz_jesup,jetE_jesup);
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
                new ( (*jet_p4_jesdn)[njets_pt30_eta4p7_jesdn-1] ) TLorentzVector(jetPx_jesdn,jetPy_jesdn,jetPz_jesdn,jetE_jesdn);
            }
        }
        
    } // loop over jets

    if(njets_pt30_eta4p7>1){
        TLorentzVector *jet1, *jet2;
        jet1 = (TLorentzVector*) jet_p4->At(0);
        jet2 = (TLorentzVector*) jet_p4->At(1);
        TLorentzVector Dijet;
        Dijet = (*jet1)+(*jet2); 
        DijetMass = Dijet.M();
        DijetDEta = fabs(jet1->Eta()-jet2->Eta());
        // OLD MORIOND --- FisherDiscrim = 0.09407*fabs(VBFDeltaEta) + 4.1581e-4*VBFDiJetMass;
        DijetFisher = 0.18*fabs(DijetDEta) + 1.92e-4*DijetMass;
    }

    // Double loop over jets, for V-jet tagging
    for (int i=0; i<njets_pt30_eta4p7; i++) {
        for (int j=i+1; j<njets_pt30_eta4p7; j++) {
            if (i==j) continue;
            TLorentzVector *ijet, *jjet;
            ijet = (TLorentzVector*) jet_p4->At(i);
            jjet = (TLorentzVector*) jet_p4->At(j);
            if (ijet->Pt()<40.0 || abs(ijet->Eta())>2.4) continue;
            if (jjet->Pt()<40.0 || abs(jjet->Eta())>2.4) continue;
            TLorentzVector Dijet;
            Dijet = (*ijet)+(*jjet);
            double mass = Dijet.M();
            if (mass > 60 && mass < 120) nvjets_pt40_eta2p4++;            
        }
    }

    // Higgs Variables
    if( RecoFourMuEvent ){ finalState = 1;}
    if( RecoFourEEvent  ){ finalState = 2;}
    if( RecoTwoETwoMuEvent ){ finalState = 3;}
    if( RecoTwoMuTwoEEvent ){ finalState = 4;}

    new ( (*H_p4)[0] ) TLorentzVector(HVec.Px(),HVec.Py(),HVec.Pz(),HVec.Energy());
    new ( (*H_p4_noFSR)[0] ) TLorentzVector(HVecNoFSR.Px(),HVecNoFSR.Py(),HVecNoFSR.Pz(),HVecNoFSR.Energy());

    mass4l = HVec.M();
    mass4l_noFSR = HVecNoFSR.M();

    if(RecoFourMuEvent){mass4mu = HVec.M();}
    else{mass4mu = -1;}
    if(RecoFourEEvent){ mass4e = HVec.M();}
    else{ mass4e = -1;}
    if(RecoTwoETwoMuEvent || RecoTwoMuTwoEEvent){ mass2e2mu = HVec.M(); }
    else{ mass2e2mu = -1;}

    pT4l = HVec.Pt();
    eta4l = HVec.Eta();
    rapidity4l = HVec.Rapidity();
    phi4l = HVec.Phi();

    pTZ1 = Z1Vec.Pt();
    pTZ2 = Z2Vec.Pt();
    massZ1 = Z1Vec.M();
    massZ2 = Z2Vec.M();

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

        TLorentzVector *Lep1FSR, *Lep2FSR, *Lep3FSR, *Lep4FSR;
        Lep1FSR = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[0]);
        Lep2FSR = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[1]);
        Lep3FSR = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[2]);
        Lep4FSR = (TLorentzVector*) lep_p4_FSR->At(lep_Hindex[3]);
        pTL1FSR = Lep1FSR->Pt();
        pTL2FSR = Lep2FSR->Pt();
        pTL3FSR = Lep3FSR->Pt();
        pTL4FSR = Lep4FSR->Pt();

        //FIXME put some protections on number of leptons
        TLorentzVector *Lep1;
        Lep1 = (TLorentzVector*) lep_p4->At(lep_Hindex[0]);
        idL1 = lep_id[lep_Hindex[0]];
        pTL1 = Lep1->Pt();
        etaL1 = Lep1->Eta();

        TLorentzVector *Lep2;
        Lep2 = (TLorentzVector*) lep_p4->At(lep_Hindex[1]);
        idL2 = lep_id[lep_Hindex[1]];
        pTL2 = Lep2->Pt();
        etaL2 = Lep2->Eta();

        TLorentzVector *Lep3;
        Lep3 = (TLorentzVector*) lep_p4->At(lep_Hindex[2]);
        idL3 = lep_id[lep_Hindex[2]];
        pTL3 = Lep3->Pt();
        etaL3 = Lep3->Eta();
        
        TLorentzVector *Lep4;
        Lep4 = (TLorentzVector*) lep_p4->At(lep_Hindex[3]);
        idL4 = lep_id[lep_Hindex[3]];
        pTL4 = Lep4->Pt();
        etaL4 = Lep4->Eta();

    }

    for(unsigned int i = 0; i < lep_pt.size(); i++) {
        if ((int)i==lep_Hindex[0] || (int)i==lep_Hindex[1] || (int)i==lep_Hindex[2] || (int)i==lep_Hindex[3]) { nisoleptons++; }
        else {
            if (abs(lep_id[i])==11 && lep_tightId[i]==1 && lep_RelIso[i]<isoCutEl) { nisoleptons++; }
            else if (abs(lep_id[i])==13 && lep_tightId[i]==1 && lep_RelIso[i]<isoCutMu) { nisoleptons++; }
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
    int nGENLeptons=-1;

    if (verbose) cout<<"begin looping on gen particles"<<endl;
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
        j++;

        if (abs(genPart->pdgId())==11  || abs(genPart->pdgId())==13 || abs(genPart->pdgId())==15) {

            if (!(genPart->status()==1)) continue;
            if ( !(genAna.MotherID(&prunedgenParticles->at(j))==23 || abs(genAna.MotherID(&prunedgenParticles->at(j)))==24) ) continue;

            nGENLeptons++;

            // Collect FSR photons within dR<0.1

            TLorentzVector lep_dressed;            
            lep_dressed.SetPtEtaPhiE(genPart->pt(),genPart->eta(),genPart->phi(),genPart->energy());
            set<int> gen_fsrset;
            for(size_t k=0; k<packedgenParticles->size();k++){

                //const Candidate * motherInPrunedCollection = (*packed)[k].mother(0) ;
                if( (*packedgenParticles)[k].status() != 1) continue; // stable particles only
                if( (*packedgenParticles)[k].pdgId() != 22) continue; // only photons
                double this_dR_lgamma = deltaR(genPart->eta(), genPart->phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
                if (verbose) cout<<"lep id: "<<genPart->pdgId()<<" lep eta: "<< genPart->eta()<<" pho pt: "<<(*packedgenParticles)[k].pt()<<" pho eta: "<<(*packedgenParticles)[k].eta() <<" mother0 id: "<<(*packedgenParticles)[j].mother(0)->pdgId()<<" dr: "<<this_dR_lgamma<<endl;
                //if( (*packedgenParticles)[j].mother(0)->pdgId() != genPart->pdgId()) continue; //only photons with mother id equal to this lepton's id
                bool idmatch=false;
                if ((*packedgenParticles)[k].mother(0)->pdgId()==genPart->pdgId() ) idmatch=true;
                const reco::Candidate * mother = (*packedgenParticles)[k].mother(0);
                for(size_t m=0;m<mother->numberOfMothers();m++) {
                    if (verbose) cout<<"pho mother "<<m<<" pho mother id: "<<(*packedgenParticles)[k].mother(m)->pdgId()<<endl;
                    if ( (*packedgenParticles)[k].mother(m)->pdgId() == genPart->pdgId() ) idmatch=true;
                }
                if (!idmatch) continue;
                if(this_dR_lgamma<genFSRDrCut) {
                    gen_fsrset.insert(k);
                    TLorentzVector gamma;
                    gamma.SetPtEtaPhiE((*packedgenParticles)[k].pt(),(*packedgenParticles)[k].eta(),(*packedgenParticles)[k].phi(),(*packedgenParticles)[k].energy());
                    lep_dressed = lep_dressed+gamma;
                }
            } // Dressed leptons loop
            if (verbose) cout<<"lep pt "<<genPart->pt()<< " dressed pt: " << lep_dressed.Pt()<<endl;
  
            new ( (*GENlep_p4)[nGENLeptons] ) TLorentzVector(lep_dressed.Px(),lep_dressed.Py(),lep_dressed.Pz(),lep_dressed.Energy());
            GENlep_id.push_back( genPart->pdgId() );
            GENlep_status.push_back(genPart->status());
            GENlep_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
            GENlep_MomMomId.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j)));

            TLorentzVector *thisLep;
            thisLep = (TLorentzVector*) GENlep_p4->At(nGENLeptons);
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
                double this_dRvL = deltaR(thisLep->Eta(), thisLep->Phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
                if(this_dRvL<0.4) {
                    if (verbose) cout<<"adding to geniso id: "<<(*packedgenParticles)[k].pdgId()<<" status: "<<(*packedgenParticles)[k].status()<<" pt: "<<(*packedgenParticles)[k].pt()<<" dR: "<<this_dRvL<<endl;
                    this_GENiso = this_GENiso + (*packedgenParticles)[k].pt();
                    if ((*packedgenParticles)[k].charge()==0) this_GENneutraliso = this_GENneutraliso + (*packedgenParticles)[k].pt();
                    if ((*packedgenParticles)[k].charge()!=0) this_GENchargediso = this_GENchargediso + (*packedgenParticles)[k].pt();
                }
            } // GEN iso loop
            this_GENiso = this_GENiso/thisLep->Pt();
            GENlep_RelIso.push_back(this_GENiso);
            // END GEN iso calculation

        } // leptons
        
        if (genPart->pdgId()==25) {
            GENMH=genPart->mass();
        }

        
        if (genPart->pdgId()==23 && (genPart->status()>=20 && genPart->status()<30) ) {
            const reco::Candidate *Zdau0=genPart->daughter(0);
            if (fabs(Zdau0->pdgId())==23) {
                int ndau = genPart->numberOfDaughters();
                for (int d=0; d<ndau; d++) {
                    const reco::Candidate *Zdau=genPart->daughter(d);
                    if (verbose) cout<<"ZDau "<<d<<" id "<<fabs(Zdau->pdgId())<<endl;
                }
            }
            if (Zdau0) GENZ_DaughtersId.push_back(fabs(Zdau0->pdgId()));
           
            GENZ_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));                
            if (verbose) cout<<"GENZ status "<<genPart->status()<<" MomId: "<<genAna.MotherID(&prunedgenParticles->at(j))<< "DauId: "<<fabs(Zdau0->pdgId())<< endl;
        }
      
    }
    
    /////// DO THE FIDUCIAL VOLUME CALCULATION //////////////
    if (verbose) cout<<"begin fiducial volume calculation"<<endl;
    passedFiducialSelection=false;
    int nFiducialLeptons = 0;
    int nFiducialPtLead=0;
    int nFiducialPtSublead=0;
    
    for (unsigned int i=0; i<GENlep_id.size(); ++i) {    

        TLorentzVector *thisLep;
        thisLep = (TLorentzVector*) GENlep_p4->At(i);
        
        if ( ( (abs(GENlep_id[i]) == 13 && thisLep->Pt() > 5.0 && abs(thisLep->Eta()) < 2.4)
               || (abs(GENlep_id[i]) == 11 && thisLep->Pt() > 7.0 && abs(thisLep->Eta()) < 2.5) )
             && GENlep_RelIso[i]<genIsoCut ) {
                nFiducialLeptons++;
                if (thisLep->Pt()>leadingPtCut) nFiducialPtLead++;
                if (thisLep->Pt()>subleadingPtCut) nFiducialPtSublead++;
        }                
    }

    if (nFiducialLeptons>=4 && nFiducialPtLead>=1 && nFiducialPtSublead>=2) {
                                
        // START FIDUCIAL EVENT TOPOLOGY CUTS
        unsigned int L1=99; unsigned int L2=99; unsigned int L3=99; unsigned int L4=99;
        GENmass4l = -1.0; GENmass4e = -1.0; GENmass4mu = -1.0; GENmass2e2mu = -1.0;
        GENmassZ1 = -1.0; GENmassZ2 = -1.0; GENpT4l = -1.0; GENeta4l = 999.; GENrapidity4l = 999.;
        passedFiducialSelection = mZ1_mZ2(L1, L2, L3, L4);      
                
        GENlep_Hindex[0] = L1; GENlep_Hindex[1] = L2; GENlep_Hindex[2] = L3; GENlep_Hindex[3] = L4;

        if (passedFiducialSelection) {
            
            TLorentzVector *LS3_Z1_1 = (TLorentzVector*) GENlep_p4->At(L1);
            TLorentzVector *LS3_Z1_2 = (TLorentzVector*) GENlep_p4->At(L2);
            TLorentzVector *LS3_Z2_1 = (TLorentzVector*) GENlep_p4->At(L3);
            TLorentzVector *LS3_Z2_2 = (TLorentzVector*) GENlep_p4->At(L4);
            
            GENmass4l = ((*LS3_Z1_1)+(*LS3_Z1_2)+(*LS3_Z2_1)+(*LS3_Z2_2)).M();
            if (abs(GENlep_id[L1])==11 && abs(GENlep_id[L3])==11) {GENmass4e = GENmass4l;};
            if (abs(GENlep_id[L1])==13 && abs(GENlep_id[L3])==13) {GENmass4mu = GENmass4l;};
            if ( (abs(GENlep_id[L1])==11 || abs(GENlep_id[L1])==13) &&
                 (abs(GENlep_id[L3])==11 || abs(GENlep_id[L3])==13) &&
                 (abs(GENlep_id[L1])!=abs(GENlep_id[L3]) ) ) {GENmass2e2mu = GENmass4l;};
            GENpT4l = ((*LS3_Z1_1)+(*LS3_Z1_2)+(*LS3_Z2_1)+(*LS3_Z2_2)).Pt();
            GENeta4l = ((*LS3_Z1_1)+(*LS3_Z1_2)+(*LS3_Z2_1)+(*LS3_Z2_2)).Eta();
            GENrapidity4l = ((*LS3_Z1_1)+(*LS3_Z1_2)+(*LS3_Z2_1)+(*LS3_Z2_2)).Rapidity();
            GENmassZ1 = ((*LS3_Z1_1)+(*LS3_Z1_2)).M();
            GENmassZ2 = ((*LS3_Z2_1)+(*LS3_Z2_2)).M();
                    
            int tmpIdL1,tmpIdL2,tmpIdL3,tmpIdL4;
            TLorentzVector GENL11P4, GENL12P4, GENL21P4, GENL22P4;
            if(GENlep_id[L1] < 0){ GENL11P4.SetPxPyPzE((*LS3_Z1_1).Px(),(*LS3_Z1_1).Py(),(*LS3_Z1_1).Pz(),(*LS3_Z1_1).E()); tmpIdL1 = GENlep_id[L1];}
            else{ GENL11P4.SetPxPyPzE((*LS3_Z1_2).Px(),(*LS3_Z1_2).Py(),(*LS3_Z1_2).Pz(),(*LS3_Z1_2).E()); tmpIdL1 = GENlep_id[L2];}
            if(GENlep_id[L2] > 0){ GENL12P4.SetPxPyPzE((*LS3_Z1_2).Px(),(*LS3_Z1_2).Py(),(*LS3_Z1_2).Pz(),(*LS3_Z1_2).E()); tmpIdL2 = GENlep_id[L2];}
            else{ GENL12P4.SetPxPyPzE((*LS3_Z1_1).Px(),(*LS3_Z1_1).Py(),(*LS3_Z1_1).Pz(),(*LS3_Z1_1).E()); tmpIdL2 = GENlep_id[L1];}
            if(GENlep_id[L3] < 0){ GENL21P4.SetPxPyPzE((*LS3_Z2_1).Px(),(*LS3_Z2_1).Py(),(*LS3_Z2_1).Pz(),(*LS3_Z2_1).E()); tmpIdL3 = GENlep_id[L3];}
            else{ GENL21P4.SetPxPyPzE((*LS3_Z2_2).Px(),(*LS3_Z2_2).Py(),(*LS3_Z2_2).Pz(),(*LS3_Z2_2).E()); tmpIdL3 = GENlep_id[L4];}
            if(GENlep_id[L4] > 0) { GENL22P4.SetPxPyPzE((*LS3_Z2_2).Px(),(*LS3_Z2_2).Py(),(*LS3_Z2_2).Pz(),(*LS3_Z2_2).E()); tmpIdL4 = GENlep_id[L4];}
            else{ GENL22P4.SetPxPyPzE((*LS3_Z2_1).Px(),(*LS3_Z2_1).Py(),(*LS3_Z2_1).Pz(),(*LS3_Z2_1).E()); tmpIdL4 = GENlep_id[L3];}
                    
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
        unsigned int N=(u_int)GENlep_p4->GetLast()+1;
        for(unsigned int i = 0; i<N; i++) {
            for(unsigned int j = i+1; j<N; j++) {
                        
                // only consider the leptons from Z1 and Z2
                if (!(i==L1 || i==L2 || i==L3 || i==L4)) continue; 
                if (!(j==L1 || j==L2 || j==L3 || j==L4)) continue;
                        
                TLorentzVector *li, *lj;
                li = (TLorentzVector*) GENlep_p4->At(i);
                lj = (TLorentzVector*) GENlep_p4->At(j);
                TLorentzVector mll = (*li)+(*lj);
                        
                if(GENlep_id[i]*GENlep_id[j]<0) {
                    if(mll.M()<=4) { passedMassOS = false; break; }
                }
                        
                if(abs(GENlep_id[i]) != abs(GENlep_id[j])) {
                    double deltaR = (*li).DeltaR((*lj));
                    if(deltaR<=0.02) { passedElMuDeltaR = false; break; }
                }
                double deltaRll = (*li).DeltaR((*lj));
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
                unsigned int N=(u_int)GENlep_p4->GetLast()+1;
                for(unsigned int i = 0; i<N; i++) {
                    //if (GENlep_status[i]!=1) continue;
                    if (!(abs(GENlep_id[i])==11 || abs(GENlep_id[i])==13)) continue;
                    TLorentzVector *genlep;
                    genlep = (TLorentzVector*) GENlep_p4->At(i);
                    double dR = deltaR(genlep->Eta(), genlep->Phi(), genjet->eta(),genjet->phi());                        
                    if(dR<0.5) {
                        inDR_pt30_eta4p7=true;
                    }                                
                }
                
                if (verbose) cout<<"check overlap of gen jet with gen leptons"<<endl;
                // count number of gen jets which no gen leptons are inside its cone
                if (!inDR_pt30_eta4p7) { 
                    GENnjets_pt30_eta4p7++;
                    new ( (*GENjet_p4)[GENnjets_pt30_eta4p7] ) TLorentzVector(genjet->px(), genjet->py(), genjet->pz(), genjet->energy());
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

bool UFHZZ4LAna::mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4)
{

    double offshell = 999.0; bool findZ1 = false; bool passZ1 = false;

    L1 = 0; L2 = 0;

    unsigned int N = GENlep_p4->GetLast()+1;

    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<N; j++){

            if((GENlep_id[i]+GENlep_id[j])!=0) continue;

            TLorentzVector *li, *lj;
            li = (TLorentzVector*) GENlep_p4->At(i);
            lj = (TLorentzVector*) GENlep_p4->At(j);

            if ( abs(GENlep_id[i]) == 13 && (li->Pt() < 5.0 || abs(li->Eta()) > 2.4)) continue;
            if ( abs(GENlep_id[i]) == 11 && (li->Pt() < 7.0 || abs(li->Eta()) > 2.5)) continue;
            if ( GENlep_RelIso[i]>genIsoCut ) continue;

            if ( abs(GENlep_id[j]) == 13 && (lj->Pt() < 5.0 || abs(lj->Eta()) > 2.4)) continue;
            if ( abs(GENlep_id[j]) == 11 && (lj->Pt() < 7.0 || abs(lj->Eta()) > 2.5)) continue;
            if ( GENlep_RelIso[j]>genIsoCut ) continue;

            TLorentzVector mll = (*li)+(*lj);
            if(abs(mll.M()-Zmass)<offshell){
                double mZ1 = mll.M();
                L1 = i; L2 = j; findZ1 = true; offshell = abs(mZ1-Zmass);          
            }
        }    
    }

    TLorentzVector *l1, *l2;
    l1 = (TLorentzVector*) GENlep_p4->At(L1);
    l2 = (TLorentzVector*) GENlep_p4->At(L2);
    TLorentzVector ml1l2 = (*l1)+(*l2);

    if(ml1l2.M()>40 && ml1l2.M()<120 && findZ1) passZ1 = true;

    double pTL34 = 0.0; bool findZ2 = false;

    for(unsigned int i=0; i<N; i++){
        if(i==L1 || i==L2) continue; // can not be the lep from Z1
        for(unsigned int j=i+1; j<N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;            

            TLorentzVector *li, *lj;
            li = (TLorentzVector*) GENlep_p4->At(i);
            lj = (TLorentzVector*) GENlep_p4->At(j);            
            TLorentzVector Z2 = (*li)+(*lj);

            if ( abs(GENlep_id[i]) == 13 && (li->Pt() < 5.0 || abs(li->Eta()) > 2.4)) continue;
            if ( abs(GENlep_id[i]) == 11 && (li->Pt() < 7.0 || abs(li->Eta()) > 2.5)) continue;
            if ( GENlep_RelIso[i]>genIsoCut ) continue;

            if ( abs(GENlep_id[j]) == 13 && (lj->Pt() < 5.0 || abs(lj->Eta()) > 2.4)) continue;
            if ( abs(GENlep_id[j]) == 11 && (lj->Pt() < 7.0 || abs(lj->Eta()) > 2.5)) continue;
            if ( GENlep_RelIso[j]>genIsoCut ) continue;

            if ( ( (*li).Pt()+(*lj).Pt() ) >=pTL34 ) { // choose high sum pT pair satisfy the following selection
                double mZ2 = Z2.M();
                if(mZ2>12 && mZ2<120){
                    L3 = i; L4 = j; findZ2 = true; pTL34 = (*li).Pt()+(*lj).Pt();
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
