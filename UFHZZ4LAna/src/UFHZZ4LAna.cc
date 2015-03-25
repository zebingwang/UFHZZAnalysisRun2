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
#include <cmath>
#include <iomanip>

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
//#include "ZZMatrixElement/MELA/src/computeAngles.h" // removed for miniAOD

//MEKD
//#include "ZZMatrixElement/MEKD/interface/MEKD.h" // removed for miniAOD
//#include "ZZMatrixElement/MEKD/interface/MEKD_MG.h" // removed for miniAOD
//#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h" // removed for miniAOD

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
  
    void findHiggsCandidate(TClonesArray *lep_p4,
                            std::vector< pat::PFParticle > fsrPhotons, std::vector<double> deltaRVec,
                            int lep_Hindex[],
                            std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons,
                            std::vector< pat::PFParticle > &selectedFsrPhotons,const edm::Event& iEvent);

    void bookResolutionHistograms();
    void fillResolutionHistograms(edm::Handle<edm::View<pat::Muon> > muons);
    bool findZ(std::vector<pat::PFParticle> photons, std::vector<double> &deltaRVec, TLorentzVector &li, unsigned int &index1, TLorentzVector &lj, unsigned int &index2, int takenZ1, int taken, int assocMuon, TLorentzVector &ZVec, TLorentzVector &photVec, bool &foundPhoton);
  
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
    void setGENVariables(edm::Handle<reco::GenParticleCollection> genParticles,
                         edm::Handle<edm::View<reco::GenJet> > genJets);

    bool mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4);

    // -------------------------
    // RECO level information
    // -------------------------

    // Event Variables
    ULong64_t Run, Event, LumiSect;
    int nVtx;
    int finalState;
    bool passedFullSelection, passedZ4lSelection, passedQCDcut;

    // Event Weights
    double pileupWeight, dataMCWeight, eventWeight;

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
    vector<double> lep_Sip;
    vector<double> lep_IP;
    vector<double> lep_dIP;
    vector<double> lep_isoNH;
    vector<double> lep_isoCH;
    vector<double> lep_isoPhot;
    vector<double> lep_isoPU;
    vector<double> lep_RelIso;
    vector<int> lep_missingHits;
    double muRho, elRho;

    // Higgs candidate variables
    TClonesArray *H_p4;
    TClonesArray *H_p4_noFSR;
    double mass4l, mass4l_noFSR, mass4e, mass4mu, mass2e2mu, pT4l, eta4l, phi4l, rapidity4l;
    float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;

    // Z candidate variables
    TClonesArray *Z_p4;
    double massZ1, massZ2, pTZ1, pTZ2;

    // MET
    TClonesArray *met_p4;
    double met;

    // Jets
    TClonesArray *jet_p4;
    TClonesArray *jet_p4_jesup;
    TClonesArray *jet_p4_jesdn;
    TClonesArray *jet_p4_jerup;
    TClonesArray *jet_p4_jerdn;
    
    int njets_pt30_eta4p7;
    int njets_pt30_eta4p7_jesup;
    int njets_pt30_eta4p7_jesdn;
    int njets_pt30_eta4p7_jerup;
    int njets_pt30_eta4p7_jerdn;

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
    double DijetFisherDiscrim;

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
    TLorentzVector HVec, HVecNoFSR, Z1Vec, Z2Vec;
    bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;
    bool eventPassedPtAndIsoCuts, foundHiggsCandidate, twoLooseIsoLeptons;
    bool isIsolated, passedSIP3D, passedPtCuts;
    bool fourLep_Cleaned;
    bool isRecord; 
    bool foundZ1;

    // hist container
    std::map<std::string,TH1F*> histContainer_;
    std::map<std::string,TH2F*> histContainer2D_;
 
    //Input tags
    edm::InputTag photonSrc_;
    edm::InputTag elecSrc_;
    edm::InputTag muonSrc_;
    edm::InputTag correctedJetSrc_;
    edm::InputTag jetSrc_;
    edm::InputTag metSrc_;
    edm::InputTag vertexSrc_;
    edm::InputTag muRhoSrc_;
    edm::InputTag elRhoSrc_;

    // Configuration
    const double Zmass;
    double mZ1Low, mZ2Low;
    double mZ1High, mZ2High;
    double m4lLowCut;
    double pt_cut, eta_cut;
    std::string elecID;
    bool isMC, isSignal;
    double mH;
    double crossSection, filterEff, Lumi;
    bool weightEvents;
    double isoCut, earlyIsoCut, sip3dCut;
    double leadingPtCut, subleadingPtCut;
    double genIsoCut;
    double _elecPtCut, _muPtCut;
    bool reweightForPU;
    bool interactiveRun;
    std::string PUVersion;
    bool doVarDump, doBlinding, doFsrRecovery;
    bool bStudyResolution;
    bool bStudyDiLeptonResolution;
    bool bStudyFourLeptonResolution;
    bool verbose;

    // register to the TFileService
    edm::Service<TFileService> fs;

    // Counters
    double nEventsTotal;
};


UFHZZ4LAna::UFHZZ4LAna(const edm::ParameterSet& iConfig) :
    histContainer_(),
    histContainer2D_(),
    photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
    elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
    muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
    correctedJetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("correctedJetSrc" )),
    jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" )),
    metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc" )),
    vertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")),
    muRhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc")),
    elRhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc")),
    Zmass(91.1876),
    mZ1Low(iConfig.getUntrackedParameter<double>("mZ1Low",40)),
    mZ2Low(iConfig.getUntrackedParameter<double>("mZ2Low",12)), // was 12
    mZ1High(iConfig.getUntrackedParameter<double>("mZ1High",120)),
    mZ2High(iConfig.getUntrackedParameter<double>("mZ2High",120)),
    m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",100)),
    pt_cut(iConfig.getUntrackedParameter<double>("pt_cut",10.0)),
    eta_cut(iConfig.getUntrackedParameter<double>("eta_cut",4.7)),
    //elecID(iConfig.getUntrackedParameter<std::string>("elecID","mvaNonTrigV0")),
    elecID(iConfig.getUntrackedParameter<std::string>("elecID","cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium")), //FIXME
    isMC(iConfig.getUntrackedParameter<bool>("isMC",true)),
    isSignal(iConfig.getUntrackedParameter<bool>("isSignal",false)),
    mH(iConfig.getUntrackedParameter<double>("mH",0.0)),
    crossSection(iConfig.getUntrackedParameter<double>("CrossSection",1.0)),
    filterEff(iConfig.getUntrackedParameter<double>("FilterEff",1.0)),
    Lumi(iConfig.getUntrackedParameter<double>("Lumi",1.0)),
    weightEvents(iConfig.getUntrackedParameter<bool>("weightEvents",false)),
    isoCut(iConfig.getUntrackedParameter<double>("isoCut",0.4)),
    earlyIsoCut(iConfig.getUntrackedParameter<double>("earlyIsoCut",0.4)),
    sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
    leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
    subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
    genIsoCut(iConfig.getUntrackedParameter<double>("genIsoCut",0.4)), // was 0.4
    _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",7)),
    _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",5)),
    reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
    interactiveRun(iConfig.getUntrackedParameter<bool>("interactiveRun",false)),
    PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Legacy53X")),
    doVarDump(iConfig.getUntrackedParameter<bool>("doVarDump",false)),
    doBlinding(iConfig.getUntrackedParameter<bool>("doBlinding",false)),
    doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",true)),
    bStudyResolution(iConfig.getUntrackedParameter<bool>("bStudyResolution",false)),
    bStudyDiLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyDiLeptonResolution",false)),
    bStudyFourLeptonResolution(iConfig.getUntrackedParameter<bool>("bStudyFourLeptonResolution",false)),
    verbose(iConfig.getUntrackedParameter<bool>("verbose",false))
{
  
    if(!isMC){reweightForPU = false;}

    nEventsTotal=0;
    histContainer_["NEVENTS"]=fs->make<TH1F>("nEvents","nEvents in Sample",2,0,2);
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
    //using namespace MEMNames;

    if (verbose) cout<<"starting to analyze"<<endl;

    Run = iEvent.id().run();
    Event = iEvent.id().event();
    LumiSect = iEvent.id().luminosityBlock();
    if (verbose) cout<<Run<<":"<<Event<<":"<<LumiSect<<endl;

    // ======= Get Collections ======= //

    // vertex collection
    edm::Handle<reco::VertexCollection> vertex;
    iEvent.getByLabel(vertexSrc_,vertex);
    const reco::Vertex *PV = 0;
    if(!vertex->empty() && vertex->size() > 0) PV = &(vertex->at(0));
    if (verbose) cout<<"got vertex collection"<<endl;

    // photon collection 
    edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByLabel(photonSrc_,photons);
    if (verbose) cout<<"got photon collection"<<endl;
  
    // electron collection
    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByLabel(elecSrc_,electrons);
    if (verbose) cout<<"got electron collection"<<endl;

    // muon collection
    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByLabel(muonSrc_,muons);
    if (verbose) cout<<"got muon collection"<<endl;

    // met collection 
    edm::Handle<edm::View<pat::MET> > mets;
    iEvent.getByLabel(metSrc_,mets);
    if (verbose) cout<<"got met collection"<<endl;
    
    // beamspot collection
    edm::Handle<reco::BeamSpot> beamspot;
    iEvent.getByLabel("offlineBeamSpot",beamspot);
    if (verbose) cout<<"got beamspot collection"<<endl;
  
    // Rho Correction
    edm::Handle<double> eventRhoMu;
    iEvent.getByLabel(muRhoSrc_,eventRhoMu);
    muRho = *eventRhoMu;
    if (verbose) cout<<"got muon rho"<<endl;

    edm::Handle<double> eventRhoE;
    iEvent.getByLabel(elRhoSrc_,eventRhoE);
    elRho = *eventRhoE;
    if (verbose) cout<<"got electron rho"<<endl;

    // Particle Flow Cands
    edm::Handle<reco::PFCandidateCollection> pfCands;
    iEvent.getByLabel("packedPFCandidates",pfCands);
    if (verbose) cout<<"got pfcand collection"<<endl;

    edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;
    iEvent.getByLabel("boostedFsrPhotons",photonsForFsr);
    if (verbose) cout<<"got photons for fsr collection"<<endl;
  
    // GEN collection
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("prunedGenParticles", genParticles);
    if (verbose) cout<<"got genParticles collection"<<endl;

    //VBF 
    edm::Handle<edm::View<pat::Jet> > correctedJets;
    iEvent.getByLabel(correctedJetSrc_,correctedJets);
    if (verbose) cout<<"got corrected jet collection"<<endl;

    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel(jetSrc_,jets);
    if (verbose) cout<<"got jet collection"<<endl;
  
    edm::Handle<edm::ValueMap<float> > puJetIdMva;
    iEvent.getByLabel("puJetMva","full53xDiscriminant",puJetIdMva);
    //iEvent.getByLabel("puJetMva","fullDiscriminant",puJetIdMva);
    if (verbose) cout<<"got pileupjet mva"<<endl;

    edm::Handle<edm::ValueMap<int> > puJetIdFlag;
    iEvent.getByLabel("puJetMva","full53xId",puJetIdFlag);  
    //iEvent.getByLabel("puJetMva","fullId",puJetIdFlag);  
    if (verbose) cout<<"got pileupjet id flags"<<endl;
  
    edm::Handle<edm::View<reco::GenJet> > genJets;
    iEvent.getByLabel("slimmedGenJets", genJets);
    if (verbose) cout<<"got gen jet collection"<<endl;

    // ============ Initialize Variables ============= //

    // Event Variables
    if (verbose) cout<<"clear event variables"<<endl;
    nVtx = -1.0;
    finalState = -1;
    passedFullSelection=false; passedZ4lSelection=false; passedQCDcut=false;

    // Event Weights
    pileupWeight=0.0; dataMCWeight=0.0; eventWeight=0.0;

    //lepton variables
    if (verbose) cout<<"clear lepton variables"<<endl;
    if (lep_p4->GetLast()!=-1) lep_p4->Clear(); 
    cout<<"test1"<<endl;
    if (lep_p4_FSR->GetLast()!=-1) lep_p4_FSR->Clear();
    cout<<"test2"<<endl;
    for (int i=0; i<4; ++i) {lep_Hindex[i]=-1;}
    cout<<"test3"<<endl;
    pTL1=-1.0; pTL2=-1.0; pTL3=-1.0; pTL4=-1.0;
    cout<<"test4"<<endl;
    etaL1=9999.0; etaL2=9999.0; etaL3=9999.0; etaL4=9999.0;
    cout<<"test5"<<endl;
    idL1=9999; idL2=9999; idL3=9999; idL4=9999;
    cout<<"test6"<<endl;
    pTL1FSR=-1.0; pTL2FSR=-1.0; pTL3FSR=-1.0; pTL4FSR=-1.0;
    cout<<"test7"<<endl;
    lep_genindex.clear(); lep_id.clear(); 
    cout<<"test8"<<endl;
    lep_mva.clear();
    cout<<"test9"<<endl;
    lep_Sip.clear(); lep_IP.clear(); lep_dIP.clear();
    cout<<"test10"<<endl;
    lep_isoNH.clear(); lep_isoCH.clear(); lep_isoPhot.clear(); lep_isoPU.clear(); lep_RelIso.clear();
    cout<<"test11"<<endl;
    lep_missingHits.clear();
 
    // Higgs candidate variables
    if (verbose) cout<<"clear higgs/Z variables"<<endl;
    if (H_p4->GetLast()!=-1) H_p4->Clear();
    if (H_p4_noFSR->GetLast()!=-1) H_p4_noFSR->Clear();
    mass4l=-1.0; mass4l_noFSR=-1.0; mass4e=-1.0; mass4mu=-1.0; mass2e2mu=-1.0; pT4l=-1.0; eta4l=9999.0; phi4l=9999.0; rapidity4l=9999.0;
    cosTheta1=9999.0; cosTheta2=9999.0; cosThetaStar=9999.0; Phi=9999.0; Phi1=9999.0;

    // Z candidate variables
    if (Z_p4->GetLast()!=-1) Z_p4->Clear();
    massZ1=-1.0; massZ2=-1.0; pTZ1=-1.0; pTZ2=-1.0;

    // MET
    if (verbose) cout<<"clear jet/met variables"<<endl;
    if (met_p4->GetLast()!=-1) met_p4->Clear();
    met=-1.0;

    // Jets
    if (jet_p4->GetLast()!=-1) jet_p4->Clear();
    if (jet_p4_jesup->GetLast()!=-1) jet_p4_jesup->Clear();
    if (jet_p4_jesdn->GetLast()!=-1) jet_p4_jesdn->Clear();
    if (jet_p4_jerup->GetLast()!=-1) jet_p4_jerup->Clear();
    if (jet_p4_jerdn->GetLast()!=-1) jet_p4_jerdn->Clear();

    njets_pt30_eta4p7=-1;
    njets_pt30_eta4p7_jesup=-1;
    njets_pt30_eta4p7_jesdn=-1;
    njets_pt30_eta4p7_jerup=-1;
    njets_pt30_eta4p7_jerdn=-1;

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
    DijetFisherDiscrim=9999.0;

    // FSR Photons
    if (verbose) cout<<"clear other reco variables"<<endl;
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
    if (verbose) cout<<"clear gen event variables"<<endl;
    GENfinalState=-1;
    passedFiducialSelection=false;

    // lepton variables
    if (verbose) cout<<"clear gen lepton variables"<<endl;
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
    if (verbose) cout<<"clear gen higgs/Z variables"<<endl;
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
    if (verbose) cout<<"clear gen jet variables"<<endl;
    if (GENjet_p4->GetLast()!=-1) GENjet_p4->Clear();
    GENnjets_pt30_eta4p7=-1;
    GENpt_leadingjet_pt30_eta4p7=-1.0;
    GENabsrapidity_leadingjet_pt30_eta4p7=-1.0;
    GENabsdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;

    if (verbose) cout<<"clear other variables"<<endl;
    // Resolution
    massErrorUCSD=-1.0; massErrorUCSDCorr=-1.0; massErrorUF=-1.0; massErrorUFCorr=-1.0; massErrorUFADCorr=-1.0;

    // Global variables not stored in tree
    lep_pt.clear(); lep_ptid.clear(); lep_ptindex.clear();
    recoMuons.clear(); recoElectrons.clear();
    RecoFourMuEvent = false; RecoFourEEvent = false;
    RecoTwoETwoMuEvent = false; RecoTwoMuTwoEEvent = false;
    eventPassedPtAndIsoCuts = false;
    foundHiggsCandidate = false;
    foundZ1 = false;
    twoLooseIsoLeptons = false;
    isIsolated = false;
    passedSIP3D = false;
    passedPtCuts = false;
    passedFullSelection = false;
    passedZ4lSelection = false;
    passedQCDcut = false;
    isRecord = false; 

    // ====================== Do Analysis ======================== //

    if (verbose) cout<<"start pileup reweighting"<<endl;
    // PU information
    if(isMC && reweightForPU) {        
        edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);      
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
    eventWeight = pileupWeight;
    if (verbose) cout<<"finished pileup reweighting"<<endl;
    
    nEventsTotal += eventWeight;

    if(isMC) {
        if (verbose) cout<<"setting gen variables"<<endl;       
        setGENVariables(genParticles,genJets); 
        if (verbose) cout<<"finshed setting gen variables"<<endl;       
    }

    //Check for duplicate events in data
    if (verbose) cout<<"checking duplicates"<<endl;       
    std::vector<ULong64_t> runVec, lumiVec, eventVec;
    bool notDuplicateEvent = true;
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
    if (verbose) cout<<"finished checking duplicates"<<endl;       

    if( notDuplicateEvent && !vertex->empty()) {

    
        //N Vertex 
        if (verbose) cout<<"fill nvtx histogram"<<endl;       
        nVtx = vertex->size();
        histContainer_["NVTX"]->Fill(nVtx);
        histContainer_["NVTX_RW"]->Fill(nVtx,pileupWeight);

        //MET
        if (verbose) cout<<"get met value"<<endl;       
        met = mets->empty() ? 0 : (*mets)[0].et();

        if (verbose) cout<<"start lepton analysis"<<endl;           
        vector<pat::Muon> AllMuons;
        vector<pat::Electron> AllElectrons;  
        AllMuons = helper.goodLooseMuons2012(muons,_muPtCut);
        AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut);

        helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);
        recoMuons = helper.goodMuons2012_noIso(AllMuons,_muPtCut,PV);
        recoElectrons = helper.goodElectrons2012_noIso(AllElectrons,_elecPtCut,elecID,PV,iEvent);

        //sort electrons and muons by pt
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
                lep_isoPU.push_back(helper.getPUIso(recoElectrons[lep_ptindex[i]],elRho));
                lep_Sip.push_back(helper.getSIP3D(recoElectrons[lep_ptindex[i]]));           
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
                lep_isoPU.push_back(helper.getPUIso(recoMuons[lep_ptindex[i]],muRho));
                lep_Sip.push_back(helper.getSIP3D(recoMuons[lep_ptindex[i]]));            
                lep_genindex.push_back(-1.0);
            }
            if (verbose) cout<<"( RelIso: "<<lep_RelIso[i]<<" isoCH: "<<lep_isoCH[i]<<" isoNH: "<<lep_isoNH[i]
                             <<" isoPhot: "<<lep_isoPhot[i]<<" isoPU: "<<lep_isoPU[i]<<" Sip: "<<lep_Sip[i]<<endl;
        }

        // GEN matching
        if(isMC) {
            if (verbose) cout<<"begin gen matching"<<endl;
            // for each reco lepton find the nearest status 1 gen lepton with same ID
            for(unsigned int i = 0; i < lep_pt.size(); i++) {

                double minDr=9999.0;

                for (unsigned int j = 0; j < GENlep_id.size(); j++) {

                    if (GENlep_status[j]!=1) continue;
                    if (GENlep_id[j]!=lep_id[i]) continue;

                    TLorentzVector *reco, *gen;
                    reco = (TLorentzVector*) lep_p4->At(i);
                    gen = (TLorentzVector*) GENlep_p4->At(j);
                    double thisDr = deltaR(reco->Eta(),reco->Phi(),gen->Eta(),gen->Phi());

                    if (thisDr<minDr && thisDr<0.5) {
                        lep_genindex[i]=j;
                        minDr=thisDr;
                    }

                } // all gen leptons 
            } // all reco leptons
            if (verbose) cout<<"finished gen matching"<<endl;
        } //isMC

        bool twoLep_ID = false;
        if( (recoMuons.size() + recoElectrons.size()) >= 2 ){twoLep_ID = true;}

        if(twoLep_ID) {
            if (verbose) cout<<"found two leptons"<<endl;
            // Mass Resolution Study
            if(bStudyResolution) {
                if (verbose) cout<<"begin 2lep mass resolution study"<<endl;
                vector<pat::Muon> recoIsoMuons;
                vector<pat::Electron> recoIsoElectrons;	    
                recoIsoMuons = helper.goodMuons2012_Iso(AllMuons,_muPtCut, muRho, isoCut, PV);
                recoIsoElectrons = helper.goodElectrons2012_Iso(AllElectrons,_elecPtCut, elRho, isoCut, elecID,PV);	   
                PerLepReso->fillHistograms(hContainer_, hContainer2D_, hContainer3D_, recoIsoElectrons, recoIsoMuons, eventWeight, !isMC);
                if(bStudyDiLeptonResolution) { DiLepReso->fillHistograms(hContainer_, hContainer2D_, recoIsoElectrons, recoIsoMuons, eventWeight, !isMC); }
                if (verbose) cout<<"finished 2lep mass resolution study"<<endl;
            }

            bool fourLep_ID = false;
            if((recoMuons.size() + recoElectrons.size()) >= 4 ) {fourLep_ID = true;}  

            if(fourLep_ID){
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
                if(properLep_ID){     
                    if (verbose) cout<<"found 2 OSSF lepton pairs"<<endl;

                    // FSR Photons
                    vector<pat::PFParticle> fsrPhotons; 
                    vector<double> deltaRVec;
                    if(doFsrRecovery) {
                        if (verbose) cout<<"filling fsr photon candidates"<<endl;                                    
                        for(edm::View<pat::PFParticle>::const_iterator phot=photonsForFsr->begin(); phot!=photonsForFsr->end(); ++phot) {

                            bool matched = false; double chosenDeltaRPh = 999;

                            for(unsigned int i = 0; i < recoElectrons.size(); i++) {
                                double tmpDeltaREPh = deltaR(recoElectrons[i].eta(), recoElectrons[i].phi(), phot->eta(), phot->phi());
                                double fsrDeltaPhi = fabs(deltaPhi(phot->phi(),recoElectrons[i].phi()));
                                double fsrDeltaEta = fabs(phot->eta()-recoElectrons[i].eta());

                                if( tmpDeltaREPh < 0.15){matched = true;}	
                                if( fsrDeltaPhi < 2 && fsrDeltaEta < 0.05 ){matched = true;}
                                if( tmpDeltaREPh < chosenDeltaRPh ){chosenDeltaRPh = tmpDeltaREPh;}
                            }

                            for(unsigned int i = 0; i < recoMuons.size(); i++) {
                                double tmpDeltaRMPh = deltaR(recoMuons[i].eta(), recoMuons[i].phi(), phot->eta(),phot->phi());
                                if( tmpDeltaRMPh < chosenDeltaRPh ){chosenDeltaRPh = tmpDeltaRMPh;}
                            }
                            
                            if( (phot->pt() > 2 && chosenDeltaRPh < 0.07) || (phot->pt() > 4 && chosenDeltaRPh < 0.5) ) {
                                if(!matched && fabs(phot->eta()) < 2.4){ fsrPhotons.push_back(*phot); deltaRVec.push_back(chosenDeltaRPh);}
                            }

                        } // loop over fsr photon candidates
                        if (verbose) cout<<"finished filling fsr photon candidates"<<endl;            
                    } // doFsrRecovery
                    	           
                    vector<pat::Muon> selectedMuons;
                    vector<pat::Electron> selectedElectrons;
                    vector<pat::PFParticle> selectedFsrPhotons;
                    
                    if (verbose) cout<<"begins looking for higgs candidate"<<endl;                    
                    findHiggsCandidate(lep_p4,fsrPhotons,deltaRVec,lep_Hindex,selectedMuons,selectedElectrons,selectedFsrPhotons,iEvent);
                    if (verbose) cout<<"found higgs candidate? "<<foundHiggsCandidate<<endl;                    
                        
                    //VBF Jets
                    vector<pat::Jet> goodJets;
                    double tempDeltaR = -999;

                    if (verbose) cout<<"begin filling jet candidates"<<endl;                                        
                    for(unsigned int i = 0; i < jets->size(); ++i) {

                        const pat::Jet & patjet = jets->at(i);
                        const pat::Jet & correctedJet = correctedJets->at(i);

                        //int  idflag = (*puJetIdFlag)[jets->refAt(i)];
                        //if( PileupJetIdentifier::passJetId(idflag, PileupJetIdentifier::kLoose) ) {
                        if (1==1) { //FIXME
  
                            //PF ID
                            if (verbose) cout<<"checking jetid..."<<endl;                                        
                            if(jetHelper.patjetID(correctedJet) == 1) {

                                if (verbose) cout<<"checking overlap with fsr photons..."<<endl;                                        
                                bool isDeltaR = true;
                                for(unsigned int phIndex = 0; phIndex < selectedFsrPhotons.size(); phIndex++) {
                                    tempDeltaR = deltaR(patjet.eta(),patjet.phi(),selectedFsrPhotons[phIndex].eta(),selectedFsrPhotons[phIndex].phi());
                                    if (tempDeltaR < 0.5) isDeltaR = false;
                                }

                                if(correctedJet.pt() > pt_cut && fabs(patjet.eta()) < eta_cut && isDeltaR) {
                                    
                                    // apply scale factor for PU Jets by demoting 1-data/MC % of jets jets in certain pt/eta range 
                                    // Configured now that SF is 1.0
                                    if (verbose) cout<<"adding pu jet scale factors..."<<endl;       
                                    bool dropit=false;
                                    if (abs(patjet.eta())>3.0) {
                                        TRandom3 rand;
                                        rand.SetSeed(abs(static_cast<int>(sin(patjet.phi())*100000)));
                                        float coin = rand.Uniform(1.); 
                                        if (correctedJet.pt()>=20.0 && correctedJet.pt()<36.0 && coin>1.0) dropit=true;
                                        if (correctedJet.pt()>=36.0 && correctedJet.pt()<50.0 && coin>1.0) dropit=true;
                                        if (correctedJet.pt()>=50.0 && coin>1.0) dropit=true;
                                    }                                        
                                    
                                    if (!dropit) {
                                        if (verbose) cout<<"adding jet candidate, pt: "<<correctedJet.pt()<<" eta: "<<correctedJet.eta()<<endl;
                                        goodJets.push_back(correctedJet);
                                    } // pu jet scale factor

                                } // pass deltaR jet/fsr photons		    
                            } // pass loose pf jet id
                        } // pass pu jet id
                    } // all jets
      
                    if( foundHiggsCandidate ){
                        
                        if (verbose) cout<<"checking lepton pt cuts"<<endl;  
                        passedPtCuts = helper.passedPtandEarlyIso(selectedMuons, selectedElectrons, leadingPtCut, subleadingPtCut, 1000, "PF", "dB","PFEffAreaRho", muRho);

                        if( passedPtCuts ){
		  

                            // ========= m2l > 4 cut ========= //
                            if (verbose) cout<<"checking m2l cuts"<<endl;  
                            bool fourLep_Cleaned = helper.passedM2lCut_OS(selectedMuons,selectedElectrons,4,10000);
                            passedQCDcut = helper.passedM2lCut_OS(selectedMuons,selectedElectrons,4,10000);

                            if(fourLep_Cleaned){

                                //M4L Error
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
	      
                                //Set All the Variables for Saved Trees --- must be done after variables are available
                                if (verbose) cout<<"begin setting tree variables"<<endl; 
                                setTreeVariables(iEvent, iSetup, selectedMuons, selectedElectrons, recoMuons, recoElectrons, goodJets);
                                if (verbose) cout<<"finshed setting tree variables"<<endl; 

                                //mela::computeAngles(P4s[0], tmpIDs[0], P4s[1], tmpIDs[1], P4s[2], tmpIDs[2], P4s[3], tmpIDs[3], cosThetaStar,cosTheta1,cosTheta2,Phi,Phi1);

                                //Z4lmaxP = helper.largestLepMomentum(L11P4,L12P4,L21P4,L22P4);
                                //vector<double> thetas = angles.angleBetweenLep(L11P4,L12P4,L21P4,L22P4);
                                //theta12 = thetas[0]; theta13 = thetas[1]; theta14 = thetas[2];
                                //thetaPhoton = angles.minAngleOfPhoton(L11P4,L12P4,L21P4,L22P4);
                                //thetaPhotonZ = angles.angleOfPhotonZframe(L11P4,L12P4,L21P4,L22P4,chargeL3);
                                //maxMass2Lep = helper.maxMass2l(recoMuons,recoElectrons); 

                                if(Z2Vec.M() > 0) {
                                    passedZ4lSelection = true;
                                    if(Z2Vec.M() > mZ2Low) passedFullSelection = true;
                                }

                                if(bStudyResolution && bStudyFourLeptonResolution){
                                    if (passedFullSelection) {
                                        FourLepReso->fillHistograms(hContainer_, hContainer2D_, selectedElectrons, selectedMuons, selectedFsrPhotons, eventWeight,!isMC);
                                    }
                                } // bStudyFourLepResolution

                            }//fourLep_Cleaned
                            //else {std::cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed fourLep_Cleaned"<<std::endl;}
                        }//passedPtCuts
                        //else {std::cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed PtCuts"<<std::endl;}
                    }//if HC
                    //else {std::cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed HC"<<std::endl;}
                }// if 4 properID
                //else {std::cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed 4 properID"<<std::endl;}
            }//if 4lepID
            //else {std::cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed 4ID"<<std::endl;}
        if (!isMC) passedEventsTree_All->Fill();		  
        }//if 2lepID
        //else {std::cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed ID"<<std::endl;}     
    }//notDuplicate
    //else {std::cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed notDuplicate"<<std::endl;}
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

}

// ------------ method called once each job just after ending the event loop  ------------
void 
UFHZZ4LAna::endJob() 
{
    histContainer_["NEVENTS"]->SetBinContent(1,nEventsTotal);
    histContainer_["NEVENTS"]->GetXaxis()->SetBinLabel(1,"N Events in Sample");
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
        nEventsTotal += numEventsCounter->value;
    }
}



// ============================ UF Functions ============================= //


//Find Z1,Z2, and Higgs candidate
//Pass good leptons for analysis as candMuons and candElectrons
//Pass empty vectors of pat leptons as selectedMuons and selectedElectrons
// these will be filled in the function and then useable for more analysis.
void
UFHZZ4LAna::findHiggsCandidate(TClonesArray *lep_p4,
                               std::vector< pat::PFParticle > fsrPhotons, std::vector<double> deltaRVec,
                               int lep_Hindex[],
                               std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons,
                               std::vector< pat::PFParticle > &selectedFsrPhotons,const edm::Event& iEvent )
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    double ZmassDiff = 1000.0;
    const double Zmass = 91.1876;
    bool foundZ1=false, foundZ2=false;

    int tmpTakenPhotonZ1 = 999;
    int takenPhotonZ1 = 999;
    int associatedPh1 = 999, associatedPh2 = 999;
    int tmpAssociatedPh = 999;
    TLorentzVector photVecZ1, photVecZ2, tmpPhotVec, tmpZVec;
    bool foundPhotZ1 = false, foundPhotZ2 = false;
    bool foundPhot = false;


    unsigned int N = lep_p4->GetLast()+1;
    if (verbose) cout<<N<<" leptons in total"<<endl;

    // Select Z1 by 2 leptons with same flavor opposite charge
    // with m(2l) closest to mZ
    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<N; j++){

            // same flavor opposite charge
            if((lep_id[i]+lep_id[j])!=0) continue;

            if (verbose) cout<<"OSSF pair: id1="<<lep_id[i]<<" id2="<<lep_id[j]<<endl;    
            TLorentzVector *li, *lj;
            li = (TLorentzVector*) lep_p4->At(i);
            lj = (TLorentzVector*) lep_p4->At(j);

            if (verbose) cout<<"begin Z1 candidate formation"<<endl;
            foundZ1 = findZ(fsrPhotons, deltaRVec, *lj, i, *lj, j, 999, tmpTakenPhotonZ1, tmpAssociatedPh, tmpZVec, tmpPhotVec, foundPhot);
            if (verbose) cout<<"found Z1 candidate?"<<foundZ1<<endl;
            if (!foundZ1) continue;
            
            if (verbose) cout<<"this m(Z1): "<<tmpZVec.M()<<" best m(Z1): "<<Z1Vec.M()<<endl;
            double dm = abs(Zmass-tmpZVec.M());
            if(dm < ZmassDiff ) {
                ZmassDiff = dm;
                Z1Vec = tmpZVec;
                if (li->Pt()>lj->Pt()) {lep_Hindex[0] = i; lep_Hindex[1] = j;}
                else {lep_Hindex[0] = j; lep_Hindex[1] = i;}                
                if (verbose) cout<<"fsr? "<<foundPhotZ1<<" new lep_Hindex[0]="<<lep_Hindex[0]<<" lepHindex[1]="<<lep_Hindex[1]<<endl;
                foundPhotZ1 = foundPhot;
                if (foundPhotZ1) {
                    photVecZ1 = tmpPhotVec;                
                    if(tmpAssociatedPh == 1) associatedPh1 = i;
                    if(tmpAssociatedPh == 2) associatedPh1 = j;
                    takenPhotonZ1=tmpTakenPhotonZ1;
                } // better Z1 candidate has FSR photon
            } // better Z1 candidate
        } // lepton i
    } // lepton j

    if (!foundZ1) return;

    if (verbose) cout<<"best Z1: lep_Hindex[0]="<<lep_Hindex[0]<<" lep_Hindex[1]="<<lep_Hindex[1]<<endl;
    TLorentzVector *lep_1, *lep_2;
    lep_1 = (TLorentzVector*) lep_p4->At(lep_Hindex[0]);
    lep_2 = (TLorentzVector*) lep_p4->At(lep_Hindex[1]);

    if( foundPhotZ1 ) {
        FSR_Z1 = true; 
        nFSRPhotons++;
        new ( (*phofsr_p4)[nFSRPhotons-1] ) TLorentzVector(photVecZ1.Px(),photVecZ1.Py(), photVecZ1.Pz(), photVecZ1.E());
        phofsr_lepindex.push_back(associatedPh1);
        TLorentzVector *pho;
        pho = (TLorentzVector*) phofsr_p4->At(nFSRPhotons-1);
        TLorentzVector lepfsr;
        if (associatedPh1==lep_Hindex[0]) lepfsr = (*lep_1)+(*pho);
        else if (associatedPh1==lep_Hindex[1]) lepfsr = (*lep_2)+(*pho);
        new ( (*lep_p4_FSR)[associatedPh1] ) TLorentzVector(lepfsr.Px(),lepfsr.Py(),lepfsr.Pz(),lepfsr.E());

    }

    // Z2 selection
    double sumPtZ2 = 0, sumPtZ2_tmp = 0;
    tmpAssociatedPh = 999;
    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<N; j++){

            // same flavor opposite charge
            if((lep_id[i]+lep_id[j])!=0) continue;
            // not Z1 leptons
            if( (int)i == lep_Hindex[0] || (int)i == lep_Hindex[1] || (int)j == lep_Hindex[0] || (int)j == lep_Hindex[1] ) continue;

            TLorentzVector *li, *lj;
            li = (TLorentzVector*) lep_p4->At(i);
            lj = (TLorentzVector*) lep_p4->At(j);

            double tmpDeltaR=999.0;
            tmpDeltaR = deltaR(li->Eta(),li->Phi(),lep_1->Eta(),lep_1->Phi());
            if(tmpDeltaR < 0.02) continue;
            tmpDeltaR = deltaR(li->Eta(),li->Phi(),lep_2->Eta(),lep_2->Phi());
            if(tmpDeltaR < 0.02) continue;
            tmpDeltaR = deltaR(lj->Eta(),lj->Phi(),lep_1->Eta(),lep_1->Phi());
            if(tmpDeltaR < 0.02) continue;
            tmpDeltaR = deltaR(lj->Eta(),lj->Phi(),lep_2->Eta(),lep_2->Phi());
            if(tmpDeltaR < 0.02) continue;
                    
            // tmpTakenPhotonZ1 now is just a dummy, needed so we can use the same function findZ() for both Z1 and Z2
            foundZ2 = findZ(fsrPhotons, deltaRVec, *li, i, *lj, j, takenPhotonZ1, tmpTakenPhotonZ1, tmpAssociatedPh, tmpZVec, tmpPhotVec, foundPhot);
            sumPtZ2_tmp = li->Pt() + lj->Pt();
            
            if(sumPtZ2_tmp > sumPtZ2 && foundZ2) {
                sumPtZ2 = sumPtZ2_tmp;
                Z2Vec = tmpZVec;
                if (li->Pt() > lj->Pt()){lep_Hindex[2] = (int)i; lep_Hindex[3] = (int)j;}
                else {lep_Hindex[2] = (int)j; lep_Hindex[3] = (int)i;}
                foundPhotZ2 = foundPhot;
                if (foundPhotZ2) {
                    photVecZ2 = tmpPhotVec;
                    if(tmpAssociatedPh == 1) associatedPh2 = i;
                    if(tmpAssociatedPh == 2) associatedPh2 = j;
                } // better Z2 has fsr photon
            } // better Z2
        } // lepton j
    } // lepton i
        
    if (!foundZ2) return;

    TLorentzVector *lep_3, *lep_4;
    lep_3 = (TLorentzVector*) lep_p4->At(lep_Hindex[2]);
    lep_4 = (TLorentzVector*) lep_p4->At(lep_Hindex[3]);

    if( foundPhotZ2 ) {
        FSR_Z2 = true;
        nFSRPhotons++;
        new ( (*phofsr_p4)[nFSRPhotons-1] ) TLorentzVector(photVecZ2.Px(),photVecZ2.Py(), photVecZ2.Pz(), photVecZ2.E());
        phofsr_lepindex.push_back(associatedPh2);
        TLorentzVector *pho;
        pho = (TLorentzVector*) phofsr_p4->At(nFSRPhotons-1);
        TLorentzVector lepfsr;
        if (associatedPh2==lep_Hindex[2]) lepfsr = (*lep_3)+(*pho);
        else if (associatedPh2==lep_Hindex[3]) lepfsr = (*lep_4)+(*pho);
        new ( (*lep_p4_FSR)[associatedPh2] ) TLorentzVector(lepfsr.Px(),lepfsr.Py(),lepfsr.Pz(),lepfsr.E());
    }

    //Determine whether a Higgs candidate was formed
    massZ1 = Z1Vec.M();
    massZ2 = Z2Vec.M();
    HVec = Z1Vec + Z2Vec;
    mass4l = HVec.M();

    if( massZ1 > mZ1Low && massZ1 < mZ1High && massZ2 > mZ2Low && massZ2 < mZ2High) {
        
        foundHiggsCandidate = true;

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
    tree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    tree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    tree->Branch("passedQCDcut",&passedQCDcut,"passedQCDcut/O");
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
    tree->Branch("lep_genindex",&lep_genindex);
    tree->Branch("lep_missingHits",&lep_missingHits);
    tree->Branch("lep_mva",&lep_mva);
    tree->Branch("lep_Sip",&lep_Sip);
    tree->Branch("lep_IP",&lep_IP);
    tree->Branch("lep_dIPL1",&lep_dIP);
    tree->Branch("lep_isoNH",&lep_isoNH);
    tree->Branch("lep_isoCH",&lep_isoCH);
    tree->Branch("lep_isoPhot",&lep_isoPhot);
    tree->Branch("lep_RelIso",&lep_RelIso);
    tree->Branch("muRho",&muRho,"muRho/D");
    tree->Branch("elRho",&elRho,"elRho/D");

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
    tree->Branch("massZ1",&massZ1,"massZ1/D");
    tree->Branch("massZ2",&massZ2,"massZ2/D");  
    tree->Branch("pTZ1",&pTZ1,"pTZ1/D");
    tree->Branch("pTZ2",&pTZ2,"pTZ2/D");

    // MET
    met_p4 = new TClonesArray("TLorentzVector", 10);
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

    tree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
    tree->Branch("njets_pt30_eta4p7_jesup",&njets_pt30_eta4p7_jesup,"njets_pt30_eta4p7_jesup/I");
    tree->Branch("njets_pt30_eta4p7_jesdn",&njets_pt30_eta4p7_jesdn,"njets_pt30_eta4p7_jesdn/I");
    tree->Branch("njets_pt30_eta4p7_jerup",&njets_pt30_eta4p7_jerup,"njets_pt30_eta4p7_jerup/I");
    tree->Branch("njets_pt30_eta4p7_jerdn",&njets_pt30_eta4p7_jerdn,"njets_pt30_eta4p7_jerdn/I");

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
    tree->Branch("DijetFisherDiscrim",&DijetFisherDiscrim,"DijetFisherDiscrim/D");

    // FSR Photons
    phofsr_p4 = new TClonesArray("TLorentzVector", 20);
    tree->Branch("phofsr_p4","TClonesArray", &phofsr_p4, 128000, 0);    
    tree->Branch("nFSRPhotons",&nFSRPhotons,"nFSRPhotons/I");
    tree->Branch("FSR_Z1",&FSR_Z1,"FSR_Z1/O");
    tree->Branch("FSR_Z2",&FSR_Z1,"FSR_Z2/O");
    tree->Branch("phofsr_lepindex",&phofsr_lepindex);

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

}



void UFHZZ4LAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                   std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, std::vector<pat::Jet> goodJets)
{

    using namespace edm;
    using namespace pat;
    using namespace std;

    // Tree Variables
    if( RecoFourMuEvent ){ finalState = 1;}
    if( RecoFourEEvent  ){ finalState = 2;}
    if( RecoTwoETwoMuEvent ){ finalState = 3;}
    if( RecoTwoMuTwoEEvent ){ finalState = 4;}

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

    double tempDeltaR = 999.0;
    vector<pat::Jet> finalVBFJets;

    for( unsigned int k = 0; k < goodJets.size(); k++) {

        bool isDeltaR_eta4p7 = true;

        // check overlap with isolated leptons
        for(unsigned int i = 0; i < recoMuons.size(); i++) {                       
            if (helper.pfIso(recoMuons[i],muRho)>0.4) continue;
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),recoMuons[i].eta(),recoMuons[i].phi());
            if (tempDeltaR<0.5) {
                isDeltaR_eta4p7 = false;
            }
        }
        for(unsigned int i = 0; i < recoElectrons.size(); i++) {
            if (helper.pfIso(recoElectrons[i],elRho)>0.4) continue;
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),recoElectrons[i].eta(),recoElectrons[i].phi());
            if (tempDeltaR<0.5) {
                isDeltaR_eta4p7 = false;
            }
        }
        // check overlap with fsr photons
        unsigned int N = phofsr_p4->GetLast()+1;
        for(unsigned int i=0; i<N; i++) {
            TLorentzVector *pho;
            pho = (TLorentzVector*) phofsr_p4->At(i);
            tempDeltaR=999.0;
            tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),pho->Eta(),pho->Phi());
            if (tempDeltaR<0.5) {
                isDeltaR_eta4p7 = false;
            }
        }


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
                finalVBFJets.push_back(goodJets[k]);
                new ( (*jet_p4)[njets_pt30_eta4p7-1] ) TLorentzVector(jet_jer->Px(), jet_jer->Py(), jet_jer->Pz(), jet_jer->Energy());
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
    
    if(njets_pt30_eta4p7>1){
        TLorentzVector *jet1, *jet2;
        jet1 = (TLorentzVector*) jet_p4->At(0);
        jet2 = (TLorentzVector*) jet_p4->At(1);
        TLorentzVector Dijet;
        Dijet = (*jet1)+(*jet2); 
        DijetMass = Dijet.M();
        DijetDEta = fabs(jet1->Eta()-jet2->Eta());
        // OLD MORIOND --- FisherDiscrim = 0.09407*fabs(VBFDeltaEta) + 4.1581e-4*VBFDiJetMass;
        DijetFisherDiscrim = 0.18*fabs(DijetDEta) + 1.92e-4*DijetMass;
    } 
 
}




void UFHZZ4LAna::setGENVariables(edm::Handle<reco::GenParticleCollection> genParticles,
                                 edm::Handle<edm::View<reco::GenJet> > genJets)
{

    reco::GenParticleCollection::const_iterator genPart;
    int j = -1;
    int nGENLeptons=-1;

    if (verbose) cout<<"begin looping on gen particles"<<endl;
    for(genPart = genParticles->begin(); genPart != genParticles->end(); genPart++) {
        j++;
        if (abs(genPart->pdgId())==11  || abs(genPart->pdgId())==13 || abs(genPart->pdgId())==15) {
            if (genPart->pt()<3.0) continue;
            if (abs(genPart->eta())>2.7) continue;
            nGENLeptons++;
            new ( (*GENlep_p4)[nGENLeptons] ) TLorentzVector(genPart->px(),genPart->py(),genPart->pz(),genPart->energy());
            GENlep_id.push_back( genPart->pdgId() );
            GENlep_status.push_back(genPart->status());
            GENlep_MomId.push_back(genAna.MotherID(&genParticles->at(j)));
            GENlep_MomMomId.push_back(genAna.MotherMotherID(&genParticles->at(j)));
            TLorentzVector *thisLep;
            thisLep = (TLorentzVector*) GENlep_p4->At(nGENLeptons);
            // GEN iso calculation
            reco::GenParticleCollection::const_iterator OtherGenPart;
            int k = -1;
            double this_GENiso=0.0;
            double this_GENneutraliso=0.0;
            double this_GENchargediso=0.0;
            for(OtherGenPart = genParticles->begin(); OtherGenPart != genParticles->end(); OtherGenPart++) {
                k++;
                if (k==j) continue;
                if( OtherGenPart->status() != 1 ) continue; // stable particles only         
                if (abs(OtherGenPart->pdgId())==12 || abs(OtherGenPart->pdgId())==14 || abs(OtherGenPart->pdgId())==16) continue;
                if ( (abs(OtherGenPart->pdgId())==11 || abs(OtherGenPart->pdgId())==13) && genAna.IsMotherStatus3EorMu(&genParticles->at(k))) continue;
                if (OtherGenPart->pdgId()==22 && genAna.IsFSR(&genParticles->at(k),GENlep_id[j])) continue;
                double this_dRvL = deltaR(thisLep->Eta(), thisLep->Phi(), OtherGenPart->eta(), OtherGenPart->phi());
                if(this_dRvL<0.4) {
                    this_GENiso = this_GENiso + OtherGenPart->pt();
                    if (genPart->charge()==0) this_GENneutraliso = this_GENneutraliso + OtherGenPart->pt();
                    if (genPart->charge()!=0) this_GENchargediso = this_GENchargediso + OtherGenPart->pt();
                }
            } // GEN iso loop
            this_GENiso = this_GENiso/genPart->pt();
            GENlep_RelIso.push_back(this_GENiso);
            // END GEN iso calculation
        } //leptons

        if (genPart->pdgId()==25) {
            GENMH=genPart->mass();
        }

        if (genPart->pdgId()==23) {
            const reco::Candidate *Zdau0=genPart->daughter(0);
            if (Zdau0) GENZ_DaughtersId.push_back(fabs(Zdau0->pdgId()));
            const reco::Candidate *Zmom0=genPart->mother(0);
            if (Zmom0) GENZ_MomId.push_back(Zmom0->pdgId());
                
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
                    if (GENlep_status[i]!=1) continue;
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
        if(GENlep_status[i]!=1) continue;
        for(unsigned int j=i+1; j<N; j++){

            if((GENlep_id[i]+GENlep_id[j])!=0) continue;
            if(GENlep_status[j]!=1) continue;

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
        if(GENlep_status[i]!=1) continue;
        for(unsigned int j=i+1; j<N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;            
            if(GENlep_status[j]!=1) continue;

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
                    //cout<<"mZ2: "<<mZ2<<endl;
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


////////// Find Z ////////////
bool UFHZZ4LAna::findZ(std::vector<pat::PFParticle> photons, std::vector<double> &deltaRVec, TLorentzVector &l1, unsigned int &index1, TLorentzVector &l2, unsigned int &index2, int takenZ1, int taken, int assocLep, TLorentzVector &ZVec, TLorentzVector &photVec, bool &foundPhoton)
{

    using namespace std;
    using namespace pat;
    using namespace reco;

    double deltaR1, deltaR2, smallestDeltaR, totalSmallestDeltaR = 999;
    double photHighestPt = 0;
    TLorentzVector mll, mllgam;
    bool foundPhot = false;
    double massDiffPhot = 999, massDiffNoPhot = 999;
    double coneSize = 0.4;
    double lep1Iso = 999, lep2Iso = 999;
    double assocLepTmp = 999;
    double photIsoCut = 1.0;
    bool foundZ = false;

    if( !photons.empty() && photons.size() > 0 && doFsrRecovery) {

        for(unsigned int i = 0; i < photons.size(); i++) {

            if( takenZ1 == (int)i ) continue;

            //pt, eta checks
            if( photons[i].pt() < 4 ) continue;
            if( photons[i].eta() > 2.4 ) continue;
            if( (photons[i].userFloat("fsrPhotonPFIsoChHad03pt02")+photons[i].userFloat("fsrPhotonPFIsoNHad03")
                 +photons[i].userFloat("fsrPhotonPFIsoPhoton03")
                 +photons[i].userFloat("fsrPhotonPFIsoChHadPU03pt02"))/photons[i].pt() > photIsoCut) continue;

            //calc both deltaRs
            deltaR1 = deltaR(l1.Eta(),l1.Phi(),photons[i].eta(),photons[i].phi());
            deltaR2 = deltaR(l2.Eta(),l2.Phi(),photons[i].eta(),photons[i].phi());

            //associate with closest lepton
            if( deltaR1 < deltaR2 ){ assocLepTmp = 1; smallestDeltaR = deltaR1;}
            else{ assocLepTmp = 2; smallestDeltaR = deltaR2;}
            if( smallestDeltaR > 0.5 || smallestDeltaR > deltaRVec[i] ) continue;
            if( photons[i].pt() < photHighestPt ) continue;

            //calc P vectors
            TLorentzVector phoi;
            phoi.SetPtEtaPhiE(photons[i].pt(),photons[i].eta(),photons[i].phi(),photons[i].energy());
            mllgam = l1+l2+phoi;
            mll = l1+l2;

            //check inv mass
            if( mllgam.M() < 4 || mllgam.M() > 100) continue;
            massDiffPhot = fabs(mllgam.M() - Zmass);
            massDiffNoPhot = fabs(mll.M() - Zmass);

            //if its smaller with phot, keep phot
            if( massDiffPhot < massDiffNoPhot )
            {
                //check iso cone
                if( deltaR1 < coneSize && deltaR1 > 1e-06){
                    lep1Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index1]+lep_isoPhot[index1]-lep_isoPU[index1]-photons[i].pt(),0.0))/l1.Pt();
                }
                else {
                    lep1Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index1]+lep_isoPhot[index1]-lep_isoPU[index1],0.0))/l1.Pt();
                }

                if( deltaR2 < coneSize && deltaR2 > 1e-06){
                    lep2Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index2]+lep_isoPhot[index2]-lep_isoPU[index2]-photons[i].pt(),0.0))/l2.Pt();
                }
                else {
                    lep2Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index2]+lep_isoPhot[index2]-lep_isoPU[index2]-photons[i].pt(),0.0))/l2.Pt();
                }

                if (lep1Iso < isoCut && lep2Iso < isoCut) {
                    foundPhot = true;
                    taken = (int)i;
                    photHighestPt = photons[i].pt();                    
                    photVec.SetPtEtaPhiE(photons[i].pt(),photons[i].eta(),photons[i].phi(),photons[i].energy());
                    ZVec = mllgam;
                    assocLep = assocLepTmp;
                    foundZ = true;
                }
            }
        }

        if(!foundPhot) {

            for(unsigned int i = 0; i < photons.size(); i++) {
                //FIXME not sure about the useDR and usePT
                bool useDR = false, usePT = false;
                if( takenZ1 == (int)i ) continue;
                //pt, eta checks
                if( photons[i].pt() < 2 ) continue;
                if( photons[i].eta() > 2.4 ) continue;
                if( photons[i].pt() < 4 ) useDR = true;
                if( photons[i].pt() > 4 ) usePT = true;
                if( usePT && photons[i].pt() < photHighestPt) continue;
                //calc both deltaRs
                deltaR1 = deltaR(l1.Eta(),l1.Phi(),photons[i].eta(),photons[i].phi());
                deltaR2 = deltaR(l2.Eta(),l2.Phi(),photons[i].eta(),photons[i].phi());
                //associate with closest lepton
                if( deltaR1 < deltaR2 ){ assocLepTmp = 1; smallestDeltaR = deltaR1;}
                else{ assocLepTmp = 2; smallestDeltaR = deltaR2;}
                if( smallestDeltaR > 0.07  || smallestDeltaR > deltaRVec[i] ) continue;
                if( smallestDeltaR > totalSmallestDeltaR && useDR ) continue;
                //calc P vectors
                TLorentzVector phoi;
                phoi.SetPtEtaPhiE(photons[i].pt(),photons[i].eta(),photons[i].phi(),photons[i].energy());
                mllgam = l1+l2+phoi;
                mll = l1+l2;

                //check inv mass
                if( mllgam.M() < 4 || mllgam.M() > 100) continue;
                massDiffPhot = fabs(mllgam.M() - Zmass);
                massDiffNoPhot = fabs(mll.M() - Zmass);

                //if its smaller with phot, keep phot
                if( massDiffPhot < massDiffNoPhot ) {
                    if( deltaR1 < coneSize && deltaR1 > 1e-06){
                        lep1Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index1]+lep_isoPhot[index1]-lep_isoPU[index1]-photons[i].pt(),0.0))/l1.Pt();
                    }
                    else {
                        lep1Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index1]+lep_isoPhot[index1]-lep_isoPU[index1],0.0))/l1.Pt();
                    }
                    if( deltaR2 < coneSize && deltaR2 > 1e-06){
                        lep2Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index2]+lep_isoPhot[index2]-lep_isoPU[index2]-photons[i].pt(),0.0))/l2.Pt();
                    }
                    else {
                        lep2Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index2]+lep_isoPhot[index2]-lep_isoPU[index2]-photons[i].pt(),0.0))/l2.Pt();
                    }
                    if(lep1Iso < isoCut && lep2Iso < isoCut) {
                        foundPhot = true;
                        taken = (int)i;
                        if(useDR) totalSmallestDeltaR = smallestDeltaR;
                        if(usePT) photHighestPt = photons[i].pt();
                        photVec.SetPtEtaPhiE(photons[i].pt(),photons[i].eta(),photons[i].phi(),photons[i].energy());
                        ZVec = mllgam;
                        assocLep = assocLepTmp;
                        foundZ = true;
                    }            
                }
            }
        }

        if(!foundPhot) {

            lep1Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index1]+lep_isoPhot[index1]-lep_isoPU[index1],0.0))/l1.Pt();
            lep1Iso = (lep_isoCH[index2]+std::max(lep_isoNH[index2]+lep_isoPhot[index2]-lep_isoPU[index2],0.0))/l2.Pt();

            mll = l1+l2;

            if( lep1Iso < isoCut && lep2Iso < isoCut ) {
                ZVec = l1+l2;
                assocLep = 999;
                foundPhot = false;
                foundZ = true;
            }

        }
    }
    else {

        lep1Iso = (lep_isoCH[index1]+std::max(lep_isoNH[index1]+lep_isoPhot[index1]-lep_isoPU[index1],0.0))/l1.Pt();
        lep1Iso = (lep_isoCH[index2]+std::max(lep_isoNH[index2]+lep_isoPhot[index2]-lep_isoPU[index2],0.0))/l2.Pt();

        mll = l1+l2;

        if( lep1Iso < isoCut && lep2Iso < isoCut ) {
            ZVec = l1+l2;
            assocLep = 999;
            foundPhot = false;
            foundZ = true;
        }

    } 
  
    foundPhoton = foundPhot;
    return foundZ;
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
