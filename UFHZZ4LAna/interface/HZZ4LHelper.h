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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

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
  
    std::vector<pat::Muon> goodLooseMuons2012(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut);
    std::vector<pat::Electron> goodLooseElectrons2012(edm::Handle<edm::View<pat::Electron> > Electrons, double elPtCut);

    std::vector<pat::Muon> goodMuons2015_noIso_noPf(std::vector<pat::Muon> Muons, double muPtCut, const reco::Vertex *&vertex, double sip3dCut);
    std::vector<pat::Electron> goodElectrons2015_noIso_noBdt(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID, const reco::Vertex *&vertex,const edm::Event& iEvent, double sip3dCut); 

    void cleanOverlappingLeptons(std::vector<pat::Muon> &Muons, std::vector<pat::Electron> &Electrons,const reco::Vertex *&vertex);

    double getSIP3D(pat::Muon muon);
    double getSIP3D(pat::Electron electron);
    double getPUIso(pat::Muon muon, double Rho);
    double getPUIso03(pat::Muon muon, double Rho);
    double getPUIso(pat::Electron elec, double Rho);
    double getPUIso03(pat::Electron elec, double Rho);
    double pfIso(pat::Muon muon, double Rho);
    double pfIso03(pat::Muon muon, double Rho);
    double pfIso(pat::Electron elec, double Rho);
    double pfIso03(pat::Electron elec, double Rho);
    double miniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands, const reco::Candidate* ptcl, double r_iso_min, double r_iso_max, double kt_scale, double rho, bool charged_only);
    double ptRatio(pat::Electron electron, edm::Handle<edm::View<pat::Jet> > jets,bool isMC);
    double ptRatio(pat::Muon muon, edm::Handle<edm::View<pat::Jet> > jets,bool isMC);
    double ptRel(pat::Electron electron, edm::Handle<edm::View<pat::Jet> > jets,bool isMC);
    double ptRel(pat::Muon muon, edm::Handle<edm::View<pat::Jet> > jets,bool isMC);

    bool passTight_BDT_Id(pat::Electron electron, std::string elecID);
    bool passTight_Id_SUS(pat::Muon muon, const reco::Vertex *&vertex);
    bool passTight_Id_SUS(pat::Electron electron, std::string elecID, const reco::Vertex *&vertex, const reco::BeamSpot BS, edm::Handle< std::vector<reco::Conversion> > theConversions);
    
    float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState);
    float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState);
    float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState);

    float dataMC(pat::Electron electron, TH2F* hElecScaleFac, TH2F* hElecScalFac_Cracks);
    float dataMCErr(pat::Electron electron, TH2F* hElecScaleFac, TH2F* hElecScalFac_Cracks);
    float dataMC(pat::Muon muon, TH2F* hElecScaleFac, TH2F* hElecScalFac_Cracks);
    float dataMCErr(pat::Muon muon, TH2F* hElecScaleFac, TH2F* hElecScalFac_Cracks);
    
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
        kElEAData2015,
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

};

#endif


#ifndef HZZ4LHELPER_CC
#define HZZ4LHELPER_CC

HZZ4LHelper::HZZ4LHelper()
{
    muEAtype = kMuGammaAndNeutralHadronIso04;
    muEAtarget = kMuEAData2012;
    elEAtype = kEGammaNeutralHadIso04;
    elEAtarget = kElEAData2015;
    //declarations
}


HZZ4LHelper::~HZZ4LHelper()
{
    //destructor ---do nothing
}


std::vector<pat::Muon> HZZ4LHelper::goodLooseMuons2012(edm::Handle<edm::View<pat::Muon> > Muons, double muPtCut) {
    using namespace pat;
    using namespace std;    
    vector<pat::Muon> bestMuons;    
    for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu != Muons->end(); ++mu) {
        if( (mu->isGlobalMuon() || mu->isTrackerMuon() || mu->isPFMuon()) && fabs(mu->eta()) < 2.4 && mu->pt() > muPtCut) {
            //cout << "dB: " << mu->userIsolation("PfPUChargedHadronIso") << endl;
            //double PUCorrVal = 0.5*mu->userIsolation("PfPUChargedHadronIso");
            //double tmpIso = (mu->chargedHadronIso()+max(mu->photonIso()+mu->neutralHadronIso()-PUCorrVal,0.0))/mu->pt();
            //cout << tmpIso << endl;          
            bestMuons.push_back(*mu);
        }
    }
    return bestMuons;    
}

std::vector<pat::Electron> HZZ4LHelper::goodLooseElectrons2012(edm::Handle<edm::View<pat::Electron> > Electrons, double elPtCut) {
    using namespace pat;
    using namespace std;    
    vector<pat::Electron> bestElectrons;    
    for(edm::View<pat::Electron>::const_iterator elec=Electrons->begin(); elec!=Electrons->end(); ++elec) {        
        if( abs(elec->eta()) < 2.5 && elec->pt() > elPtCut) {
            bestElectrons.push_back(*elec);
        }
    }
    return bestElectrons;    
}

std::vector<pat::Muon> HZZ4LHelper::goodMuons2015_noIso_noPf(std::vector<pat::Muon> Muons, double muPtCut, const reco::Vertex *&vertex, double sip3dCut)
{
    //using namespace edm;
    using namespace pat;
    using namespace std;
    vector<pat::Muon> bestMuons;
    /********** M U O N  C U T S **********/
    double muEtaCut = 2.4;
    double dxyCut = 0.5;
    double dzCut = 1;
    /**************************************/
    for(unsigned int i = 0; i < Muons.size(); i++)
        if( Muons[i].pt() > muPtCut && abs(Muons[i].eta()) < muEtaCut &&
            (Muons[i].isGlobalMuon() || (Muons[i].isTrackerMuon() && Muons[i].numberOfMatches() > 0 ) ) &&
            Muons[i].muonBestTrackType() != 2 ) {

            if( abs(getSIP3D(Muons[i])) < sip3dCut ) {
                if( fabs(Muons[i].muonBestTrack()->dxy(vertex->position())) < dxyCut ) { //miniAOD 
                    if( fabs(Muons[i].muonBestTrack()->dz(vertex->position())) < dzCut ) {// miniAOD       

                        bestMuons.push_back(Muons[i]);
                    } //else {cout<<"muon "<<i<<" failed dz cut, |dz|="<<fabs(Muons[i].muonBestTrack()->dxy(vertex->position()))<<endl;}
                } //else {cout<<"muon "<<i<<" failed dxy cut, |dxz|="<<fabs(Muons[i].muonBestTrack()->dz(vertex->position()))<<endl;}
            } //else {cout<<"muon "<<i<<" failed sip cut, |sip|="<<abs(getSIP3D(Muons[i]))<<endl;}
        } //else {cout<<"muon "<<i<<" failed pt, eta, or (isGlobal || isTracker)"<<endl;}

    return bestMuons;

}

std::vector<pat::Electron> HZZ4LHelper::goodElectrons2015_noIso_noBdt(std::vector<pat::Electron> Electrons, double elecPtCut, std::string elecID, const reco::Vertex *&vertex,const edm::Event& iEvent, double sip3dCut)
{
    using namespace edm;
    using namespace pat;
    using namespace std;

    vector<pat::Electron> bestElectrons;
  
    /****** E L E C T R O N  C U T S ******/
    int missingHitsCuts = 99999;
    double dxyCut = 0.5;
    double dzCut = 1;

    for(unsigned int i = 0; i < Electrons.size(); i++) {

        if( abs(getSIP3D(Electrons[i])) < sip3dCut) {
            if (fabs(Electrons[i].gsfTrack()->dxy(vertex->position())) < dxyCut) {
                if (fabs(Electrons[i].gsfTrack()->dz(vertex->position())) < dzCut ) {                  
                    int misHits = Electrons[i].gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS); // for miniAOD
                    if (misHits < missingHitsCuts) { bestElectrons.push_back(Electrons[i]);}                  
                }
            }
        }
    }
  
    return bestElectrons;
}


double HZZ4LHelper::miniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands, const reco::Candidate* ptcl, double r_iso_min, double r_iso_max, double kt_scale, double rho, bool charged_only) {

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    double EA=0.0;
    if(ptcl->isElectron()) {
        if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
        if (fabs(ptcl->eta())>=0.0   && fabs(ptcl->eta())<1.0)   EA=0.1752;
        if (fabs(ptcl->eta())>=1.0   && fabs(ptcl->eta())<1.479) EA=0.1862;
        if (fabs(ptcl->eta())>=1.479 && fabs(ptcl->eta())<2.0)   EA=0.1411;
        if (fabs(ptcl->eta())>=2.0   && fabs(ptcl->eta())<2.2)   EA=0.1534;
        if (fabs(ptcl->eta())>=2.2   && fabs(ptcl->eta())<2.3)   EA=0.1903;
        if (fabs(ptcl->eta())>=2.3   && fabs(ptcl->eta())<2.4)   EA=0.2243;
        if (fabs(ptcl->eta())>=2.4   && fabs(ptcl->eta())<2.5)   EA=0.2687;
    } else if(ptcl->isMuon()) {
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
        if (fabs(ptcl->eta())>=0.0 && fabs(ptcl->eta())<0.8) EA=0.0735;
        if (fabs(ptcl->eta())>=0.8 && fabs(ptcl->eta())<1.3) EA=0.0619;
        if (fabs(ptcl->eta())>=1.3 && fabs(ptcl->eta())<2.0) EA=0.0465;
        if (fabs(ptcl->eta())>=2.0 && fabs(ptcl->eta())<2.2) EA=0.0433;
        if (fabs(ptcl->eta())>=2.2 && fabs(ptcl->eta())<2.5) EA=0.0577;
    }

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if (ptcl->isElectron()) ptThresh = 0;
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));

    for (const pat::PackedCandidate &pfc : *pfcands) {
        if (abs(pfc.pdgId())<7) continue;

        double dr = deltaR(pfc, *ptcl);
        if (dr > r_iso) continue;

        //////////////////  NEUTRALS  /////////////////////////
        if (pfc.charge()==0){
            if (pfc.pt()>ptThresh) {
                /////////// PHOTONS ////////////
                if (abs(pfc.pdgId())==22) {
                    if(dr < deadcone_ph) continue;
                    iso_ph += pfc.pt();
                    /////////// NEUTRAL HADRONS ////////////
                } else if (abs(pfc.pdgId())==130) {
                    if(dr < deadcone_nh) continue;
                    iso_nh += pfc.pt();
                }
            }
            //////////////////  CHARGED from PV  /////////////////////////
        } else if (pfc.fromPV()>1){
            if (abs(pfc.pdgId())==211) {
                if(dr < deadcone_ch) continue;
                iso_ch += pfc.pt();
            }
            //////////////////  CHARGED from PU  /////////////////////////
        } else {
            if (pfc.pt()>ptThresh){
                if(dr < deadcone_pu) continue;
                iso_pu += pfc.pt();
            }
        }
    }
    double iso(0.);
    if (charged_only){
        iso = iso_ch;
    } else {
        iso = iso_ph + iso_nh;
        iso -= rho*EA*(r_iso/0.3)*(r_iso/0.3);
        if (iso>0) iso += iso_ch;
        else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;

}


void HZZ4LHelper::cleanOverlappingLeptons(std::vector<pat::Muon> &Muons, std::vector<pat::Electron> &Electrons,const reco::Vertex *&vertex) {
  
    using namespace pat;
    using namespace std;
    /********** M U O N  C U T S **********/
    double muEtaCut = 2.4;
    double dxyCut = 0.5;
    double dzCut = 1;
    double tmpDeltR =100;
    double sip3dCut = 4;
    double muPtCut = 5;
    /**************************************/
    // changed above cuts by Ahmad

    for( unsigned int i = 0; i < Muons.size(); i++ ) {
        if( Muons[i].pt() > muPtCut && abs(Muons[i].eta()) < muEtaCut &&
            Muons[i].isPFMuon() == 1 &&
            (Muons[i].isGlobalMuon() || (Muons[i].isTrackerMuon() && Muons[i].numberOfMatches(/*reco::Muon::SegmentArbitration*/) > 0 /*numberOfMatchedStations() > 0*/ ) ) &&
            Muons[i].muonBestTrackType() != 2 ) {
            if( abs(getSIP3D(Muons[i])) < sip3dCut ) {	       
                if( fabs(Muons[i].muonBestTrack()->dxy(vertex->position())) < dxyCut ) { //miniAOD                     
                    if( fabs(Muons[i].muonBestTrack()->dz(vertex->position())) < dzCut ) {// miniAOD                        
                        for( unsigned int j = 0; j < Electrons.size(); j++ ) {
                            tmpDeltR = deltaR(Muons[i].eta(),Muons[i].phi(),Electrons[j].eta(),Electrons[j].phi());
                            if( tmpDeltR < 0.05 ) {		
                                Electrons.erase(Electrons.begin()+j);
                            }
                        } 
                    }
                }
            }
        }   
    }    
}

double HZZ4LHelper::getSIP3D(pat::Muon muon) {
    using namespace pat;
    using namespace std;
    
    double ip = fabs(muon.dB(pat::Muon::PV3D));
    double ipError = muon.edB(pat::Muon::PV3D);
    
    double sip = ip/ipError;
    
    return sip;
    
}

double HZZ4LHelper::getSIP3D(pat::Electron electron) {
    using namespace pat;
    using namespace std;
    
    double ip = fabs(electron.dB(pat::Electron::PV3D));
    double ipError = electron.edB(pat::Electron::PV3D);
    
    double sip = ip/ipError;
    
    return sip;    
}

bool HZZ4LHelper::passTight_BDT_Id(pat::Electron electron, std::string elecID) {
    float cutVal=1000;
    float fSCeta = fabs(electron.superCluster()->eta());
    if(electron.pt()<=10){ 
        if(fSCeta < 0.8) cutVal = -0.265; 
        if(fSCeta >= 0.8 && fSCeta < 1.479) cutVal = -0.556;
        if(fSCeta >= 1.479) cutVal = -0.551;  
    }
    else {
        if(fSCeta < 0.8) cutVal = -0.072; 
        if(fSCeta >= 0.8 && fSCeta < 1.479) cutVal = -0.286;
        if(fSCeta >= 1.479) cutVal = -0.267;
    }
    if (electron.electronID(elecID) > cutVal ) { return true;}
    return false;
}

bool HZZ4LHelper::passTight_Id_SUS(pat::Muon muon, const reco::Vertex *&vertex) {
    if (muon.isPFMuon() != 1) return false;
    if (muon.isMediumMuon() != 1) return false;
    double dxyCut = 0.05;
    double dzCut = 0.1;
    if( fabs(muon.muonBestTrack()->dxy(vertex->position())) >= dxyCut ) return false;
    if( fabs(muon.muonBestTrack()->dz(vertex->position())) >= dzCut ) return false;
    return true;
}

bool HZZ4LHelper::passTight_Id_SUS(pat::Electron electron, std::string elecID, const reco::Vertex *&vertex, const reco::BeamSpot BS, edm::Handle< std::vector<reco::Conversion> > theConversions) {

    double dxyCut = 0.05;
    double dzCut = 0.1;
    //std::cout<<"dxy: "<<fabs(electron.gsfTrack()->dxy(vertex->position()))<<" dz: "<<fabs(electron.gsfTrack()->dz(vertex->position()))<<std::endl;
    if( fabs(electron.gsfTrack()->dxy(vertex->position())) >= dxyCut ) return false;
    if( fabs(electron.gsfTrack()->dz(vertex->position())) >= dzCut ) return false;

    float cutVal=1000;
    float fSCeta = fabs(electron.eta());
    if(electron.pt()<=10){
        if(fSCeta < 0.8) cutVal = 0.87;
        if(fSCeta >= 0.8 && fSCeta < 1.479) cutVal = 0.60;
        if(fSCeta >= 1.479) cutVal = 0.17;
    }
    else {
        if(fSCeta < 0.8) cutVal = 0.87;
        if(fSCeta >= 0.8 && fSCeta < 1.479) cutVal = 0.60;
        if(fSCeta >= 1.479) cutVal = 0.17;
    }
    //std::cout<<"el |eta|: "<<fabs(electron.eta())<<" el |SCEta|: "<< fabs(electron.superCluster()->eta())<<" el mva: "<<electron.electronID(elecID)<<std::endl;
    if (electron.electronID(elecID) <= cutVal ) return false;

    bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron(electron), theConversions, BS.position());
    //std::cout<<"vtxFitConverstion: "<<vtxFitConversion<<std::endl;
    if( vtxFitConversion )  return false;

    int missingHitsCuts = 1;
    int misHits = electron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    //std::cout<<"misHits: "<<misHits<<std::endl;
    if (misHits >= missingHitsCuts) return false;
    
    //std::cout<<"check id emulation"<<std::endl;
    // Trigger ID Emulation
    //if (fabs(electron.eta())<1.479) {
    //std::cout<<"isEB: "<<electron.isEB()<<" iEtaieta: "<<electron.full5x5_sigmaIetaIeta()<<" hOe: "<<electron.hcalOverEcal()<<" dEta: "<<electron.deltaEtaSuperClusterTrackAtVtx()<<" dPhi: "<<electron.deltaPhiSuperClusterTrackAtVtx()<<" |1/e - 1/p| "<<fabs(1.0/electron.correctedEcalEnergy() - electron.eSuperClusterOverP()/electron.correctedEcalEnergy())<<std::endl;

    if (electron.isEB()) {
        if (electron.full5x5_sigmaIetaIeta()>=0.011) return false;
        if (electron.hcalOverEcal()>=0.08) return false;
        if (fabs(electron.deltaEtaSuperClusterTrackAtVtx())>=0.01) return false;
        if (fabs(electron.deltaPhiSuperClusterTrackAtVtx())>=0.04) return false;
        if (fabs(1.0/electron.correctedEcalEnergy() - electron.eSuperClusterOverP()/electron.correctedEcalEnergy())>=0.01) return false;
    }
    if (electron.isEE()) {
        if (electron.full5x5_sigmaIetaIeta()>=0.031) return false;
        if (electron.hcalOverEcal()>=0.08) return false;
        if (fabs(electron.deltaEtaSuperClusterTrackAtVtx())>=0.01) return false;
        if (fabs(electron.deltaPhiSuperClusterTrackAtVtx())>=0.08) return false;
        if (fabs(1.0/electron.correctedEcalEnergy() - electron.eSuperClusterOverP()/electron.correctedEcalEnergy())>=0.01) return false;
    }

    //std::cout<<"passed id emulation"<<std::endl;

    return true;
}

double HZZ4LHelper::ptRatio(pat::Electron electron, edm::Handle<edm::View<pat::Jet> > jets, bool isMC) {

    double ptratio = 9999.0;
    double mindR=9999.0;
    
    for(unsigned int j = 0; j < jets->size(); ++j) {        
        const pat::Jet & patjet = jets->at(j);
        
        TLorentzVector jet_lepawareJEC(patjet.correctedP4(0).px(),patjet.correctedP4(0).py(),patjet.correctedP4(0).pz(),patjet.correctedP4(0).energy());
        TLorentzVector lep(electron.px(),electron.py(),electron.pz(),electron.energy());
        double corr_L1 = patjet.jecFactor("L1FastJet")/patjet.jecFactor("Uncorrected");
        double corr_L2L3 = patjet.jecFactor("L3Absolute")/patjet.jecFactor("Uncorrected");
        if (!isMC) corr_L2L3=patjet.jecFactor("L2L3Residual")/patjet.jecFactor("Uncorrected");
        jet_lepawareJEC *= corr_L1;
        jet_lepawareJEC -= lep;
        jet_lepawareJEC *= corr_L2L3/corr_L1;
        jet_lepawareJEC += lep;

        double tmpdR = deltaR(electron.eta(),electron.phi(),jet_lepawareJEC.Eta(),jet_lepawareJEC.Phi());        
        if (tmpdR<mindR && tmpdR<0.4) {
            mindR = tmpdR;

            ptratio = electron.pt()/jet_lepawareJEC.Pt();
        }        
    }    
    return ptratio;
}

double HZZ4LHelper::ptRatio(pat::Muon muon, edm::Handle<edm::View<pat::Jet> > jets,bool isMC) {

    double ptratio = 9999.0;
    double mindR=9999.0;
    
    for(unsigned int j = 0; j < jets->size(); ++j) {        
        const pat::Jet & patjet = jets->at(j);

        TLorentzVector jet_lepawareJEC(patjet.correctedP4(0).px(),patjet.correctedP4(0).py(),patjet.correctedP4(0).pz(),patjet.correctedP4(0).energy());
        TLorentzVector lep(muon.px(),muon.py(),muon.pz(),muon.energy());
        double corr_L1 = patjet.jecFactor("L1FastJet")/patjet.jecFactor("Uncorrected");
        double corr_L2L3 = patjet.jecFactor("L3Absolute")/patjet.jecFactor("Uncorrected");
        if (!isMC)corr_L2L3=patjet.jecFactor("L2L3Residual")/patjet.jecFactor("Uncorrected");
        jet_lepawareJEC *= corr_L1;
        jet_lepawareJEC -= lep;
        jet_lepawareJEC *= corr_L2L3/corr_L1;
        jet_lepawareJEC += lep;

        double tmpdR = deltaR(muon.eta(),muon.phi(),jet_lepawareJEC.Eta(),jet_lepawareJEC.Phi());
        if (tmpdR<mindR && tmpdR<0.4) {
            mindR = tmpdR;
            ptratio = muon.pt()/jet_lepawareJEC.Pt();
        }
    }    
    return ptratio;
}

double HZZ4LHelper::ptRel(pat::Electron electron, edm::Handle<edm::View<pat::Jet> > jets,bool isMC) {

    double ptrel = 9999.0;
    double mindR=9999.0;
    
    for(unsigned int j = 0; j < jets->size(); ++j) {        
        const pat::Jet & patjet = jets->at(j);

        TLorentzVector jet_lepawareJEC(patjet.correctedP4(0).px(),patjet.correctedP4(0).py(),patjet.correctedP4(0).pz(),patjet.correctedP4(0).energy());
        TLorentzVector lep(electron.px(),electron.py(),electron.pz(),electron.energy());
        double corr_L1 = patjet.jecFactor("L1FastJet")/patjet.jecFactor("Uncorrected");
        double corr_L2L3 = patjet.jecFactor("L3Absolute")/patjet.jecFactor("Uncorrected");
        if (!isMC) corr_L2L3=patjet.jecFactor("L2L3Residual")/patjet.jecFactor("Uncorrected");
        jet_lepawareJEC *= corr_L1;
        jet_lepawareJEC -= lep;
        jet_lepawareJEC *= corr_L2L3/corr_L1;
        jet_lepawareJEC += lep;

        double tmpdR = deltaR(electron.eta(),electron.phi(),jet_lepawareJEC.Eta(),jet_lepawareJEC.Phi());
        if (tmpdR<mindR && tmpdR<0.4) {
            mindR = tmpdR;
            TVector3 pl(electron.px(),electron.py(),electron.pz());
            TLorentzVector p4jminusl;
            p4jminusl = lep-jet_lepawareJEC;
            TVector3 pjminusl(p4jminusl.Px(),p4jminusl.Py(),p4jminusl.Pz());
            double sin_alpha = (pl.Cross(pjminusl)).Mag()/pl.Mag()/pjminusl.Mag();
            ptrel = pl.Mag()*sin_alpha;     
        }        
    }    
    return ptrel;
}

double HZZ4LHelper::ptRel(pat::Muon muon, edm::Handle<edm::View<pat::Jet> > jets,bool isMC) {

    double ptrel = 9999.0;
    double mindR=9999.0;

    for(unsigned int j = 0; j < jets->size(); ++j) {
        const pat::Jet & patjet = jets->at(j);

        TLorentzVector jet_lepawareJEC(patjet.correctedP4(0).px(),patjet.correctedP4(0).py(),patjet.correctedP4(0).pz(),patjet.correctedP4(0).energy());
        TLorentzVector lep(muon.px(),muon.py(),muon.pz(),muon.energy());
        double corr_L1 = patjet.jecFactor("L1FastJet")/patjet.jecFactor("Uncorrected");
        double corr_L2L3 = patjet.jecFactor("L3Absolute")/patjet.jecFactor("Uncorrected");
        if (!isMC) corr_L2L3=patjet.jecFactor("L2L3Residual")/patjet.jecFactor("Uncorrected");
        jet_lepawareJEC *= corr_L1;
        jet_lepawareJEC -= lep;
        jet_lepawareJEC *= corr_L2L3/corr_L1;
        jet_lepawareJEC += lep;

        double tmpdR = deltaR(muon.eta(),muon.phi(),jet_lepawareJEC.Eta(),jet_lepawareJEC.Phi());
        if (tmpdR<mindR && tmpdR<0.4) {
            mindR = tmpdR;
            TVector3 pl(muon.px(),muon.py(),muon.pz());
            TLorentzVector p4jminusl;
            p4jminusl = lep-jet_lepawareJEC;
            TVector3 pjminusl(p4jminusl.Px(),p4jminusl.Py(),p4jminusl.Pz());
            double sin_alpha = (pl.Cross(pjminusl)).Mag()/pl.Mag()/pjminusl.Mag();
            ptrel = pl.Mag()*sin_alpha;
        }
    }

    return ptrel;
}

double HZZ4LHelper::pfIso(pat::Muon muon, double Rho) {
    double PUCorr = 0.5*muon.puChargedHadronIso();
    double iso = (muon.chargedHadronIso()+std::max(muon.photonIso()+muon.neutralHadronIso()-PUCorr,0.0))/muon.pt();
    return iso;
}

double HZZ4LHelper::pfIso03(pat::Muon muon, double Rho) {
    double PUCorr = 0.5*muon.pfIsolationR03().sumPUPt;
    double isoCH = muon.pfIsolationR03().sumChargedHadronPt;
    double isoNH = muon.pfIsolationR03().sumNeutralHadronEt;
    double isoPhot = muon.pfIsolationR03().sumPhotonEt;
    double iso = (isoCH+std::max(isoPhot+isoNH-PUCorr,0.0))/muon.pt();    
    return iso;
}


double HZZ4LHelper::pfIso(pat::Electron elec, double Rho) {
    double PUCorr = Rho*ElecEffArea(elEAtype,elec.eta(),elEAtarget);
    double iso = (elec.chargedHadronIso()+std::max(elec.photonIso()+elec.neutralHadronIso()-PUCorr,0.0))/elec.pt();
    return iso;
}

double HZZ4LHelper::pfIso03(pat::Electron elec, double Rho) {
    double PUCorr = Rho*ElecEffArea(kEGammaNeutralHadIso03,elec.superCluster()->eta(),elEAtarget);
    double isoCH = elec.pfIsolationVariables().sumChargedHadronPt;
    double isoNH = elec.pfIsolationVariables().sumNeutralHadronEt;
    double isoPhot = elec.pfIsolationVariables().sumPhotonEt;
    double iso = (isoCH+std::max(isoPhot+isoNH-PUCorr,0.0))/elec.pt();    
    return iso;
}

double HZZ4LHelper::getPUIso(pat::Muon muon, double Rho) {
    double puiso = 0.5*muon.puChargedHadronIso();
    if (Rho<0.0) puiso = 0.;
    return puiso;
}

double HZZ4LHelper::getPUIso03(pat::Muon muon, double Rho) {
    double puiso = 0.5*muon.pfIsolationR03().sumPUPt;
    if (Rho<0.0) puiso = 0.;
    return puiso;
}

double HZZ4LHelper::getPUIso(pat::Electron elec, double Rho) {
    double puiso = Rho*ElecEffArea(elEAtype,elec.eta(),elEAtarget); 
    return puiso;
}

double HZZ4LHelper::getPUIso03(pat::Electron elec, double Rho) {
    double puiso = Rho*ElecEffArea(kEGammaNeutralHadIso03,elec.superCluster()->eta(),elEAtarget); 
    return puiso;
}

double HZZ4LHelper::MuonEffArea(MuonEffectiveAreaType type, double SCEta, MuonEffectiveAreaTarget EffectiveAreaTarget) {

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

    if(EffectiveAreaTarget == kElEAData2015)
    {

        if (type == kEGammaNeutralHadIso04){
            if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 0.8 )   EffectiveArea = 0.1830;
            if (fabs(SCEta) >= 0.8 && fabs(SCEta) < 1.3 )   EffectiveArea = 0.1734;
            if (fabs(SCEta) >= 1.3 && fabs(SCEta) < 2.0 )   EffectiveArea = 0.1077;
            if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )   EffectiveArea = 0.1565;
            if (fabs(SCEta) >= 2.2 )                        EffectiveArea = 0.2680;
        }
        if (type == kEGammaNeutralHadIso03){
            if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 )    EffectiveArea = 0.1752;
            if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 )  EffectiveArea = 0.1862;
            if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 )  EffectiveArea = 0.1411;
            if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 )    EffectiveArea = 0.1534;
            if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 )    EffectiveArea = 0.1903;
            if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 )    EffectiveArea = 0.2243;
            if (fabs(SCEta) >= 2.4 && fabs(SCEta) < 2.5 )    EffectiveArea = 0.2687;
        }

        return EffectiveArea;
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

 
float HZZ4LHelper::kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {        
        k+=1.515838921760*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
        k+=1.496256665410*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
        k+=1.495522061910*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
        k+=1.483273154250*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
        k+=1.465589701130*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
        k+=1.491500887510*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
        k+=1.441183580450*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
        k+=1.440830603990*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
        k+=1.414339019120*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
        k+=1.422534218560*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
        k+=1.401037066000*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
        k+=1.408539428810*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
        k+=1.381247744080*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
        k+=1.370553357430*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
        k+=1.347323316000*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
        k+=1.340113437450*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
        k+=1.312661036510*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
        k+=1.290055062010*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
        k+=1.255322614790*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
        k+=1.254455642450*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
        k+=1.224047664420*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
        k+=1.178816782670*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
        k+=1.162624827140*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
        k+=1.105401140940*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
        k+=1.074749265690*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
        k+=1.021864599380*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
        k+=0.946334793286*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
        k+=0.857458082628*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
        k+=0.716607670482*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
        k+=1.132841784840*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);
    }

    if (finalState==2) {
       k+=1.513834489150*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
       k+=1.541738780180*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
       k+=1.497829632510*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
       k+=1.534956782920*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
       k+=1.478217033060*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
       k+=1.504330859290*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
       k+=1.520626246850*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
       k+=1.507013090030*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
       k+=1.494243156250*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
       k+=1.450536096150*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
       k+=1.460812521660*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
       k+=1.471603622200*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
       k+=1.467700038200*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
       k+=1.422408690640*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
       k+=1.397184022730*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
       k+=1.375593447520*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
       k+=1.391901318370*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
       k+=1.368564350560*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
       k+=1.317884804290*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
       k+=1.314019950800*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
       k+=1.274641749910*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
       k+=1.242346606820*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
       k+=1.244727403840*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
       k+=1.146259351670*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
       k+=1.107804993520*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
       k+=1.042053646740*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
       k+=0.973608545141*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
       k+=0.872169942668*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
       k+=0.734505279177*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
       k+=1.163152837230*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);       
    }
    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;

}

float HZZ4LHelper::dataMC(pat::Electron electron, TH2F* hElecScaleFac, TH2F* hElecScaleFac_Cracks)
{
    float pt = std::min(electron.pt(),199.0);
    float eta = electron.superCluster()->eta();
    if (electron.isGap()) {
        return hElecScaleFac_Cracks->GetBinContent(hElecScaleFac_Cracks->FindBin(pt,eta));        
    } else {
        return hElecScaleFac->GetBinContent(hElecScaleFac->FindBin(pt,eta));        
    }
    return 1.0;
}

float HZZ4LHelper::dataMCErr(pat::Electron electron, TH2F* hElecScaleFac, TH2F* hElecScaleFac_Cracks)
{
    float pt = std::min(electron.pt(),199.0);
    float eta = electron.superCluster()->eta();
    if (electron.isGap()) {
        return hElecScaleFac_Cracks->GetBinError(hElecScaleFac_Cracks->FindBin(pt,eta));        
    } else {
        return hElecScaleFac->GetBinError(hElecScaleFac->FindBin(pt,eta));        
    }
    return 1.0;
}

float HZZ4LHelper::dataMC(pat::Muon muon, TH2F* hElecScaleFac, TH2F* hElecScaleFac_Cracks)
{
    return 1.0;
}

float HZZ4LHelper::dataMCErr(pat::Muon muon, TH2F* hElecScaleFac, TH2F* hElecScaleFac_Cracks)
{
    return 0.0;
}

#endif
