#ifndef HZZ4LJETS_H
#define HZZ4LJETS_H

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
 #include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/ServiceRegistry/interface/Service.h"
 #include "CommonTools/UtilAlgos/interface/TFileService.h"


class HZZ4LJets
{

 public:

  HZZ4LJets();
  ~HZZ4LJets();

  int patjetID(const pat::Jet& jet);

 private:

};


#endif


#ifndef HZZ4LJETS_CC
#define HZZ4LJETS_CC

HZZ4LJets::HZZ4LJets()
{

  //declarations

}


HZZ4LJets::~HZZ4LJets()
{

  //destructor ---do nothing

}


// 0 is fail 1 is loose, 2 is medium, 3 is tight
int HZZ4LJets::patjetID(const pat::Jet& jet)
{

  double NHF = jet.neutralHadronEnergyFraction();
  double NEMF = jet.neutralEmEnergyFraction();
  double CHF = jet.chargedHadronEnergyFraction();
  //double MUF = jet.muonEnergyFraction();
  double CEMF = jet.chargedEmEnergyFraction();
  double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
  double NumNeutralParticles =jet.neutralMultiplicity();
  double CHM = jet.chargedMultiplicity(); 

  bool tightestJetID_CEN = ((NHF<0.70 && NEMF<0.90 && NumConst>1) && ((fabs(jet.eta())<=2.4 && CHF>0.2 && CHM>0 && CEMF<0.99) || fabs(jet.eta())>2.4) && fabs(jet.eta())<=3.0);
  bool tightestJetID_FWD = (NEMF<0.90 && NumNeutralParticles>10 && fabs(jet.eta())>3.0);
  if (tightestJetID_CEN || tightestJetID_FWD) return 3;

  bool tightJetID_CEN = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(jet.eta())>2.4) && fabs(jet.eta())<=3.0);
  bool tightJetID_FWD = (NEMF<0.90 && NumNeutralParticles>10 && fabs(jet.eta())>3.0);
  if (tightJetID_CEN || tightJetID_FWD) return 2;
 
  bool looseJetID_CEN = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(jet.eta())>2.4) && fabs(jet.eta())<=3.0);
  bool looseJetID_FWD = (NEMF<0.90 && NumNeutralParticles>10 && fabs(jet.eta())>3.0);
  if (looseJetID_CEN || looseJetID_FWD) return 1;

  return 0;

}

#endif
