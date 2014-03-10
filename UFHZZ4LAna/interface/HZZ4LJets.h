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

  bool looseID = jet.neutralHadronEnergyFraction() < 0.99
    && jet.neutralEmEnergyFraction() < 0.99
    && jet.nConstituents() > 1;

      
  bool loEtaID = jet.chargedMultiplicity() > 0
    && jet.chargedHadronEnergyFraction() > 0.0
    && jet.chargedEmEnergyFraction() < 0.99;


  if(fabs(jet.eta()) < 2.4 && looseID && loEtaID) return 1;
  if(fabs(jet.eta()) >= 2.4 && looseID) return 1;
  

  return 0;
}

#endif
