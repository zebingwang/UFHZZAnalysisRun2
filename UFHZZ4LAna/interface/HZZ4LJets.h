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

  int patjetID(const pat::Jet& jet, int year);

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
int HZZ4LJets::patjetID(const pat::Jet& jet, int year)
{

  double NHF = jet.neutralHadronEnergyFraction();
  double NEMF = jet.neutralEmEnergyFraction();
  double CHF = jet.chargedHadronEnergyFraction();
  double CHM = jet.chargedMultiplicity(); 
  double CEMF = jet.chargedEmEnergyFraction();
  double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
  double NumNeutralParticle =jet.neutralMultiplicity();


  bool looseJetID=false;
  bool tightJetID=false;

  double eta = fabs(jet.eta());


  // 2017 Jet ID
  /*
  if (eta<=2.7) {

      looseJetID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0) || eta>2.4) && eta<=2.7);
      tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0) || eta>2.4) && eta<=2.7);

  } else if (eta>2.7 && eta<=3.0) {

      looseJetID = ( NEMF>0.02 && NEMF<0.99 && eta>2.7 && eta<=3.0 );
      tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 && eta>2.7 && eta<=3.0 );

  } else if (eta>3.0) {

      looseJetID = ( NEMF<0.90 && NumNeutralParticle>10 && eta>3.0 );
      tightJetID = ( NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10 && eta>3.0 );

  }
  */
 
  /*
  // 2018 Jet ID
  if (eta<=2.6) {

      looseJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
      tightJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );

  } else if (eta>2.6 && eta<=2.7) {

      looseJetID = ( CHM>0 && NEMF<0.99 && NHF < 0.9 );
      tightJetID = ( CHM>0 && NEMF<0.99 && NHF < 0.9 );

  } else if (eta>2.7 && eta<=3.0) {

      looseJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 );
      tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 );

  } else if (eta>3.0) {

      looseJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10 ); 
      tightJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10 ); 

  }
  */
  if(year==2018)
  {
      if (eta<=2.6) {

          looseJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
          tightJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );

      } else if (eta>2.6 && eta<=2.7) {
  
          looseJetID = ( CHM>0 && NEMF<0.99 && NHF < 0.9 );
          tightJetID = ( CHM>0 && NEMF<0.99 && NHF < 0.9 );
    
      } else if (eta>2.7 && eta<=3.0) {
    
          looseJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 );
          tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 );
    
      } else if (eta>3.0) {
    
          looseJetID = (NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 );
          tightJetID = (NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 );
    
      }
  //tightJetID = ( abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9);
  //tightJetID = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 );
  //tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
  //tightJetID = ( NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 && abs(eta)>3.0 )
 
  }

  if(year==2017)
  {
      if (eta<=2.7) {
    
          looseJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0) || eta>2.4) && eta<=2.7);
          tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0) || eta>2.4) && eta<=2.7);
    
      } else if (eta>2.7 && eta<=3.0) {
    
          looseJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 && eta>2.7 && eta<=3.0 );
          tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2 && eta>2.7 && eta<=3.0 );
    
      } else if (eta>3.0) {
    
          looseJetID = ( NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10 && eta>3.0 );
          tightJetID = ( NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10 && eta>3.0 );
    
      }
  }

  if(year==2016)
  {
      if (eta<=2.7) {

          looseJetID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0 && CEMF <0.99) || eta>2.4) && eta<=2.7);
          tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0 && CHM>0 && CEMF <0.99) || eta>2.4) && eta<=2.7);

      } else if (eta>2.7 && eta<=3.0) {

          looseJetID = ( NEMF>0.01 && NHF<0.98 && NumNeutralParticle>2 && eta>2.7 && eta<=3.0 );
          tightJetID = ( NEMF>0.01 && NHF<0.98 && NumNeutralParticle>2 && eta>2.7 && eta<=3.0 );

      } else if (eta>3.0) {

          looseJetID = ( NEMF<0.90 && NumNeutralParticle>10 && eta>3.0 );
          tightJetID = ( NEMF<0.90 && NumNeutralParticle>10 && eta>3.0 );

      }
  }


  if (tightJetID) {return 2;}
  else if (looseJetID) {return 1;}
  else {return 0;}

}

#endif
