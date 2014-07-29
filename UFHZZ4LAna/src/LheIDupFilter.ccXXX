#include "UFHZZAnalysis8TeV/UFHZZ4LAna/interface/LheIDupFilter.h"

using namespace edm;
using namespace std;

LheIDupFilter::LheIDupFilter(const edm::ParameterSet& iPSet):  
  skipLheIDup_(iPSet.getUntrackedParameter<int>("skipLheIDup",0)),
  numLheIDup_(iPSet.getUntrackedParameter<int>("numLheIDup",1))
{    
	nSkippedLheIDup = 0;
	nScannedLhe = 0;
	nScanned = 0;
}

LheIDupFilter::~LheIDupFilter() {
	cout << "[LHE IDUP Filter]" << endl << "Total scanned events: " << nScanned << /*", Events with LHE: " << nScannedLhe << */", Skipped LHE ("<<numLheIDup_<<"x)IDUP="<<skipLheIDup_<<" events: " << nSkippedLheIDup << endl << endl;
}


bool LheIDupFilter::filter(edm::Event& iEvent,const edm::EventSetup& iSetup){ 
  // skipp only the LHE event with "numLheIDup_" particles with IDUP=="skipLheIDup_", all the other events should pass
  bool event_passed = true;
	
  Handle<LHEEventProduct> lheProd;
  bool LHEEventProduct_found= iEvent.getByType( lheProd );
	
  if(LHEEventProduct_found){ 
    const lhef::HEPEUP hepeup = lheProd->hepeup();
	const std::vector<lhef::HEPEUP::FiveVector> pup = hepeup.PUP;
	int foundLheIDup = 0;
	for(unsigned int i=0; i<pup.size(); ++i){
	  if(hepeup.ISTUP[i]!=1) continue;
	  if( abs(hepeup.IDUP[i]) != skipLheIDup_ ) continue; // pdgID to skip
	  foundLheIDup++;
	  if (foundLheIDup >= numLheIDup_) {
        event_passed = false;
		nSkippedLheIDup++;
		//cout << "hepeup.IDUP="<<abs(hepeup.IDUP[i])<<" found "<<foundLheIDup<<" time(s)." << endl;		
		break;  
	  }
	}
	nScannedLhe++;
  }else {
     //cout << "LHEEventProduct not found." << endl;		
  }
  nScanned++;
	
  return event_passed;
}

DEFINE_FWK_MODULE(LheIDupFilter);

