#include "UFHZZAnalysis8TeV/UFHZZ4LAna/interface/LheNup6Nup7Filter.h"

using namespace edm;
using namespace std;

LheNup6Nup7Filter::LheNup6Nup7Filter(const edm::ParameterSet& iPSet):  
  skipLheNup6Nup7_(iPSet.getUntrackedParameter<bool>("skipLheNup6Nup7",false))
{    
	nSkippedLheNup6Nup7 = 0;
	nScannedLhe = 0;
	nScanned = 0;
}

LheNup6Nup7Filter::~LheNup6Nup7Filter() {
	cout << "[LHE NUP6 & NUP7 Filter]" << endl
	     << "Total scanned events: " << nScanned << /*", Events with LHE: " << nScannedLhe << */", Skipped LHE NUP6 & NUP7 events: " << nSkippedLheNup6Nup7 << endl << endl;
}


bool LheNup6Nup7Filter::filter(edm::Event& iEvent,const edm::EventSetup& iSetup){ 
  // skipp only the LHE event with NUP==6 || NUP==7, all the other events should pass
  bool event_passed = true;
	
  Handle<LHEEventProduct> lheProd;
  bool LHEEventProduct_found= iEvent.getByType( lheProd );
	
  if(LHEEventProduct_found){ 
    const lhef::HEPEUP hepeup = lheProd->hepeup();
	if (skipLheNup6Nup7_ && (hepeup.NUP==6 || hepeup.NUP==7)) {
      event_passed = false;
	  nSkippedLheNup6Nup7++;
	}
	nScannedLhe++;
  }else {
    //cout << "LHEEventProduct not found." << endl;		
  }
	nScanned++;
	
  return event_passed;
}

DEFINE_FWK_MODULE(LheNup6Nup7Filter);

