/*
 *  $Date: 2008/08/29 12:16 $
 *  
 * 
 *  26/10/05
 *
 */

#include "GeneratorInterface/PompytInterface/interface/PompytProducer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "SimDataFormats/GeneratorProducts/interface/GenInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "GeneratorInterface/PompytInterface/interface/PYR.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Utilities/General/interface/FileInPath.h"
#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include <iostream>
#include "time.h"

using namespace edm;
using namespace std;

#include "HepMC/PythiaWrapper6_4.h"
#include "GeneratorInterface/PompytInterface/interface/Pompyt.h"
#include "HepMC/IO_HEPEVT.h"


#define PYGIVE pygive_
extern "C" {
  void PYGIVE(const char*,int length);
}


namespace {
 CLHEP::HepRandomEngine& getEngineReference()
  {

   Service<RandomNumberGenerator> rng;
   if(!rng.isAvailable()) {
    throw cms::Exception("Configuration")
       << "The RandomNumberProducer module requires the RandomNumberGeneratorService\n"
          "which appears to be absent.  Please add that service to your configuration\n"
          "or remove the modules that require it.";
   }

// The Service has already instantiated an engine.  Make contact with it.
   return (rng->getEngine());
  }
}




HepMC::IO_HEPEVT conv2;

//used for defaults
  static const unsigned long kNanoSecPerSec = 1000000000;
  static const unsigned long kAveEventPerSec = 200;

PompytProducer::PompytProducer( const ParameterSet & pset) :
  EDProducer(), evt(0), 
  pythiaPylistVerbosity_ (pset.getUntrackedParameter<int>("pythiaPylistVerbosity",0)),
  pythiaHepMCVerbosity_ (pset.getUntrackedParameter<bool>("pythiaHepMCVerbosity",false)),
  maxEventsToPrint_ (pset.getUntrackedParameter<int>("maxEventsToPrint",1)),
  comenergy(pset.getUntrackedParameter<double>("comEnergy",14000.)),
  extCrossSect(pset.getUntrackedParameter<double>("crossSection", -1.)),
  extFilterEff(pset.getUntrackedParameter<double>("filterEfficiency", -1.)),
  pomz(pset.getUntrackedParameter<double>("pomZ", -1.)),
  pompdiss(pset.getUntrackedParameter<int>("pomPdiss",1)),
  xi_max(pset.getUntrackedParameter<double>("xiMax",0.1)),
  pdiss_b0(pset.getUntrackedParameter<double>("pdiss_B0",1.5)),
  pdiss_ap(pset.getUntrackedParameter<double>("pdiss_AP",0.016)),
  pdiss_glu(pset.getUntrackedParameter<double>("pdiss_GLU",0.5)),
  pdiss_mmax(pset.getUntrackedParameter<double>("pdiss_MMAX",2000.0)),
  fRandomEngine(getEngineReference()),
  eventNumber_(0)
{
  // PYLIST Verbosity Level
  // Valid PYLIST arguments are: 1, 2, 3, 5, 7, 11, 12, 13
  pythiaPylistVerbosity_ = pset.getUntrackedParameter<int>("pythiaPylistVerbosity",0);
  
  // HepMC event verbosity Level
  pythiaHepMCVerbosity_ = pset.getUntrackedParameter<bool>("pythiaHepMCVerbosity",false);

  //Max number of events printed on verbosity level 
  maxEventsToPrint_ = pset.getUntrackedParameter<int>("maxEventsToPrint",0);
  // CMS energy
  //comenergy = pset.getUntrackedParameter<double>("comEnergy",14000.);

 // POMPYT
//  pomz = pset.getUntrackedParameter<double>("pomZ",1.);
//  pompdiss = pset.getUntrackedParameter<int>("pomPdiss",1);
//  xi_max = pset.getUntrackedParameter<double>("xiMax",0.1);


 // Initialize the random engine unconditionally!
  randomEngine = &fRandomEngine;
  fRandomGenerator = new CLHEP::RandFlat(fRandomEngine) ;



  // Set PYTHIA parameters in a single ParameterSet
  ParameterSet pythia_params = 
    pset.getParameter<ParameterSet>("PythiaParameters") ;
  
  // The parameter sets to be read (default, min bias, user ...) in the
  // proper order.
  vector<string> setNames = 
    pythia_params.getParameter<vector<string> >("parameterSets");
  

 
  // Loop over the sets
  for ( unsigned i=0; i<setNames.size(); ++i ) {
    
    string mySet = setNames[i];

    // Read the PYTHIA parameters for each set of parameters
    vector<string> pars = 
      pythia_params.getParameter<vector<string> >(mySet);
    
    if (mySet != "SLHAParameters" && mySet != "CSAParameters"){
    
    // Loop over all parameters and stop in case of mistake
    for( vector<string>::const_iterator  
	   itPar = pars.begin(); itPar != pars.end(); ++itPar ) {

      static string sRandomValueSetting("MRPY(1)");
      if( 0 == itPar->compare(0,sRandomValueSetting.size(),sRandomValueSetting) ) {
	throw edm::Exception(edm::errors::Configuration,"PythiaError")
	  <<" attempted to set random number using pythia command 'MRPY(1)' this is not allowed.\n  Please use the RandomNumberGeneratorService to set the random number seed.";
      }
      if( ! call_pygive(*itPar) ) {
	throw edm::Exception(edm::errors::Configuration,"PythiaError") 
	  <<" pythia did not accept the following \""<<*itPar<<"\"";
      }
    }
  } 
  }

 
  //////////////////////////////////////////////////////////////////////////////////////////////

  pompyt1.pomi[2-1]=pomz;      //pomeron direction

  pomcut1.pomcut[1-1]=0.0001;  //xi min
  pomcut1.pomcut[2-1]=xi_max;    //xi max
  pomcut1.pomcut[3-1]=-3.;       // t max
  pomcut1.pomcut[4-1]=0.;        //t min

  pompyt1.ipom[1-1]=pompdiss;    // pp->Xp (1) or pp->XY (2)
  pompyt1.pomi[19-1]=pdiss_b0;   // B0 for pdiss events
  pompyt1.pomi[20-1]=pdiss_ap;   //alpha' for pdiss events

  pomcut1.pomcut[10-1]=pdiss_glu;
  pomcut1.pomcut[11-1]=1.1;
  pomcut1.pomcut[12-1]=pdiss_mmax;


  call_pomini(comenergy);

  cout << endl; 
  
  produces<HepMCProduct>();
   //  produces<GenInfoProduct, edm::InRun>();
  produces<GenEventInfoProduct>();
  produces<GenRunInfoProduct, edm::InRun>();

}



PompytProducer::~PompytProducer(){
  call_pomout();
  clear(); 
}

void PompytProducer::clear() {
 
}

void PompytProducer::endRun(Run & r) {

  double cs = 0;
  auto_ptr<GenRunInfoProduct> genRunInfoProd (new GenRunInfoProduct());
  genRunInfoProd->setInternalXSec(cs);
  r.put(genRunInfoProd);

}

void PompytProducer::produce(Event & e, const EventSetup& es ) {

    auto_ptr<HepMCProduct> bare_product(new HepMCProduct());  

    //  call_pyevnt();      // generate one event with Pythia
     call_pompyt();

    // convert to stdhep format
     call_pyhepc( 1 );
    
    // convert stdhep (hepevt) to hepmc
     HepMC::GenEvent* evt = conv2.read_next_event();

     evt->set_signal_process_id(pypars.msti[0]);
     evt->set_event_scale(pypars.pari[16]);
     ++eventNumber_;
     evt->set_event_number(eventNumber_);

     evt->weights().push_back( pomcut1.pomcut[9-1] );



     /*

     //   evt->set_signal_process_id(pypars.msti[0]);
     //   evt->set_event_scale(pypars.pari[16]);
     //   evt->set_event_number(numberEventsInRun() - remainingEvents() - 1);

    // int id1 = pypars.msti[14];
    // int id2 = pypars.msti[15];
    int id1 = pyint1.mint[14];
    int id2 = pyint1.mint[15];
    if ( id1 == 21 ) id1 = 0;
    if ( id2 == 21 ) id2 = 0; 
    double x1 = pyint1.vint[40];
    double x2 = pyint1.vint[41];  
    double Q  = pyint1.vint[50];
    double pdf1 = pyint1.vint[38];
    pdf1 /= x1 ;
    double pdf2 = pyint1.vint[39];
    pdf2 /= x2 ;
    evt->set_pdf_info( HepMC::PdfInfo(id1,id2,x1,x2,Q,pdf1,pdf2) ) ;
    
    evt->weights().push_back( pyint1.vint[96] );

    
     */

    //******** Verbosity ********

       if(e.id().event() <= maxEventsToPrint_ &&
        (pythiaPylistVerbosity_ || pythiaHepMCVerbosity_)) {
 
       // Prints PYLIST info
       if(pythiaPylistVerbosity_) {
         call_pylist(pythiaPylistVerbosity_);
       }
       
       // Prints HepMC event
       if(pythiaHepMCVerbosity_) {
         cout << "Event process = " << pypars.msti[0] << endl 
         << "----------------------" << endl;
         evt->print();
       }
     }

    


    if(evt)  bare_product->addHepMCData(evt );
    e.put(bare_product);

    std::auto_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(evt));
    e.put(genEventInfo);

  return;

}


bool 
PompytProducer::call_pygive(const std::string& iParm ) {

  int numWarn = pydat1.mstu[26]; //# warnings
  int numErr = pydat1.mstu[22];// # errors
  
//call the fortran routine pygive with a fortran string
  PYGIVE( iParm.c_str(), iParm.length() );  
  //  PYGIVE( iParm );  
//if an error or warning happens it is problem
  return pydat1.mstu[26] == numWarn && pydat1.mstu[22] == numErr;   
 
}

