#ifndef Pompyt_Producer_h
#define Pompyt_Producer_h

/** \class PompytProducer
 *
 * Generates Pompyt HepMC events
 * 
 ***************************************/

#define PYCOMP pycomp_

#include "FWCore/Framework/interface/EDProducer.h"
//#include "FWCore/Framework/interface/GeneratedInputSource.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <map>
#include <string>
#include "HepMC/GenEvent.h"

//#include "CLHEP/Random/JamesRandom.h"
//#include "CLHEP/Random/RandFlat.h"


class Run;
namespace CLHEP {
class HepRandomEngine;
class RandFlat;
}



namespace edm
{
  class PompytProducer : public EDProducer {
  public:

    /// Constructor
    PompytProducer(const ParameterSet &);
    /// Destructor
    virtual ~PompytProducer();

    void endRun( Run& r);

  private:

    /// Interface to the PYGIVE pythia routine, with add'l protections
        bool call_pygive(const std::string& iParm );

  private:
    
    virtual void produce(Event & e,const EventSetup& es);
    void clear();
    
    HepMC::GenEvent  *evt;
    
    /// Pythia PYLIST Verbosity flag
    unsigned int pythiaPylistVerbosity_;
    /// HepMC verbosity flag
    bool pythiaHepMCVerbosity_;
    /// Events to print if verbosity
    unsigned int maxEventsToPrint_;    
   
    // external cross section and filter efficiency
    double comenergy;
    double extCrossSect;
    double extFilterEff;    
 
    //pompyt
    double pomz;       //+-1 pomeron goes in +-z direction
    int pompdiss;      //1- elastic,2-pdiss
    double xi_max;     //max xi
    double pdiss_b0;   // B0 for pdiss events
    double pdiss_ap;   // alpha' for pdiss events
    double pdiss_glu;  // fraction of gluon in p-diss fragmentation 
    double pdiss_mmax; // max MY mass


    CLHEP::HepRandomEngine& fRandomEngine;
    CLHEP::RandFlat*        fRandomGenerator; 
    int                     eventNumber_;

  };
} 

#endif
