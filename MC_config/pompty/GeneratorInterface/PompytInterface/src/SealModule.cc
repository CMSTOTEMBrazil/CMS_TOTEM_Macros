#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//#include "PompytSource.h"
#include "../interface/PompytProducer.h"

//using edm::PompytSource;
using edm::PompytProducer;
DEFINE_FWK_MODULE(PompytProducer);


//  DEFINE_SEAL_MODULE();
//  DEFINE_ANOTHER_FWK_MODULE(PompytProducer);
//  DEFINE_ANOTHER_FWK_INPUT_SOURCE(PompytSource);
//  DEFINE_ANOTHER_FWK_MODULE(HepMCProductAnalyzer);

