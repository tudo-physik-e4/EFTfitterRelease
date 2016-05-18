#include "EFTfitter.h"

#include "BAT/BCMath.h"
#include "BAT/BCLog.h"
#include "BAT/BCParameter.h"
#include "../../src/mvclocal/BCMVCObservable.h"

#include <cmath>

// ---------------------------------------------------------
EFTfitter::EFTfitter() : BCMVCPhysicsModel_EFT()
{ //constructor
}

// ---------------------------------------------------------
EFTfitter::~EFTfitter()
{ //destructor
}

// ---------------------------------------------------------
// the theoretical model needs to be implemented here
double EFTfitter::CalculateObservable(int index, const std::vector<double> &parameters){

  //delete the following line when implementing the model:
  throw std::logic_error("\n\nPlease define your model in \"CalculateObservable()\" in file \"EFTfitter.cxx\" \n");
  
  //return something;

}
