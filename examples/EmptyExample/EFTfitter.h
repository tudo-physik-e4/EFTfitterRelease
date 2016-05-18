#ifndef __EFTFITTER__H
#define __EFTFITTER__H

#include "../../src/mvclocal/BCMVCPhysicsModel_EFT.h"

// ---------------------------------------------------------
class EFTfitter : public BCMVCPhysicsModel_EFT
{
 public:

  // Constructor
  EFTfitter();

  // Destructor
  ~EFTfitter();

	// return a value for an observable
	// index: the index of the observable
	// parameters: the physics parameters
  double CalculateObservable(int index, const std::vector<double> &parameters);

 private:
  
};

#endif

