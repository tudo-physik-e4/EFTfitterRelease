#ifndef __XSECCOEFFS__H
#define __XSECCOEFFS__H

#include "../../src/mvclocal/BCMVCPhysicsModel_EFT.h"

// ---------------------------------------------------------
class XsecCoeffs
{
 public:

  // Constructor
  XsecCoeffs();

  // Destructor
  ~XsecCoeffs();

  // read in the coefficients from a file
  // FileName: name of the file the coefficients should be read from
  // factor: the overall factor the cross section should be multiplied with
  // factorError: the corresponding error
  void ReadIn(char FileName[80], double factor, double factorError);
   
  // returns a vector of the coefficients 
  std::vector<double> GetCoeffsValues();

  // returns a vector of the errors of the coefficients 
  std::vector<double> GetCoeffsErrors();

  // returns the factor (slot 0) and the corresponding error (slot 1) the 
  // cross section should be multiplied with
  std::vector<double> GetFactor();  

 private:

  std::vector<double> value, error;
  std::vector<double> XsecFactor;
};
// ---------------------------------------------------------
#endif
