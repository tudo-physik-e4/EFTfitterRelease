#include "XsecCoeffs.h"

#include <fstream>

// ---------------------------------------------------------
XsecCoeffs::XsecCoeffs()
{
}

// ---------------------------------------------------------
XsecCoeffs::~XsecCoeffs()
{
}

// ---------------------------------------------------------
void XsecCoeffs::ReadIn(char FileName[80], double factor, double factorError)
{
  ifstream ifile(FileName);
  double val, err;
  
  XsecFactor.push_back(factor);
  XsecFactor.push_back(factorError);

  Int_t i = 0;
  while(ifile >> val >> err)
    {
      if(i<9)
	{
	  value.push_back(val);
	  error.push_back(err);
	
	  i++;
	}
    }

  i=0;
  ifile.close();
}

// ---------------------------------------------------------
std::vector<double> XsecCoeffs::GetCoeffsValues()
{
  return value;
}

// ---------------------------------------------------------
std::vector<double> XsecCoeffs::GetCoeffsErrors()
{
  return error;
}

// ---------------------------------------------------------
std::vector<double> XsecCoeffs::GetFactor()
{
  return XsecFactor;
}
