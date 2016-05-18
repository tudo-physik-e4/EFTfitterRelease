/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#ifndef __BCMVCPhysicsModel_EFT__H
#define __BCMVCPhysicsModel_EFT__H

#include "BCMVCombination_EFT.h"
#include "BAT/BCParameter.h"

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "BCMVCObservable.h"

#include "../tinyxml2/tinyxml2.h"
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>


// ---------------------------------------------------------
class BCMVCPhysicsModel_EFT : public BCMVCombination_EFT
{
 public:

  // Constructor
  BCMVCPhysicsModel_EFT();

  // Destructor
  ~BCMVCPhysicsModel_EFT();

  // Add a parameter
  // name: the name of the parameter
  // min:  the minimum value of the parameter
  // max:  the maximum value of the parameter
  void AddObservable(std::string name, double min, double max);

  // return a value for an observable
  // index: the index of the variable
  // parameters: the physics parameters
  virtual double CalculateObservable(int index, const std::vector<double> &parameters)
  { (void) index; // suppress compiler warning about unused parameters
    (void) parameters; // suppress compiler warning about unused parameters
    return 0; };
    

  // the log of the likelihood
  double LogLikelihood(const std::vector<double> &parameters);

  void ReadParameters(std::string file);

// plot all observables as a function of all parameters
  void PlotObservables(std::string filename);

// read input file
  void ReadMeasurements(std::string filename);

  void AddNuisanceParameter(std::string filename);

// ranking of measurements and uncertainties
  void GetMeasurementRanking(double p, std::string file); 
  void GetUncertaintyRanking(double p, std::string file); 
  void GetRanking(std::string flag, double p, std::string file); // flag: uncertainty/measurement
  void PrintRanking1D(std::string flag, std::ofstream& out, double p,std::vector<double> orig_areas1D, std::vector< std::vector<std::pair<int,double> > > areas1D);
  void PrintRanking2D(std::string flag, std::ofstream& out, double p,std::vector<double> orig_areas2D, std::vector< std::vector<std::pair<int,double> > > areas2D);
  
  double GetArea1D(int index, double p);
  double GetArea2D(int i1, int i2,double p);

  // functions for sorting the measurements/uncertainties ranks
  static bool pairCompare1(const std::pair<int,double>& a, const std::pair<int,double>& b);
  static bool pairCompare2(const std::pair<int,double>& a, const std::pair<int,double>& b);

//----------------------------------------------------------

 private:
    std::vector<double> Parameters_SM;

};
// ---------------------------------------------------------

#endif

