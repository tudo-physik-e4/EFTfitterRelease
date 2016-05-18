#ifndef __ANOMCOUPLINGS__H
#define __ANOMCOUPLINGS__H

#include "model_specific/AccClass.h"
#include "../../src/mvclocal/BCMVCPhysicsModel_EFT.h"

// ---------------------------------------------------------
class AnomCouplings : public BCMVCPhysicsModel_EFT
{
 public:

  // Constructor
  // are there two cross sections for top and antitop? 
  // else: they are merged
  AnomCouplings();

  // Destructor
  ~AnomCouplings();

  // initialize physics parameters
  void InitializePhysicsParameters();
  
  // for an input that has the cross section
  // for top and antitop merged
  void MergedXsecs()
  {TwoXsecs = false;}

  // no cross section in the input file
  void SuppressXsecs()
  {NoXsecs = true;};

  // no W helicity in the input file
  void SuppressWhelicity()
  {NoWhel = true;};

	// return a value for an observable
	// index: the index of the observable
	// parameters: the physics parameters
  double CalculateObservable(int index, const std::vector<double> &parameters);

  // calculates the cross section 
  double CalculateCrossection(double vl, double vr, double gl, double gr, double vo, double sqrts);
  double CalculateCrossectionTop(double vl, double vr, double gl, double gr, double vo, double sqrts);
  double CalculateCrossectionAntitop(double vl, double vr, double gl, double gr, double vo, double sqrts);

  // Returns a corresponding parameter value for fixed observables. 
  // This doesn't take into account the "real" input observables
  // and is for cross checking purposes.
  // name = "F0" returns F0
  // name = "FL" returns FL
  // name = "XS" returns the cross section
  double PrintObservable(TString name, double vl, double vr, double gl, double gr, double sqrts);  

  // calculate the W-boson helicity fractions as a function of
  // the anomalous couplings
  // parameters: the fractions and the couplings
  void HelicityFractions(double &f0, double &fl, double &fr, const std::vector<double> &parameters);

  // calculate the phase space for the single top (ST) cross-section
  double STPhaseSpace (double m, double sqrts); 

  // calculate the correction of the cross section for the change of 
  // efficiency/acceptance due to anomalous couplings
  double CorrEfficiency(double vl, double vr, double gl, double gr, double vo, double sqrts);

  // disables (flag = true) or enables (false) the measurements
  void DisableWhel(bool flag = true);
  void DisableXsec(bool flag = true);


 private:

  // this is for the value of the acceptance
  // in dependence on the anomalous couplings
  AccClass acceptance;

  // method for determining the acceptance effect
  // the methods are "quadratic", "linear" and "const"
  string AccMethod;
  
  // physics constants

  // top-quark pole mass
  double fMtop;
  double fMtop2;

  // bottom-quark pole mass
  double fMb;
  double fMb2;

  // W-boson pole mass
  double fMW;
  double fMW2;

  // weak coupling constant
  double fGWeak;
  double fGWeak4;

  // center-of-mass energy in GeV
  double fSqrtS;
  double fSqrtS2;

  // single top cross-sections
  double fXS_ud_tb;
  double fXS_bu_td;
  double fXS_bd_tu;

  // are there two cross sections for top and antitop?
  // else: they are merged
  bool TwoXsecs;

  // no cross section in the input file
  bool NoXsecs;

  // no W helicity in the input file
  bool NoWhel;

  // is the acceptance symmetrized at VL=VR=gL=gR=0 or not?
  bool IsSymmetrized;

};

#endif

