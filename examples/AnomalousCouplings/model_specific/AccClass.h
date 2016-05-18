#ifndef __ACCCLASS__H
#define __ACCCLASS__H

#include<iostream> 
#include<TAxis.h>
#include<TH1.h>
#include<TH2.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TFunction.h>
#include<TF1.h>
#include<TFitResultPtr.h>

class TGraph;

using namespace std;

// ---------------------------------------------------------
class AccClass
{
 public:

  // Constructor
  AccClass();

  // Destructor
  ~AccClass();

  // draws control plots
  void DrawControlPlots(TString folder = "");

  // returns value of acceptance with the couplings and the defined method
  // the methods are "quadratic", "linear" and "const"
  double GetValue(double vl = 1., double vr = 0., double gl = 0., double gr = 0., string method = "quadratic", bool IsSymmetrized = true);

 private:

  //Double_t VRfunc(Double_t *x, Double_t *par);

  // histograms with the center at the actual value
  TH2D *VLVR;
  TH2D *gLgR;
  TH2D *VLgR;
  TH2D *VLgL;
  TH2D *VRgR;
  TH2D *VRgL;

  // histograms with the center at the actual value
  TH1D *hVL;
  TH1D *hVR;
  TH1D *hgL;
  TH1D *hgR;

  // graphs with values from the histograms
  TGraph *gVL;
  TGraph *gVR;
  TGraph *ggL;
  TGraph *ggR;
  
  // functions resulting from the fit
  TF1 *fVL;
  TF1 *fVR;
  TF1 *fgL;
  TF1 *fgR;

  // parameters of the fit result
  Double_t parVL[2];
  Double_t parVR[3];
  Double_t pargL[3];
  Double_t pargR[3];

};
// ---------------------------------------------------------
#endif
