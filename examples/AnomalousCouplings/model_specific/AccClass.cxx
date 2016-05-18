#include"AccClass.h"
#include<stdlib.h>
#include<TGraph.h>
#include<math.h>

// for plotting and scaling the fitted functions
Double_t VLfunc(Double_t *x, Double_t *par);
Double_t VRfunc(Double_t *x, Double_t *par);
Double_t gLfunc(Double_t *x, Double_t *par);
Double_t gRfunc(Double_t *x, Double_t *par);

// ---------------------------------------------------------
AccClass::AccClass()
{
  double Ntotal = 300000.;
  vector<double> VL, VR, gL, gR;
  vector<double> Nselected;

  VL.push_back(0.); 
  VR.push_back(0.);
  gL.push_back(0.4);
  gR.push_back(0.);
  
  Nselected.push_back(34414.);
  
  VL.push_back(0.); 
  VR.push_back(0.);
  gL.push_back(0.4);
  gR.push_back(0.4);
  
  Nselected.push_back(30224.);
  
  VL.push_back(0.); 
  VR.push_back(0.);
  gL.push_back(0.);
  gR.push_back(0.4);
  
  Nselected.push_back(32301.);
  
  VL.push_back(0.); 
  VR.push_back(1.);
  gL.push_back(0.4);
  gR.push_back(0.);
  
  Nselected.push_back(32223.);
  
  VL.push_back(0.); 
  VR.push_back(1.);
  gL.push_back(0.4);
  gR.push_back(0.4);
  
  Nselected.push_back(30721.);
  
  VL.push_back(0.); 
  VR.push_back(1.);
  gL.push_back(0.);
  gR.push_back(0.);
  
  Nselected.push_back(26893.);
  
  VL.push_back(0.); 
  VR.push_back(1.);
  gL.push_back(0.);
  gR.push_back(0.4);
  
  Nselected.push_back(27964.);
  
  VL.push_back(1.); 
  VR.push_back(0.2);
  gL.push_back(0.);
  gR.push_back(0.);
  
  Nselected.push_back(28362.);
  
  VL.push_back(1.); 
  VR.push_back(0.4);
  gL.push_back(0.);
  gR.push_back(0.);
  
  Nselected.push_back(28972.);
  
  VL.push_back(1.); 
  VR.push_back(0.6);
  gL.push_back(0.);
  gR.push_back(0.);
  
  Nselected.push_back(28749.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.2);
  gR.push_back(0.);
  
  Nselected.push_back(28409.);

  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.4);
  gR.push_back(0.);
  
  Nselected.push_back(29043.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.6);
  gR.push_back(0.);
  
  Nselected.push_back(29720.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.4);
  gR.push_back(0.4);

  Nselected.push_back(29872.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.);
  gR.push_back(0.);

  Nselected.push_back(28232.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.);
  gR.push_back(0.2);
  
  Nselected.push_back(30861.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.);
  gR.push_back(0.4);
  
  Nselected.push_back(32269.);
 
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.);
  gR.push_back(0.6);
  
  Nselected.push_back(31683.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(0.);
  gR.push_back(1.);
  
  Nselected.push_back(30598.);
  
  VL.push_back(1.); 
  VR.push_back(0.);
  gL.push_back(1.);
  gR.push_back(0.);
  
  Nselected.push_back(31427.);
  
  VL.push_back(1.); 
  VR.push_back(1.);
  gL.push_back(0.4);
  gR.push_back(0.);
  
  Nselected.push_back(30002.);
  
  VL.push_back(1.); 
  VR.push_back(1.);
  gL.push_back(0.4);
  gR.push_back(0.4);
  
  Nselected.push_back(31744.);
  
  VL.push_back(1.); 
  VR.push_back(1.);
  gL.push_back(0.);
  gR.push_back(0.);
  
  Nselected.push_back(29073.);
  
  VL.push_back(1.); 
  VR.push_back(1.);
  gL.push_back(0.);
  gR.push_back(0.4);
  
  Nselected.push_back(29120.);

  // histograms with the center at the actual value
  VLVR = new TH2D("VLVR","V_{L}V_{R}",2,-0.5,1.5,6,-0.1,1.1);
  gLgR = new TH2D("gLgR","g_{L}g_{R}",6,-0.1,1.1,6,-0.1,1.1);
  VLgR = new TH2D("VLgR","V_{L}g_{R}",2,-0.5,1.5,6,-0.1,1.1);
  VLgL = new TH2D("VLgL","V_{L}g_{L}",2,-0.5,1.5,6,-0.1,1.1);
  VRgR = new TH2D("VRgR","V_{R}g_{R}",6,-0.1,1.1,6,-0.1,1.1);
  VRgL = new TH2D("VRgL","V_{R}g_{L}",6,-0.1,1.1,6,-0.1,1.1);
  
  // histograms with the center at the actual value
  hVL = new TH1D("hVL","V_{L}",2,-0.5,1.5);
  hVR = new TH1D("hVR","V_{R}",6,-0.1,1.1);
  hgL = new TH1D("hgL","g_{L}",6,-0.1,1.1);
  hgR = new TH1D("hgR","g_{R}",6,-0.1,1.1);

  for(unsigned int i = 0 ; i < Nselected.size() ; i++)
    {
      // fill the histograms
      if(gL.at(i) == 0. && gR.at(i) == 0.)
	VLVR -> Fill(VL.at(i),VR.at(i),Nselected.at(i)/Ntotal);
      if(VL.at(i) == 1. && VR.at(i) == 0.)
	gLgR -> Fill(gL.at(i),gR.at(i),Nselected.at(i)/Ntotal);
      if(gL.at(i) == 0. && VR.at(i) == 0.)
	VLgR -> Fill(VL.at(i),gR.at(i),Nselected.at(i)/Ntotal);
      if(gR.at(i) == 0. && VR.at(i) == 0.)
	VLgL -> Fill(VL.at(i),gL.at(i),Nselected.at(i)/Ntotal);
      if(VL.at(i) == 1. && gL.at(i) == 0.)
	VRgR -> Fill(VR.at(i),gR.at(i),Nselected.at(i)/Ntotal);
      if(VL.at(i) == 1. && gR.at(i) == 0.)
	VRgL -> Fill(VR.at(i),gL.at(i),Nselected.at(i)/Ntotal);

      if(VR.at(i) == 0. && gL.at(i) == 0. && gR.at(i) == 0.)
	hVL -> Fill(VL.at(i),Nselected.at(i)/Ntotal);
      if(VL.at(i) == 1. && gL.at(i) == 0. && gR.at(i) == 0.)
	hVR -> Fill(VR.at(i),Nselected.at(i)/Ntotal);
      if(VL.at(i) == 1. && VR.at(i) == 0. && gR.at(i) == 0.)
	hgL -> Fill(gL.at(i),Nselected.at(i)/Ntotal);
      if(VL.at(i) == 1. && gL.at(i) == 0. && VR.at(i) == 0.)
	hgR -> Fill(gR.at(i),Nselected.at(i)/Ntotal);
    }


  // option for fitting
  TString option = "Q";

  // put the histograms into graphs to avoid problems with zeroes
  // number of bins in the corresponding histograms
  const int NbinsVL = (hVL -> GetNbinsX());
  // fill array according to the histograms
  double aVL[NbinsVL], acceptanceVL[NbinsVL];
  int counterVL = 0; // number of valid entries to avoid problems with the gap
  for(int i = 1 ; i <= NbinsVL ; i++)
    {
      // keep the zero to make the Xsec vanish at VL=0
      aVL[counterVL]          = hVL->GetBinCenter(i);
      acceptanceVL[counterVL] = hVL->GetBinContent(i);
      counterVL ++;
    }

  gVL = new TGraph(counterVL,aVL,acceptanceVL);
  fVL = new TF1("fVL","[0] + [1]*x",-0.1,1.1);
  gVL -> Fit(fVL,option);

  // number of bins in the corresponding histograms
  const int NbinsVR = (hVR -> GetNbinsX());
  // fill array according to the histograms
  double aVR[NbinsVR], acceptanceVR[NbinsVR];
  int counterVR = 0; // number of valid entries to avoid problems with the gap
  for(int i = 1 ; i <= NbinsVR ; i++)
    {
      if(hVR->GetBinContent(i) != 0)
	{
	  aVR[counterVR]          = hVR->GetBinCenter(i);
	  acceptanceVR[counterVR] = hVR->GetBinContent(i);
	  counterVR ++;
	}     	
    }
  gVR = new TGraph(counterVR,aVR,acceptanceVR);
  fVR = new TF1("fVR","[0] + [1]*x + [2]*x**2",-0.1,1.1);
  gVR -> Fit(fVR,option);

  // number of bins in the corresponding histograms
  const int NbinsgL = (hgL -> GetNbinsX());
  // fill array according to the histograms
  double agL[NbinsgL], acceptancegL[NbinsgL];
  int countergL = 0; // number of valid entries to avoid problems with the gap
  for(int i = 1 ; i <= NbinsgL ; i++)
    {
      if(hgL->GetBinContent(i) != 0)
	{
	  agL[countergL]          = hgL->GetBinCenter(i);
	  acceptancegL[countergL] = hgL->GetBinContent(i);
	  countergL ++;
	}     	
    }
  ggL = new TGraph(countergL,agL,acceptancegL);
  fgL = new TF1("fgL","[0] + [1]*x + [2]*x**2",-0.1,1.1);
  ggL -> Fit(fgL,option);


  // number of bins in the corresponding histograms
  const int NbinsgR = (hgR -> GetNbinsX());
  // fill array according to the histograms
  double agR[NbinsgR], acceptancegR[NbinsgR];
  int countergR = 0; // number of valid entries to avoid problems with the gap
  for(int i = 1 ; i <= NbinsgR ; i++)
    {
      if(hgR->GetBinContent(i) != 0)
	{
	  agR[countergR]          = hgR->GetBinCenter(i);
	  acceptancegR[countergR] = hgR->GetBinContent(i);
	  countergR ++;
	}    	
    }
  ggR = new TGraph(countergR,agR,acceptancegR);
  fgR = new TF1("fgR","[0] + [1]*x + [2]*x**2",-0.1,1.1);
  ggR -> Fit(fgR,option);

  parVL[0] = fVL->GetParameter(0);
  parVL[1] = fVL->GetParameter(1);

  parVR[0] = fVR->GetParameter(0);
  parVR[1] = fVR->GetParameter(1);
  parVR[2] = fVR->GetParameter(2);

  pargL[0] = fgL->GetParameter(0);
  pargL[1] = fgL->GetParameter(1);
  pargL[2] = fgL->GetParameter(2);

  pargR[0] = fgR->GetParameter(0);
  pargR[1] = fgR->GetParameter(1);
  pargR[2] = fgR->GetParameter(2);
}

// ---------------------------------------------------------
AccClass::~AccClass()
{
}

// ---------------------------------------------------------
void AccClass::DrawControlPlots(TString folder)
{
  gStyle->SetCanvasColor(0);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.4f");  // What precision to put numbers if plotted with "TEXT"
    
  VLVR->SetFillColor(0);
  gLgR->SetFillColor(0);
  VLgR->SetFillColor(0);
  VLgL->SetFillColor(0);
  VRgR->SetFillColor(0);
  VRgL->SetFillColor(0);
  
  hVL->SetFillColor(0);
  hgL->SetFillColor(0);
  hgR->SetFillColor(0);
  hVR->SetFillColor(0);

  gVL->SetFillColor(0);
  ggL->SetFillColor(0);
  ggR->SetFillColor(0);
  gVR->SetFillColor(0);

  fVL->SetFillColor(0);
  fgL->SetFillColor(0);
  fgR->SetFillColor(0);
  fVR->SetFillColor(0);
  
  if(folder != "")
    folder += "/";

  TCanvas *c1 = new TCanvas("c1","2d histos",1200, 900);
  c1 -> Divide(3,2);

  TCanvas *c2 = new TCanvas("c2","1d histos",900, 900);
  c2 -> Divide(2,2);
   
  Double_t levels[] = {0.080,0.082,0.084,0.086,0.088,0.090,0.092,0.094,0.096,0.098,0.100,0.102,0.104,0.106,0.108,0.110,0.112,0.114,0.116,0.120};
 
  c1 -> cd(1);
  VLVR -> GetXaxis() -> SetTitle("V_{L}");
  VLVR -> GetYaxis() -> SetTitle("V_{R}");
  VLVR -> SetContour(20, levels);
  VLVR -> SetMarkerSize(2); //size of the numbers in the histogram
  VLVR -> Draw("COL,TEXT");
  
  c1 -> cd(2);
  gLgR -> GetXaxis() -> SetTitle("g_{L}");
  gLgR -> GetYaxis() -> SetTitle("g_{R}");
  gLgR -> SetContour(20, levels);
  gLgR -> SetMarkerSize(2); //size of the numbers in the histogram
  gLgR -> Draw("COL,TEXT");

  c1 -> cd(3);
  VLgR -> GetXaxis() -> SetTitle("V_{L}");
  VLgR -> GetYaxis() -> SetTitle("g_{R}");
  VLgR -> SetContour(20, levels);
  VLgR -> SetMarkerSize(2); //size of the numbers in the histogram
  VLgR -> Draw("COL,TEXT");

  c1 -> cd(4);
  VLgL -> GetXaxis() -> SetTitle("V_{L}");
  VLgL -> GetYaxis() -> SetTitle("g_{L}");
  VLgL -> SetContour(20, levels);
  VLgL -> SetMarkerSize(2); //size of the numbers in the histogram
  VLgL -> Draw("COL,TEXT");

  c1 -> cd(5);
  VRgR -> GetXaxis() -> SetTitle("V_{R}");
  VRgR -> GetYaxis() -> SetTitle("g_{R}");
  VRgR -> SetContour(20, levels);
  VRgR -> SetMarkerSize(2); //size of the numbers in the histogram
  VRgR -> Draw("COL,TEXT");

  c1 -> cd(6);
  VRgL -> GetXaxis() -> SetTitle("V_{R}");
  VRgL -> GetYaxis() -> SetTitle("g_{L}");
  VRgL -> SetContour(20, levels);
  VRgL -> SetMarkerSize(2); //size of the numbers in the histogram
  VRgL -> Draw("COL,TEXT");
      
  c1 -> Print(folder + "2d.pdf");
  c1 -> Print(folder + "2d.png");

  hVL  -> GetXaxis() -> SetTitle("V_{L}");
  hVL  -> GetYaxis() -> SetTitle("acceptance");
  hVR  -> GetXaxis() -> SetTitle("V_{R}");
  hVR  -> GetYaxis() -> SetTitle("acceptance");
  hgL  -> GetXaxis() -> SetTitle("g_{L}");
  hgL  -> GetYaxis() -> SetTitle("acceptance");
  hgR  -> GetXaxis() -> SetTitle("g_{R}");
  hgR  -> GetYaxis() -> SetTitle("acceptance");

  Float_t TitleSize = 0.08;
  hVL  -> GetXaxis() -> SetTitleSize(TitleSize);
  hVL  -> GetYaxis() -> SetTitleSize(TitleSize);
  hVR  -> GetXaxis() -> SetTitleSize(TitleSize);
  hVR  -> GetYaxis() -> SetTitleSize(TitleSize);
  hgL  -> GetXaxis() -> SetTitleSize(TitleSize);
  hgL  -> GetYaxis() -> SetTitleSize(TitleSize);
  hgR  -> GetXaxis() -> SetTitleSize(TitleSize);
  hgR  -> GetYaxis() -> SetTitleSize(TitleSize);

  gVL -> SetMarkerStyle(21);
  gVL -> SetMarkerSize(2);
  c2  -> cd(1);
  gVL -> SetTitle("");
  gVL -> GetXaxis() -> SetTitle("V_{L}");
  // gVL -> GetXaxis() -> SetLimits(-1.,1.);
  gVL -> GetYaxis() -> SetTitle("acceptance   ");
  gVL -> GetYaxis() -> SetTitleOffset(-0.9);
  gVL -> Draw("AP");
  fVL -> Draw("same");

  gVR -> SetMarkerStyle(21);
  gVR -> SetMarkerSize(2);
  c2  -> cd(2);
  gVR -> SetTitle("");
  gVR -> GetXaxis() -> SetTitle("V_{R}");
  // gVR -> GetXaxis() -> SetLimits(-1.,1.);
  gVR -> GetYaxis() -> SetTitle("acceptance   ");
  gVR -> GetYaxis() -> SetTitleOffset(-0.9);
  gVR -> Draw("AP");
  fVR -> Draw("same");

  ggL -> SetMarkerStyle(21);
  ggL -> SetMarkerSize(2);
  c2  -> cd(3);
  ggL -> SetTitle("");
  ggL -> GetXaxis() -> SetTitle("g_{L}");
  // ggL -> GetXaxis() -> SetLimits(-1.,1.);
  ggL -> GetYaxis() -> SetTitle("acceptance   ");
  ggL -> GetYaxis() -> SetTitleOffset(-0.9);
  ggL -> Draw("AP");
  fgL -> Draw("same");

  ggR -> SetMarkerStyle(21);
  ggR -> SetMarkerSize(2);
  c2  -> cd(4);
  ggR -> SetTitle("");
  ggR -> GetXaxis() -> SetTitle("g_{R}");
  // ggR -> GetXaxis() -> SetLimits(-1.,1.);
  ggR -> GetYaxis() -> SetTitle("acceptance   ");
  ggR -> GetYaxis() -> SetTitleOffset(-0.9);
  ggR -> Draw("AP");
  fgR -> Draw("same");

  TF1 *fVL2 = new TF1("fVL2",VLfunc,-1.1,1.1,0);
  TF1 *fVR2 = new TF1("fVR2",VRfunc,-1.1,1.1,0);
  TF1 *fgL2 = new TF1("fgL2",gLfunc,-1.1,1.1,0);
  TF1 *fgR2 = new TF1("fgR2",gRfunc,-1.1,1.1,0);

  fVL2 -> SetTitle("");
  fVL2 -> GetXaxis() -> SetTitle("V_{L}");
  fVL2 -> GetYaxis() -> SetTitle("acceptance normed to SM  ");
  fVL2 -> GetYaxis() -> SetTitleOffset(-0.9);

  fVR2 -> SetTitle("");
  fVR2 -> GetXaxis() -> SetTitle("V_{R}");
  fVR2 -> GetYaxis() -> SetTitle("acceptance normed to SM  ");
  fVR2 -> GetYaxis() -> SetTitleOffset(-0.9);

  fgL2 -> SetTitle("");
  fgL2 -> GetXaxis() -> SetTitle("g_{L}");
  fgL2 -> GetYaxis() -> SetTitle("acceptance normed to SM  ");
  fgL2 -> GetYaxis() -> SetTitleOffset(-0.9);
  
  fgR2 -> SetTitle("");
  fgR2 -> GetXaxis() -> SetTitle("g_{R}");
  fgR2 -> GetYaxis() -> SetTitle("acceptance normed to SM  ");
  fgR2 -> GetYaxis() -> SetTitleOffset(-0.9);

  TCanvas *c10 = new TCanvas("c10","fitted functions",900, 900);
  c10 -> Divide(2,2);
  // fVL -> GetYaxis() -> SetLimits(-0.1,0.1);
  // // fVR2 -> GetYaxis() -> SetLimits(-0.1,0.1);
  // fgL -> GetYaxis() -> SetLimits(-0.1,0.1);
  // fgR -> GetYaxis() -> SetLimits(-0.1,0.1);
  c10 -> cd(1);
  fVL2 -> Draw();
  c10 -> cd(2);
  fVR2 -> Draw();
  c10 -> cd(3);
  fgL2 -> Draw();
  c10 -> cd(4);
  fgR2 -> Draw();

  c2  -> Print(folder + "1d.pdf");
  c2  -> Print(folder + "1d.png");
  c10 -> Print(folder + "FittedFunctions.png"); 
  c10 -> Print(folder + "FittedFunctions.pdf"); 
}

// ---------------------------------------------------------
Double_t VLfunc(Double_t *x, Double_t *par)
{
  double xx =x[0];
  //Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
  AccClass dummy;
  Double_t f = dummy.GetValue(xx,0.,0.,0.,"quadratic",false)/dummy.GetValue(1.,0.,0.,0.,"quadratic",false);
  return f;
}
Double_t VRfunc(Double_t *x, Double_t *par)
{
  double xx =x[0];
  //Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
  AccClass dummy;
  Double_t f = dummy.GetValue(1.,xx,0.,0.,"quadratic",false)/dummy.GetValue(1.,0.,0.,0.,"quadratic",false);
  return f;
}
Double_t gLfunc(Double_t *x, Double_t *par)
{
  double xx =x[0];
  //Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
  AccClass dummy;
  Double_t f = dummy.GetValue(1.,0.,xx,0.,"quadratic",false)/dummy.GetValue(1.,0.,0.,0.,"quadratic",false);
  return f;
}
Double_t gRfunc(Double_t *x, Double_t *par)
{
  double xx =x[0];
  //Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
  AccClass dummy;
  Double_t f = dummy.GetValue(1.,0.,0.,xx,"quadratic",false)/dummy.GetValue(1.,0.,0.,0.,"quadratic",false);
  return f;
}

// ---------------------------------------------------------
double AccClass::GetValue(double vl, double vr, double gl, double gr, string method, bool IsSymmetrized)
{
  if(IsSymmetrized)
    {
      if(!method.compare("quadratic"))
	{
	  // assume that the acceptance at VL=1 and VR=gL=gR=0 is equal
	  return parVL[1]*fabs(vl) + parVR[1]*fabs(vr) + pargL[1]*fabs(gl) + pargR[1]*fabs(gr) 
	    + parVR[2]*vr*vr + pargL[2]*gl*gl + pargR[2]*gr*gr;
	}
      else if(!method.compare("linear"))
	{
	  double slopeVL =  hVR->GetBinContent(1);
	  double slopeVR = ((hVR->GetBinContent(2)) - (hVR->GetBinContent(1)))/0.2;    
	  double slopegL = ((hgL->GetBinContent(2)) - (hgL->GetBinContent(1)))/0.2;    
	  double slopegR = ((hgR->GetBinContent(2)) - (hgR->GetBinContent(1)))/0.2;
	  
	  return slopeVL*fabs(vl) + slopeVR*fabs(vr) + slopegL*fabs(gl) + slopegR*fabs(gr);
	}
      else if(!method.compare("const"))
	return hVR->GetBinContent(1);
      else
	{
	  std::cout << "Error:\t\"" << method << "\" is not a method for approximating the acceptance. Use \"quadratic\" \"linear\" or \"const\"." 
		    << std::endl;
	  exit(1);
	}
    }
  else
    {
      if(!method.compare("quadratic"))
	{
	  // assume that the acceptance at VL=1 and VR=gL=gR=0 is equal
	  // return parVL[1]*vl + parVR[1]*vr + pargL[1]*gl + pargR[1]*gr  // linear acceptance for VL
	  //   + parVR[2]*vr*vr + pargL[2]*gl*gl + pargR[2]*gr*gr;
	  return parVL[1] + parVR[1]*vr + pargL[1]*gl + pargR[1]*gr   // constant acceptance for VL
	    + parVR[2]*vr*vr + pargL[2]*gl*gl + pargR[2]*gr*gr;
	}
      else if(!method.compare("linear"))
	{
	  double slopeVL =  hVR->GetBinContent(1);
	  double slopeVR = ((hVR->GetBinContent(2)) - (hVR->GetBinContent(1)))/0.2;    
	  double slopegL = ((hgL->GetBinContent(2)) - (hgL->GetBinContent(1)))/0.2;    
	  double slopegR = ((hgR->GetBinContent(2)) - (hgR->GetBinContent(1)))/0.2;
	  
	  return slopeVL + slopeVR*vr + slopegL*gl + slopegR*gr;// constant acceptance for VL
	}
      else if(!method.compare("const"))
	return hVR->GetBinContent(1);
      else
	{
	  std::cout << "Error:\t\"" << method << "\" is not a method for approximating the acceptance. Use \"quadratic\" \"linear\" or \"const\"." 
		    << std::endl;
	  exit(1);
	}
    }
}
