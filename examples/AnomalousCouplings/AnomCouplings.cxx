#include "AnomCouplings.h"

#include "model_specific/XsecCoeffs.h"
#include "model_specific/Observables.h"

#include "BAT/BCMath.h"
#include "BAT/BCLog.h"
#include "BAT/BCParameter.h"
#include "../../src/mvclocal/BCMVCObservable.h"

#include <cmath>

// ---------------------------------------------------------
AnomCouplings::AnomCouplings() : BCMVCPhysicsModel_EFT()
{
  TwoXsecs  = true;  // top- and antitop cross sections are splitted
  NoXsecs   = false; // there are cross sections in the input file
  NoWhel    = false; // there are cross sections in the input file 

  // method for determining the acceptance effect
  // the methods are "quadratic", "linear" and "const"
  AccMethod = "const";

  // is the acceptance symmetrized at VL=VR=gL=gR=0 or not?
  IsSymmetrized = true;
}

// ---------------------------------------------------------
AnomCouplings::~AnomCouplings()
{
}

// ---------------------------------------------------------
void AnomCouplings::InitializePhysicsParameters()
{
  fMtop = 172.5;
  fMtop2 = fMtop * fMtop;

  fMb = 4.8; 
  fMb2 = fMb * fMb;

  fMW = 80.4; 
  fMW2 = fMW * fMW;

  fGWeak = 0.65; //sqrt(8/sqrt(2)*FermiConstant)*WBosonMass
  fGWeak4 = fGWeak * fGWeak * fGWeak * fGWeak;

  fSqrtS = 7000.;
  fSqrtS2 = fSqrtS * fSqrtS;

  fXS_ud_tb = 1;
  fXS_bu_td = 1;
  fXS_bd_tu = 1;
}

// ---------------------------------------------------------
double AnomCouplings::CalculateObservable(int index, const std::vector<double> &parameters)
{
  // calculate the W-boson helicity fractions
  double f0 = 0;
  double fl = 0;
  double fr = 0;
  double xs_tchannel = 0;
  double xs_tchannelTop = 0;
  double xs_tchannelAntitop = 0;

  if ((index == 0 || index == 1) && !NoWhel)
    HelicityFractions(f0, fl, fr, parameters);

  else if ((index == 2 && !TwoXsecs && !NoXsecs && !NoWhel) || (index == 0 && !TwoXsecs && !NoXsecs && NoWhel)) {
    double vl = parameters[0];
    double vr = parameters[1];
    double gl = parameters[2];
    double gr = parameters[3];;
    
    // acceptance effect (= acceptance(VL,VR,gL,gR)/acceptance(1,0,0,0))
    double AcceptCorr = acceptance.GetValue(vl,vr,gl,gr,AccMethod,IsSymmetrized) / acceptance.GetValue(1.,0.,0.,0.,AccMethod,IsSymmetrized);

    xs_tchannel = CalculateCrossection(vl, vr, gl, gr, 0., fSqrtS)*AcceptCorr;
  }

  else if ((index == 2 && TwoXsecs && !NoXsecs && !NoWhel) || (index == 0 && TwoXsecs && !NoXsecs && NoWhel)) {
    double vl = parameters[0];
    double vr = parameters[1];
    double gl = parameters[2];
    double gr = parameters[3];
    
    // acceptance effect (= acceptance(VL,VR,gL,gR)/acceptance(1,0,0,0))
    double AcceptCorr = acceptance.GetValue(vl,vr,gl,gr,AccMethod,IsSymmetrized) / acceptance.GetValue(1.,0.,0.,0.,AccMethod,IsSymmetrized);

    xs_tchannelTop = CalculateCrossectionTop(vl, vr, gl, gr, 0., fSqrtS)*AcceptCorr;
  }

  else if ((index == 3 && TwoXsecs && !NoXsecs && !NoWhel) || (index == 1 && TwoXsecs && !NoXsecs && NoWhel)) {
    double vl = parameters[0];
    double vr = parameters[1];
    double gl = parameters[2];
    double gr = parameters[3];
    
    // acceptance effect (= acceptance(VL,VR,gL,gR)/acceptance(1,0,0,0))
    double AcceptCorr = acceptance.GetValue(vl,vr,gl,gr,AccMethod,IsSymmetrized) / acceptance.GetValue(1.,0.,0.,0.,AccMethod,IsSymmetrized);

    xs_tchannelAntitop = CalculateCrossectionAntitop(vl, vr, gl, gr, 0., fSqrtS)*AcceptCorr;
  }

  // return the index of the observable under study
  if (index == 0 && !NoWhel) {
    return f0;
  }
  else if (index == 1 && !NoWhel) {
    return fl;
  }
  else if ((index == 2 && !TwoXsecs && !NoXsecs && !NoWhel) || (index == 0 && !TwoXsecs && !NoXsecs && NoWhel)) {
    return xs_tchannel;
  }
  else if ((index == 2 && TwoXsecs && !NoXsecs && !NoWhel) || (index == 0 && TwoXsecs && !NoXsecs && NoWhel)) {
    return xs_tchannelTop;
  }
  else if ((index == 3 && TwoXsecs && !NoXsecs && !NoWhel) || (index == 1 && TwoXsecs && !NoXsecs && NoWhel)) {
    return xs_tchannelAntitop;
  }
  else {
    std::cout << "ERROR: index #" << index << " does not correspond to an Observable!" << std::endl;
    BCLog::OutWarning("CalculateObservable. Index of observable out of range. Return 0");
    return 0;
  }
}

// ---------------------------------------------------------
double AnomCouplings::CalculateCrossection(double vl, double vr, double gl, double gr, double vo, double sqrts)
{
  // top+antitop 2->2
  return SingleTXsec(11,vl,vr,gl,gr,vo)[0] + SingleTXsec(21,vl,vr,gl,gr,vo)[0]; 
}
// ---------------------------------------------------------
double AnomCouplings::CalculateCrossectionTop(double vl, double vr, double gl, double gr, double vo, double sqrts)
{
  // top 2->2
  return SingleTXsec(11,vl,vr,gl,gr,vo)[0]; 
}
// ---------------------------------------------------------
double AnomCouplings::CalculateCrossectionAntitop(double vl, double vr, double gl, double gr, double vo, double sqrts)
{

  // antitop 2->2
  return SingleTXsec(21,vl,vr,gl,gr,vo)[0]; 
}
// ---------------------------------------------------------
double AnomCouplings::PrintObservable(TString name, double vl, double vr, double gl, double gr, double sqrts) 
{
  // acceptance effect (= acceptance(VL,VR,gL,gR)/acceptance(1,0,0,0))
  double AcceptCorr = 1.;//acceptance.GetValue(vl,vr,gl,gr,AccMethod,IsSymmetrized) / acceptance.GetValue(1.,0.,0.,0.,AccMethod,IsSymmetrized);

  double f0 = 0;
  double fl = 0;
  double fr = 0;
  double xs_tchannel = 0;
  double xs_tchannelTop = 0;
  double xs_tchannelAntitop = 0;
  double Br = 0;
  std::vector<double> SMCouplings;

  if(name == "F0" || name == "FL" || name == "FR")
    {
      SMCouplings.push_back(vl); 
      SMCouplings.push_back(vr); 
      SMCouplings.push_back(gl); 
      SMCouplings.push_back(gr); 
      
      HelicityFractions(f0,fl,fr,SMCouplings);
    }
  else if(name == "XS") // top+antitop Xsec
    xs_tchannel = CalculateCrossection(vl, vr, gl, gr, 0., sqrts)*AcceptCorr;
  else if(name == "XT") // only top Xsec
    xs_tchannelTop = CalculateCrossectionTop(vl, vr, gl, gr, 0., sqrts)*AcceptCorr;
  else if(name == "XA") // only antitop Xsec
    xs_tchannelAntitop = CalculateCrossectionAntitop(vl, vr, gl, gr, 0., sqrts)*AcceptCorr;

  if(name == "F0") return f0;  // in %
  else if(name == "FL") return fl;  // in %
  else if(name == "FR") return fr;  // in %
  else if(name == "XS") return xs_tchannel;
  else if(name == "XT") return xs_tchannelTop;  // in pb
  else if(name == "XA") return xs_tchannelAntitop;  // in pb
  else if(name == "BR") return Br;  // 10^-5
 
  else 
    {
      std::cout << "Error: Observable with the name " << name << " doesn't exist!" << std::endl;
      return 0;
    }
}

// ---------------------------------------------------------
void AnomCouplings::HelicityFractions(double &f0, double &fl, double &fr, const std::vector<double> &parameters)  // in %
{
  f0 = HelicityFracs( 0,parameters[0], parameters[1], parameters[2], parameters[3],fMtop,fMb,fMW);
  fl = HelicityFracs(-1,parameters[0], parameters[1], parameters[2], parameters[3],fMtop,fMb,fMW);
  fr = HelicityFracs( 1,parameters[0], parameters[1], parameters[2], parameters[3],fMtop,fMb,fMW);

}

// ---------------------------------------------------------
double AnomCouplings::STPhaseSpace(double m, double sqrts)
{
  double s  = pow(sqrts,2.);
  double m2 = pow(m,2.);
  double pi = 3.14159265359;
  double gev2fb = 3.89379338*1E11;
  return gev2fb/(s*32.*pi)*(1-m2/s);
}


// ---------------------------------------------------------
void AnomCouplings::DisableWhel(bool flag)
{
  NoWhel = flag;
}

// ---------------------------------------------------------
void AnomCouplings::DisableXsec(bool flag)
{
  NoXsecs = flag;
}



