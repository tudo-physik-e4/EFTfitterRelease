// this code calculates the observables according to
// F. Bach, T. Ohl, Phys. Rev. D 90, 113007 (2014)
// changed the normalization in that way that the plots are consistent with the approximate NNLO values

#include <iostream>

#include<TMath.h>


// -----------------------------------------------------------------------------------------------
// calculates the ST cross-section in fb
// first entry ([0]) is the value, second entry ([1]) is the error/uncertainty
// type = 11 calculates from single_t_t-channel_onshell_2to2
// type = 12 calculates from single_t_t-channel_onshell_2to3
// type = 13 calculates from single_t_t-channel_onshell_matched
// type = 21 calculates from single_tbar_t-channel_onshell_2to2
// type = 22 calculates from single_tbar_t-channel_onshell_2to3
// type = 23 calculates from single_tbar_t-channel_onshell_matched


std::vector<double> SingleTXsec(int type, double vl, double vr, double gl, double gr, double vo) // in pb
{
  double NNLO22top = 41.9;        // reference values (approximate NNLO) from arxiv:1103.2792v1 at 7 TeV and mtop = 172.5 GeV
  double NNLO22antitop = 22.7;
  double value[15], error[15];
  /* double constant, constantError; // this one forces the result to be consistent with the NNLO prediction.  */
  /*                                 // ONLY FOR 11 AND 21. THE LATTER ONE IS NOT YET IMPLEMENTED. */
  double factor, factorError;
  std::vector<double> result;
  double Xsec = 0, XsecError2 = 0;


  if(type == 11) 
    {
      /* factor      = 43.4558; // in pb */
      factor      = NNLO22top; // in pb
      factorError = 0.0024;

      /* constant = NNLO22top - factor; */

      value[0] = 1.;			
      error[0] = 0.;
      value[1] = 0.;
      error[1] = 0.;
      value[2] = 0.; 
      error[2] = 0.;
      value[3] = -0.5371918132907458;	
      error[3] = 0.00017495550233783498;
      value[4] = -5.6024006001500375;	
      error[4] = 0.0016412869341733355;
      value[5] = 0.8633047832510273;	
      error[5] = 0.00004917592196278635;
      value[6] = 0.04966425655493634;	
      error[6] = 0.00014645050182936022;
      value[7] = 0.; 
      error[7] = 0.;
      value[8] = 0.; 
      error[8] = 0.;
      value[9] = 1.1990735413914826;	
      error[9] = 0.00007195139209137117;
      value[10] = 0.; 
      error[10] = 0.;
      value[11] = 0.; 
      error[11] = 0.;
      value[12] = 1.8290170702184745;	
      error[12] = 0.00011010041778328444;
      value[13] = 2.8805959158501295;	
      error[13] = 0.001969966769041145;
      value[14] = 24.004850905978028;	
      error[14] = 0.001257648846757457;
    }
  else if(type == 12) 
    {
      factor      = 23.3786; // in pb
      factorError = 0.0050;

      /* constant = 0; */

      value[0]  = 1.;			
      error[0]  = 0.;	
      value[1]  = -0.10798336940620923;	
      error[1]  = 0.000489648525159936;
      value[2]  = -0.021233093512870393;	
      error[2]  = 0.0006342982967334109;
      value[3]  = -0.6704678637728523;
      error[3]  = 0.0008317792571288818;
      value[4]  = -5.824497617479231;	
      error[4]  = 0.007051475460122934;
      value[5]  = 0.8771953838125465;	
      error[5]  = 0.00019492027567152853;
      value[6]  = -0.030455202621200606;	
      error[6]  = 0.0006054265609166766;
      value[7]  = -0.03128929876040498;	
      error[7]  = 0.000910959422443907;
      value[8]  = 0.221159521955979;	
      error[8]  = 0.007896444171478792;
      value[9]  = 1.3494306759172918;	
      error[9]  = 0.00030560134144931936;
      value[10] = -0.005932776128596107;	
      error[10] = 0.0010421121714775574;
      value[11] = 0.09098919524693727;	
      error[11] = 0.00790741078734306;
      value[12] = 2.4677653922818306;	
      error[12] = 0.0005314807679703893;
      value[13] = 3.7649431531400523;	
      error[13] = 0.008609977670272937;
      value[14] = 26.605913100014543;	
      error[14] = 0.00545639165776817;
    } 
  else if(type == 13) 
    {
      factor      = 36.9825; // in pb
      factorError = 0.0051;

      /* constant = 0; */

      value[0] = 1.;			
      error[0] = 0.;
      value[1] = -0.06840803082539026;	
      error[1] = 0.0003159605571308092;
      value[2] = -0.0129818157236532;	
      error[2] = 0.0004056191979096772;
      value[3] = -0.6044534577164877;	
      error[3] = 0.000529462371522793;
      value[4] = -5.808544581896843;	
      error[4] = 0.004530872180183161;
      value[5] = 0.8768065977151356;	
      error[5] = 0.0001254800997290676;
      value[6] = -0.005067261542621715;	
      error[6] = 0.0003864770606398027;
      value[7] = -0.0184708983978914;	
      error[7] = 0.0005773793145661192;
      value[8] = 0.13594267558980633;	
      error[8] = 0.005022097191674445;
      value[9] = 1.3299912120597581;	
      error[9] = 0.00019635851155110746;
      value[10] = -0.0026769417967957843;	
      error[10] = 0.0006671727315287779;
      value[11] = 0.0607746907320994;	
      error[11] = 0.005062564650450958;
      value[12] = 2.246213749746502;	
      error[12] = 0.00033724348697110166;
      value[13] = 3.426036638950851;	
      error[13] = 0.005523547233300082;
      value[14] = 26.952369363888327;	
      error[14] = 0.003468985570064029;
    } 
  else if(type == 21) 
    {
      /* factor      = 22.5673; // in pb */
      factor      = NNLO22antitop; // in pb
      factorError = 0.0012;

      /* constant = NNLO22antitop - factor; */

      value[0] = 1.;			
      error[0] = 0.;
      value[1] = 0.; 
      error[1] = 0.;
      value[2] = 0.; 
      error[2] = 0.;
      value[3] = 0.001985173237383453;	
      error[3] = 0.00016599864602315903;
      value[4] = -4.702525335330325;	
      error[4] = 0.0008606795574999549;
      value[5] = 1.1394938694482726;	
      error[5] = 0.00006373656670438958;
      value[6] = -0.596712943063636;	
      error[6] = 0.0001963348439603702;
      value[7] = 0.; 
      error[7] = 0.;
      value[8] = 0.; 
      error[8] = 0.;
      value[9] = 2.008202133174992;	
      error[9] = 0.00011983271951638762;
      value[10] = 0.; 
      error[10] = 0.;
      value[11] = 0.; 
      error[11] = 0.;
      value[12] = 1.3662733246777417;	 
      error[12] = 0.00008238164325497241;
      value[13] = -0.0711250348956245;	
      error[13] = 0.000994532756274122;
      value[14] = 12.742153469843537;	
      error[14] = 0.0006709209633598152;
    } 
  else if(type == 22)
    {
      factor      = 11.9556; // in pb
      factorError = 0.0025;

      /* constant = 0; */

      value[0] = 1.;			
      error[0] = 0.;
      value[1] = -0.12142426979825327;	
      error[1] = 0.0005546738212415445;
      value[2] = -0.032570510890293836;	
      error[2] = 0.0009872888235332793;
      value[3] = -0.08991602261701637;	
      error[3] = 0.0006679785596718598;
      value[4] = -4.997122687276255;	
      error[4] = 0.004145582113661473;
      value[5] = 1.1231974974070726;	
      error[5] = 0.00024723558591206507;
      value[6] = -0.7256515775034291;	
      error[6] = 0.0008932967222997929;
      value[7] = -0.02243300210779875;	
      error[7] = 0.0007139016450470215;
      value[8] = 0.2403476195255756;	
      error[8] = 0.004834693431336954;
      value[9] = 2.617200307805547;	
      error[9] = 0.0005692097643990187;
      value[10] = -0.00482619023721087;	
      error[10] = 0.0011188900389560708;
      value[11] = 0.16074475559570445;	
      error[11] = 0.0050404683264037055;
      value[12] = 1.5389273645822876;	
      error[12] = 0.0003414366740912507;
      value[13] = 0.3760748101308167;	
      error[13] = 0.004934405164332686;
      value[14] = 16.025209943457458;	
      error[14] = 0.0032495544740053766;
    } 
  else if(type == 23)
    {
      factor      = 19.3051; // in pb
      factorError = 0.0026;

      /* constant = 0; */

      value[0] = 1.;			
      error[0] = 0.;
      value[1] = -0.07565876374636771;	
      error[1] = 0.00035043598041387814;
      value[2] = -0.021507270099611286;	
      error[2] = 0.0006127060717558157;
      value[3] = -0.057202500893546615;	
      error[3] = 0.0004206862656440399;
      value[4] = -4.9392181340682;	
      error[4] = 0.002581128712214851;
      value[5] = 1.1242417806693568;	
      error[5] = 0.00015533899700344924;
      value[6] = -0.6584633076233741;	
      error[6] = 0.000556933990054195;
      value[7] = -0.014457319568404126;	
      error[7] = 0.00045019762094692273;
      value[8] = 0.1503436915633678;	
      error[8] = 0.0030394656314059773;
      value[9] = 2.3905496475024735;	
      error[9] = 0.0003538810465629399;
      value[10] = -0.004775940036570558;	
      error[10] = 0.0006996925334646549;
      value[11] = 0.09563276025506084;	
      error[11] = 0.003164475124760582;
      value[12] = 1.5002201490797769;	
      error[12] = 0.00021658811156858207;
      value[13] = 0.20839570890593784;
      error[13] = 0.003080699758858034;
      value[14] = 15.329783321505717;	
      error[14] = 0.0020364080899829555;
    }
  else
    {
      factor      = 0;
      factorError = 0;

      /* constant = 0; */

      value[0] = 0.;			
      error[0] = 0.;
      value[1] = 0.;	
      error[1] = 0.;
      value[2] = 0.;	
      error[2] = 0.;
      value[3] = 0.;	
      error[3] = 0.;
      value[4] = 0.;	
      error[4] = 0.;
      value[5] = 0.;	
      error[5] = 0.;
      value[6] = 0.;
      error[6] = 0.;
      value[7] = 0.;	
      error[7] = 0.;
      value[8] = 0.;	
      error[8] = 0.;
      value[9] = 0.;	
      error[9] = 0.;
      value[10] = 0.;	
      error[10] = 0.;
      value[11] = 0.;	
      error[11] = 0.;
      value[12] = 0.;	
      error[12] = 0.;
      value[13] = 0.;
      error[13] = 0.;
      value[14] = 0.;	
      error[14] = 0.;
    }


  // changed the normalization in that way that the plot is consistent with the NNLO values
  //calculates the value of the cross section in fb
  Xsec = factor*(value[0]  *vl*vl
		 +value[1] *vl*vr
		 +value[2] *vl*gl
		 +value[3] *vl*gr
		 +value[4] *vl*vo
		 +value[5] *vr*vr
		 +value[6] *vr*gl
		 +value[7] *vr*gr
		 +value[8] *vr*vo
		 +value[9] *gl*gl
		 +value[10]*gl*gr
		 +value[11]*gl*vo
		 +value[12]*gr*gr
		 +value[13]*gr*vo
		 +value[14]*vo*vo)
    /* +constant; */;


  result.push_back(Xsec);

  //calculates the error of the cross section in fb
  XsecError2 = 
    pow(Xsec*factorError/factor,2.0)
    +pow(factor,2.0)
   *(pow(error[0] *vl*vl,2.0)
    +pow(error[1] *vl*vr,2.0)
    +pow(error[2] *vl*gl,2.0)
    +pow(error[3] *vl*gr,2.0)
    +pow(error[4] *vl*vo,2.0)
    +pow(error[5] *vr*vr,2.0)
    +pow(error[6] *vr*gl,2.0)
    +pow(error[7] *vr*gr,2.0)
    +pow(error[8] *vr*vo,2.0)
    +pow(error[9] *gl*gl,2.0)
    +pow(error[10]*gl*gr,2.0)
    +pow(error[11]*gl*vo,2.0)
    +pow(error[12]*gr*gr,2.0)
    +pow(error[13]*gr*vo,2.0)
    +pow(error[14]*vo*vo,2.0)
	     );

  result.push_back(sqrt(XsecError2));

  return result;
}


// -----------------------------------------------------------------------------------------------
// helicity fractions F0, FL, FR
double HelicityFracs(int type, double vl, double vr, double gl, double gr, double fMtop, double fMb, double fMW)  // in %
                                                                          // type = -1: FL
                                                                          // type = 0:  F0
                                                                          // type = 1:  FR                            
{
  double NNLOF0 = 68.7; // theory values from arXiv:1005:2625
  double NNLOFL = 31.1;
  // double fMtop = 172.5; 
  // double fMb = 4.8;    
  // double fMW = 80.4;

  double xw = fMW/fMtop;
  double xb = fMb/fMtop;

  double q = sqrt(fMtop*fMtop*fMtop*fMtop + fMW*fMW*fMW*fMW +
                  fMb*fMb*fMb*fMb - 2.*(fMtop*fMW)*(fMtop*fMW) -
                  2.*(fMtop*fMb)*(fMtop*fMb) - 2.*(fMb*fMW)*(fMb*fMW))/(2.*fMtop);

  double Gamma_0 = q*((vl*vl+vr*vr)*(1.-xw*xw-2.*xb*xb-
                                     (xw*xb)*(xw*xb)+xb*xb)/(xw*xw)+
                      (gl*gl+gr*gr)*(1.-xw*xw+xb*xb)-
                      4.*xb*(vl*vr+gl*gr)-
                      2.*vl*(gr-xb*gl)*(1.-xb*xb)/xw-
                      2.*vr*(gl-xb*gr)*(1.-xb*xb)/xw+
                      2.*xw*vl*(gr+xb*gl)+2.*xw*vr*(gl+xb*gr));

  double g1 = q*((vl*vl+vr*vr)*(1.-xw*xw+xb*xb)-
                 4.*xb*(vl*vr+gl*gr)+
                 (gl*gl+gr*gr)*(1.-xw*xw-2.*xb*xb-
                                (xw*xb)*(xw*xb)+xb*xb)/(xw*xw)-
                 2.*vl*(gr-xb*gl)*(1.-xb*xb)/xw-
                 2.*vr*(gl-xb*gr)*(1.-xb*xb)/xw+
                 2.*xw*vl*(gr+xb*gl)+2.*xw*vr*(gl+xb*gr));

  double g2 = (fMtop*fMtop*fMtop/(2.*fMW*fMW))*(
                                                (gl*gl-gr*gr)*(1.-xb*xb)+2.*xw*vl*(gr+xb*gl)-
                                                (xw*xw)*(vl*vl-vr*vr)-
                                                2.*xw*vr*(gl+xb*gr))*
    (1.-2.*xw*xw-2.*xb*xb-2.*(xb*xw)*(xb*xw)+
     xw*xw*xw*xw+xb*xb*xb*xb);

  double Gamma_p = g1 + g2;
  double Gamma_m = g1 - g2;

  // values at VL=1, VR=gL=gR=0;
  double Gamma_0_0 = q*((1.-xw*xw-2.*xb*xb-(xw*xb)*(xw*xb)+xb*xb)/(xw*xw));
  
  double g1_0 = q*((1.*1.+0.*0.)*(1.-xw*xw+xb*xb)-
                 4.*xb*(1.*0.+0.*0.)+
                 (0.*0.+0.*0.)*(1.-xw*xw-2.*xb*xb-
                                (xw*xb)*(xw*xb)+xb*xb)/(xw*xw)-
                 2.*1.*(0.-xb*0.)*(1.-xb*xb)/xw-
                 2.*0.*(0.-xb*0.)*(1.-xb*xb)/xw+
                 2.*xw*1.*(0.+xb*0.)+2.*xw*0.*(0.+xb*0.));
  
  double g2_0 = (fMtop*fMtop*fMtop/(2.*fMW*fMW))*(
                                                (0.*0.-0.*0.)*(1.-xb*xb)+2.*xw*1.*(0.+xb*0.)-
                                                (xw*xw)*(1.*1.-0.*0.)-
                                                2.*xw*0.*(0.+xb*0.))*
    (1.-2.*xw*xw-2.*xb*xb-2.*(xb*xw)*(xb*xw)+
     xw*xw*xw*xw+xb*xb*xb*xb);

  double Gamma_p_0 = g1_0 + g2_0;
  double Gamma_m_0 = g1_0 - g2_0;

  double aux = Gamma_0 + Gamma_m + Gamma_p;
  double aux_0 = Gamma_0_0 + Gamma_m_0 + Gamma_p_0;
      


  // changed the normalization in that way that the plot is consistent with the approximate NNLO values
  if(type ==  0)
    return Gamma_0*aux_0/(Gamma_0_0*aux)*NNLOF0;  // F0 in %
  else if(type == -1)
    return Gamma_m*aux_0/(Gamma_m_0*aux)*NNLOFL;  // FL in %
  else if(type ==  1) 
    return 100*Gamma_p/aux;  // FR in %
  else 
    {
      std::cout << "Error! \"type = " << type << "\" is not valid because it doesn't describe F0, FL or FR. Return 0." << std::endl; 
      return 0;
    }
}
