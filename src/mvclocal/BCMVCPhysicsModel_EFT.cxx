/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMVCPhysicsModel_EFT.h"
#include "BCMVCObservable.h"

#include "BAT/BCModel.h"
#include "BAT/BCMath.h"
#include "BAT/BCLog.h"
#include "BAT/BCH1D.h"
#include "BAT/BCH2D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <iomanip>


// ---------------------------------------------------------
BCMVCPhysicsModel_EFT::BCMVCPhysicsModel_EFT() : BCMVCombination_EFT()
{
}

// ---------------------------------------------------------
BCMVCPhysicsModel_EFT::~BCMVCPhysicsModel_EFT()
{
}

// ---------------------------------------------------------
double BCMVCPhysicsModel_EFT::LogLikelihood(const std::vector<double> &parameters)
{
  std::vector<double> observables;

  for (int i = 0; i < fNObservables; ++i)
    observables.push_back(CalculateObservable(i, parameters));

  // add nuisance parameters
  for(int j=0; j<fNNuisanceCorrelation; j++)
    observables.push_back(parameters.at(fNObservables+j));

  return BCMVCombination_EFT::LogLikelihood(observables);
}

// ---------------------------------------------------------
void BCMVCPhysicsModel_EFT::AddObservable(std::string name, double min, double max)
{
  // check if observable exists already
  int index = GetIndexObservable(name);
 
  if (index >= 0)
    return;

  BCMVCObservable* observable = new BCMVCObservable();
  observable->SetName(name);
  observable->SetMinMax(min, max);
  fObservables.push_back(observable);

  fNObservables++;
}

// ---------------------------------------------------------

void BCMVCPhysicsModel_EFT::ReadParameters(std::string file)
{
  //open file:
  tinyxml2::XMLDocument doc;
  doc.LoadFile(file.c_str());

  // set number of parameters
  int NParams=0;
  if(doc.FirstChildElement("options")->QueryIntAttribute("NParams",&NParams)){
    throw std::runtime_error("Please check attribute \"NParams\" in file \""+file+"\"\n");
  }

  int NBins=0;
  if(doc.FirstChildElement("options")->QueryIntAttribute("NBins",&NBins)){
    throw std::runtime_error("Please check attribute \"NBins\" in file \""+file+"\"\n");
  }

  // check if NParams is smaller than the given parameters
  std::stringstream par;
  par<<"parameter"<<NParams;
  if(doc.FirstChildElement(par.str().c_str())!=NULL){
    std::cerr<<"WARNING: NParams="<<NParams<<" is smaller than the given set of parameters in file \""<<file.c_str()<<"\"."<<std::endl;
    std::cerr<<"Press ENTER to continue with only "<<NParams<<" parameters."<<std::endl;
    std::cin.get();
    std::cout<<"continuing..."<<std::endl;
  }


  // read parameter names, min and max values from settings and add new parameter
  for(int i=0; i<NParams; i++){ 
    std::stringstream param;
    param<<"parameter"<<i;

    if(doc.FirstChildElement(param.str().c_str())==NULL){
      throw std::runtime_error("Could not find \""+param.str()+"\"! Please check attribute \"NParams\" in file \""+file+"\"\n");}
    
    std::string name=doc.FirstChildElement(param.str().c_str())->Attribute("name");
   
    double min=0;
    if(doc.FirstChildElement(param.str().c_str())->QueryDoubleAttribute("min",&min)){
    throw std::runtime_error("Please check attribute \"min\" in \""+param.str()+"\" in file \""+file+"\"\n");}

    double max=0;
    if(doc.FirstChildElement(param.str().c_str())->QueryDoubleAttribute("max",&max)){
    throw std::runtime_error("Please check attribute \"max\" in \""+param.str()+"\" in file \""+file+"\"\n");}

    AddParameter(name.c_str(), min,max);
    // set number of bins 
    GetParameter(name.c_str())->SetNbins(NBins);

    // get SM value
    double sm_value=0;
    if(doc.FirstChildElement(param.str().c_str())->QueryDoubleAttribute("SM_value",&sm_value)){
    throw std::runtime_error("Please check attribute \"SM_value\" in \""+param.str()+"\" in file \""+file+"\"\n");}
    Parameters_SM.push_back(sm_value);


    // set prior according to runcard:
    std::string prior=doc.FirstChildElement(param.str().c_str())->Attribute("prior");

    //const prior
    if(prior=="const"){
      SetPriorConstant(name.c_str());}

    //delta prior
    else if (prior=="delta"){
      // throw error if delta_value attribute does not exist
    double delta_value=0;
    if(doc.FirstChildElement(param.str().c_str())->QueryDoubleAttribute("delta_value",&delta_value)){
    throw std::runtime_error("Please check attribute \"delta_value\" in \""+param.str()+"\" in file \""+file+"\"\n");}
    SetPriorDelta(name.c_str(), delta_value);
    }

    // gaussian prior
   else if (prior=="gauss"){
      // throw error if attribute for value does not exist
    double mean=0;
    if(doc.FirstChildElement(param.str().c_str())->QueryDoubleAttribute("gauss_mean",&mean)){
    throw std::runtime_error("Please check attribute \"gauss_mean\" in \""+param.str()+"\" in file \""+file+"\"\n");}

    double sigma=0;
    if(doc.FirstChildElement(param.str().c_str())->QueryDoubleAttribute("gauss_sigma",&sigma)){
    throw std::runtime_error("Please check attribute \"gauss_sigma\" in \""+param.str()+"\" in file \""+file+"\"\n");}

    SetPriorGauss(name.c_str(),mean,sigma);
    }

  else{
    throw std::runtime_error("Please check attribute \"prior\" in \""+param.str()+"\" in file \""+file+"\"\n"); }

  }//end of parameter loop
}

// ---------------------------------------------------------
void BCMVCPhysicsModel_EFT::PlotObservables(std::string filename){

  TFile *file = new TFile((filename + ".root").c_str(), "RECREATE"); 
  std::vector<double> Parameters_new (GetNParameters()-fNNuisanceCorrelation);
  TGraph *gr1; // TGraph for plots
  TGraph *gr2; //Tgraph for SM values
  unsigned int NParams = GetNParameters()-fNNuisanceCorrelation; // number of parameters without nuisance parameters

  // for all Observables:
  for(unsigned int j=0; j<fObservables.size(); j++){

    // for all (not nuisance) parameters p:
    for(unsigned int p=0; p<NParams; p++){
      std::stringstream jp;
      jp<<j<<p;

       // set all Parameters_new on SM value:
      for(unsigned int i=0; i<NParams; i++){
          Parameters_new.at(i)=Parameters_SM.at(i);
      }

      // get Parameter limits
      double min=GetParameter(p)->GetLowerLimit();
      double max=GetParameter(p)->GetUpperLimit();

      // Number of steps for plotting
      int NSteps = 1e4;
      double x[NSteps], y[NSteps];

      // calculate observable for each step
      for(int i=0; i<NSteps; i++){
        x[i]=min+i*((max-min)/NSteps);
        Parameters_new[p]=x[i];
        y[i]=CalculateObservable(j, Parameters_new);
      }

      //fill Graph with line
      gr1 = new TGraph(NSteps,x,y); 

      // Graph for SM value
      gr2 = new TGraph(1); 
      gr2->SetPoint(0,Parameters_SM.at(p),CalculateObservable(j, Parameters_SM)); // add SM value to Graph 2

      // combine line and SM-point in one TMultiGraph:
      TMultiGraph *mg = new TMultiGraph();
      mg->Add(gr1,"C");
      mg->Add(gr2,"*");

      // Set Axis titles and legend
      TCanvas *c1 = new TCanvas(jp.str().c_str(),"Graph",200,10,600,400);
      mg->Draw("a");
      mg->GetXaxis()->SetTitle(GetParameter(p)->GetName().c_str());
      mg->GetYaxis()->SetTitle(fObservables.at(j)->GetName().c_str());

      gPad->Modified();
      file->cd();
      mg->Write();

      if(p==0 && j==0){c1->Print((filename+".pdf(").c_str());}
      else if(p==NParams-1 && j==fObservables.size()-1){c1->Print((filename+".pdf)").c_str());}
      else{c1->Print((filename+".pdf").c_str());}

      delete gr1;
      delete gr2;
      delete c1;
      delete mg;
    }
  }  
}

// ---------------------------------------------------------
void BCMVCPhysicsModel_EFT::ReadMeasurements(std::string filename){

  if(GetNParameters()==0){ // check that ReadParameters was called first
     BCLog::OutWarning("BCMVCPhysicsModel_EFT::ReadMeasurements: No parameters defined!");
     throw std::runtime_error("BCMVCPhysicsModel_EFT::ReadMeasurements: No parameters defined!");
  }

   tinyxml2::XMLDocument doc;
   doc.LoadFile(filename.c_str());

   // get number of observables, measurements and uncertainties
   int nobservables=0;
   if(doc.FirstChildElement("options")->QueryIntAttribute("NObservables",&nobservables)){
    throw std::runtime_error("Please check attribute \"NObservables\" in file \""+filename+"\"\n"); }
 
   int nmeasurements=0;
   if(doc.FirstChildElement("options")->QueryIntAttribute("NMeasurements",&nmeasurements)){
    throw std::runtime_error("Please check attribute \"NMeasurements\" in file \""+filename+"\"\n"); }

   int nuncertainties=0;
   if(doc.FirstChildElement("options")->QueryIntAttribute("NUncertainties",&nuncertainties)){
    throw std::runtime_error("Please check attribute \"NUncertainties\" in file \""+filename+"\"\n"); }

    // check if NObservables is smaller than the given set of observables
    std::stringstream o;
    o<<"observable"<<nobservables;
    if(doc.FirstChildElement(o.str().c_str())!=NULL){
      std::cerr<<"WARNING: NObservables="<<nobservables<<" is smaller than the given set of observables in file \""<<filename.c_str()<<"\"."<<std::endl;
      std::cerr<<"Press ENTER to continue with only "<<nobservables<<" observables."<<std::endl;
      std::cin.get();
      std::cout<<"continuing..."<<std::endl;
    }

   // add Observables
   for (int i = 0; i < nobservables; ++i) {
      std::stringstream obs;
      obs<<"observable"<<i;

      // check if obs exists
      if(doc.FirstChildElement(obs.str().c_str())==NULL){
      throw std::runtime_error("Could not find \""+obs.str()+"\"! Please check attribute \"NObservables\" in file \""+filename+"\"\n");}

      if(doc.FirstChildElement(obs.str().c_str())->Attribute("name")==NULL){
        throw std::runtime_error("Please check attribute \"name\" in \""+obs.str()+"\" in file \""+filename+"\"\n");}
      std::string name=doc.FirstChildElement(obs.str().c_str())->Attribute("name");

      double min=0;
      if(doc.FirstChildElement(obs.str().c_str())->QueryDoubleAttribute("min",&min)){
        throw std::runtime_error("Please check attribute \"min\" in \""+obs.str()+"\" in file \""+filename+"\"\n");}
     
     double max=0;
      if(doc.FirstChildElement(obs.str().c_str())->QueryDoubleAttribute("max",&max)){
        throw std::runtime_error("Please check attribute \"max\" in \""+obs.str()+"\" in file \""+filename+"\"\n");}

      // add observable
      AddObservable(name.c_str(), min, max);
   }


   // add uncertainties

   // check if NUncertainties is smaller than the given set of uncertainties
    std::stringstream un;
    un<<"uncertainty"<<nuncertainties;
    if(doc.FirstChildElement("uncertainties")->Attribute(un.str().c_str())!=NULL || doc.FirstChildElement("measurement0")->Attribute(un.str().c_str())!=NULL){
      std::cerr<<"WARNING: NUncertainties="<<nuncertainties<<" is smaller than the given set of uncertainties in file \""<<filename.c_str()<<"\"."<<std::endl;
      std::cerr<<"Press ENTER to continue with only "<<nuncertainties<<" uncertainties."<<std::endl;
      std::cin.get();
      std::cout<<"continuing..."<<std::endl;
    }

   for (int i = 0; i < nuncertainties; ++i) {
      std::stringstream unc;
      unc<<"uncertainty"<<i;
      //check if unc exists
      if(doc.FirstChildElement("uncertainties")->Attribute(unc.str().c_str())==NULL){
      throw std::runtime_error("Could not find \""+unc.str()+"\"! Please check attribute \"NUncertainties\" in file \""+filename+"\"\n");}
      
      std::string name=doc.FirstChildElement("uncertainties")->Attribute(unc.str().c_str());
      AddUncertainty(name);
   }

   // add measurements

   // check if NMeasurements is smaller than the given set of measurements
    std::stringstream measurm;
    measurm<<"measurement"<<nmeasurements;
    if(doc.FirstChildElement(measurm.str().c_str())!=NULL){
      std::cerr<<"WARNING: NMeasurements="<<nmeasurements<<" is smaller than the given set of measurements in file \""<<filename.c_str()<<"\"."<<std::endl;
      std::cerr<<"Press ENTER to continue with only "<<nmeasurements<<" measurements."<<std::endl;
      std::cin.get();
      std::cout<<"continuing..."<<std::endl;
    }

   for (int i = 0; i < nmeasurements; ++i) {
      std::stringstream meas;
      meas<<"measurement"<<i;

      // check if meas exists
      if(doc.FirstChildElement(meas.str().c_str())==NULL){
      throw std::runtime_error("Could not find \""+meas.str()+"\"! Please check attribute \"NMeasurements\" in file \""+filename+"\"\n");}

      if(doc.FirstChildElement(meas.str().c_str())->Attribute("name")==NULL){
        throw std::runtime_error("Please check attribute \"name\" in \""+meas.str()+"\" in file \""+filename+"\"\n");}
      std::string name=doc.FirstChildElement(meas.str().c_str())->Attribute("name");

      if(doc.FirstChildElement(meas.str().c_str())->Attribute("observable")==NULL){
        throw std::runtime_error("Please check attribute \"observable\" in \""+meas.str()+"\" in file \""+filename+"\"\n");}
      std::string observable=doc.FirstChildElement(meas.str().c_str())->Attribute("observable");

      double central = 0;
      if(doc.FirstChildElement(meas.str().c_str())->QueryDoubleAttribute("value",&central)){
        throw std::runtime_error("Please check attribute \"value\" in \""+meas.str()+"\" in file \""+filename+"\"\n");}

      std::vector<double> uncertainties(0);
      for (int j = 0; j < nuncertainties; ++j) {
        std::stringstream unc;
        unc<<"uncertainty"<<j;
        double uncertainty=0;
        if(doc.FirstChildElement(meas.str().c_str())->QueryDoubleAttribute(unc.str().c_str(),&uncertainty)){
          throw std::runtime_error("Please check attribute \""+unc.str()+"\" in \""+meas.str()+"\" in file \""+filename+"\"\n");}

        uncertainties.push_back(uncertainty);
      }

      // add measurement
      AddMeasurement(name, observable, central, uncertainties);

      //check if measurements active
      bool active = 0;
      if(doc.FirstChildElement(meas.str().c_str())->QueryBoolAttribute("active",&active)){
        throw std::runtime_error("Please check attribute \"active\" in \""+meas.str()+"\" in file \""+filename+"\"\n");}
      GetMeasurement(i)->SetFlagActive(active);
   }

   // read covariance matrices
   for(int i=0; i<nuncertainties; i++){
     std::stringstream matrix;
     matrix<<"CorrelationMatrix"<<i;
     if(doc.FirstChildElement(matrix.str().c_str())==NULL){
      throw std::runtime_error("Please check attribute \""+matrix.str()+"\" in file \""+filename+"\"\n");}
     std::string correlation = doc.FirstChildElement(matrix.str().c_str())->GetText();
     std::stringstream corr;
     corr<<correlation;
   
      TMatrixD mat(nmeasurements, nmeasurements);

      for (int j = 0; j < nmeasurements; ++j){
        for (int k = 0; k < nmeasurements; ++k) {
          double c;
          corr >> c;
          mat[j][k] = c;
        }
      }
        
      GetUncertainty(i)->SetCorrelationMatrix(mat);     
   }

  AddNuisanceParameter(filename); // call function for adding nuisance parameters
  
 PrepareAnalysis();
}

// ---------------------------------------------------------
void BCMVCPhysicsModel_EFT::AddNuisanceParameter(std::string filename){
  tinyxml2::XMLDocument doc;
  doc.LoadFile(filename.c_str());

  int nnuisance=0;
  if(doc.FirstChildElement("options")->QueryIntAttribute("NNuisance",&nnuisance)){
    throw std::runtime_error("Please check attribute \"NNuisance\" in file \""+filename+"\"\n"); }


  // check if NNuisance is smaller than the given parameters
  std::stringstream nu;
  nu<<"nuisance_parameter"<<nnuisance;
  if(doc.FirstChildElement(nu.str().c_str())!=NULL){
    std::cerr<<"WARNING: NNuisance="<<nnuisance<<" is smaller than the given set of nuisance parameters in file \""<<filename.c_str()<<"\"."<<std::endl;
    std::cerr<<"Press ENTER to continue with only "<<nnuisance<<" nuisance parameters."<<std::endl;
    std::cin.get();
    std::cout<<"continuing..."<<std::endl;
  }


   for (int i = 0; i < nnuisance; ++i) {

      std::stringstream nuis;
      nuis<<"nuisance_parameter"<<i;

      // check ig nuis exists
      if(doc.FirstChildElement(nuis.str().c_str())==NULL){
      throw std::runtime_error("Could not find \""+nuis.str()+"\"! Please check attribute \"NNuisance\" in file \""+filename+"\"\n");}


      if(doc.FirstChildElement(nuis.str().c_str())->Attribute("uncertainty")==NULL){
        throw std::runtime_error("Please check attribute \"uncertainty\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
      std::string uncertainty=doc.FirstChildElement(nuis.str().c_str())->Attribute("uncertainty"); 

      if(doc.FirstChildElement(nuis.str().c_str())->Attribute("measurement1")==NULL){
        throw std::runtime_error("Please check attribute \"measurement1\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
      std::string measurement1=doc.FirstChildElement(nuis.str().c_str())->Attribute("measurement1");  

      if(doc.FirstChildElement(nuis.str().c_str())->Attribute("observable1")==NULL){
        throw std::runtime_error("Please check attribute \"observable1\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
      std::string observable1=doc.FirstChildElement(nuis.str().c_str())->Attribute("observable1");  

      if(doc.FirstChildElement(nuis.str().c_str())->Attribute("measurement2")==NULL){
        throw std::runtime_error("Please check attribute \"measurement2\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
      std::string measurement2=doc.FirstChildElement(nuis.str().c_str())->Attribute("measurement2");  

      if(doc.FirstChildElement(nuis.str().c_str())->Attribute("observable2")==NULL){
        throw std::runtime_error("Please check attribute \"observable2\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
      std::string observable2=doc.FirstChildElement(nuis.str().c_str())->Attribute("observable2");  

      if(doc.FirstChildElement(nuis.str().c_str())->Attribute("parameter_name")==NULL){
        throw std::runtime_error("Please check attribute \"parameter_name\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
      std::string parname=doc.FirstChildElement(nuis.str().c_str())->Attribute("parameter_name");  

       if(doc.FirstChildElement(nuis.str().c_str())->Attribute("prior")==NULL){
        throw std::runtime_error("Please check attribute \"prior\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
      std::string priorshape=doc.FirstChildElement(nuis.str().c_str())->Attribute("prior");  


    // check if parameter exists already
    for (unsigned int j = 0; j < GetNParameters(); j++){
      if (parname.c_str() == GetParameter(j)->GetName())
      throw std::runtime_error("Nuisance parameter \""+parname+"\" already exists as a parameter.");
    }


     // read properties of parameter
    double min=0;
    if(doc.FirstChildElement(nuis.str().c_str())->QueryDoubleAttribute("min",&min)){
    throw std::runtime_error("Please check attribute \"min\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}

    double max=0;
    if(doc.FirstChildElement(nuis.str().c_str())->QueryDoubleAttribute("max",&max)){
    throw std::runtime_error("Please check attribute \"max\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}

     // add nuisance parameter
     AddParameter(parname.c_str(), min, max);

    // set pre-factor
    double pre = 1;

    // set index
   int index = GetNObservables()-1;

     if (priorshape == "const") {
      SetPriorConstant(parname.c_str());
     }

    else if (priorshape=="delta"){
      // throw error if delta_value attribute does not exist
    double delta_value=0;
    if(doc.FirstChildElement(nuis.str().c_str())->QueryDoubleAttribute("delta_value",&delta_value)){
    throw std::runtime_error("Please check attribute \"delta_value\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}
    SetPriorDelta(parname.c_str(), delta_value);
    }

    else if (priorshape == "gauss") {
     double mean=0;
     if(doc.FirstChildElement(nuis.str().c_str())->QueryDoubleAttribute("gauss_mean",&mean)){
      throw std::runtime_error("Please check attribute \"gauss_mean\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}

     double sigma=0;
     if(doc.FirstChildElement(nuis.str().c_str())->QueryDoubleAttribute("gauss_sigma",&sigma)){
      throw std::runtime_error("Please check attribute \"gauss_sigma\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");}

     SetPriorGauss(parname.c_str(),mean,sigma);
    }

    else {
     throw std::runtime_error("Please check attribute \"prior\" in \""+nuis.str()+"\" in file \""+filename+"\"\n");
    }
    

    // increase counter of nuisance parameters
    fNNuisanceCorrelation++;

    NuisanceParameter p;
    p.index_uncertainty  = GetIndexUncertainty(uncertainty);
    p.index_measurement1 = GetIndexMeasurement(measurement1, observable1);
    p.index_measurement2 = GetIndexMeasurement(measurement2, observable2);
    p.index_rhoparameter = index;
    p.pre = pre;
    fNuisanceCorrelation.push_back(p);
  }
}

// ---------------------------------------------------------
void BCMVCPhysicsModel_EFT::GetRanking(std::string flag, double p, std::string file){
  std::ofstream out(file.c_str());

  int nmeas=GetNMeasurements();
  int nmeas_active=GetNActiveMeasurements();
  int NUncertainties=GetNUncertainties();
  unsigned int NParams = GetNParameters()-fNNuisanceCorrelation; // number of parameters without nuisance parameters

  // return if only one measuements/uncertainty is active
  if(flag=="measurement"){
    if(nmeas_active <= 1){
      std::cout<<"WARNING : Can't do measurement ranking. Only 1 measurement active!"<<std::endl;
      return;
    }
  }

  else if(flag=="uncertainty"){
    if(NUncertainties <= 1){
      std::cout<<"WARNING : Can't do uncertainty ranking. Only 1 uncertainty active!"<<std::endl;
      return;
    }
  }

  std::string outstring="Calculate "+flag+" ranking";
  BCLog::Out(outstring.c_str());

  // get original activitiy status of measurements
  std::vector<bool> orig_status;
  if(flag=="measurement"){
    for(int i=0; i<nmeas; i++){
      orig_status.push_back(GetMeasurement(i)->GetFlagActive());
    }
  }

  // get original activitiy status of uncertainties
  if(flag=="uncertainty"){
    for(int i=0; i<NUncertainties; i++){
      orig_status.push_back(GetUncertainty(i)->GetFlagActive());
    }
  }

  // get original areas 1D
  PrepareAnalysis();
  MarginalizeAll();    

  std::vector<double> orig_areas1D(NParams);
  int NFreeParams=0;
  for(unsigned int i=0; i<NParams; i++){
    if(GetParameter(i)->Fixed()!=0){continue;}
    orig_areas1D[i]=GetArea1D(i,p);
    NFreeParams++;
  }

  int NCombinations=0.5*NFreeParams*(NFreeParams-1)+0.5; // number of different combinations of two parameters (+0.5 for correct casting to int)

 // get original areas 2D
 std::vector<double> orig_areas2D(NCombinations);
 int Counter_2D=0;
  for(int index1=0; index1<=NParams/2.+0.5; index1++){
    for(unsigned int index2=index1+1; index2<NParams; index2++){
      if(GetParameter(index1)->Fixed()==0 && GetParameter(index2)->Fixed()==0){
        orig_areas2D[Counter_2D]=GetArea2D(index1,index2,p);
        Counter_2D++;
      }
    }
  }


  // vectors for saving areas
  std::vector< std::vector<std::pair<int,double> > > areas1D(NParams, std::vector<std::pair<int,double> >(nmeas));
  std::vector< std::vector<std::pair<int,double> > > areas2D(NCombinations, std::vector<std::pair<int,double> >(nmeas));

  // loop over all active measurements/uncertainties
  int N=0;
  if(flag=="measurement"){N=nmeas;}
  else if(flag=="uncertainty"){N=NUncertainties;}

  for(int i=0; i<N; i++){

    // if measurement/uncertainty already deactivated, continue
    std::vector<bool> active(orig_status);
    if(active.at(i)==0){ 
      continue;
    }

    // deactivate one measurement/uncertainty
    active.at(i)=0; // active=false
    if(flag=="measurement"){
      for(int j=0; j<N;j++){
        GetMeasurement(j)->SetFlagActive(active.at(j));
      }
    }

    else if(flag=="uncertainty"){
      for(int j=0; j<N;j++){
        GetUncertainty(j)->SetFlagActive(active.at(j));
      }
    }

    PrepareAnalysis();
    MarginalizeAll();    

    // 1D ranking
    for(unsigned int j=0; j<NParams; j++){
      if(GetParameter(j)->Fixed()==0){ // if parameter.fixed==false
        // get area
        areas1D[j][i]=std::make_pair(i,GetArea1D(j,p));
      }
    }

    // 2D ranking
    Counter_2D=0;
    for(int index1=0; index1<=NParams/2.+0.5; index1++){
      for(unsigned int index2=index1+1; index2<NParams; index2++){
        if(GetParameter(index1)->Fixed()==0 && GetParameter(index2)->Fixed()==0){
        areas2D[Counter_2D][i]=std::make_pair(i,GetArea2D(index1,index2,p));
        Counter_2D++;
        }
      }
    }

    //check for error
    if(Counter_2D!=NCombinations){
      throw std::runtime_error("Error: Range of Counter_2D in GetRanking wrong");
    }
 
  } // end of measurement/uncertainty loop

  // set measurements/uncertainties status back to original
  if(flag=="measurement"){
    for(int j=0; j<nmeas;j++){
      GetMeasurement(j)->SetFlagActive(orig_status.at(j));
    }

    //print Ranking
    PrintRanking1D("measurement", out, p,orig_areas1D, areas1D);
    PrintRanking2D("measurement",out, p,orig_areas2D, areas2D);
  }

  else if(flag=="uncertainty"){
    for(int j=0; j<NUncertainties;j++){
      GetUncertainty(j)->SetFlagActive(orig_status.at(j));
    }

    //print Ranking
    PrintRanking1D("uncertainty", out, p,orig_areas1D, areas1D);
    PrintRanking2D("uncertainty",out, p,orig_areas2D, areas2D);
  }

  // get original results
  PrepareAnalysis();
  MarginalizeAll(); 
}


// ---------------------------------------------------------
void BCMVCPhysicsModel_EFT::PrintRanking2D(std::string flag, std::ofstream& out,double p, std::vector<double> orig_areas2D,std::vector< std::vector<std::pair<int,double> > > areas2D){
 unsigned int NParams = GetNParameters()-fNNuisanceCorrelation; // number of parameters without nuisance parameters

  out<<"\n\n\nResults of the two-dimensional "<<flag<<" ranking based on the "<<p*100<<" % interval: \n"<<std::endl;
    int Counter_2D=0;
    for(int index1=0; index1<=NParams/2.+0.5; index1++){
      for(unsigned int index2=index1+1; index2<NParams; index2++){
        if(GetParameter(index1)->Fixed()==1 || GetParameter(index2)->Fixed()==1){continue;}

      out<<"Ranking based on parameters \""<<GetParameter(index1)->GetName()<<"\" and \""<< GetParameter(index2)->GetName()<<"\" (parameter "<<index1<<" and parameter "<<index2<<"):"<<std::endl;
      if(flag=="measurement"){out<< std::left << std::setw(5)<<"Rank"<<"        "<<"Measurement"<<"                 "<<"Increase (%)"<<std::endl;
        out<<"------------------------------------------------------"<<std::endl;}
      else if(flag=="uncertainty"){out<< std::left << std::setw(5)<<"Rank"<<"        "<<"Uncertainty"<<"          "<<"Decrease (%)"<<std::endl;
        out<<"------------------------------------------------"<<std::endl;}
      
      // calculate ranking
      std::vector<std::pair<int,double> > rank; // first: measurement_index, second: rank

      // for measurements
      if(flag=="measurement"){
        for(int i=0; i<GetNMeasurements();i++){
          if(GetMeasurement(i)->GetFlagActive()==0){continue;}
          double x = (areas2D[Counter_2D][i].second-orig_areas2D.at(Counter_2D))/orig_areas2D.at(Counter_2D);
          rank.push_back(std::make_pair(i,x));
        }
      }

      // for uncertainties
      if(flag=="uncertainty"){
        for(int i=0; i<GetNUncertainties();i++){
          if(GetUncertainty(i)->GetFlagActive()==0){continue;}
          double x = (areas2D[Counter_2D][i].second-orig_areas2D.at(Counter_2D))/orig_areas2D.at(Counter_2D);
          rank.push_back(std::make_pair(i,x));
        }
      }

  
     if(flag=="measurement"){std::sort(rank.begin(), rank.end(), BCMVCPhysicsModel_EFT::pairCompare1);}
     else if(flag=="uncertainty"){std::sort(rank.begin(), rank.end(), BCMVCPhysicsModel_EFT::pairCompare2);}


      // print ranking
      for(unsigned int i=0; i<rank.size();i++){
        int index=rank.at(i).first;
      if(flag=="measurement"){out<<" "<<std::setw(5)<<i+1<< std::setw(12)<<std::right<<GetMeasurement(index)->GetName()<<"  "<<std::left <<std::setw(23)<<fObservables.at(index)->GetName()<< std::setw(8)<<rank.at(i).second*100<<std::endl;}
      else if(flag=="uncertainty"){out<<" "<<std::setw(5)<<i+1<< std::setw(12)<<std::right<<GetUncertainty(index)->GetName()<<"                  "<<std::left << std::setw(8)<<-1*rank.at(i).second*100<<std::endl;}
      }

    out<<"\n";
    Counter_2D++;
  }
 }
}



// ---------------------------------------------------------
void BCMVCPhysicsModel_EFT::PrintRanking1D(std::string flag, std::ofstream& out,double p, std::vector<double> orig_areas1D,std::vector< std::vector<std::pair<int,double> > > areas1D){

 unsigned int NParams = GetNParameters()-fNNuisanceCorrelation; // number of parameters without nuisance parameters

  out<<"Results of the one-dimensional "<<flag<<" ranking based on the "<<p*100<<" % interval: \n"<<std::endl;

  for(unsigned int j=0; j<NParams; j++){ 
    if(GetParameter(j)->Fixed()!=0){continue;}

    out<<"Ranking based on parameter \""<<GetParameter(j)->GetName()<<"\" (parameter "<<j<<"):"<<std::endl;
   if(flag=="measurement"){out<< std::left << std::setw(5)<<"Rank"<<"        "<<"Measurement"<<"                 "<<"Increase (%)"<<std::endl;
    out<<"------------------------------------------------------"<<std::endl;}
   else if (flag=="uncertainty"){out<< std::left << std::setw(5)<<"Rank"<<"        "<<"Uncertainty"<<"          "<<"Decrease (%)"<<std::endl;
    out<<"------------------------------------------------"<<std::endl;}

    std::vector<std::pair<int,double> > rank; // first: measurement_index, second: rank

    // for Measurements
    if(flag=="measurement"){
      // calculate ranking
      for(int i=0; i<GetNMeasurements();i++){
        if(GetMeasurement(i)->GetFlagActive()==0){continue;}
        double x = (areas1D[j][i].second-orig_areas1D.at(j))/orig_areas1D.at(j);
        rank.push_back(std::make_pair(i,x));
      }
    }

    // for uncertainties
    else if(flag=="uncertainty"){
      // calculate ranking
      for(int i=0; i<GetNUncertainties();i++){
        if(GetUncertainty(i)->GetFlagActive()==0){continue;}
        double x = (areas1D[j][i].second-orig_areas1D.at(j))/orig_areas1D.at(j);
        rank.push_back(std::make_pair(i,x));
      }
    }

    if(flag=="measurement"){std::sort(rank.begin(), rank.end(), BCMVCPhysicsModel_EFT::pairCompare1);}
    else if(flag=="uncertainty"){std::sort(rank.begin(), rank.end(), BCMVCPhysicsModel_EFT::pairCompare2);}
    
    // print ranking
    for(unsigned int i=0; i<rank.size();i++){
      int index=rank.at(i).first;
      if(flag=="measurement"){out<<" "<<std::setw(5)<<i+1<< std::setw(12)<<std::right<<GetMeasurement(index)->GetName()<<"  "<<std::left <<std::setw(23)<<fObservables.at(index)->GetName()<< std::setw(8)<<rank.at(i).second*100<<std::endl;}
      else if(flag=="uncertainty"){out<<" "<<std::setw(5)<<i+1<< std::setw(12)<<std::right<<GetUncertainty(index)->GetName()<<"                  "<<std::left << std::setw(8)<<-1*rank.at(i).second*100<<std::endl;}
    }
    out<<"\n";
  }
}

// ---------------------------------------------------------

void BCMVCPhysicsModel_EFT::GetUncertaintyRanking(double p, std::string file){
  GetRanking("uncertainty",p,file);
}

void BCMVCPhysicsModel_EFT::GetMeasurementRanking(double p, std::string file){
  GetRanking("measurement",p,file);
}

// ---------------------------------------------------------

double BCMVCPhysicsModel_EFT::GetArea1D(int index, double p){
  BCH1D* bchist = GetMarginalized(GetParameter(index));
  
  std::vector< double > s;
  s=bchist->GetSmallestIntervals(p); // s.size()=5*i (0: min, 1: max, 2: relative height, 3: local mode, 4: relative area)
  double sum=0; 
  // add area of all smallest intervals
  for (int i = 0; i <(s.size()/5.); i++){
      double x=0;
      x=s.at(i*5+1)-s.at(i*5);
      sum+=x;
    }  
    
  return sum;
}

// --------------------------------------------------------- 

double BCMVCPhysicsModel_EFT::GetArea2D(int i1, int i2, double p){
  BCH2D* bchist = GetMarginalized(GetParameter(i1),GetParameter(i2));

  // copy histograms for bands
   TH2D hist_temp(*bchist->GetHistogram());

  //integral
   double integral = hist_temp.Integral();

   // calculate total number of bins
   int nbins = hist_temp.GetNbinsX() * hist_temp.GetNbinsY();

   // area
   double area = 0.;

   // integrated probability
   double intp = 0.;

   // a counter
   int counter = 0;

   // loop over bins
   while ( (intp < p) && (counter < nbins) ) {

      // find maximum bin
      int binx;
      int biny;
      int binz;
      hist_temp.GetBinXYZ(hist_temp.GetMaximumBin(), binx, biny, binz);

      // increase probability
      double dp = hist_temp.GetBinContent(binx, biny)/integral;
      intp += dp ;

      // reset maximum bin
      hist_temp.SetBinContent(binx, biny, 0.);

      // increase area
      area += hist_temp.GetXaxis()->GetBinWidth(binx) * hist_temp.GetYaxis()->GetBinWidth(biny);

      // increase counter
      counter++;
   }

  return area;
}



// --------------------------------------------------------- 
 
bool BCMVCPhysicsModel_EFT::pairCompare1(const std::pair<int,double>& a, const std::pair<int,double>& b) {
  return a.second > b.second;

}

bool BCMVCPhysicsModel_EFT::pairCompare2(const std::pair<int,double>& a, const std::pair<int,double>& b) {
  return a.second < b.second;

}