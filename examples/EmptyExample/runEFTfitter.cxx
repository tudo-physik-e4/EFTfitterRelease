#include <string>

#include "BAT/BCLog.h"
#include "BAT/BCAux.h"
#include "BAT/BCParameter.h"
#include "BAT/BCSummaryTool.h"
#include "BAT/BCModelOutput.h"

#include "EFTfitter.h"

//------------------------------------------------------------------------------------------

int main(){

// set nicer style for drawing than the ROOT default
BCAux::SetStyle();

// open log file
BCLog::OpenLog("results/log.txt");
BCLog::SetLogLevelFile(BCLog::detail);
BCLog::SetLogLevelScreen(BCLog::summary);

// create new instance EFT model
EFTfitter * m  = new EFTfitter();

  // read parameter-settings from file
  std::string input_parameters("input/parameters.xml");
  m -> ReadParameters(input_parameters);
  
  // read measurements from file
  std::string input_measurements("input/measurements.xml");
  m -> ReadMeasurements(input_measurements);

  // set precision
  m -> MCMCSetPrecision(BCIntegrate::kLow); // options: kLow / kMedium / kHigh / kVeryHigh      

  // write model output to ROOT file    
  BCModelOutput *mout = NULL;
  mout = new BCModelOutput(m, "results/result_plots.root"); 
  mout->WriteMarkovChain(true);
   
  // plot all observables as a function of all free parameters as .pdf and .root file
  m->PlotObservables("results/observables");

  // do analysis    
  m -> MarginalizeAll();    
  m -> FindMode( m ->GetBestFitParameters() );

  // calculate ranking of measurements based on the given probability
  //m->GetMeasurementRanking(0.683,"results/measurement_ranking.txt");
  //m->GetUncertaintyRanking(0.683,"results/uncertainty_ranking.txt");

  // make ḱnowledge-update-plot and plot correlation between parameters
  BCSummaryTool * summary = new BCSummaryTool(m);
  summary -> PrintKnowledgeUpdatePlots("results/knowledge_update.pdf");
  summary -> PrintCorrelationMatrix("results/parameter_correlation.pdf");
    
  // write marginalized distributions to output ROOT-file                                                                                                         
  mout->WriteMarginalizedDistributions();                                                                                                                                   
  mout->Close();

   // print results to file and make plots    
  m -> PrintResults("results/results.txt");    
  m -> PrintAllMarginalized("results/result_plots.pdf");
  
  // show results on screen
  m->PrintSummary();
  
  delete m;

  // close log file
  BCLog::CloseLog();

return 0;
}
//------------------------------------------------------------------------------------------