# EFTfitter version 1.03

## About
EFTfitter is a tool for interpreting measurements in the context of effective field theories. It allows to combine measurements and take their correlations into account.

## Authors
Nuno Castro, Johannes Erdmann, Cornelius Grunwald, Kevin Kroeninger, Nils-Arne Rosien

## Versions
* 1.03: fix nuisance parameters
* 1.02: fixes for the out-of-the-box example
* 1.01: one minor fix
* 1.0: first version

## Required software
* BAT (Tested with version 0.9.4.1, BAT 1.0 will not work at the moment), https://www.mppmu.mpg.de/bat/
* ROOT (Tested with version 5.34)

## Content
* The folder "examples/EmptyExample" contains an empty example for implementing a new analysis.
* The folder "examples/AnomalousCouplings" contains an out-of-the-box physics example for anomalous couplings at the Wtb-vertex.
* The folder "src" contains the source files for EFTfitter, which normally do not need to be modified.

## Structure
* Model: The theoretical model that provides the functional relation between the measurable quantity (the observables) and the model parameters has to be defined in the function "CalculateObservable".
* Input: The parameters have to be defined in "parameters.xml" (name, range, prior). The observables have to be defined in "measurements.xml" (name, range) and the measured values, uncertainties and correlations of the observables have to be provided.

## Features
The following additional functions can be called in the run file.
* PlotObservables(string outputfile): Plots the functional relation between the model parameters and the observables defined in "CalculateObservable".
* GetMeasurementRanking(double p, string outputfile): Calculates the ranking of all activated measurements based on the impact of each measurement on the 1D and 2D marginalized distributions of the free parameters, i.e. one measurement is deactivated at a time and the resulting increase of the smallest area containing the given probability p is used as the ranking criterion.
* GetUncertaintyRanking(double p, string outputfile): Analogously to GetMeasurementRanking, the uncertainties are deactiveted one-by-one and the decrease in the probability area is used to calculate the ranking.

## Results
The resulting plots and values of the analysis will be found in the folder "results".
 1. "result_plots.pdf" and "result_plots.root" contain histogramms of all marginalized 1D and 2D distributions of the parameters.
 2. "results.txt" contains the numerical results of the analysis.
 3. "parameter_correlation.pdf" contains the estimated correlation matrix of the parameters.
 4. "knowledge_update.pdf" contains a comparison between the prior and posterior distributions of the parameters.
 5. "observables.pdf" and "observables.root" contain plots of the functional relations between the observables and the parameters.
 6. "measurement_ranking.txt" contains the results of the measurement ranking.
 7. "uncertainty_ranking.txt" contains the results of the uncertainty ranking.

## Run
1. compile using "make"
2. run: "./runEFTfitter" (or: "./runAnomalousCouplings")

After changing the input xml-files it is sufficient to run "./runEFTfitter" and not necessary to compile again.

## License
EFTfitter is distributed under the GNU lesser general public license. A copy of the license is included in the file LICENSE.

## Credit
The folder "src/tinyxml2" contains the (unmodified) xml parser "TinyXML-2", which has been released on: "http://www.grinninglizard.com/tinyxml2docs/index.html" under the zlib license. The folder also contains a copy of the zlib license.
