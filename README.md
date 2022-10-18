# _Escherichia coli_ (str. K12) growing on different carbon sources (lactose and D-glucose)
![Flyer](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/notes/4Md/Flyer_Sum_G-L_model.png)

## Petri Net-based model of bacteria growing with different combinations of carbon sources
**Epimod** [[1](#references)] is the first analysis tool that can integrate chemical reaction networks analysis and Petri Net-based dynamical models (based on Ordinary Differential Equations) into a general modelling framework. This repository contains a working example that introduces the framework's ability to simulate multi-state systems. We performed this analysis showing the effect of glucose and lactose concentrations on the metabolic response of _E. coli_. This can readily identify different metabolic behaviours depending on environmental status providing the translation of quantitative values into discrete intervals and allowing flux analyses. Flux Balance Analysis (FBA) become an established constraint-based method approach to studying metabolic networks. However, mechanism-based modelling approaches are the best optionto to model with finer-grained biological detail systems up perturbations.  By the way, this method is not only based on the steady-state analysis of the metabolic network but on the dynamics of external metabolites concentration over time given an initial state. The so-implemented time-dependent analysis is based on an input metabolic network. At each time, the values of external metabolites concentrations are updated by solving ODEs and values are then translated into constraints for FBA. In turn, flux analòysis solution can be upcycled as kinetic parameters for the ODEs.

Simulations were performed using the genome-scale metabolic model of _E. coli_ K-12 MG1655 [[2](#references)] while the inclusion of a specific environment is achieved by tuning constraints of the target metabolic model's boundary reactions. 

## Repository files
The model describes metabolism for _E. coli_ str. K-12 MG1655. Culture conditions can be simulated by applying the respective boundaries set, which is implemented as fluxes values elaborated from nutritional data available at [Virtual Metabolic Human](https://www.vmh.life/#nutrition) database. The paradigm and computational environment applied to working with Petri net is described in the [GreatSPN home page](http://www.di.unito.it/~amparore/mc4cslta/editor.html)

The following scenarios of carbon sources culture administration are simulated:

  * [AlternedPulses](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv/AlternedPulses.csv)
  * [PairedPulses](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv/PairedPulses.csv)
  * [StepwisePulses](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv/StepwisePulses.csv)

To simulate the experiments, standard application **epimod** protocols are implemented within containerized environments.

* Model parametrization can be made by modifying the values of the selected parameters. Sets of parameters can be imported given a textual file in which the parameters to be studied are listed and associated with their range of variability. An example of such a file is provided as “**CarbonAdmin.csv**”. The file includes a sheet that is utilized either to manage parameters deposited into additional sets or to recall methods of parametrization.

* Simulations can be performed in the R environment. The script file “**Main.R**” shows how to simulate the above-mentioned protocols from a single generic simulation file. Also, the scripts generate images comparing simulation results given different conditions dataset as in “**CarbonSourcesUtilization.pdf**”. The metabolic model used is a plain text file provided as “**FBAmodelK12**” is located in the folder [CompiledModels](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/CompiledModels). This file was obtained by compilation of the corresponding S4 object provided as “**Ec_K12Diet.RData**”) located in the folder [K12](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/Models/K12). 

* The script applies the parameter values defined in the folder [csv](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv). By changing values in the CSV file (or addressing alternative files), new parameter sets can be run out without altering the model itself. Petri Net model definition is provided by two files located in the folder [Net](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/Models/Net). The file “**Ec_corePN.PNPRO**” stores the Petri Net model. File “**transitions.cpp**” is a C++ file defining the functions managing the behaviour of Petri Net's general transition.


## Harmonization of metabolic network and dynamic model

scaricare dal sito il file matlab, usare le funzioni R che offriamo (to do) per generare 

### Download and generate the file to pass to GreatMod

### Net description (transitions etc)

## Model Analysis

### How to read the trace and fluxes generated

### How to do the FVA
(TO DO)

### Simulation 

The sample images generated are the following:

![ModelSim](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Results/CarbonSourcesUtilization.png)

MANCA UNA DESCRIZIONE
## References
[1] (https://github.com/qBioTurin/epimod)

[2] [Feist AM, Henry CS, Reed JL, Krummenacker M, Joyce AR, Karp PD, et al. A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. Mol Sys Biol. 2007;3:121](https://www.embopress.org/doi/full/10.1038/msb4100155)

[3] [Amparore, E.G., Balbo, G., Beccuti, M., Donatelli, S., Franceschinis, G. 30 Years of GreatSPN.  In: Fiondella, L., Puliafito, A. (eds) Principles of Performance and Reliability Modeling and Evaluation. Springer Series in Reliability Engineering. Springer, Cham. 2016](https://link.springer.com/chapter/10.1007/978-3-319-30599-8_9)
