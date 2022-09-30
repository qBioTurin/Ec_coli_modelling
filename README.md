# _Escherichia coli_ (str. K12) growing on different carbon sources (lactose and D-glucose)
![Flyer](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/notes/4Md/Flyer_Sum_G-L_model.png)

## Petri Net-based model of bacteria growing with different combinations of carbon sources
**Epimod** [[1](#references)] is the first analysis tool that can integrate chemical reaction networks analysis and Petri Net based dynamical models (based on Ordinary Differential Equations) into a general modeling framework. In this repository a working example where the framework's ability to simulate multi-state systems is introduced. We performed this analysis showing the effect of glucose and lactose concentrations on the metabolic response of _E. coli_. This can readily identify different metabolic behaviours depending on environmental status providing transaltion of quantitative values into discrete intervals and allowing flux analyses. Flux Balance Analysis (FBA) become an established constraint-based method approach to study metabolic networks. However, to model of finer-grained biological detail and of systems up perturbations, mechanism-based modeling approaches are the best option.  By the way, this method is not only based on the steady-state analysis of the metabolic network, but also on the dynamics of external metabolites concentration over time given and initial state. The so implemented time-dependent analyis is based on a input metabolic network. At each time, the values of external metabolites concentrations are updated by solving ODEs and values are then transalted into constraints for FBA. In turn, flux analòysis solution can be upcycled as kinetic parameters for the ODEs.

Simulations were performed using the genome-scale metabolic model of _E. coli_ K-12 MG1655 [[2](#references)] while the nclusion of specific environment is achieved by tuning constraints of target metabolic model's boundary reactions. 

## Repository files
The model describes metabolism for _E. coli_ str. K-12 MG1655. Culture conditions can be simulated by applying the respective boundaries set, which is implemented as fluxes values elaborated from nutritional data available at [Virtual Metabolic Human](https://www.vmh.life/#nutrition) database. The paradigm and computational environemnt applyed to working with Petri net is described in the [GreatSPN home page](http://www.di.unito.it/~amparore/mc4cslta/editor.html)

Following scenarios of carbon sources culture administration are simulated:

  * [AlternedPulses](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv/AlternedPulses.csv)
  * [PairedPulses](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv/PairedPulses.csv)
  * [StepwisePulses](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv/StepwisePulses.csv)

To simulate the experiments, a standard application **epimod** protocols are implemented within containerized environments.

* Model parametrization can be made by modifing the values of the selected parameters. Sets of parameters can be imported into a textual file in which the parameters to be studied are listed associated with their range of variability. Example of such a file is provided as “**CarbonAdmin.csv**”. The file includes a sheet that is utilized either to manage parameters deposited into additional sets or to recall methods of paramterization.

* Simulations can be performed in the R environemnt. The script file “**Main.R**” shows how to simulate the above mentioned protocols from a single generic simulation file. Also, the scripts generate images comparing simulation results given differnt conditions dataset as in “**CarbonSourcesUtilization.pdf**”. Metabolic model used is a plain test file provided as “**FBAmodelK12**” is located in the folder [CompiledModels](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/CompiledModels). This file was obtained by compilation of the corresponding S4 object of the class modelorg provided as “**Ec_K12Diet.RData**”) located in the folder [K12](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/Models/K12). 

* The script applies the parameter values defined in the folder [csv](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/csv). By changing values in csv file (or addressing alternative files), new parameter sets can be run out without altering the model itself. Petri Net model definition is provided by two file located in the folder [Net](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/Models/Net). The file “**Ec_corePN.PNPRO**” stores the Petri Net model. File “**transitions.cpp**” is a C++ file defining the functions managing the behaviour Petri Net general transition.

The sample images generated are the following:

![ModelSim](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Results/CarbonSourcesUtilization.pdf)


## References
[1] (https://github.com/qBioTurin/epimod)

[2] [Feist AM, Henry CS, Reed JL, Krummenacker M, Joyce AR, Karp PD, et al. A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. Mol Sys Biol. 2007;3:121](https://www.embopress.org/doi/full/10.1038/msb4100155)

[3] [Amparore, E.G., Balbo, G., Beccuti, M., Donatelli, S., Franceschinis, G. 30 Years of GreatSPN.  In: Fiondella, L., Puliafito, A. (eds) Principles of Performance and Reliability Modeling and Evaluation. Springer Series in Reliability Engineering. Springer, Cham. 2016](https://link.springer.com/chapter/10.1007/978-3-319-30599-8_9)
