# _Escherichia coli_ (str. K12) growing on different carbon sources (lactose and D-glucose)
![Flyer](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/notes/4Md/Flyer_Sum_G-L_model.png)

## Petri Net-based model of bacteria growing with different combinations of carbon sources
**Epimod** [[1](#references)] is the first analysis tool that can integrate chemical reaction networks analysis and Petri Net-based dynamical models (based on Ordinary Differential Equations) into a general modelling framework. This repository contains a working example that introduces the framework's ability to simulate multi-state systems. We performed this analysis showing the effect of glucose and lactose concentrations on the metabolic response of _E. coli_. 

The model can readily identify different metabolic behaviours depending on environmental status providing the translation of quantitative values into discrete intervals and allowing flux analyses. Flux Balance Analysis (FBA) become an established constraint-based method approach to studying metabolic networks. However, mechanism-based modelling approaches are the best optionto to model with finer-grained biological detail systems up perturbations. By the way, this method is not only based on the steady-state analysis of the metabolic network but on the dynamics of external metabolites concentration over time given an initial state. The so-implemented time-dependent analysis is based on an input metabolic network. At each time, the values of external metabolites concentrations are updated by solving ODEs and values are then translated into constraints for FBA. In turn, flux analòysis solution can be upcycled as kinetic parameters for the ODEs.

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

## Harmonization of metabolic networks and dynamic models

The GreatMod is a framework for quantitative modelling of biological systems and biochemical networks with constraints-based modelling. The software includes a broad set of modelling methods, among which generation of model and Petri Net-based analysis techniques. GreatMod is presently used for model analysis and simulations of metabolic systems using chemical reaction networks.

### Download and generate the file to pass to GreatMod

Metabolic models can be loaded into GreatMod, which supports models in MAT-file (.mat) format. This commonly used model format is a MATLAB structure containing fields as noted in The COBRA Toolbox [Documentation](https://github.com/opencobra/cobratoolbox/blob/master/docs/source/notes/COBRAModelFields.md).

**Databases**

Before designing a model using GreatMod, it might be informative to browse for some of the main platforms integrating open-source genome-scale metabolic models (GEMs). At the current state-of-art and applications of GEMs, we report three projects focused on high-throughput generation of GEMs: (i) [Path2Models](https://www.ebi.ac.uk/biomodels/path2models, (ii) [AGORA](https://github.com/VirtualMetabolicHuman/AGORA/tree/master/CurrentVersion), (iii) [CarveMe](https://github.com/cdanielmachado/carveme). These studies allowed the generation of the biggest number of GEMs (especially for bacteria) [[4](#references)]. The largest part of GEMs is now available and collected into high-quality databases.


| Repository      | About GEMs | Website     |
| :---        |    :----:   |          ---: |
| BiGG Models      | 84 high quality manually-curated GEMs         | [http://bigg.ucsd.edu](http://bigg.ucsd.edu)|
| Human Metabolic Atlas   |   Human, microbial species, fruitfly, mouse, rat, zebrafish, C. elegans GEMs      | [https://metabolicatlas.org/](https://metabolicatlas.org/)      |
| Virtual Metabolic Human   | 818 GEMs for gut microbes (AGORA models) and Recon3D         | [http://www.vmh.life](http://www.vmh.life)    |
| BioModels   | Over 1000 GEMs and over 640 small-scale metabolic models        | [http://www.ebi.ac.uk/biomodels](http://www.ebi.ac.uk/biomodels)      |
| MEMOSys 2.0    | GEMs in a standardized format        | [http://memosys.i-med.ac.at](http://memosys.i-med.ac.at)      |
| Model SEED   | Plant GEMs        | [http://modelseed.org](http://modelseed.org)      |
| MetaNetX    | GEMs collected from GEM-relevant databases        | [http://www.metanetx.org](http://www.metanetx.org)      |

Human and microbial models and maps are available at the [Virtual Metabolic Human](http://www.vmh.life) website, which collects human and gut metabolism data and links this information with nutrition. [BiGG Models ](http://bigg.ucsd.edu/) is another database that integrates several GEMs. Notably, BiGG allocates only networks with a set of standardized identifiers called BiGG IDs. [Human Metabolic Atlas](https://metabolicatlas.org/) goal is to collect curated GEMs and to bring these models closer to FAIR principles.

**Compiling procedure**

We will use a MAT-file formated AGORA model _E. coli_ downloaded from [Virtual Metabolic Human](http://www.vmh.life). First, we will load the metabolic model directory [Models](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/Input/Models). To compile model files for GreatMod, we will use the (?) function. A file with the given filename containing the model in the specified text format will be generated and saved to the directory [CompiledModels](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/Input/CompiledModels).

The procedure consists of the following steps:

* Select a metabolic model to convert into a complied model and enter parameters
* Identification of boundary reactions and metabolite exclusively involved in these reactions
* Summary of pivotal properties of the metabolic model

The time required to compile a model depends on the model's complexity and size. The compiling procedure could take a few minutes to conclude. In general, the following fields (mandatory fields) should always be present after compiling :

* Identifiers of the reactions
* Stoichiometric matrix
* Lower and upper bounds of the reactions
* Definition of the objective function 

This format only uses the essential metabolic modelling fields and will neglect other optional fields. Optional fields include:

* Identifiers of the metabolites
* List of associated-reaction genes
* Gene-protein-reaction rules

### Petri Net model description

This working example aim is to map the carbon sources available in bacteria environments to their representations in the metabolic model. Of primary importance is that the list of compounds unequivocally matches metabolites coded on the model. The results of model simulations rely on the accurate mapping of metabolites to specific exchange reactions. Exchange reactions are pseudo-reactions that allow direct coupling of metabolites with the compartment external to the metabolic model (i.e. the hypothetical biochemical adjacent environment of the model). The biochemical exchanges interested in the culture media are described by the Petri Model model shown in [Fig. 1](#Flyer). The model consists of two places and six transitions. Transitions modelling different kinetics whose parameters are defined in a dedicated CSV file. The C++ file coordinates the Linear Programming Problems (LPP) solving and general transitions behaviour. Metabolic reactions were converted into irreversible format to render the exchange reactions in Petri Net formalism. Reversible reactions split into two reactions, and the corresponding reaction identifiers are prolonged by "in" or "out" to avoid losing original direction information. The split reactions will have a lower bound  < 0 ("in") and upper bound > 0 ("out"), and lower bounds equal 0 and upper bounds as maximum uptake flux).

This model is organized into two modules corresponding to the environments characterizing the batch culture: Metabolites and FBAenvironment.

**Metabolites**

Metabolites are characterized by two places: the glc_e and lcts_e. The transitions Tglc and Tlcts represent the inflow of fresh nutrients (glucose and lactose, respectively) from the external. Its rate is defined by the general function CarbonFlux, to govern the nutrient delivery in the system. The function CarbonFlux is parametrized utilizing textual files (see _Repository files_ section) modelling different carbon sources administration regimes. In the case of the metabolic modelling approach, which goal is to predict the bacterial's metabolic reactions rates (h<sup>-1</sup>), the concentrations defined after modelling steps will be transformed to the maximum uptake fluxes for metabolic reactions in the unit *mmol/(gDWh)*.

**FBAenvironment**

Our practical example is built on the _E. coli_ K-12 MG1655 model downloaded from the database of models [Virtual Metabolic Human](http://www.vmh.life). In all, 1347 compounds and 1795 reactions are present in the model metabolite list. In metabolic modelling, exchange reactions are conventionally associated with identifiers with the prefix “EX” followed by the identifier of the metabolite and the suffix identifying the external compartment (“(e)” or "e" or "[e]"). Variations are according to the nomenclature of the underlying reaction database of reference).

Importantly, flux estimations for the exchange reactions ought to be constrained to represent the availability of the medium nutrients. Exchange reactions were constrained by upper and lower bounds, where the uptake of a compound is represented by a negative and excretion by a positive rate. This operation focuses on adjusting the exchange reactions linked to the metabolite identifiers. Furthermore, this control is essential for the accuracy of the flux solutions, as the FBA resolution depends on the availability of each metabolite according to the actual concentration. For two (i.e. glucose and lactose), we identified the appropriate exchange reactions in the model, so they are added in Petri Net as FBA-related general transitions. Concentrations of a given metabolite modelled in the Petri Net model constrained the matched exchange reactions in the metabolic model. Unmatched exchange reactions are constrained according to the concentration of diet formulation. Lower bounds were set to a default value for those compounds unlisted in the computational diet. Computational diet data locates in the folder [diets](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Input/diets). The file “**diets.RData**” is an R compressed file containing dietary data for six different nutritional information needed for the _in-silico_ representation of the chemical environment.

## Model Analysis

### How to read the trace and fluxes generated

The _E. coli_ growth in modelled glucose-lactose medium results are collected in the folder [Results](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/Results). The ODEs solver output filename can be specified with "analysis", it receives automatically the extension ".trace". The filename is optional, the default name is attributable to the name of file storing the Petri Net model, e.g. "PetriNetName-analysis.trace". Similarly, FBA problem solutions output filename boils down to "<#>.flux" (where "#" is the indeitfier number of the FBA problem). Also this file contains the value for the objective function.

### How to do the FVA

Flux variability analysis (FVA) allows addressing new questions to study the flexibility of chemical reaction networks in different environmental and internal conditions. Using GreatMod is now possible to conceive computational experiments to explore the range of multiple solutions through FVA.

### Simulation 

The sample images generated are the following:

![ModelSim](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/Results/CarbonSourcesUtilization.png)

We used the dynamical model to represent the direct effects of glucose and lactose supplementation on bacterial metabolic output. The ODE model couples to the metabolic model via the concentration of the two carbon sources, capturing the variation in sugar content of the batch reactor over time (bottom-right and left). The resulting unified model predicted the carbon source consumption given different administration scenarios. The prediction of sugar levels shows a variation in exchange reaction activity following the nutrient administration (upper-left). Furthermore, we predicted the effects of sugar batch content on enzymatic metabolic pathways key to mechanisms of catabolite repression in the bacterium during growth on differential carbon sources. Indeed, in addition to the bacterial sugar transport system, represented by reaction lactose transport via proton symport (LCTSt), lactose also influences b-galactosidase (LACZ) activity (upper-right).

## References
[1] (https://github.com/qBioTurin/epimod)

[2] [Feist AM, Henry CS, Reed JL, Krummenacker M, Joyce AR, Karp PD, et al. A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. Mol Sys Biol. 2007;3:121](https://www.embopress.org/doi/full/10.1038/msb4100155)

[3] [Amparore, E.G., Balbo, G., Beccuti, M., Donatelli, S., Franceschinis, G. 30 Years of GreatSPN.  In: Fiondella, L., Puliafito, A. (eds) Principles of Performance and Reliability Modeling and Evaluation. Springer Series in Reliability Engineering. Springer, Cham. 2016](https://link.springer.com/chapter/10.1007/978-3-319-30599-8_9)

[4] [Gu, Changdai & Kim, Gi & Kim, Won Jun & Kim, Hyun & Lee, Sang Yup. Current status and applications of genome-scale metabolic models. Genome Biology. 20. 10.1186/s13059-019-1730-3. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1730-3)
