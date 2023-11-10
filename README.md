# Using _UnifiedGreatMod_ to integrate transcriptional data onto _Escherichia coli_ genome-scale metabolic model (iML1515) growing on different regime of carbon feeding

![flyer](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/4md/Flyer_Sum_G-L_model.png)

## Petri Net-based model of _E. coli_ strain K-12 growing with different combinations of carbon sources

### Background
**Epimod** [ [1](#references) ] is an analysis tool that can integrate Flux Balance Analysis (FBA) and Petri Net dynamical models based on Ordinary Differential Equations (ODEs) into a general and unified modelling framework. This repository contains a working example introducing the framework's ability to simulate multi-state systems. The FBA optimization problem is solved along with constraints on the rate of change of metabolite levels at specific time instants. We have shown that **Epimod** can mechanistacally simulate the metabolic output of _E. coli_ growing under various transcriptional programs and environmental nutrient perturbations. Therefore **Epimod** provides a framework for analyzing the transience of metabolism due to metabolic reprogramming and obtaining insights for the design of metabolic networks.

### Methods
FBA become an established constraint-based method approach to studying metabolic networks. However, mechanism-based modelling approaches are the best option to model biological systems with finer-grained detail-up perturbations. Indeed, systems of ODEs can not only account for the steady-state analysis of the metabolic network but even the dynamics of medium culture metabolites concentration over time given an initial state. The so-implemented time-dependent analysis depends on an input metabolic network. At each time, the values of external metabolite concentrations are updated by solving ODEs and values are then translated into constraints for FBA. In turn, FBA solutions upcycled as kinetic parameters for the ODEs. In doing so, the flow of metabolites through the _E. coli_  metabolic network was calculated  given a specific media enriched with various carbon regimes supporting growth.

By imposing experimentally derived constraints, a metabolic model can be transformed into a condition-specific model. These constraints can be defined by setting upper and lower flux bounds for each reaction. A metabolic model can have several types of constraints that represent specific intra- and extracellular conditions, such as biomass maintenance requirements, environmental constraints, or enzyme capacities.

Simulations were performed using the genome-scale metabolic COBRA model of _E. coli_ K-12 strain MG1655 [ [2](#references) ] while the inclusion of a specific environment is achieved by tuning constraints of the target metabolic model's boundary reactions. We have implemented a strategy for integrating metabolomics and transcriptomics data. This approach uses constraint-based stoichiometric metabolic models as a framework. Our computational pipeline computes differential reaction activities from transcriptomics data [ [3](#references) ] to predict whether the differential expression of metabolic enzymes directly originates differences in metabolic fluxes.

### Results
The model can readily identify different metabolic behaviours depending on environmental and transcriptional status, and providing the translation of quantitative values into discrete intervals and allowing flux analyses. It is important to note that FBA makes several assumptions, e.g. steady-state metabolite concentrations (meaning there is no change in metabolite concentrations over time). These assumptions might not hold in the real-world batch culture scenario. To harmonize the FBA, which presumes a steady state for all metabolites, we exploited the Petri Net formalism to model the time dynamics of a select group of medium metabolites. 

Our results point to the steady-state constraint parametrization of fluxes as a determinant of model inaccuracy. To overcome these limitations, expanding the Petri Net model to represent other medium constituents could be a promising route toward more accurate models of metabolic physiology.

## Repository files

The model describes metabolism for _E. coli_ str. K-12. Cultural conditions are simulated by applying the respective boundaries set implemented as fluxes values. The paradigm and computational environment applied to working with Petri Net are described in the [GreatSPN home page](http://www.di.unito.it/~amparore/mc4cslta/editor.html). 

The following carbon source supplementation scenarios for a fed-batch were simulated:

* _Constant_feeding_:
The feed rate of glucose and lactose is constant throughout the culture.
* _Linear_feeding_:
The feed rate of glucose and lactose are increased linearly over time to match the growth rate of the cells
* _Pulsed_feeding_60_: Glucose and lactose are added in repeated cycles of short duration, followed by a period of no feeding; the same total amount of glucose was fed in repeated 300s (5 min) cycles for 60s
* _Pulsed_feeding_150_: Glucose and lactose are added in repeated cycles of short duration, followed by a period of no feeding; the same total amount of glucose was fed in repeated 300s (5 min) cycles for 150s
* _blank_:
The culture does not have glucose and lactose added to it.

The model was parameterized by adjusting the values of selected parameters. These sets of parameters were imported from a text file, where each parameter is listed and associated with its range of variability. An example of such a file is provided, named **CarbonAdmin.csv**. This file contains a sheet that can be used either to manage parameters stored in additional sets or to invoke parameterization methods.

Simulations were performed in the R environment. The script file **main.R** shows how to perform simulations from a single generic simulation file. Also, the scripts generate visualizations comparing simulation results given different conditions dataset as in [ [plots](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/results/plots) ]. The metabolic model used is a plain text file provided as **iML1515.txt** is located in the folder [ [compiled_models](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/compiled_models) ]. This file was obtained by compilation of the corresponding S4 object provided as **iML1515.RData** located in the folder [ [iML1515](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/models/ecoli/iML1515) ].

The script applies the parameter values defined in the folder [ [carbon_regimes](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/csv/carbon_regimes) ]. By changing values in the CSV file (or addressing alternative files), new parameter sets can be run out without altering the model itself. Petri Net model definition is provided by two files located in the folder [ [net](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/net) ]. The file **iML1515_PN.PNPRO** stores the Petri Net model. File **transitions.cpp** is a C++ file defining the functions managing the behaviour of Petri Net's general transition.

To simulate the experiments, standard application **Epimod** protocols are implemented within containerized environments.

## Harmonization of FBA and ODEs-based models

The *GreatMod* includes a broad set of modelling methods, among which Petri Net-based model generation and analysis techniques. Now, *UnifiedGreatMod* is presently used for model analysis and simulations of metabolic systems.

### Download and generate the file to pass to *GreatMod*
Metabolic models are loaded into *GreatMod*, which supports models in MAT-file (".mat") format. This commonly used model format is a MATLAB structure containing fields as noted in The COBRA Toolbox [Documentation](https://github.com/opencobra/cobratoolbox/blob/master/docs/source/notes/COBRAModelFields.md).

Before designing a model using *GreatMod*, it might be informative to browse for some of the platforms integrating open-source genome-scale metabolic models (GEMs). At the current state-of-art and applications of GEMs, we report three projects focused on high-throughput generation of GEMs: (i) [Path2Models](https://www.ebi.ac.uk/biomodels/path2models), (ii) [AGORA](https://github.com/VirtualMetabolicHuman/AGORA/tree/master/CurrentVersion), (iii) [CarveMe](https://github.com/cdanielmachado/carveme), (iv)  [BiGG Models](http://bigg.ucsd.edu/). These studies allowed the generation of high number of GEMs (especially for bacteria) [ [5](#references) ]. The most extensive part of GEMs is now available and collected into high-quality databases.
<br />

| Repository      | About GEMs | Website     |
| :---        |    :----:   |          ---: |
| BiGG Models      | 84 high quality manually-curated GEMs         | [http://bigg.ucsd.edu](http://bigg.ucsd.edu)|
| Human Metabolic Atlas   |   Human, microbial species, fruitfly, mouse, rat, zebrafish, C. elegans GEMs      | [https://metabolicatlas.org/](https://metabolicatlas.org/)      |
| Virtual Metabolic Human   | 818 GEMs for gut microbes (AGORA models) and Recon3D         | [http://www.vmh.life](http://www.vmh.life)    |
| BioModels   | Over 1000 GEMs and over 640 small-scale metabolic models        | [http://www.ebi.ac.uk/biomodels](http://www.ebi.ac.uk/biomodels)      |
| MEMOSys 2.0    | GEMs in a standardized format        | [http://memosys.i-med.ac.at](http://memosys.i-med.ac.at)      |
| Model SEED   | Plant GEMs        | [http://modelseed.org](http://modelseed.org)      |
| MetaNetX    | GEMs collected from GEM-relevant databases        | [http://www.metanetx.org](http://www.metanetx.org)      |

Human and microbial models and maps are available at the [Virtual Metabolic Human](http://www.vmh.life) website, which collects human and gut metabolism data and links this information with nutrition. [BiGG Models](http://bigg.ucsd.edu/) is another database integrating several GEMs. Notably, BiGG allocates only networks with standardized identifiers called _BiGG IDs_. [Human Metabolic Atlas](https://metabolicatlas.org/) goal is to collect curated GEMs and to bring these models closer to FAIR principles. 

We will use a MAT-file Escherichia coli str. K-12 substr (model: iML1515) downloaded from [BiGG Models](http://bigg.ucsd.edu/). The model has the following metrics: 
<br />

| Components | Counts |
| :--- | ---: |
| Metabolites| 1877 |
| Reactions | 2712 |
| Genes | 1516 |

First, we load the metabolic model directory [ [models](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/models/) ]. To compile model files for GreatMod, we use the functions collected [ [processing_functions](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/code/processing_function) ]. A file with the given filename containing the model in the specified text format is generated and saved to the directory [ [compiled_models](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/models) ]. The time required to compile a model depends on the model's complexity and size. The compiling procedure could take a few minutes to conclude. In general, the following fields would always be present after compiling:

* Identifiers of the reactions
* Stoichiometric matrix
* Lower and upper bounds of the reactions
* Definition of the objective function


The process can be broken down into these steps: (i) choose a metabolic model to compile, (ii) identify the reactions that need to be constrained, and (iii) introduce FBA constraints, which include determining the range of metabolites included in silico medium. Directory [ [iML1515](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/models/ecoli/iML1515) ] contains the following files:

* **iML1515.mat**: MATLAB structure containing one or more fields defined in The COBRA Toolbox Documentation.
* **iML1515.RData** R S4 structure containing the fields required for compiling
* **all_react.rds**: R dataframe structure containing the identifiers, and the lower and upper bounds of the reactions before contextualization (by convention, the bounds set on reaction rates in a metabolic model range from -1000 to 1000 for reversible, and from 0 to 1000 for irreversible reactions)
* **param_all_react.rds**: R dataframe structure containing the identifiers, and the lower and upper bounds of the reactions after contextualization
* **BIGGdata_iML1515_react.tsv**: TSV file storing identifiers and stochiometric formula of reactions in the COBRA model
* **equation.rds**: R vector containing the stochiometric formula of reactions in the COBRA model
* **genesfromGRrules.rds**: R vector storing the genes included in the model involved in gene-protein-reaction association relationships.
* **geni.rds**: R vector containing genes from the model (represented by their gene locus number (e.g. b2097)). The gene name can be achieved with a quick search of the EcoCyc website (https://www.ecocyc.org/).
* **GRrules.rds**: R vector containing the gene-protein-reaction association (i.e. the Boolean relationship between the genes that are required to produce a specific reaction
* **officialName.rds**: R vector containing the official
name of reactions in the COBRA model
* **subsystem.rds**: R vector referring to the subsystems included in the COBRA model.

### Petri Net model description
The model aims to map the carbon sources' availability variation in batch to their representations in the metabolic model. Concentrations of a given metabolite modelled in the Petri Net constrained the matched exchange reactions in the metabolic model. Unmatched exchange reactions are constrained according to the concentration of the diet input, constrained through the integration of further information, or kept unconstrained.

The biochemical exchanges interested in the batch media are described by the Petri Net model as shown in ( [Fig. 1](#flyer) ). The model consists of two places and six transitions. Transitions modelling different kinetics whose parameters are defined in a dedicated CSV file. The C++ file **transition.cpp** coordinates the Linear Programming Problems (LPP) solving and general transitions behaviour. Metabolic reactions were converted into irreversible format to render the exchange reactions in Petri Net formalism. Reversible reactions split into two components; the corresponding reaction identifiers are prolonged by "in" or "out" to avoid losing original direction information. The split reactions will have a lower bound < 0 ("in") and upper bound > 0 ("out"), and lower bounds equal 0 and upper bounds as maximum uptake flux.

Petri Net model is organized into two modules corresponding to the compartments characterizing the batch culture: **_Metabolites_** and **_FBAenvironment_** ( [Fig. 1](#flyer) ).

**_Metabolites_**
Metabolites are indicated by two places: the "glc_e" and "lcts_e". The transitions "Tglc" and "Tlcts" represent the inflow of fresh nutrients (glucose and lactose, respectively) from the external. Its rate is defined by the general function **CarbonFlux** governing the nutrient delivery in the system. The function **CarbonFlux** is parametrized utilizing textual files (see _Repository files_ section) modelling different carbon sources administration regimes. The "glc_e" and "lcts_e" amounts are transformed to the maximum uptake fluxes for metabolic reactions in the unit *mmol/(gDWh)* after each modelling steps.

**_FBAenvironment_**
Our practical example assembles on the _E. coli_ K-12 MG1655 model downloaded from the database of models [VMH](http://www.vmh.life). Here, exchange reactions are conventionally associated with identifiers with the prefix "EX" followed by the identifier of the metabolite and the suffix indicating the external compartment ("_e"). Variations are according to the nomenclature of the underlying reaction database of reference.

### Set dietary inputs and others constraints
Importantly, flux estimations for the exchange reactions ought to be constrained to represent the availability of the medium nutrients. Exchange reactions are pseudo-reactions that allow direct coupling of metabolites with the compartment external to the metabolic model (i.e. the hypothetical biochemical adjacent environment of the model). Exchange reactions are constrained by upper and lower bounds, where the uptake of a compound is represented by a negative and efflux by a positive rate. This operation focuses on adjusting the exchange reactions linked to the metabolite identifiers. For two glucose and lactose, we identified the appropriate exchange reactions in the model, so they are added in Petri Net as FBA-related general transitions (the transitions called **FBA**). Of primary importance is that the list of compounds unequivocally matches metabolites coded on the model. The results of model simulations rely on the accurate mapping of metabolites to specific exchange reactions.

The expression levels were obtained from data of _E. coli_ growing into a glucose-enriched minimal medium obtained via the commercial microarray platform Affymetrix GeneChip system. Expression data integrated into this modelling coming from MG1655 cells cultivated aerobically in MOPS minimal medium supplemented with 0.1% glucose [ [7](#references) ]. Gene expression data integration allows the metabolic model to reflect not just the internal state of the cell but also the external environment in which the cell is growing. In Genome-Scale Metabolic Modelling, integrating nutrient concentrations through the adjustment of bounds associated with exchange reactions is a way to simulate the availability of different nutrients in the environment. This can lead to more accurate predictions of cellular behaviour under a given environmental conditions. While gene expression data and media composition formulation were utilized to constrain a significant portion of the metabolic model, approximately 330 reactions remain unconstrained.

**Setting in silicoi media** 

We are going to adjust these bounds based on nutrient concentrations to simulate aerobic in MOPS minimal medium supplemented with 0.1% glucose environmental conditions:

- **Aerobic growth**: This refers to the growth of _E. coli_ strain K-12 MG1655 in the presence of oxygen. We consider the presence of environmental oxygen to efficiently generate energy through aerobic respiration. The reaction and "FHL" were constrained at _lb_ -1e-12 and _ub_ 1e-12 _mmol/gDW*h_ ( genes encoding FHL are known to be active under anaerobic conditions this reaction is constrained almost to zero to avoid unrealistic aerobic hydrogen production )

- **MOPS minimal medium**: This is a type of culture medium used for growing bacteria. It is "minimal" because it contains only the essential nutrients for bacterial growing. MOPS is a buffering agent that helps maintain the pH of the medium.

- **Supplemented with 0.1% glucose**: This means that glucose has been added to the medium as a carbon source, which the bacteria can use for energy and growth. The concentration of glucose is 0.1%, or 1 gram per liter.

- **MG1655 growth**: The reactions "BIOMASS_Ec_iML1515_core_75p37M" is constrained at _lb_ -1e-12 (_mmol/gDW*h_) to avoid biomass consuption. The doubling time of _E. coli_ strain K-12 MG1655 cultured in MOPS minimal medium supplemented with 0.1% glucose during the mid-log phase is approximately 63.5 minutes [ [8](#references) ]. Therefore, The reactions "BIOMASS_Ec_iML1515_core_75p37M" _ub_ is constrained to 0.655 (_mmol/gDW*h_) as:

        # Given data
        doubling_time <- 1.0583 # (h) 63.5 (minutes)
        
        # Define the molecular weight of biomass (MW)
        MW_biomass <- 1 # g mmol^(-1)

        # Calculate specific growth rate
        growth_rate <- log(2) / doubling_time

        # Calculate biomass flux (mmol/gDW*h)
        biomass_flux <- specific_growth_rate * MW_biomass

        # Print the result
        cat("Biomass flux (mmol/gDW*h):", biomass_flux, "\n")

The log(2) is used because the biomass doubles during each doubling time. The division by the doubling time (in hours) then gives the rate of biomass production per hour. This is the biomass flux corresponding to the upper bounds for the objective. The result will be in the unit _mmol/gDWh_, as required.

### Setting modelling conditions
* We could define a dry mass weight of a _E. coli_ cell cellular population equal to 1 (g)
* We could model the bacterial culture growing in a bioreactor with a given working volume
* We could set 1 (g/L) of glc_e initial concentration

**Setting bacterial biomass**

To estimate the number of cells that can constitute 1 g of bacterial dry mass weight, we can divide 1 g by the dry mass weight of a single _E. coli_ cell to estimate the dry mass weight of a _E. coli_ cellular population in batch [ [10](#references) ]: 

        gDW_ec = 3e-13 # (g/cell)
        n_ec_cell = 1/gDW_ec # (cells)

**Setting batch culture volume**

A saturated cell culture of _E. coli_ contains about 1e+09 cell/mL.

        # assuming number of cells during the log-mid phase is half of a saturated
        n_ec_cell_log_mid <- 0.5*1e+09 # (cells/mL)

        # Minimum volume that can contain these cells
        v_batch = n_ec_cell / n_ec_cell_log_mid # (mL)

**Setting initial glucose concentration**

To convert 1g/L of glucose to mmol/L, we need to know the molar mass of glucose:

        conc_glc_e_gL = 1 # (g/L)
        mm_glc = 180.16 # (g/mol)
        
        # then can then apply the following relation:
        conc_glc_e_mmolL = (conc_glc_e_gL/mm_glc)*1000 # (mmol/L)
        
        # setting glucose amount given a batch culture volume
        iG = (conc_glc_e_mmolL/1000)*v_batch # (mmol)

**Setting minimal medium**

From [ [9](#references) ] MOPS minimal medium as defined by (_Bergholz et al, 2007_) and rescale by the batch volume considered: 
<br />

| compound | amount_mM | amount_batch | bigg_id | exchanges |
| :--- | :----: | :----: | :----: | ---: |    
| D-Glucose | 5.55 | 18.5 | glc_D_e | EX_glc_D_e  |
| Magnesium chloride | 0.523 | 3.48 | mg2_e | EX_mg2_e    |
| Calcium chloride anhydrous | 0.001 | 0.00167 | ca2_e | EX_ca2_e    |
| Potassium sulfate | 0.276 | 1.84 | k_e | EX_k_e      |
| Sodium chloride | 50 | 334 | na1_e | EX_na1_e    |
| Ammonium chloride | 9.52 | 63.4 | nh4_e | EX_nh4_e    |
| Ferrous sulfate | 0.01 | 0.0666 | fe3_e | EX_fe3_e  |  
| Cobalt chloride | 0.00006 | 0.0004 | cobalt2_e | EX_cobalt2_e |
| Zinc sulfate | 0.0002 | 0.001334 | zn2_e | EX_zn2_e    |
| Manganese chloride | 0.0016 | 0.01066 | mn2_e | EX_mn2_e    |
| Cupric sulfate | 0.0002 | 0.001334 | cu2_e | EX_cu2_e  |  
| Molybdic acid ammonium salt tetrahydrate | 0.00006 | 0.0004 | mobd_e | EX_mobd_e   |

**Altenative constainting strategies**

Pre-made diets are also available at [_nutrition_](https://www.vmh.life/#nutrition) and in [_resources_](https://github.com/opencobra/COBRA.papers/tree/64f997ffb2ecac660381c6f4a401792ed7c165e0/2018_microbiomeModelingToolbox/input). We could simulate growth on a _Average European_ dietary formulation (this diet covers 58 out of 331 exchange reactions that were constraints on the dietary input). Computational diet data locates in the folder [ [diets](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/diets) ]. The file **diets.RData** is an R compressed file containing dietary data for six different nutritional information for the _in silico_ representation of the chemical environment.

### Gene Expression Preprocessing
The critical aspect of integrating transcriptomic data lies in the selection of the gene expression data integration method. In this practical demonstration, we have outlined the data integration process by the approach used in reference [ [6](#references) ]. The typical method for deriving reaction expression values from gene expression values involves evaluating the Gene-Protein-Reaction rules associated with each reaction.

The mapping utilized data from Affymetrix high-density oligonucleotide expression arrays [ [7](#references) ]. The study examined the impact of different carbon sources on the global gene expression of _E. coli_ MG1655 cells using the *GPL199* - Ecoli_ASv2 Affymetrix _E. coli_ Antisense Genome Array platform. The authors cultivated MG1655 cells aerobically in MOPS minimal medium with one of the following:

- Glucose
- Glycerol
- Succinate
- L-alanine
- Acetate
- L-proline

Samples were collected from each culture during the mid-log phase. The dataset, which is available in the GEO database under the identifier GSE2037, includes a total of 15 samples. We used a generic reconstruction and the set of expression data, specifically **GSM37063 | glucose replicate number 1**, to automatically derive a context-specific network. This particular data was gathered under conditions where E. coli cells were exposed to glucose as their carbon source.

**_Gene expression data integration_**

We will adjust the lower and upper bounds of a metabolic model based on gene expression data for _E. coli_ growing in a glucose-enriched minimal medium obtained through the commercial microarray platform Affymetrix GeneChip system. Gene expressions evaluated with the Gene-Protein-Reaction relations are saved in the folder [ [gene_expression](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/gene_expression) ] under the name **GSE_iML1515.csv**.

Each gene in the metabolic model is associated with one or more reactions through Gene-Protein-Reaction associations. These associations can be represented as Boolean expressions, where **AND** typically represents the interaction of multiple subunits of a protein complex, and **OR** represents isozymes for each reaction in the metabolic model, the gene expression level of the associated gene(s) is obtained from the Affymetrix data. If a reaction associates with multiple genes through an **AND** relationship, the minimum expression level among the associated genes could be used. If the relationship is **OR**, the maximum expression level could be used. The result is a vector where each element represents the gene expression level associated with a specific reaction in the metabolic model [ [6](#references) ]:

* **If the gene expression level is greater than or equal to 1:x** The lower and upper bounds of the corresponding reaction in the metabolic model are multiplied by ( _1 + log(GSE[reaction_i])_ ). When a gene is highly expressed the bounds of the corresponding reaction will be increases, allowing for more flux through that reaction. The use of the logarithm ensures that the adjustment is more sensitive to changes in gene expression when the expression level is low and less sensitive when the expression level is high.

* **_ If the gene expression level is less than 1:** The lower and upper bounds of the corresponding reaction in the metabolic model are divided by ( _1 + abs(log(GSE[reaction_i]))_ ). When a gene is lowly expressed, the bounds of the corresponding reaction will be decreased, restricting the flux through that reaction. The use of the absolute value of the logarithm ensures that the adjustment is always positive, even when GSE[reaction_i] is less than 1.

## Model Analysis

![results](https://github.com/qBioTurin/Ec_coli_modelling/blob/main/4md/results_figures.png)

### How to read the trace and fluxes generated
The simulation results are collected in the folder [ [results](https://github.com/qBioTurin/Ec_coli_modelling/tree/main/input/results) ]. The ODEs solver output filename can be specified with "-analysis", then it receives automatically the extension ".trace". The default name is attributable to the name of the file storing the Petri Net model, e.g. "PetriNetName-analysis.trace". Similarly, the FBA problem solutions output filename boils down to "<#>.flux" (where "#" is the identifier number of the FBA problem). Also, this file contains the value for the objective function.

### Simulation Results
We used the dynamical model to represent the direct effects of glucose and lactose supplementation on bacterial metabolic output. The ODE-based model couples to the FBA via the concentration of the two carbon sources, capturing the variation in the sugar content of the batch reactor over time. The resulting unified model predicted the carbon source consumption given different administration scenarios. The prediction of sugar levels shows a variation in exchange reaction activity following the nutrient administration ( [Fig. 2](#results) ).

FBA resolution estimates biomass flux ( i.e. the growth rate of the organism) under the given conditions. A biomass flux of about 0.65 mmol/gDW*h means that each gram of dry weight of the _E. coli_ strain K-12 MG1655 is producing 0.65 mmol of biomass per hour. This value is a measure of the metabolic activity estimating as the cell's ability to convert the nutrients in the medium into new cell material. In a culture at the mid-log phase, the value might resemble that of the exponential growth phase. After this point, the cells have already utilized a substantial portion of the available nutrients, leading to a deceleration in their growth rate (refer to 1.2 hours). Keep in mind that FBA is a mathematical model that makes certain assumptions about the metabolic network of the organism and the distribution of metabolic fluxes.

### Timing
We measured the simulation execution time from Genome-scale metabolic modelling compiling initiation at gene expression data integration to the workflow of analysis ( CPU; Intel i7-3520M (4) @ 3.600GHz ).

* time compiling phase: ~ 10 minutes
* time analysis phase: ~ 5 minutes

## References
[1] (https://github.com/qBioTurin/epimod)

[2] [Monk, Jonathan & Lloyd, Colton & Brunk, Elizabeth & Mih, Nathan & Sastry, Anand & King, Zachary & Takeuchi, Rikiya & Nomura, Wataru & Zhang, Zhen & Mori, Hirotada & Feist, Adam & Palsson, Bernhard. (2017). iML1515, a knowledgebase that computes Escherichia coli traits. Nature Biotechnology. 35. 904-908](https://www.nature.com/articles/nbt.3956)

[3] [Liu, Mingzhu & Durfee, Tim & Cabrera, Julio & Zhao, Kai & Jin, Ding & Blattner, Frederick. (2005). Global Transcriptional Programs Reveal a Carbon Source Foraging Strategy by. The Journal of biological chemistry. 280. 15921-7. 10.1074/jbc.M414050200](https://www.sciencedirect.com/science/article/pii/S0021925820692699?via%3Dihub)

[4] [Amparore, E.G., Balbo, G., Beccuti, M., Donatelli, S., Franceschinis, G. 30 Years of GreatSPN.  In: Fiondella, L., Puliafito, A. (eds) Principles of Performance and Reliability Modeling and Evaluation. Springer Series in Reliability Engineering. Springer, Cham. 2016](https://link.springer.com/chapter/10.1007/978-3-319-30599-8_9)

[5] [Gu, Changdai & Kim, Gi & Kim, Won Jun & Kim, Hyun & Lee, Sang Yup. Current status and applications of genome-scale metabolic models. Genome Biology. 20. 10.1186/s13059-019-1730-3. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1730-3)

[6] [Kashaf, Sara & Angione, Claudio & Lio, Pietro. (2017). Making life difficult for Clostridium difficile: Augmenting the pathogen's metabolic model with transcriptomic and codon usage data for better therapeutic target characterization. BMC Systems Biology](https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0395-3)

[7] [Liu, Mingzhu & Durfee, Tim & Cabrera, Julio & Zhao, Kai & Jin, Ding & Blattner, Frederick. (2005). Global Transcriptional Programs Reveal a Carbon Source Foraging Strategy by. The Journal of biological chemistry. 280. 15921-7.](https://www.sciencedirect.com/science/article/pii/S0021925820692699?via%3Dihub)

[8] (https://www.genome.wisc.edu/resources/k12growth/mopscurve.htm)

[9] (https://mediadb.systemsbiology.net/defined_media/media/102/)

[10] (https://ecmdb.ca/e_coli_stats)
