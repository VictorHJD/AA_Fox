# AA_Fox

## Setup 
For a fully reproducible analysis starting from raw sequencing files R
script have to be executed in numerical order (starting with a "\\d_";
e.g. 0_, 1_, ..).

We follow one simple convention for executing the scripts: Your R
session should run in the repository folder directly (in the folder
that you cloned, not in 'input_data/', '/intermediate_data' or 'R/').

At the beginning of each script we test for the availability of the
required data (in a potentially interactive R session) and either
recompute it executing ('sourcing') the previous scripts (recompute =
TRUE) or read it from 'intermediate_data'. 

All neccesary intermediate data to run the script are contained within
this repository. You can run each script seperately to scrutinze the
methods and code an to reproduce the results. 

## Matching of files from this repository with references in the manuscript

| File reference | Created in                      | File in ropository                                        |
|:---------------|---------------------------------|:----------------------------------------------------------|
| Figure 1       | outside, manually               | In this readme (ADD!)                                     |
| Figure 2       | R/0_Extract_Einvir_Covariates.R | figures/map_study_overview_multi.png                      |
| suppl. table 1 | outside, manually               | input_data/primer_file_foxes.csv                          |
| Table 1        | outside, manually               | text only                                                 |
| Suppl. text F1 | R/1_Fox_general_MA.R            | figures/suppl/AmpSampleHeatmap.png                        |
| Suppl. text T1 | R/1_Fox_general_MA.R            | tables/suppl/quant.html                                   |
| Suppl. text T2 | R/1_Fox_general_MA.R            | tables/suppl/qual1.html                                   |
| Suppl. text T3 | R/1_Fox_general_MA.R            | tables/suppl/qual2.html                                   |
| suppl. table 2 | outside, manually               | input_data/helminth_traits.csv                            |
| Figure 3       | R/2_iNEXT_fox.R                 | figures/Div_Model.png                                     |
| Suppl. text F2 | R/2_iNEXT_fox.R                 | figures/suppl/DiversityHelminth.png                       |
| Table 2        | R/2_iNEXT_fox.R                 | tables/HelmDiversityArea.html                             |
| Suppl. text T4 | R/2_iNEXT_fox.R                 | tables/HelmDiversityConti.html                            |
| Table 3        | R/3_compositionHelmDiet.R       | tables/Permanova.csv                                      |
| suppl.table 4  | R/3_compositionHelmDiet.R       | tables/PermanovaConti.csv                                 |
| Figure 4       | R/3_compositionHelmDiet.R       | figures/CompositionEnvHelm.png                            |
| Table 4        | R/3_compositionHelmDiet.R       | tables/suppl/EnvFitnMDS.csv                               |
| Table 5        | R/4_abu_models.R                | tables/IndHelmAbu.html                                    |
| Figure 5       | R/5b_JSDM_helmAnalysis.R        | figures/PAModel_area_varpart.png                          |
| Figure 6       | R/5b_JSDM_helmAnalysis.R        | figures/PAModel_area_BetaCoefs.png                        |
| Figure 7       | R/5b_JSDM_helmAnalysis.R        | JSDM_models/figures_PA/PAModel_area_GammaCoefs_traits.png |
| Figure S2      | R/5b_JSDM_helmAnalysis.R        | JSDM_models/figures_PA/PAModel_area_sp_assoc.png          |
| Figure S3      | R/5b_JSDM_helmAnalysis.R        | /JSDM_models/figures_PA/VarPart_PAModel_grad.png          |
| Figure S4      | R/5b_JSDM_helmAnalysis.R        | JSDM_models/figures_PA/PAModel_grad_BetaCoefs.png         |
| Figure S5      | R/5b_JSDM_helmAnalysis.R        | JSDM_models/figures_PA/PAModel_grad_GammaCoefs_traits.png |
| Figure S6      | R/5b_JSDM_helmAnalysis.R        | JSDM_models/figures_PA/PAModel_grad_sp_assoc.png          |


## 0) Environmental variables

-> R/0_Extract_Einvir_Covariates.R (by Cedric Scherer and Aimara Planillo)
 
The input raw (layer) files this is based on are in:
input_data/tifs/*.tif

As this contains only impervious surface, tree cover and human
footprint index, unlike data for previous (pre-mid-2021) versions of
the manuscript, now allows us analysis of landscape variables for both
Berlin and Brandenburg.

The script reads data on the sampled foxes from
"input_data/Fox_data.csv', appends environmental variables for each
fox and writes them (together with the 'basic data') to
'intermediate_data/Fox_data_envir.RDS'.

 
## 1) Sequencing data, sample data and helminith trait data

-> R/1_Fox_general_MA.R (by Victor Victor Jarquin-Diaz and Emanuel Heitlinger)

The script reads the raw sequencing data (and stors and reads
intermediate files) from the compute server of the Heitlinger group
(harriet@biologie.hu-berlin.de) and is in the present form not
executable on other systems. In oder to execute it you'll need to
download the data from NCBI-SRA [BioProject
PRJNA386767](https://www.ncbi.nlm.nih.gov/sra/PRJNA386767), as
sequecing data is to large for storage on github. Then adapt the path
for reading the data in the script accoringly.

We process the raw sequencing data based on matching of the primer
sequences in 'input_data/primer_file_foxes.csv'. We are using the
MultiAmplicon wrapper of the dada2 package to produce amplified
sequence variant (ASV) abundances for each fox.

The script also adds the environmenal covariates (from
intermediate_data/Fox_data_envir.RDS produced in
0_Extract_Einvir_Covariates.R) to the central "phyloseq" object of the
pipeline. The environmental covariates are stored as "sample_data" in
this object.

Similarly, the helminth trait data (found in
input_data/helminth_traits.csv; compiled by Carolin Scholz and Emanuel
Heitlinger) is storted in the "phloseq" object's "taxon_table"). 
 
We store output as a phyloseq object in
'intermediate_data/PhyloSeqCombi.Rds'

## 2) Diversity analysis


-> R/2_iNEXT_fox.R (by Carolin Scholz, Aimara Planillo, Cedric Scherer
and Emanuel Heitlinger)

We look at three differen masures for alpha diverstiy: species
richness [q=0], Shannon diversity [q=1] and Simpson diversity [q=2]),
based on the iNext package. We compare beta diversity (jaccard
centroid distances) using the package vegan and gamma diversity
(rarefied diversity by fox individual) again using iNext.

Diversity data for each fox are appended to the "sample_data" in the
phyloseq object.

## 3) Comunity composition analyses

-> 3_CommunityComposition.R (by Emanuel Heitlinger)

Performs community analyses using the vegan package


## 4) Generalized linear models for helminth abundance

-> R/4_abu_models.R (by Emanuel Heitlinger)

Runs generalized linear models on the abundance of individual helminth species

## 5) jSDM analysis 


-> 5a_JSDM_helminths.R (by Aimara Planillo)
Runs the JSDM models 

-> 5b_JSDM_helmAnalysis.R
Evaluates the models and extract posterior effects and plot in 

 
