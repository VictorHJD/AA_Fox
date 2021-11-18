# AA_Fox

## Setup 
For a fully reproducible analysis starting from raw sequencing files R
script have to be executed in numerical order (starting with a "\\d_";
e.g. 0_, 1_, ..), scripts prefixed with X_ are currently not used.

We follow one simple convention for executing the scripts: Your R
session should run in the repository folder directly (in the folder
that you cloned, not in 'input_data/', 'intermediate_data' or 'R/').

At the beginning of each script we test for the availability of the
required data (in a potentially interactive R session) and either
recompute it executing ('sourcing') the previous scripts (recompute =
TRUE) or read it from 'intermediate_data'. 

For coauthors, here is how this was migrated from the Dropbox folder
and how the data flows through the pipeline and analyses:

## 0) Environmental variables

Dropbox/Project_Canid_Metabarcoding/5_scripts/Extract_EnvirCovariates_BE_BB_20200903.Rmd
(by Cedric Scherer)
 
-> 0_Extract_Einvir_Covariates.R
 
The input raw (layer) files this is based on are in:
input_data/tifs/*.tif

This unlike data for previous (pre-mid-2021) versions of the
manuscript now allows us analysis of landscape variables for both
Berlin and Brandenburg. 

The script reads data on the sampled foxes from
"input_data/Fox_data.csv', appends environmental variables (x, y, v)
for each fox and writes them (together with the 'basic data') to
'intermediate_data/Fox_data_envir.RDS'.
 
## 1) Sequencing data 
We (Victor Jarquin-Diaz and Emanuel Heitlinger) process the raw
sequencing data in 1_Fox_general_MA.R based on matching of the primer
sequences in 'input_data/primer_file_foxes.csv'. We are using the
MultiAmplicon wrapper of the dada2 package to produce amplified
sequence variant (ASV) abundances for each fox.

The sequencing data is most easiely available on the compute server of
the Heitlinger group (harriet@biologie.hu-berlin.de). The script could
alternatively be run on the sequencing data after download from
NCBI-SRA/ or ENA (ADD LINK HERE) (for full reproducibilty; sequencing
data is to large for storage on github).

The script also adds the environmenal covariates (from
intermediate_data/Fox_data_envir.RDS produced in
0_Extract_Einvir_Covariates.R) to the central "phyloseq" object of the
pipeline. The environmental covariates are stored as "sample_data" in
this object.

The script stores intermediate data on our server and is in the
present not executable anywhere else. You can continue with it's
output object. 
 
We store output as a phyloseq object in
'intermediate_data/PhyloSeqCombi.Rds'

## 2) Helminth trait data
Caro Scholz compiled helminth traits in 
Dropbox/Project_Canid_Metabarcoding/6_processed_data/traits_grouped.RDS
 
-> input_data/traits_grouped.csv
 
This manual input data can be edited directly. It still contains many
NAs we might want to add missing information as we progress with
analysis/publication.
 
AT THIS POINT WE HAVE ALL THE DATA collected and processed (in the
phloseq object and the traits table) and the remaining two scripts are
only analysing this.
 
## 3) JSDM analysis 
JSDM_parasites_helminths_AP_20200828.Rmd (by Aimara
Planillo). 
 
-> 3_JSDM_helminths.R

Runs the model, the model check and shows the results
 
## 4) Comunity composition analyses

Dropbox/Project_Canid_Metabarcoding/5_scripts/PCA_CS_20200903_1946 
 
-> 4_CommunityComposition.R
 
