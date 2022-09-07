library(phyloseq)
library(tidyverse)
library(Hmsc)
library(abind)
library(MCMCvis)
library(corrplot)
library(reshape2)

## # JSDM model for fox parasites

## Explanatory variables:
## - impervious surface (%): buffer 1000m around fox. <- IS LEAFT OUT
## AS HIGHLY CORRELATED WITH human footprint index
## - tree cover (%): buffer 1000m around fox.
## - human footprint index: buffer 1000m around fox
## (higher means more human influenc)

## - sex fox: male, female
## - weight fox (kg)
### We need to use weight as it has a BIG influence on diversity

## Traits of parasites:
## - transmission type
## - human related: yes/not affects humans
## - pet related: yes/no might be transmitted by pets

## Randon structure:
## - spatial

## Fox and helminth data from the central phyloseq object
recomputeBioinfo <- FALSE
if(!exists("PSG")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSG <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
    }
}

## We need:
## 1. a matrix of species: sites in rows (foxes) and helminth genera
## in columns

PSGHelm <- subset_taxa(PSG, category%in%"Helminth")

## Also reomve foxes with with NA for weight or na for environmental
## variables
PSGHelmR <- phyloseq::prune_samples(
                         sample_sums(PSGHelm)>0 &
                         !is.na(as.numeric(sample_data(PSGHelm)$weight_kg)) &
                         !is.na(sample_data(PSGHelm)$human_fpi_1000m) &
                         !is.na(sample_data(PSGHelm)$tree_cover_1000m), 
                         PSGHelm)

HelmCounts <- as.data.frame(otu_table(PSGHelmR))
colnames(HelmCounts) <- phyloseq::tax_table(PSGHelmR)[, "genus"]
    
HelmCounts[, colSums(HelmCounts)>0] %>%
    as.matrix() ->
    response_data 

## 2. a trait data frame with the genus in rows, same name and order
## as in species matrix (columns)

## Helminth traits
traits <- read.csv("input_data/helminth_traits.csv")

traits %>%
    column_to_rownames("t.genus")  -> traits 


## making sure the two are aligned
traits <- traits[colnames(response_data),]

## make every column a factor
traits[] <- as.data.frame(lapply(traits[], factor))

## 3.  a data frame with the environmental covariates for sites: sites
## (foxes) in rows
foxes <- phyloseq::sample_data(PSGHelmR)
class(foxes) <- "data.frame"

## some foxes are found at the same coordinates but for the
## spatial random structure we use only unique geolocation
foxes %>%
    mutate(coords.x1 =ifelse(duplicated(coords.x1, coords.x2),
                             coords.x1, coords.x1+1)) ->
    foxes


### now we need to know how the environmental data for the foxes is
### correlated
foxes %>%
    dplyr::select(weight_kg, 
                  tree_cover_1000m, imperv_1000m, human_fpi_1000m) %>%
    mutate_all(as.numeric) %>%
    cor(x=., use = "pairwise.complete.obs") -> envcov_cor

### and then create a predictor dataset (without non-predictor
### variables) and we leave out imperv_1000m as it's highly correlated
### with human_fpi_1000m, which is more relevant (for diversity at least)

### We should consider using also diet ("Diet_Species_richness") and
### maybe bacterial ("BacM_Species_richness") and fungal
### ("FunM_Species_richness") microbiome species richness as
### predictors!

foxes %>%
    dplyr::select(IZW_ID, sex, weight_kg, human_fpi_1000m, 
                  tree_cover_1000m)  %>%
    mutate_at(c("IZW_ID", "sex"), as.factor) %>%
    mutate_at(c("weight_kg",
                "human_fpi_1000m", "tree_cover_1000m"), as.numeric) ->
    envcov_data

#### now the coordinates (in the coordinate system Cedric used) for
#### the random structure
foxes %>%
    transmute(x.coord = coords.x1, y.coord = coords.x2) ->
    xyData 

## Maybe remove the spatial random effects for computational efficiency?!

## FOR NOW also removing the foxes from the same sites here
response_data <- response_data[rownames(envcov_data), ]


studyDesign <- data.frame(site = rownames(envcov_data))
studyDesign[] <- as.data.frame(lapply(studyDesign[], factor))

rL <- HmscRandomLevel(sData = xyData)

### Define MCMC parameters
thin <- 10
samples <- 10000
transient <- 1000
nChains <- 3
verbose <- 1000


## Regression formula for environmental covariates
XFormula.Genera = ~ sex + weight_kg + human_fpi_1000m + tree_cover_1000m 
## (previously: weight not included because NAs), now: rather removed
## the NAs as weight is so important for diversity (see

# Regression formula for traits
TrFormula.Genera = ~ zoonotic + transmission.fox + lifecycle + host.range

## *BINOMIAL DISTRIBUTION* ~> PROBIT MODEL
## Fit models for PRESENCE/ABSENCE  data 

PAModel <- Hmsc(Y = response_data>0, XData = envcov_data, XFormula = XFormula.Genera,
                studyDesign=studyDesign, ranLevels=list(site=rL),
                TrFormula = TrFormula.Genera, TrData = traits,
                distr = "probit")

PAModel <- sampleMcmc(PAModel, thin = 10, samples = 20, verbose=TRUE)

# the real model
PAModel <- sampleMcmc(PAModel, thin = thin, samples = samples, transient = transient, 
                      nChains = nChains, verbose = verbose, nParallel = nChains)

## save this as it takes very long to compute!
## we can't put it in the repos as it's to big
saveRDS(PAModel, "/SAN/Metabarcoding/AA_Fox/PAModel_jSDM.rds")

## *POISSON  (or negative binomial?) DISTRIBUTION* ~> POISSON MODEL

## Fit models for COUNT data THIS IS THE PRIORITY
COModel <- Hmsc(Y = response_data, XData = envcov_data, XFormula = XFormula.Genera,
                studyDesign=studyDesign, ranLevels=list(site=rL),
                TrFormula = TrFormula.Genera, TrData = traits,
                distr = "poisson")
### also try "lognormal poisson"!!??

COModel <- sampleMcmc(COModel, thin = 10, samples = 20, verbose=TRUE)

# the real model
COModel <- sampleMcmc(COModel, thin = thin, samples = samples, transient = transient, 
                      nChains = nChains, verbose = verbose, nParallel = nChains)


## save this as it takes very long to compute!
saveRDS(COModel, "/SAN/Metabarcoding/AA_Fox/COModel_jSDM.rds")

