library(phyloseq)
library(tidyverse)
library(Hmsc)
library(abind)
library(MCMCvis)
library(corrplot)
library(reshape2)

## # JSDM model for fox parasites

## Explanatory variables:
## - impervious surface (%): buffer 1000m around fox.
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
if(!exists("PS")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PS <- readRDS(file="intermediate_data/PhyloSeqCombi.Rds")
    }
}

## We need:
## 1. a matrix of species: sites in rows (foxes) and helminth genera
## in columns

PSHelm <- phyloseq::subset_taxa(PS, phylum%in%c("Nematoda", "Platyhelminthes"))
PSHelmG <- phyloseq::tax_glom(PSHelm, "genus")

## Remove: genera with no observations 
PSHelmG <- phyloseq::prune_taxa(taxa_sums(PSHelmG)>0, PSHelmG)

## Also reomve foxes with with NA for weight or na for environmental
## variables
PSHelmG <- phyloseq::prune_samples(
                         sample_sums(PSHelmG)>0 &
                         !is.na(as.numeric(sample_data(PSHelmG)$weight_kg))&
                         !is.na(sample_data(PSHelmG)$imperv_1000m),
                         PSHelmG)

HelmCounts <- as.data.frame(otu_table(PSHelmG))
colnames(HelmCounts) <- phyloseq::tax_table(PSHelmG)[, "genus"]
    
HelmCounts %>%
    as.matrix()  -> 
    response_data 

## 2. a trait data frame with the genus in rows, same name and order
## as in species matrix (columns)

## Helminth traits
traits <- read.csv("input_data/helminth_traits_.csv")

traits %>%
    column_to_rownames("t.genus")  -> traits 

## making sure the two are aligned
traits <- traits[colnames(response_data),]
rownames(traits) <- colnames(response_data)

traits %>%
  mutate(human.rel = case_when(
    human.related == 0 ~ "No",
    human.related == 1 ~ "Yes",
    TRUE ~ "Unknown")) %>%
  mutate(pet.rel = case_when(
    pet.related == 0 ~ "No",
    pet.related == 1 ~ "Yes",
    TRUE ~ "Unknown")) %>%
  mutate(livestock.rel = case_when(
    livestock.related == 0 ~ "No",
    livestock.related == 1 ~ "Yes",
    TRUE ~ "Unknown")) %>%
  mutate(transmission.3class = case_when(
    transmission == "diet" ~ "diet", 
    is.na(transmission) ~ "Unknown", 
    TRUE ~ "Others")) %>%
    dplyr::select(human.rel, pet.rel, transmission.3class) ->
  traits

all(rownames(traits)==colnames(response_data))

## 3.  a data frame with the environmental covariates for sites: sites
## (foxes) in rows
foxes <- phyloseq::sample_data(PSHelmG)
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
### variables)
foxes %>%
    ## FOR NOW also removing the foxes from the same sites here
    dplyr::select(IZW_ID, sex, weight_kg, imperv_1000m, human_fpi_1000m, 
                  tree_cover_1000m)  %>%
    mutate_at(c("IZW_ID", "sex"), as.factor) %>%
    mutate_at(c("weight_kg", "imperv_1000m",
                "human_fpi_1000m", "tree_cover_1000m"), as.numeric) ->
    envcov_data

#### now the coordinates (in the coordinate system Cedric used) for
#### the random structure
foxes %>%
    transmute(x.coord = coords.x1, y.coord = coords.x2) ->
    xyData 

## FOR NOW also removing the foxes from the same sites here
response_data <- response_data[rownames(envcov_data), ]


studyDesign <- data.frame(site = envcov_data$IZW_ID)
rL <- HmscRandomLevel(sData = xyData)

### Define MCMC parameters
thin <- 10
samples <- 10000
transient <- 1000
nChains <- 3
verbose <- 1000


# Regression formula for environmental covariates
XFormula.Genera = ~ sex + weight_kg + imperv_1000m + human_fpi_1000m + tree_cover_1000m 
#weight not included because NAs

# Regression formula for traits
TrFormula.Genera = ~ human.rel + pet.rel + transmission.3class

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


## Fit models for COUNT data
COModel <- Hmsc(Y = response_data, XData = envcov_data, XFormula = XFormula.Genera,
                studyDesign=studyDesign, ranLevels=list(site=rL),
                TrFormula = TrFormula.Genera, TrData = traits,
                distr = "poisson")
### also try "lognormal poisson"!!??

COModel <- sampleMcmc(COModel, thin = 10, samples = 20, verbose=TRUE)

