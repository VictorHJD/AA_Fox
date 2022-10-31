library(phyloseq)
#library(tidyverse)
library(dplyr)
library(tibble)
library(Hmsc)
library(abind)
library(MCMCvis)
library(corrplot)
library(reshape2)
library(ggplot2)
library(patchwork)
library(ggcorrplot)
library(tidyr)
library(forcats)

## # JSDM model for fox parasites

## Explanatory variables:
## - impervious surface (%): buffer 1000m around fox. <- IS LEFT OUT
## AS HIGHLY CORRELATED WITH human footprint index
## - tree cover (%): buffer 1000m around fox.
## - human footprint index: buffer 1000m around fox
## (higher means more human influence)

## - sex fox: male, female
## - weight fox (kg)
### We need to use weight as it has a BIG influence on diversity

## - Species richness in diet
## - Species richness in FunM
## - Species richness in BacM

## Traits of parasites:
## - host range: moderate vs wide
## - zoonotic: yes/no
## - lifecycle: one/two/three hosts

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


table(sample_data(PSG)[ , "season"])

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
    
## keep those at least in five percent of the samples
ntokeep <- nrow(HelmCounts)*0.05
response_data <- HelmCounts[, colSums(HelmCounts>0)>ntokeep] %>%
    as.matrix() 
nrow(response_data) # [1] 140

## 2. a trait data frame with the genus in rows, same name and order
## as in species matrix (columns)

## Helminth traits
traits <- read.csv("input_data/helminth_traits.csv")

traits <- traits %>%
    column_to_rownames("t.genus")  


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
foxes <- foxes %>%
  mutate(coords.x1 = as.numeric(coords.x1)) %>% 
    mutate(coords.x1 = ifelse(duplicated(coords.x1, coords.x2),
                             coords.x1+1, coords.x1)) 
nrow(foxes) # [1] 140

### now we need to know how the environmental data for the foxes is
### correlated
envcov_cor <- foxes %>%
    dplyr::select(weight_kg, 
                  tree_cover_1000m, imperv_1000m, human_fpi_1000m, 
                  Diet_Species_richness, BacM_Species_richness, FunM_Species_richness) %>%
    mutate_all(as.numeric) %>%
    cor(x=., use = "pairwise.complete.obs") 

envcov_cor
ggcorrplot(envcov_cor, type = "lower", lab = TRUE)

### and then create a predictor dataset (without non-predictor
### variables) and we leave out imperv_1000m as it's highly correlated
### with human_fpi_1000m, which is more relevant (for diversity at least)

### We should consider using also diet ("Diet_Species_richness") and
### maybe bacterial ("BacM_Species_richness") and fungal
### ("FunM_Species_richness") microbiome species richness as
### predictors!

envcov_data <- foxes %>%
    dplyr::select(IZW_ID, area, sex, age, weight_kg, season, area, human_fpi_1000m, tree_cover_1000m,
                  Diet_Species_richness, BacM_Species_richness, FunM_Species_richness)  %>%
    mutate_at(c("IZW_ID", "area", "sex", "age", "season", "area"), as.factor) %>%
    mutate_at(c("weight_kg", "human_fpi_1000m", "tree_cover_1000m",
                "Diet_Species_richness", "BacM_Species_richness", "FunM_Species_richness"), as.numeric) %>% 
  filter(!is.na(season))
nrow(envcov_data) # [1] 140

#### now the coordinates (in the coordinate system Cedric used) for
#### the random structure
xyData <- foxes %>%
  transmute(x.coord = coords.x1, y.coord = coords.x2) %>% 
  mutate(x.coord = as.numeric(x.coord), 
         y.coord = as.numeric(y.coord))
nrow(xyData)  # [1] 140

## Maybe remove the spatial random effects for computational efficiency?!

## FOR NOW also removing the foxes from the same sites here
response_data <- response_data[rownames(envcov_data), ]
nrow(response_data)

studyDesign <- data.frame(site = rownames(envcov_data))
studyDesign[] <- as.data.frame(lapply(studyDesign[], factor))

rL <- HmscRandomLevel(sData = xyData)

### Define MCMC parameters
thin <- 10
samples <- 20000
transient <- 1000
nChains <- 3
verbose <- 1000


## Regression formula for environmental covariates
XFormula.area = ~ sex + weight_kg + season + area +
  Diet_Species_richness + BacM_Species_richness + FunM_Species_richness
  
XFormula.grad = ~ sex + weight_kg + season + human_fpi_1000m + tree_cover_1000m +
  Diet_Species_richness + BacM_Species_richness + FunM_Species_richness


# Regression formula for traits
TrFormula.Genera = ~ zoonotic + lifecycle + host.range

## *BINOMIAL DISTRIBUTION* ~> PROBIT MODEL
## Fit models for PRESENCE/ABSENCE  data 

## area model
PAModel_fitarea <- Hmsc(Y = response_data>0, XData = envcov_data, XFormula = XFormula.area,
                studyDesign=studyDesign, ranLevels=list(site=rL),
                TrFormula = TrFormula.Genera, TrData = traits,
                distr = "probit")

## trial to see if model runs
PAModel_area <- sampleMcmc(PAModel_fitarea, thin = 5, samples = 20, verbose=TRUE)

# the real model
PAModel_area <- sampleMcmc(PAModel_fitarea, thin = thin, samples = samples, transient = transient, 
                      nChains = nChains, verbose = verbose, nParallel = nChains)

## save this as it takes very long to compute!
## we can't put it in the repos as it's to big
saveRDS(PAModel_area, "./JSDM_models/PAModel_area_jSDM.rds")


## gradient model
PAModel_fitgrad <- Hmsc(Y = response_data>0, XData = envcov_data, XFormula = XFormula.grad,
                        studyDesign=studyDesign, ranLevels=list(site=rL),
                        TrFormula = TrFormula.Genera, TrData = traits,
                        distr = "probit")

PAModel_grad <- sampleMcmc(PAModel_fitgrad, thin = 5, samples = 20, verbose=TRUE)

# the real model
PAModel_grad <- sampleMcmc(PAModel_fitgrad, thin = thin, samples = samples, transient = transient, 
                           nChains = nChains, verbose = verbose, nParallel = nChains)

## save this as it takes very long to compute!
## we can't put it in the repos as it's to big
saveRDS(PAModel_grad, "./JSDM_models/PAModel_grad_jSDM.rds")

## quick look at the models
PAModel_area
# Hmsc object with 140 sampling units, 11 species, 9 covariates, 5 traits and 1 random levels
PAModel_area$XFormula

PAModel_grad
# Hmsc object with 140 sampling units, 11 species, 10 covariates, 5 traits and 1 random levels
PAModel_grad$XFormula

