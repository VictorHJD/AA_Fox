library(phyloseq)
#library(tidyverse)
library(dplyr)
library(tibble)
library(Hmsc)
library(abind)
library(MCMCvis)
library(corrplot)
library(reshape2)
library(ggcorrplot)
library(tidyr)
library(forcats)
library(GGally)

source("./R/plot_setup.R")
extrafont::loadfonts(device = "all") ## run every time

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

## Traits of parasites:
## - host range: moderate vs wide
## - zoonotic: yes/no
## - lifecycle: one/two/three hosts

## Randon structure:
## - spatial

## Fox and helminth data from the central phyloseq object
recomputeBioinfo <- FALSE
if(recomputeBioinfo){
    source("R/1_Fox_general_MA.R")
} else {
    PSGHelm <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
}



## We need:
## 1. a matrix of species: sites in rows (foxes) and helminth genera
## in columns


## Also reomve foxes with with NA for weight or na for environmental
## variables
PSGHelmR <- phyloseq::prune_samples(sample_sums(PSGHelm)>0, 
                                    PSGHelm)
## the good old 131 wormy foxes

HelmCounts <- as.data.frame(otu_table(PSGHelmR))
colnames(HelmCounts) <- phyloseq::tax_table(PSGHelmR)[, "genus"]

## ## THIS IS IMPORTANT TO REPORT!!!
## keep those at least in five percent of the samples
ntokeep <- nrow(HelmCounts)*0.05
response_data <- HelmCounts[, colSums(HelmCounts>0)>ntokeep] %>%
    as.matrix() 
nrow(response_data) # [1] 131

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
nrow(foxes) # [1] 131

### now we need to know how the environmental data for the foxes is
### correlated

### THIS IS SLIGHTLY MISPLACED HERE, SHOULD BE BEFORE ANY OTHER
### STATISTICAL ANALYSIS, but I'll leave it here now! 

foxes %>%
    dplyr::select(weight_kg, sex, season,
                  area, tree_cover_1000m, imperv_1000m, human_fpi_1000m) %>%
    mutate_all(as.numeric) %>%
    cbind(foxes$area) %>%
    ggpairs() -> Cor_Plot

ggsave("figures/suppl/CorrelatPedictors.png", Cor_Plot,
       width = 12.5, height = 10.5, units = "in")

### and then create a predictor dataset (without non-predictor
### variables) and we leave out imperv_1000m as it's highly correlated
### with human_fpi_1000m, which is more relevant (for diversity at least)

envcov_data <- foxes %>%
    dplyr::select(IZW_ID, area, sex, weight_kg, season, area, 
                  human_fpi_1000m, tree_cover_1000m)  %>%
    mutate_at(c("IZW_ID", "area", "sex", "season"), as.factor) %>%
    mutate_at(c("weight_kg", "human_fpi_1000m", "tree_cover_1000m"), as.numeric) %>% 
  filter(!is.na(season))
nrow(envcov_data) # [1] 131

#### now the coordinates (in the coordinate system Cedric used) for
#### the random structure
xyData <- foxes %>%
  transmute(x.coord = coords.x1, y.coord = coords.x2) %>% 
  mutate(x.coord = as.numeric(x.coord), 
         y.coord = as.numeric(y.coord))
nrow(xyData)  # [1] 131

## Maybe remove the spatial random effects for computational efficiency?!

## FOR NOW also removing the foxes from the same sites here
response_data <- response_data[rownames(envcov_data), ]
nrow(response_data) # 131

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
XFormula.area = ~ sex + weight_kg + season + area 
  
XFormula.grad = ~ sex + weight_kg + season + human_fpi_1000m + tree_cover_1000m 

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
PAModel_area <- sampleMcmc(PAModel_fitarea, thin = thin, samples = samples,
                           transient = transient, 
                           nChains = nChains, verbose = verbose, nParallel = nChains)

## save this as it takes very long to compute!
## we can't put it in the repos as it's to big
saveRDS(PAModel_area, "./JSDM_models/PAModel_area_jSDM_DNA.rds")


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
saveRDS(PAModel_grad, "./JSDM_models/PAModel_grad_jSDM_DNA.rds")

## quick look at the models
PAModel_area
## Hmsc object with 131 sampling units, 10 species, 6 covariates, 5
## traits and 1 random levels Posterior MCMC sampling with 3 chains
## each with 20000 samples, thin 10 and transient 1000
PAModel_area$XFormula

PAModel_grad
## Hmsc object with 131 sampling units, 10 species, 7 covariates, 5
## traits and 1 random levels Posterior MCMC sampling with 3 chains
## each with 20000 samples, thin 10 and transient 1000
PAModel_grad$XFormula

