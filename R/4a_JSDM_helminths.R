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

#### now the coordinates (in the coordinate system Cedric used) for
#### the random structure
xyData <- foxes %>%
  transmute(x.coord = coords.x1, y.coord = coords.x2) %>% 
  mutate(x.coord = as.numeric(x.coord), 
         y.coord = as.numeric(y.coord))
     

## Maybe remove the spatial random effects for computational efficiency?!

## FOR NOW also removing the foxes from the same sites here
response_data <- response_data[rownames(envcov_data), ]


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



#####################################
############ check model ############
#####################################

## Model convergence 

## We evaluate MCMC convergence in terms of two kinds of parameters that
## we are especially interested in: the species niches Beta, influence of
## traits on species niches Gamma, and the residual species associations
## Omega.  The strength of phylogenetic signal rho was not included in
## this model

## Evaluate convergence: Effective sample size and gelman-rubin
## diagnostic (potencial reduction factor)


## get everything in one nice table
getConvergenceStats <- function (Mpost) {
  ## get effective size 
  cl <- list(beta = cbind(effectiveSize(Mpost$Beta),
                          gelman.diag(Mpost$Beta,
                                      multivariate = FALSE)$psrf,
                          "beta"), 
             gamma = cbind(effectiveSize(Mpost$Gamma),
                           gelman.diag(Mpost$Gamma,
                                       multivariate = FALSE)$psrf,
                           "gamma"),
             omega = cbind(effectiveSize(Mpost$Omega[[1]]),
                           gelman.diag(Mpost$Omega[[1]],
                                       multivariate = FALSE)$psrf,
                           "omega"))
  ## name the columns
  ncl <- lapply(cl, function(x) {
    df <- as.data.frame(x)
    colnames(df) <- c("ESS", "GELMAN.est", "GELMAN.CI", "variable")
    df
  })
  Reduce(rbind, ncl)
}

###########################
## Model convergence

## area model
PAMpost_area <- convertToCodaObject(PAModel_area)
PAConv_area <- getConvergenceStats(PAMpost_area)

PAConv_area[, c("ESS", "GELMAN.est", "GELMAN.CI")] <-
  apply(PAConv_area[, c("ESS", "GELMAN.est", "GELMAN.CI")], 2, as.numeric)

PAConv_area %>% 
  group_by(variable) %>% 
  summarise(mean.gd = mean(GELMAN.est))
# variable mean.gd
# <chr>      <dbl>
#   1 beta        1.01
# 2 gamma       1.00
# 3 omega       1.03

ESSPAhist <- ggplot(PAConv_area, aes(ESS)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")


GELMPAhist <- ggplot(PAConv_area, aes(GELMAN.est)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## Graphical output (which I assume is a supplementary file)
png("JSDM_models/figures_PA/jSDM_PA_area_convergence.png", width = 800, height = 1000,
    pointsize = 20)
ESSPAhist / GELMPAhist
dev.off()

MCMCtrace(PAMpost_area$Beta, 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_beta_area",
          wd= "JSDM_models/figures_PA/")

MCMCtrace(PAMpost_area$Gamma, 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_gamma_area",
          wd = "JSDM_models/figures_PA/")

MCMCtrace(PAMpost_area$Omega[[1]], 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_omega_area",
          wd = "JSDM_models/figures_PA/")


## gradient model
PAMpost_grad <- convertToCodaObject(PAModel_grad)
PAConv_grad <- getConvergenceStats(PAMpost_grad)

PAConv_grad[, c("ESS", "GELMAN.est", "GELMAN.CI")] <-
  apply(PAConv_grad[, c("ESS", "GELMAN.est", "GELMAN.CI")], 2, as.numeric)

PAConv_grad %>% 
  group_by(variable) %>% 
  summarise(mean.gd = mean(GELMAN.est))
# variable mean.gd
# <chr>      <dbl>
#   1 beta        1.00
# 2 gamma       1.00
# 3 omega       1.01

ESSPAhist <- ggplot(PAConv_grad, aes(ESS)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")


GELMPAhist <- ggplot(PAConv_grad, aes(GELMAN.est)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## Graphical output (which I assume is a supplementary file)
png("JSDM_models/figures_PA/jSDM_PA_grad_convergence.png", width = 800, height = 1000,
    pointsize = 20)
ESSPAhist / GELMPAhist
dev.off()

MCMCtrace(PAMpost_grad$Beta, 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_beta_grad",
          wd= "JSDM_models/figures_PA/")

MCMCtrace(PAMpost_grad$Gamma, 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_gamma_grad",
          wd = "JSDM_models/figures_PA/")

MCMCtrace(PAMpost_grad$Omega[[1]], 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_omega_grad",
          wd = "JSDM_models/figures_PA/")



###########################
### Model predictions

## area model

PApreds_area <- computePredictedValues(PAModel_area, expected = TRUE)
saveRDS(PApreds_area, "./JSDM_models/PAModel_jSDM_area_preds.rds")

## Median of the predictions
# PApreds_area_values <- apply(abind(PApreds_area, along=3), c(1,2), median)
# Mean of the predictions
# PApreds_areavalues_mean <- apply(abind(PApreds_area, along = 3), c (1,2), mean)


# R2 with the built in function
modelr2.explanatory <- evaluateModelFit(hM = PAModel_area, predY = PApreds_area)
# explanatory power of the model
mean(modelr2.explanatory$TjurR2, na.rm=TRUE) # [1] 0.2673145
png("./JSDM_models/figures_PA/PAmodel_area_hist_r2.png")
hist(modelr2.explanatory$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(modelr2.explanatory$TjurR2, na.rm = TRUE),2)))
dev.off()

# AUC of the model
mean(modelr2.explanatory$AUC, na.rm=TRUE) # [1] 0.8701379


## gradient model
PApreds_grad <- computePredictedValues(PAModel_grad, expected = TRUE)
saveRDS(PApreds_grad, "./JSDM_models/PAModel_jSDM_grad_preds.rds")

# R2 with the built in function
modelr2.explanatory <- evaluateModelFit(hM = PAModel_grad, predY = PApreds_grad)
mean(modelr2.explanatory$TjurR2, na.rm=TRUE) # [1] 0.2801554
png("./JSDM_models/figures_PA/PAmodel_grad_hist_r2.png")
hist(modelr2.explanatory$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(modelr2.explanatory$TjurR2, na.rm = TRUE),2)))
dev.off()

# AUC of the model
mean(modelr2.explanatory$AUC, na.rm=TRUE) # [1] 0.8732902


######################
#### PLOT RESULTS ####
######################


###################
## beta values

## Area model
Beta_area <- as.data.frame(MCMCsummary(PAMpost_area$Beta))
postBeta_area <- getPostEstimate(PAModel_area, parName = "Beta")

png("./JSDM_models/figures_PA/Betaplot_area_default_support95.png")
plotBeta(PAModel_area, post = postBeta_area, param = "Support", supportLevel = 0.95)
dev.off()

png("./JSDM_models/figures_PA/Betaplot_area_default_mean95.png")
plotBeta(PAModel_area, post = postBeta_area, param = "Mean", supportLevel = 0.95)
dev.off()

# Print a plot for each predictor
n_cov <- length(PAModel_area$covNames) # Number of covariates without the intercept
var_code <- vector()
for (i in 1:n_cov){
  var_code[i] <- paste0("C", i)
}

var_name <- as.vector(PAModel_area$covNames[1:n_cov])
predictors <- as.data.frame(cbind(var_code, var_name))

for (i in 1:nrow(predictors)){
  png(paste0("./JSDM_models/figures_PA/betas_area_covariates_coef_plot_", 
             var_name[i], ".png"), width = 5, 
      height = 8, units = "in", res = 300, pointsize = 16)
  MCMCplot(PAMpost$Beta,
           params = predictors[i,1],
           ISB = FALSE,
           exact = FALSE,
           ref_ovl = TRUE,
           rank = FALSE,
           xlab = 'ESTIMATE',
           main = predictors[i,2],
           sz_labels = 0.5,
           sz_med = 1,
           sz_thick = 1,
           sz_thin = 1,
           sz_ax = 1,
           sz_main_txt = 1)
  dev.off()
}


## Gradient model
Beta_grad <- as.data.frame(MCMCsummary(PAMpost_grad$Beta))
postBeta_grad <- getPostEstimate(PAModel_grad, parName = "Beta")

png("./JSDM_models/figures_PA/Betaplot_grad_default_support95.png")
plotBeta(PAModel_grad, post = postBeta_grad, param = "Support", supportLevel = 0.95)
dev.off()

png("./JSDM_models/figures_PA/Betaplot_grad_default_mean95.png")
plotBeta(PAModel_grad, post = postBeta_grad, param = "Mean", supportLevel = 0.95)
dev.off()

# Print a plot for each predictor
n_cov <- length(PAModel_grad$covNames) # Number of covariates without the intercept
var_code <- vector()
for (i in 1:n_cov){
  var_code[i] <- paste0("C", i)
}

var_name <- as.vector(PAModel_grad$covNames[1:n_cov])
predictors <- as.data.frame(cbind(var_code, var_name))

for (i in 1:nrow(predictors)){
  png(paste0("./JSDM_models/figures_PA/betas_grad_covariates_coef_plot_", 
             var_name[i], ".png"), width = 5, 
      height = 8, units = "in", res = 300, pointsize = 16)
  MCMCplot(PAMpost$Beta,
           params = predictors[i,1],
           ISB = FALSE,
           exact = FALSE,
           ref_ovl = TRUE,
           rank = FALSE,
           xlab = 'ESTIMATE',
           main = predictors[i,2],
           sz_labels = 0.5,
           sz_med = 1,
           sz_thick = 1,
           sz_thin = 1,
           sz_ax = 1,
           sz_main_txt = 1)
  dev.off()
}


##################
## manual plotting of beta values


## area model

# get species in model
m1_species <- colnames(PAModel_area$Y)

# get names for explanatory variables
exp_variables <- colnames(PAModel_area$X)
exp_variables
# rename exp variables
my_variables <- c("(Intercept)", "sex[male]", "weight_kg", "season[spring]", "season[winter]",
                  "area[Brandenburg]", "Diet_Species_richness", "BacM_Species_richness", 
                  "FunM_Species_richness")

ModelFrame_area <- data.frame()
# betas (coefficients) for each species
for (i in 1:length(m1_species)){
  # get the variables for each species
  mpost_beta_tmp <- PAMpost_area$Beta[,grep(m1_species[i], colnames(PAMpost_area$Beta[[1]]))]
  # rename variables
  for (j in 1:length(mpost_beta_tmp)){
    colnames(mpost_beta_tmp[[j]]) <- my_variables
  }
  # Put model estimates into temporary data.frames. Add variable "Species" for plotting
  modelFrame_tmp <- data.frame(Variable = my_variables,
                               Coefficient = summary(mpost_beta_tmp)$statistics[,1],
                               CI_low = summary(mpost_beta_tmp)$quantiles[,1],
                               Q_25 = summary(mpost_beta_tmp)$quantiles[, 2],
                               Q_50 = summary(mpost_beta_tmp)$quantiles[,3],
                               Q_75 = summary(mpost_beta_tmp)$quantiles[, 4],
                               CI_high = summary(mpost_beta_tmp)$quantiles[,5],
                               Species = m1_species[i]) 
  
  # Combine these data.frames
  ModelFrame_area <- data.frame(rbind(ModelFrame_area, modelFrame_tmp))
}

ModelFrame_area

levels(as.factor(ModelFrame_area$Variable))
# Relevel factors so they are plotted in the desired order
ModelFrame_area <- ModelFrame_area %>% 
  mutate(Variable = as.factor(Variable)) %>% 
  mutate(Variable = fct_relevel(Variable, c("(Intercept)", "sex[male]", "weight_kg", "season[spring]", "season[winter]",
                                            "area[Brandenburg]", "Diet_Species_richness", "BacM_Species_richness", 
                                            "FunM_Species_richness"))) %>% 
  mutate(Variable = fct_rev(Variable))
summary(ModelFrame_area)

write.csv(ModelFrame_area, "./JSDM_models/ModelFrame_PAModel_area.csv", row.names = FALSE)

# variables with CRI not overlapping 0
toplot_ModelFrame_area <- ModelFrame_area %>%
  mutate(significant = case_when( 
    CI_low < 0 & CI_high < 0 ~ "Yes", #both extremes of CI are negative
    CI_low > 0 & CI_high > 0 ~ "Yes", #both extremes of CI are positive
    TRUE ~ "No")) 

toplot_ModelFrame_area[toplot_ModelFrame_area$significant == "Yes",]

# Plot Effects
plot_1 <- toplot_ModelFrame_area %>% 
  filter(Variable %in% c("sex[male]","weight_kg")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high, fill = significant),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75, fill = significant),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  # scale_y_continuous(limits = c(-5, 3)) +
  # coord_flip(ylim=c(-2, 1)) + 
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold")) +
  ggtitle("Helminth presence - Area model")
plot_1

ggsave(plot = plot_1, "./JSDM_models/figures_PA/PAModel_area_BetaCoefs_plot1.png", 
      width = 9, height = 8, dpi = 600)

plot_2 <- toplot_ModelFrame_area %>% 
  filter(Variable %in% c("Diet_Species_richness", "BacM_Species_richness",
                         "FunM_Species_richness")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high, fill = significant),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75, fill = significant),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  # scale_y_continuous(limits = c(-5, 3)) +
  # coord_flip(ylim=c(-2, 1)) + 
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold")) +
  ggtitle("Helminth presence - Area model")
plot_2

ggsave(plot = plot_2, "./JSDM_models/figures_PA/PAModel_area_BetaCoefs_plot2.png", 
       width = 9, height = 8, dpi = 600)

plot_3 <- toplot_ModelFrame_area %>% 
  filter(Variable %in% c("area[Brandenburg]", "season[spring]", "season[winter]")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high, fill = significant),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75, fill = significant),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold")) +
  ggtitle("Helminth presence - Area model")
plot_3

ggsave(plot = plot_3, "./JSDM_models/figures_PA/PAModel_area_BetaCoefs_plot3.png", 
      width = 9, height = 8, dpi = 600)


## gradient model

# get species in model
m1_species <- colnames(PAModel_grad$Y)

# get names for explanatory variables
exp_variables <- colnames(PAModel_grad$X)
exp_variables
# rename exp variables
my_variables <- c("(Intercept)", "sex[male]", "weight_kg", "season[spring]", "season[winter]",          
                  "human_fpi_1000m", "tree_cover_1000m", 
                  "Diet_Species_richness", "BacM_Species_richness", "FunM_Species_richness")

ModelFrame_grad <- data.frame()
# betas (coefficients) for each species
for (i in 1:length(m1_species)){
  # get the variables for each species
  mpost_beta_tmp <- PAMpost_grad$Beta[,grep(m1_species[i], colnames(PAMpost_grad$Beta[[1]]))]
  # rename variables
  for (j in 1:length(mpost_beta_tmp)){
    colnames(mpost_beta_tmp[[j]]) <- my_variables
  }
  # Put model estimates into temporary data.frames. Add variable "Species" for plotting
  modelFrame_tmp <- data.frame(Variable = my_variables,
                               Coefficient = summary(mpost_beta_tmp)$statistics[,1],
                               CI_low = summary(mpost_beta_tmp)$quantiles[,1],
                               Q_25 = summary(mpost_beta_tmp)$quantiles[, 2],
                               Q_50 = summary(mpost_beta_tmp)$quantiles[,3],
                               Q_75 = summary(mpost_beta_tmp)$quantiles[, 4],
                               CI_high = summary(mpost_beta_tmp)$quantiles[,5],
                               Species = m1_species[i]) 
  
  # Combine these data.frames
  ModelFrame_grad <- data.frame(rbind(ModelFrame_grad, modelFrame_tmp))
}

ModelFrame_grad

levels(as.factor(ModelFrame_grad$Variable))
# Relevel factors so they are plotted in the desired order
ModelFrame_grad <- ModelFrame_grad %>% 
  mutate(Variable = as.factor(Variable)) %>% 
  mutate(Variable = fct_relevel(Variable, c("(Intercept)", "sex[male]", "weight_kg", "season[spring]", "season[winter]",          
                                            "human_fpi_1000m", "tree_cover_1000m", 
                                            "Diet_Species_richness", "BacM_Species_richness", "FunM_Species_richness"))) %>% 
  mutate(Variable = fct_rev(Variable))
summary(ModelFrame_grad)

write.csv(ModelFrame_grad, "./JSDM_models/ModelFrame_PAModel_grad.csv", row.names = FALSE)

# variables with CRI not overlapping 0
toplot_ModelFrame_grad <- ModelFrame_grad %>%
  mutate(significant = case_when( 
    CI_low < 0 & CI_high < 0 ~ "Yes", #both extremes of CI are negative
    CI_low > 0 & CI_high > 0 ~ "Yes", #both extremes of CI are positive
    TRUE ~ "No")) 

toplot_ModelFrame_grad[toplot_ModelFrame_grad$significant == "Yes",]

# Plot Effects
plot_1 <- toplot_ModelFrame_grad %>% 
  filter(Variable %in% c("sex[male]","weight_kg")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high, fill = significant),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75, fill = significant),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  # scale_y_continuous(limits = c(-5, 3)) +
  # coord_flip(ylim=c(-2, 1)) + 
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold")) +
  ggtitle("Helminth presence - Gradient model")
plot_1

ggsave(plot = plot_1, "./JSDM_models/figures_PA/PAModel_grad_BetaCoefs_plot1.png", 
      width = 9, height = 8, dpi = 600)

plot_2 <- toplot_ModelFrame_grad %>% 
  filter(Variable %in% c("Diet_Species_richness", "BacM_Species_richness",
                         "FunM_Species_richness")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high, fill = significant),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75, fill = significant),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  # scale_y_continuous(limits = c(-5, 3)) +
  # coord_flip(ylim=c(-2, 1)) + 
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold")) +
  ggtitle("Helminth presence - Gradient model")
plot_2

ggsave(plot = plot_2, "./JSDM_models/figures_PA/PAModel_grad_BetaCoefs_plot2.png", 
       width = 9, height = 8, dpi = 600)

# Plot Effects
plot_3 <- toplot_ModelFrame_grad %>% 
  filter(Variable %in% c("human_fpi_1000m",  "tree_cover_1000m", "season[spring]", "season[winter]")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high, fill = significant),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75, fill = significant),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold")) +
  ggtitle("Helminth presence - Gradient model")
plot_3

ggsave(plot = plot_3, "./JSDM_models/figures_PA/PAModel_grad_BetaCoefs_plot3.png", 
      width = 9, height = 8, dpi = 600)


##################
## plot traits

## area model
postGamma_area <- getPostEstimate(PAModel_area, parName = "Gamma")

pdf("JSDM_models/figures_PA/jSDM_PAModel_area_GammaEffects_support.pdf", width=7, height=8)
plotGamma(hM = PAModel_area, post = postGamma_area, param = "Support", supportLevel = 0.95,
          newplot = TRUE, covNamesNumbers = c(TRUE, TRUE), mar =c (10,10,1,1))
dev.off()

pdf("JSDM_models/figures_PA/jSDM_PAModel_area_GammaEffects_mean.pdf", width=7, height=8)
plotGamma(hM = PAModel_area, post = postGamma_area, param = "Mean", supportLevel = 0.95,
          newplot = TRUE, covNamesNumbers = c(TRUE, TRUE), mar =c (10,10,1,1))
dev.off()

## gradient model
postGamma_grad <- getPostEstimate(PAModel_grad, parName = "Gamma")

pdf("JSDM_models/figures_PA/jSDM_PAModel_grad_GammaEffects_support.pdf", width=7, height=8)
plotGamma(hM = PAModel_grad, post = postGamma_grad, param = "Support", supportLevel = 0.95,
          newplot = TRUE, covNamesNumbers = c(TRUE, TRUE), mar =c (10,10,1,1))
dev.off()

pdf("JSDM_models/figures_PA/jSDM_PAModel_grad_GammaEffects_mean.pdf", width=7, height=8)
plotGamma(hM = PAModel_grad, post = postGamma_grad, param = "Mean", supportLevel = 0.95,
          newplot = TRUE, covNamesNumbers = c(TRUE, TRUE), mar =c (10,10,1,1))
dev.off()

#PAMpost$Gamma



############################
## Plot sp associations

## area model
OmegaCor_area <- computeAssociations(PAModel_area)
saveRDS(OmegaCor_area, file = "./JSDM_models/PAModel_area_omegaCor.rds")

OmegaCor_area[[1]]$mean
OmegaCor_area[[1]]$support

# Default plot in Hmsc package
supportLevel <- 0.95

toPlot_PAModel_area <- ((OmegaCor_area[[1]]$support > supportLevel)
                        + (OmegaCor_area[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor_area[[1]]$mean

png("./JSDM_models/figures_PA/PAModel_area_Omegaplot_default_support95.png")
par(mfrow = c(1,1))
corrplot(toPlot_PAModel_area, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = paste0("random effect level: ", PAModel_area$rLNames[1]), 
         mar = c(0,0,1,0))
dev.off()

# get associations for manual plotting - replicate
assoc_mean <- melt(OmegaCor_area[[1]]$mean)
assoc_support <- melt(OmegaCor_area[[1]]$support)

associations_area <- cbind.data.frame(assoc_mean, support = assoc_support$value)
colnames(associations_area) <- c("species1", "species2", "mean", "support")
associations_area

write.csv(associations_area, "./JSDM_models/PAModel_area_table_sp_associations.csv",
          row.names = FALSE)


## area model
OmegaCor_grad <- computeAssociations(PAModel_grad)
saveRDS(OmegaCor_grad, file = "./JSDM_models/PAModel_grad_omegaCor.rds")

OmegaCor_grad[[1]]$mean
OmegaCor_grad[[1]]$support

# Default plot in Hmsc package
supportLevel <- 0.95

toPlot_PAModel_grad <- ((OmegaCor_grad[[1]]$support > supportLevel)
                        + (OmegaCor_grad[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor_grad[[1]]$mean

png("./JSDM_models/figures_PA/PAModel_grad_Omegaplot_default_support95.png")
par(mfrow = c(1,1))
corrplot(toPlot_PAModel_grad, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = paste0("random effect level: ", PAModel_grad$rLNames[1]), 
         mar = c(0,0,1,0))
dev.off()

# get associations for manual plotting - replicate
assoc_mean <- melt(OmegaCor_grad[[1]]$mean)
assoc_support <- melt(OmegaCor_grad[[1]]$support)

associations_grad <- cbind.data.frame(assoc_mean, support = assoc_support$value)
colnames(associations_grad) <- c("species1", "species2", "mean", "support")
associations_grad

write.csv(associations_grad, "./JSDM_models/PAModel_grad_table_sp_associations.csv",
          row.names = FALSE)


############################
## Plot Variance partitioning

## area model
head(PAModel_area$X)

VP_area <- computeVariancePartitioning(PAModel_area, group = c(1,1,1, 2,2, 3, 4,4,4), 
                                       groupnames = c("Individual", "Season", "Environment", "Diet_Microbiomes"))
plotVariancePartitioning(PAModel_area, VP_area)
saveRDS(VP_area, "./JSDM_models/PAModel_area_varpart.rds")

# Extract the values for the manual plot
VP_vals_area <- as.data.frame(VP_area$vals) 
VP_vals_area

# get mean variance explained
mean_vp_area <- as.data.frame(rowSums(VP_vals_area)/ncol(VP_vals_area))
colnames(mean_vp_area) <- "mean"
mean_vp_area <- mean_vp_area %>% 
  mutate(percent = round(mean * 100, 2), 
         Variable = factor(rownames(mean_vp_area), 
                           levels = c("Random: site", "Diet_Microbiomes", "Environment", "Season", "Individual"))
  )  
mean_vp_area

# set species names
my_species <- colnames(PAModel_area$Y) 

colnames(VP_vals_area) <- my_species

# Melt the data for plotting and add a column with the variable group 
VP_toplot_area <- VP_vals_area %>% 
  pivot_longer(everything(), names_to = "Species") %>% 
  mutate(Variable = rep(rownames(VP_vals_area), each = length(my_species))) %>% 
  mutate(Variable = factor(Variable, levels = c("Random: site", "Diet_Microbiomes", "Environment", "Season", "Individual")))

head(VP_toplot_area)
tail(VP_toplot_area)

plot(VP_toplot_area$value ~ VP_toplot_area$Variable)

vp_plot_area <- ggplot(VP_toplot_area, aes(x = Species, y = value, fill = Variable)) +
  geom_bar(stat = 'identity', colour = "grey40", alpha = 0.3) +
  scale_fill_manual(values = alpha(c("lightyellow", "darkorange3", "darkgreen", "darkblue", "firebrick4"), 0.7),
                    name = "Variable group", 
                    labels=c(paste0("Random: site\n(mean = ", 
                                    mean_vp_area$percent[mean_vp_area$Variable == "Random: site"], ")"),
                             paste0("Diet_Microbiomes\n(mean = ", 
                                    mean_vp_area$percent[mean_vp_area$Variable == "Diet_Microbiomes"], ")"), 
                             paste0("Environment\n(mean = ", 
                                    mean_vp_area$percent[mean_vp_area$Variable == "Environment"], ")"),
                             paste0("Season\n(mean = ", 
                                    mean_vp_area$percent[mean_vp_area$Variable == "Season"], ")"),
                             paste0("Individual\n(mean = ", 
                                    mean_vp_area$percent[mean_vp_area$Variable == "Individual"], ")"))) +
  labs(#title = "Variance Partitioning", 
    title = "Response = PA, Environment = area",
    x = "\nHelminths", 
    y = "Variance partitioning (%)\n", col = "black") +
  scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, 
                                   size=10, colour = "black", face = "italic"),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(hjust = 0.5, vjust = 1.5),
        legend.key.height = unit(1.5, "lines"),
  )

vp_plot_area 

ggsave(plot = vp_plot_area, "./JSDM_models/figures_PA/VarPart_PAModel_area.png",  
       dpi = 600, width = 6, height = 5)

## gradient model
head(PAModel_grad$X)

VP_grad <- computeVariancePartitioning(PAModel_grad, group = c(1,1,1, 2,2, 3,3, 4,4,4), 
                                       groupnames = c("Individual", "Season", "Environment", "Diet_Microbiomes"))
plotVariancePartitioning(PAModel_grad, VP_grad)
saveRDS(VP_grad, "./JSDM_models/PAModel_grad_varpart.rds")

# Extract the values for the manual plot
VP_vals_grad <- as.data.frame(VP_grad$vals) 
VP_vals_grad

# get mean variance explained
mean_vp_grad <- as.data.frame(rowSums(VP_vals_grad)/ncol(VP_vals_grad))
colnames(mean_vp_grad) <- "mean"
mean_vp_grad <- mean_vp_grad %>% 
  mutate(percent = round(mean * 100, 2), 
         Variable = factor(rownames(mean_vp_grad), 
                           levels = c("Random: site", "Diet_Microbiomes", "Environment", "Season", "Individual"))
  )  
mean_vp_grad

# set species names
my_species <- colnames(PAModel_grad$Y) 

colnames(VP_vals_grad) <- my_species

# Melt the data for plotting and add a column with the variable group 
VP_toplot_grad <- VP_vals_grad %>% 
  pivot_longer(everything(), names_to = "Species") %>% 
  mutate(Variable = rep(rownames(VP_vals_grad), each = length(my_species))) %>% 
  mutate(Variable = factor(Variable, levels = c("Random: site", "Diet_Microbiomes", "Environment", "Season", "Individual")))

head(VP_toplot_grad)
tail(VP_toplot_grad)

plot(VP_toplot_grad$value ~ VP_toplot_grad$Variable)

vp_plot_grad <- ggplot(VP_toplot_grad, aes(x = Species, y = value, fill = Variable)) +
  geom_bar(stat = 'identity', colour = "grey40", alpha = 0.3) +
  scale_fill_manual(values = alpha(c("lightyellow", "darkorange3", "darkgreen", "darkblue", "firebrick4"), 0.7),
                    name = "Variable group", 
                    labels=c(paste0("Random: site\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable == "Random: site"], ")"),
                             paste0("Diet_Microbiomes\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable == "Diet_Microbiomes"], ")"), 
                             paste0("Environment\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable == "Environment"], ")"), 
                             paste0("Season\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable == "Season"], ")"), 
                             paste0("Individual\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable == "Individual"], ")"))) +
  labs(#title = "Variance Partitioning", 
    title = "Response = PA, Environment = gradient",
    x = "\nHelminths", 
    y = "Variance partitioning (%)\n", col = "black") +
  scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, 
                                   size=10, colour = "black", face = "italic"),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(hjust = 0.5, vjust = 1.5),
        legend.key.height = unit(1.5, "lines"),
  )

vp_plot_grad 

ggsave(plot = vp_plot_grad, "./JSDM_models/figures_PA/VarPart_PAModel_grad.png",  
       dpi = 600, width = 6, height = 5)
