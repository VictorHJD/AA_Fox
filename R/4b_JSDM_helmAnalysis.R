library(ggplot2)
library(patchwork)
library(Hmsc)
library(dplyr)
library(forcats)
library(MCMCvis)
library(ggpubr)
library(tidyr)

recomputejSDMModels <- FALSE

if(!exists("PAModel_area") | !exists("PAModel_grad")){
    if(recomputejSDMModels){
        source("R/4a_JSDM_helminths.R")
    } else {
        PAModel_area <- readRDS(file="JSDM_models/PAModel_area_jSDM.rds")
        PAModel_grad <- readRDS(file="JSDM_models/PAModel_grad_jSDM.rds")
    }
}

########### First look at the models

## summary of models
PAModel_area
PAModel_grad

## check formula
PAModel_area$XFormula
PAModel_grad$XFormula


######################################
############ check models ############
######################################


#######################
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

## area model
PAMpost_area <- convertToCodaObject(PAModel_area)
PAConv_area <- getConvergenceStats(PAMpost_area)

PAConv_area[, c("ESS", "GELMAN.est", "GELMAN.CI")] <-
  apply(PAConv_area[, c("ESS", "GELMAN.est", "GELMAN.CI")], 2, as.numeric)

## numerical output
PAConv_area %>% 
  group_by(variable) %>% 
  summarise(mean.ess = mean(ESS, na.rm = TRUE),
            mean.gd = mean(GELMAN.est, na.rm = TRUE))
# A tibble: 3 x 3
# variable mean.ess mean.gd
# <chr>       <dbl>   <dbl>
# 1 beta       49073.    1.00
# 2 gamma      57321.    1.00
# 3 omega      32808.    1.00

## visually inspecting convergence statistics

## effective sample size
ESSPAhist <- ggplot(PAConv_area, aes(ESS)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## gelman diagnostic
GELMPAhist <- ggplot(PAConv_area, aes(GELMAN.est)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

ESSPAhist / GELMPAhist

## chain convergence
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

## numerical output
PAConv_grad %>% 
  group_by(variable) %>% 
  summarise(mean.ess = mean(ESS, na.rm = TRUE),
            mean.gd = mean(GELMAN.est, na.rm = TRUE))
# A tibble: 3 x 3
# variable mean.ess mean.gd
# <chr>       <dbl>   <dbl>
# 1 beta       47245.    1.00
# 2 gamma      56997.    1.00
# 3 omega      27881.    1.01

## visually inspecting convergence statistics

## effective sample size
ESSPAhist <- ggplot(PAConv_grad, aes(ESS)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## gelman diagnostic
GELMPAhist <- ggplot(PAConv_grad, aes(GELMAN.est)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

ESSPAhist / GELMPAhist

## chain convergence
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
### Evaluate model predictive power

## area model

PApreds_area <- computePredictedValues(PAModel_area, expected = TRUE)

# Obtain R2 of area model 
modelr2.explanatory <- evaluateModelFit(hM = PAModel_area, predY = PApreds_area)
# OBtain explanatory power of the model based on R2
mean(modelr2.explanatory$TjurR2, na.rm=TRUE) # [1] 0.2615152
# Obtain AUC of the model
mean(modelr2.explanatory$AUC, na.rm=TRUE) # [1] 0.8687716


## gradient model
PApreds_grad <- computePredictedValues(PAModel_grad, expected = TRUE)

# Obtain R2 of gradient model 
modelr2.explanatory <- evaluateModelFit(hM = PAModel_grad, predY = PApreds_grad)
# OBtain explanatory power of the model based on R2
mean(modelr2.explanatory$TjurR2, na.rm=TRUE) # [1] 0.2723732
# Obtain AUC of the model
mean(modelr2.explanatory$AUC, na.rm=TRUE) # [1] 0.8745384




############################
#### PLOT MODEL RESULTS ####
############################


###################
## beta values 

## Area model
Beta_area <- as.data.frame(MCMCsummary(PAMpost_area$Beta))
postBeta_area <- getPostEstimate(PAModel_area, parName = "Beta")

## quick exploration of variables that are significant 
plotBeta(PAModel_area, post = postBeta_area, param = "Support", supportLevel = 0.95)

## Forest plot of beta effects 
# get species in model
m1_species <- colnames(PAModel_area$Y)

# get names for explanatory variables
exp_variables <- colnames(PAModel_area$X)
exp_variables
# rename exp variables for nicer plots
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

## Save for maybe later use
# write.csv(ModelFrame_area, "./JSDM_models/ModelFrame_PAModel_area.csv", row.names = FALSE)

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
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  ylab("\n") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

# ggsave(plot = plot_1, "./JSDM_models/figures_PA/PAModel_area_BetaCoefs_plot1.png", 
       # width = 9, height = 8, dpi = 600)

plot_2 <- toplot_ModelFrame_area %>% 
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
  xlab("") +
  ylab("\nCoefficient") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))
 

# ggsave(plot = plot_2, "./JSDM_models/figures_PA/PAModel_area_BetaCoefs_plot2.png", 
#        width = 9, height = 8, dpi = 600)


plot_3 <- toplot_ModelFrame_area %>% 
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
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  xlab("") +
  ylab("\n") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)) #+
# ggtitle("Helminth presence - Area model")

# ggsave(plot = plot_3, "./JSDM_models/figures_PA/PAModel_area_BetaCoefs_plot3.png", 
#        width = 9, height = 8, dpi = 600)

(plot_betas <- ggarrange(plot_1, plot_2, plot_3, 
          ncol = 3, nrow = 1, common.legend = TRUE, legend="right", 
          labels = c("A", "B", "C"), 
          widths = c(1,1.2,1.4)))

ggsave(plot = plot_betas, "./figures/PAModel_area_BetaCoefs.png", 
       width = 12, height = 6, dpi = 600)

# Detail plots for each prediction 
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
  MCMCplot(PAMpost_area$Beta,
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
  MCMCplot(PAMpost_grad$Beta,
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

## Save for maybe later use
# write.csv(ModelFrame_grad, "./JSDM_models/ModelFrame_PAModel_grad.csv", row.names = FALSE)

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
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  ylab("\n") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

plot_2 <- toplot_ModelFrame_grad %>% 
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
  xlab("") +
  ylab("\nCoefficient") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)) 

plot_3 <- toplot_ModelFrame_grad %>% 
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
  coord_flip() +
  scale_colour_viridis_d(option = "viridis", begin = 0, end = 1, 
                         guide = guide_legend(reverse = TRUE)) +
  xlab("") +
  ylab("\n") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)) 


(plot_betas <- ggarrange(plot_1, plot_2, plot_3, 
                         ncol = 3, nrow = 1, common.legend = TRUE, legend="right", 
                         labels = c("A", "B", "C"), 
                         widths = c(1,1.2,1.4)))

ggsave(plot = plot_betas, "./JSDM_models/figures_PA/PAModel_grad_BetaCoefs.png", 
       width = 12, height = 6, dpi = 600)


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



############################
## Plot sp associations

## area model
OmegaCor_area <- computeAssociations(PAModel_area)

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

# write.csv(associations_area, "./JSDM_models/PAModel_area_table_sp_associations.csv",
#           row.names = FALSE)


## area model
OmegaCor_grad <- computeAssociations(PAModel_grad)

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
                                       groupnames = c("Host-intrinsic", "Season", "Natural envir", "Other Microbiomes"))
plotVariancePartitioning(PAModel_area, VP_area)

# Extract the values for the manual plot
VP_vals_area <- as.data.frame(VP_area$vals) 
VP_vals_area

# get mean variance explained
mean_vp_area <- as.data.frame(rowSums(VP_vals_area)/ncol(VP_vals_area))
colnames(mean_vp_area) <- "mean"
mean_vp_area <- mean_vp_area %>% 
  mutate(percent = round(mean * 100, 2), 
         Variable = factor(rownames(mean_vp_area), 
                           levels = c("Random: site", "Other Microbiomes", "Natural envir", "Season", "Host-intrinsic"))
  )  
mean_vp_area

# set species names
my_species <- colnames(PAModel_area$Y) 

colnames(VP_vals_area) <- my_species

# Melt the data for plotting and add a column with the variable group 
VP_toplot_area <- VP_vals_area %>% 
  pivot_longer(everything(), names_to = "Species") %>% 
  mutate(Variable = rep(rownames(VP_vals_area), each = length(my_species))) %>% 
  mutate(Variable = factor(Variable, levels = c("Random: site", "Other Microbiomes", "Natural envir", "Season", "Host-intrinsic")))

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

ggsave(plot = vp_plot_area, "./figures/PAModel_area_varpart.png",  
       dpi = 600, width = 6, height = 5)


## gradient model
head(PAModel_grad$X)

VP_grad <- computeVariancePartitioning(PAModel_grad, group = c(1,1,1, 2,2, 3,3, 4,4,4), 
                                       groupnames = c("Individual", "Season", "Environment", "Diet_Microbiomes"))
plotVariancePartitioning(PAModel_grad, VP_grad)

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


## OLD SCRIPT FROM HERE (theme still needed???)

## ## set theme for plots
## theme_set(theme_minimal(base_family = "Roboto", base_size = 12))
## theme_update(
##     axis.title.x = element_text(margin = margin(t = 12)),
##     axis.title.y = element_text(margin = margin(r = 12)),
##     strip.text = element_text(face = "bold", color = "black", size = 12, margin = margin(b = 10)),
##     legend.title = element_text(size = 12, face = "bold"),
##     legend.text = element_text(size = 12),
##     panel.spacing.x = unit(2, "lines"),
##     panel.grid.minor = element_blank(),
##     plot.margin = margin(rep(12, 4))
## )

## ## font for numeric label
## font_num <- "Roboto Condensed"



## ## Model convergence 

## ## We evaluate MCMC convergence in terms of two kinds of parameters that
## ## we are especially interested in: the species niches Beta, influence of
## ## traits on species niches Gamma, and the residual species associations
## ## Omega.  The strength of phylogenetic signal rho was not included in
## ## this model

## ## Evaluate convergence: Effective sample size and gelman-rubin
## ## diagnostic (potencial reduction factor)


## ## get everything in one nice table
## getConvergenceStats <- function (Mpost) {
##     ## get effective size 
##     cl <- list(beta = cbind(effectiveSize(Mpost$Beta),
##                             gelman.diag(Mpost$Beta,
##                                         multivariate = FALSE)$psrf,
##                             "beta"), 
##                gamma = cbind(effectiveSize(Mpost$Gamma),
##                              gelman.diag(Mpost$Gamma,
##                                          multivariate = FALSE)$psrf,
##                              "gamma"),
##                omega = cbind(effectiveSize(Mpost$Omega[[1]]),
##                              gelman.diag(Mpost$Omega[[1]],
##                                          multivariate = FALSE)$psrf,
##                              "omega"))
##     ## name the columns
##     ncl <- lapply(cl, function(x) {
##         df <- as.data.frame(x)
##         colnames(df) <- c("ESS", "GELMAN.est", "GELMAN.CI", "variable")
##         df
##     })
##     Reduce(rbind, ncl)
## }

## PAMpost <- convertToCodaObject(PAModel)
## ##COMpost <- convertToCodaObject(COModel)

## PAConv <- getConvergenceStats(PAMpost)
## ## COConv <- getConvergenceStats(COMpost)

## PAConv[, c("ESS", "GELMAN.est", "GELMAN.CI")] <-
##     apply(PAConv[, c("ESS", "GELMAN.est", "GELMAN.CI")], 2, as.numeric)

## ## COConv[, c("ESS", "GELMAN.est", "GELMAN.CI")] <-
## ##     apply(COConv[, c("ESS", "GELMAN.est", "GELMAN.CI")], 2, as.numeric)


## ESShist <- ggplot(PAConv, aes(ESS)) +
##     geom_histogram() +
##     facet_wrap(~variable, scales = "free_y")


## GELMhist <- ggplot(PAConv, aes(GELMAN.est)) +
##     geom_histogram() +
##     facet_wrap(~variable, scales = "free_y")


## ## Graphical output (which I assume is a supplementary file)
## png("figures/suppl/jSDM_PA_convergence.png", width = 800, height = 1000,
##     pointsize = 20)
## ESShist / GELMhist
## dev.off()


## ## ESSCOhist <- ggplot(COConv, aes(ESS)) +
## ##     geom_histogram() +
## ##     facet_wrap(~variable, scales = "free_y")


## ## GELMCOhist <- ggplot(COConv, aes(GELMAN.est)) +
## ##     geom_histogram() +
## ##     facet_wrap(~variable, scales = "free_y")



## ## ## Graphical output (which I assume is a supplementary file)
## ## png("figures/suppl/jSDM_CO_convergence.png", width = 800, height = 1000,
## ##     pointsize = 20)
## ## ESSCOhist / GELMCOhist
## ## dev.off()


## MCMCtrace(PAMpost$Beta, 
##           pdf = FALSE,
##           plot = TRUE,
##           open_pdf = FALSE,
##           filename = "jSDM_PA_MCMCtrace_beta",
##           wd= "figures/suppl/"
##           )

## MCMCtrace(PAMpost$Gamma, 
##           pdf = TRUE, 
##           open_pdf = FALSE,
##           filename = "jSDM_PA_MCMCtrace_gamma",
##           wd= "figures/suppl/")


## MCMCtrace(PAMpost$Omega[[1]], 
##           pdf = TRUE, o
##           open_pdf = FALSE,
##           filename = "jSDM_PA_MCMCtrace_omega",
##           wd = "figures/suppl/")

## ## MCMCtrace(COMpost$Beta, 
## ##           pdf = TRUE, 
## ##           open_pdf = FALSE,
## ##           filename = "jSDM_CO_MCMCtrace_beta",
## ##           wd= "figures/suppl/")

## ## MCMCtrace(COMpost$Gamma, 
## ##           pdf = TRUE, 
## ##           open_pdf = FALSE,
## ##           filename = "jSDM_CO_MCMCtrace_gamma",
## ##           wd = "figures/suppl/")

## ## MCMCtrace(COMpost$Omega[[1]], 
## ##           pdf = TRUE, 
## ##           open_pdf = FALSE,
## ##           filename = "jSDM_CO_MCMCtrace_omega",
## ##           wd = "figures/suppl/")

## ### -> No proper convergence for the Count (CO) models will continue
## ### -> only with the presence/absence (PA models) for now!!!


## PApreds <- computePredictedValues(PAModel, expected = TRUE)

## ## Median of the predictions
## PApreds.values <- apply(abind(PApreds, along=3), c(1,2), median)
## # Mean of the predictions
## PApreds.values.mean <- apply(abind(PApreds, along = 3), c (1,2), mean)


## # R2 with the built in function
## modelr2.explanatory <- evaluateModelFit(hM = PAModel, predY = PApreds)
## modelr2.explanatory

## # AUC of the model
## mean(modelr2.explanatory$AUC, na.rm=TRUE)
## ## 0.8826433 [was 0.8524357 in previous versions)



## postGamma <- getPostEstimate(PAModel, parName = "Gamma")


## pdf("figures/suppl/jSDM_PA_GammaEffects.pdf", width=10, height=20)
## ## ## don't understand this strange "heatmap"
## ##plotGamma(hM = PAModel, post = postGamma, param = "Support", supportLevel = 0.95)
## MCMCplot(PAMpost$Gamma, ref_ovl = TRUE)
## dev.off()





## pdf("figures/suppl/jSDM_PA_BetaEffects.pdf", width=10, height=20)
## ## ## don't understand this strange "heatmap"
## MCMCplot(PAMpost$Beta, ref_ovl = TRUE)
## dev.off()


## ## MCMCplot(PAMpost$Omega, ref_ovl = TRUE)

## ## MCMCplot(PAMpost$Beta, ref_ovl = TRUE)






