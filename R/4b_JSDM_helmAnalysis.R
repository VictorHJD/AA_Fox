library(ggplot2)
library(patchwork)
library(Hmsc)
library(dplyr)
library(forcats)
library(MCMCvis)
library(ggpubr)
library(tidyr)
library(igraph)
library(ggraph)
library(tibble)
library(reshape2)
library(cowplot)
# read custom theme for plotting
source("./R/plot_setup.R")
extrafont::loadfonts(device = "all") ## run every time

recomputejSDMModels <- FALSE

## to avoid doing now convergence tests (and pdfs of those) every time
## the script is run
newConvergenceTest <- FALSE

if(!exists("PAModel_area") | !exists("PAModel_grad")){
    if(recomputejSDMModels){
        source("R/4a_JSDM_helminths.R")
    } else {
        PAModel_area <- readRDS(file="JSDM_models/PAModel_area_jSDM_DNA.rds")
        PAModel_grad <- readRDS(file="JSDM_models/PAModel_grad_jSDM_DNA.rds")
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

## # A tibble: 3 × 3
##   variable mean.ess mean.gd
##   <chr>       <dbl>   <dbl>
## 1 beta       49303.    1.00
## 2 gamma      57686.    1.00
## 3 omega      24539.    1.01

## visually inspecting convergence statistics

## effective sample size
ESSPAhist <- ggplot(PAConv_area, aes(ESS)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## gelman diagnostic
GELMPAhist <- ggplot(PAConv_area, aes(GELMAN.est)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## ### To inspect graphics 
## ESSPAhist / GELMPAhist

if(newConvergenceTest){
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
}
    
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

## ## A tibble: 3 × 3
##   variable mean.ess mean.gd
##   <chr>       <dbl>   <dbl>
## 1 beta       45996.    1.00
## 2 gamma      57624.    1.00
## 3 omega      17592.    1.01
## visually inspecting convergence statistics

## effective sample size
ESSPAhist <- ggplot(PAConv_grad, aes(ESS)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## gelman diagnostic
GELMPAhist <- ggplot(PAConv_grad, aes(GELMAN.est)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free_y")

## ## To inspect the graphics
## ESSPAhist / GELMPAhist

## chain convergence
if(newConvergenceTest) {
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
}


###########################
### Evaluate model predictive power

## area model

PApreds_area <- computePredictedValues(PAModel_area, expected = TRUE)

# Obtain R2 of area model 
modelr2.explanatory <- evaluateModelFit(hM = PAModel_area, predY = PApreds_area)
# OBtain explanatory power of the model based on R2
mean(modelr2.explanatory$TjurR2, na.rm=TRUE) ## [1] 0.1779459
# Obtain AUC of the model
mean(modelr2.explanatory$AUC, na.rm=TRUE) # [1] 0.8260046


## gradient model
PApreds_grad <- computePredictedValues(PAModel_grad, expected = TRUE)

# Obtain R2 of gradient model 
modelr2.explanatory <- evaluateModelFit(hM = PAModel_grad, predY = PApreds_grad)
# OBtain explanatory power of the model based on R2
mean(modelr2.explanatory$TjurR2, na.rm=TRUE) # [1] 0.1839139
# Obtain AUC of the model
mean(modelr2.explanatory$AUC, na.rm=TRUE) # [1] 0.8326632




############################
#### PLOT MODEL RESULTS ####
############################


###################
## beta values 

## Area model
Beta_area <- as.data.frame(MCMCsummary(PAMpost_area$Beta))
postBeta_area <- getPostEstimate(PAModel_area, parName = "Beta")

## ## quick exploration of variables that are significant 
## plotBeta(PAModel_area, post = postBeta_area, param = "Support", supportLevel = 0.95)

## Forest plot of beta effects 
# get species in model
m1_species <- colnames(PAModel_area$Y)

# get names for explanatory variables
exp_variables <- colnames(PAModel_area$X)
exp_variables
# rename exp variables for nicer plots

my_variables <- c("(Intercept)", "sex[male]", "weight_kg",
                  "season[spring]", "season[winter]",
                  "area[Brandenburg]"
                  )

ModelFrame_area <- data.frame()

## betas (coefficients) for each species
for (i in 1:length(m1_species)){
  # get the variables for each species
  mpost_beta_tmp <- PAMpost_area$Beta[,grep(m1_species[i], colnames(PAMpost_area$Beta[[1]]))]
##    print(head(mpost_beta_tmp))
    ## # rename variables
    for (j in 1:length(mpost_beta_tmp)){
  ##      print(colnames(mpost_beta_tmp[[j]]))
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

## Relevel factors so they are plotted in the desired order
ModelFrame_area <- ModelFrame_area %>% 
  mutate(Variable = as.factor(Variable)) %>% 
    mutate(Variable = fct_relevel(Variable, c("(Intercept)", "sex[male]", "weight_kg",
                                              "season[spring]", "season[winter]",
                                              "area[Brandenburg]"
                                              ))) %>% 
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

# prepare colour scale for plotting
my_palette <- rcartocolor::carto_pal(10, "Prism")
my_palette <- colorspace::desaturate(my_palette, .22)

# Plot Effects
plot_beta_weight <- toplot_ModelFrame_area %>% 
  filter(Variable %in% c("sex[male]","weight_kg")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  coord_flip() +
    scale_colour_manual("Genus",
        values = my_palette, 
        guide = guide_legend(reverse = TRUE)) +
    ylab("") +
    xlab("") +
  theme_custom() +
  theme(    
    plot.margin = margin(1,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12,
                               face = "italic"))

plot_beta_area <- toplot_ModelFrame_area %>% 
  filter(Variable %in% c("area[Brandenburg]", "season[spring]", "season[winter]")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  coord_flip() +
  scale_colour_manual(values = my_palette, 
                      guide = "none") +
  xlab("Variable") +
  ylab("") +
  theme_custom() +
  theme(
    plot.margin = margin(1,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

beta_plot <- plot_grid(plot_beta_area, plot_beta_weight,
          ncol = 2, 
          rel_widths = c(1,1.4))

beta_plot <- ggdraw(add_sub(beta_plot,
                            "Coefficient",
                            vpadding = grid::unit(0, "lines"),
                            y = 0.1,
                            x = 0.5,
                            vjust = 0,
                            size = 14,
                            fontface = "bold"))

ggsave(plot = beta_plot, "./figures/PAModel_area_BetaCoefs.png", 
       width = 12, height = 6, dpi = 600, 
       bg = "white")


### Detail plots for each prediction 
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

## ## quick exploration of variables that are significant 
## plotBeta(PAModel_grad, post = postBeta_grad, param = "Support", supportLevel = 0.95)

# get species in model
m1_species <- colnames(PAModel_grad$Y)

# get names for explanatory variables
exp_variables <- colnames(PAModel_grad$X)
exp_variables
# rename exp variables
my_variables <- c("(Intercept)", "sex[male]", "weight_kg",
                  "season[spring]", "season[winter]",          
                  "human_fpi_1000m", "tree_cover_1000m"
                  )

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

ModelFrame_grad <- ModelFrame_grad %>% 
  mutate(Variable = as.factor(Variable)) %>% 
    mutate(Variable = fct_relevel(Variable, c("(Intercept)", "sex[male]", "weight_kg",
                                              "season[spring]", "season[winter]",          
                                              "human_fpi_1000m", "tree_cover_1000m"
                                              ))) %>% 
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

# define colour palette
my_palette <- rcartocolor::carto_pal(10, "Prism")
my_palette <- colorspace::desaturate(my_palette, .22)

# Plot Effects
plot_beta_grad_1 <- toplot_ModelFrame_grad %>% 
  filter(Variable %in% c("sex[male]","weight_kg")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none") +
  coord_flip() +
    scale_colour_manual(values = my_palette,
                        guide = "none") +
    ylab("") +
    xlab("Variables") +
  theme_custom() +
  theme(
    plot.margin = margin(1,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

plot_beta_grad_2 <- toplot_ModelFrame_grad %>% 
    filter(Variable %in% c("season[spring]", "season[winter]")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none") +
  coord_flip() +
    scale_colour_manual(values = my_palette,
                        guide = "none") +
  xlab("") +
  ylab("Coefficient") +
  theme_custom() +
  theme(
    plot.margin = margin(1,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)) 

plot_beta_grad_3 <- toplot_ModelFrame_grad %>% 
    filter(Variable %in% c("human_fpi_1000m",  "tree_cover_1000m")) %>%
  ggplot(aes(group = Species, colour = Species)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none") +
  coord_flip() +
    scale_colour_manual("Genus",
                        values = my_palette, 
                        guide = guide_legend(reverse = TRUE)) +
  xlab("") +
  ylab("") +
  theme_custom() +
  theme(
    plot.margin = margin(1,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12,
                               face = "italic")) 


beta_plot_grad <- plot_grid(plot_beta_grad_1, plot_beta_grad_2, plot_beta_grad_3, 
                       ncol = 3, 
                       nrow = 1,
                       rel_widths = c(0.7, 0.7, 1.3))


ggsave(plot = beta_plot_grad, "./figures/suppl/PAModel_grad_BetaCoefs.png", 
       width = 12, height = 6, dpi = 600, 
       bg = "white")


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




##################
## plot traits

## area model
postGamma_area <- getPostEstimate(PAModel_area, parName = "Gamma")

## ## quick exploration of variables that are significant?!
## plotGamma(PAModel_area, post = postGamma_area, param = "Support", supportLevel = 0.95)

## Forest plot of gamma effects 
# get traits in model
m1_trait_area <- colnames(PAModel_area$Tr)[-1] # remove intercept in traits to avoid future problems

# get names for explanatory variables
exp_variables_area <- colnames(PAModel_area$X)

# rename exp variables for nicer plots
my_variables_area <- c("(Intercept)", "sex[male]", "weight_kg",
                       "season[spring]", "season[winter]",
                       "area[Brandenburg]")

ModelFrame_Traits_area <- data.frame()
# betas (coefficients) for each species
for (i in 1:length(m1_trait_area)){
  # get the variables for each trait
  mpost_gamma_tmp <- PAMpost_area$Gamma[,grep(m1_trait_area[i], colnames(PAMpost_area$Gamma[[1]]))]
  # rename variables
  for (j in 1:length(mpost_gamma_tmp)){
    colnames(mpost_gamma_tmp[[j]]) <- my_variables_area
  }
  # Put model estimates into temporary data.frames. Add variable "Species" for plotting
  modelFrame_tmp <- data.frame(Variable = my_variables_area,
                               Coefficient = summary(mpost_gamma_tmp)$statistics[,1],
                               CI_low = summary(mpost_gamma_tmp)$quantiles[,1],
                               Q_25 = summary(mpost_gamma_tmp)$quantiles[, 2],
                               Q_50 = summary(mpost_gamma_tmp)$quantiles[,3],
                               Q_75 = summary(mpost_gamma_tmp)$quantiles[, 4],
                               CI_high = summary(mpost_gamma_tmp)$quantiles[,5],
                               Traits = m1_trait_area[i]) 
    # Combine these data.frames
  ModelFrame_Traits_area <- data.frame(rbind(ModelFrame_Traits_area, modelFrame_tmp))
}

ModelFrame_Traits_area

# Relevel factors so they are plotted in the desired order
ModelFrame_Traits_area <- ModelFrame_Traits_area %>% 
  mutate(Variable = as.factor(Variable)) %>% 
    mutate(Variable = fct_relevel(Variable, c("(Intercept)", "sex[male]", "weight_kg",
                                              "season[spring]", "season[winter]",
                                              "area[Brandenburg]"))) %>% 
  mutate(Variable = fct_rev(Variable)) %>% 
  mutate(traits_legend = case_when(
    Traits == "zoonoticYes" ~ "Zoonotic",
    Traits == "lifecycletwo.host" ~ "LC - two hosts",
    Traits == "lifecyclethree.host" ~ "LC - three hosts",
    Traits == "host.rangewide" ~ "Wide host range"))
summary(ModelFrame_Traits_area)

# variables with CRI not overlapping 0
toplot_ModelFrame_Traits_area <- ModelFrame_Traits_area %>%
  mutate(significant = case_when( 
    CI_low < 0 & CI_high < 0 ~ "Yes", #both extremes of CI are negative
    CI_low > 0 & CI_high > 0 ~ "Yes", #both extremes of CI are positive
    TRUE ~ "No")) 

toplot_ModelFrame_Traits_area[toplot_ModelFrame_Traits_area$significant == "Yes",]

# define colour palette
my_palette <- rcartocolor::carto_pal(4, "Antique")

# Plot Effects
plot_traits_area <- toplot_ModelFrame_Traits_area %>% 
    filter(Variable %in% c("sex[male]","weight_kg", "season[spring]",
                           "season[winter]", "area[Brandenburg]")) %>%
  ggplot(aes(group = traits_legend, colour = traits_legend)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none") +
  coord_flip() +
  scale_colour_manual(values = my_palette, 
                      name = "Traits",
                      guide = guide_legend(reverse = TRUE)) +
  xlab("") +
  ylab("Coefficient") +
  theme_custom() +
  theme(
    plot.margin = margin(1,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

ggsave(plot = plot_traits_area, "./figures/PAModel_area_GammaCoefs_traits.png", 
       width = 9, height = 5, dpi = 600, 
       bg = "white")


## gradient model
postGamma_grad <- getPostEstimate(PAModel_grad, parName = "Gamma")

## ## quick exploration of variables that are significant 
## plotGamma(PAModel_grad, post = postGamma_grad, param = "Support", supportLevel = 0.95)

## Forest plot of gamma effects 
# get traits in model
m1_trait_grad <- colnames(PAModel_grad$Tr)[-1] # remove intercept in traits to avoid future problems

# get names for explanatory variables
exp_variables_grad <- colnames(PAModel_grad$X)

# rename exp variables for nicer plots
my_variables_grad <- c("(Intercept)", "sex[male]", "weight_kg",
                       "season[spring]", "season[winter]",
                       "human_fpi_1000m", "tree_cover_1000m")

### ULTRA ugly, this variable is created and overwritten so often here!

ModelFrame_Traits_grad <- data.frame()
# betas (coefficients) for each species
for (i in 1:length(m1_trait_grad)){
  # get the variables for each trait
  mpost_gamma_tmp <- PAMpost_grad$Gamma[,grep(m1_trait_grad[i], colnames(PAMpost_grad$Gamma[[1]]))]
  # rename variables
  for (j in 1:length(mpost_gamma_tmp)){
    colnames(mpost_gamma_tmp[[j]]) <- my_variables_grad
  }
  # Put model estimates into temporary data.frames. Add variable "Species" for plotting
  modelFrame_tmp <- data.frame(Variable = my_variables_grad,
                               Coefficient = summary(mpost_gamma_tmp)$statistics[,1],
                               CI_low = summary(mpost_gamma_tmp)$quantiles[,1],
                               Q_25 = summary(mpost_gamma_tmp)$quantiles[, 2],
                               Q_50 = summary(mpost_gamma_tmp)$quantiles[,3],
                               Q_75 = summary(mpost_gamma_tmp)$quantiles[, 4],
                               CI_high = summary(mpost_gamma_tmp)$quantiles[,5],
                               Traits = m1_trait_grad[i]) 
  # Combine these data.frames
  ModelFrame_Traits_grad <- data.frame(rbind(ModelFrame_Traits_grad, modelFrame_tmp))
}

ModelFrame_Traits_grad

# Relevel factors so they are plotted in the desired order
ModelFrame_Traits_grad <- ModelFrame_Traits_grad %>% 
  mutate(Variable = as.factor(Variable)) %>% 
    mutate(Variable = fct_relevel(Variable, c("(Intercept)", "sex[male]", "weight_kg",
                                              "season[spring]", "season[winter]",
                                              "human_fpi_1000m", "tree_cover_1000m"))) %>% 
  mutate(Variable = fct_rev(Variable)) %>% 
  mutate(traits_legend = case_when(
    Traits == "zoonoticYes" ~ "Zoonotic",
    Traits == "lifecycletwo.host" ~ "LC - two hosts",
    Traits == "lifecyclethree.host" ~ "LC - three hosts",
    Traits == "host.rangewide" ~ "Wide host range"))
summary(ModelFrame_Traits_grad)

# variables with CRI not overlapping 0
toplot_ModelFrame_Traits_grad <- ModelFrame_Traits_grad %>%
  mutate(significant = case_when( 
    CI_low < 0 & CI_high < 0 ~ "Yes", #both extremes of CI are negative
    CI_low > 0 & CI_high > 0 ~ "Yes", #both extremes of CI are positive
    TRUE ~ "No")) 

toplot_ModelFrame_Traits_grad[toplot_ModelFrame_Traits_grad$significant == "Yes",]

# define colour palette
my_palette <- rcartocolor::carto_pal(4, "Antique")

# Plot Effects
plot_traits_grad1 <- toplot_ModelFrame_Traits_grad %>% 
    filter(Variable %in% c("sex[male]","weight_kg",
                           "season[spring]", "season[winter]")) %>%
  ggplot(aes(group = traits_legend, colour = traits_legend)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  coord_flip() +
    scale_colour_manual(
        values = my_palette, 
        guide = "none") +
  xlab("") +
  ylab("") +
  theme_custom() +
  theme(
    plot.margin = margin(1,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

plot_traits_grad2 <- toplot_ModelFrame_Traits_grad %>% 
  filter(Variable %in% c("human_fpi_1000m", "tree_cover_1000m")) %>%
  ggplot(aes(group = Traits, colour = Traits)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
  geom_linerange(aes(x = Variable, ymin = CI_low,
                     ymax = CI_high),
                 lwd = 0.8, position = position_dodge(width = 1.5/2)) + 
  geom_linerange(aes(x = Variable, ymin = Q_25,
                     ymax = Q_75),
                 lwd = 1.5, position = position_dodge(width = 1.5/2)) + 
  geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Q_25,
                      ymax = Q_75, fill = significant),
                  lwd = 1/2, shape = 21, position = position_dodge(width = 1.5/2)) +
  scale_fill_manual(values = c("White", "black"), 
                    guide = "none")+
  coord_flip() +
  scale_colour_manual(values = my_palette, 
                      guide = guide_legend(reverse = TRUE)) +
  xlab("") +
  ylab("") +
  theme_custom() +
  theme(
    plot.margin = margin(2,1,3,1),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(), 
    axis.line = element_line(colour = "black"), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)) 


traits_plot_grad <- plot_grid(plot_traits_grad1, plot_traits_grad2,
                              rel_widths = c(0.8, 1.4))

traits_plot_grad <- ggdraw(add_sub(traits_plot_grad,
                                   "Coefficient",
                                   vpadding = grid::unit(0, "lines"),
                                   y = 0.1,
                                   x = 0.5,
                                   vjust = 0,
                                   size = 14,
                                   fontface = "bold"))


ggsave(plot = traits_plot_grad, "./figures/suppl/PAModel_grad_GammaCoefs_traits.png", 
       width = 9, height = 5, dpi = 600, 
       bg = "white")


############################
## Plot sp associations

## area model
OmegaCor_area <- computeAssociations(PAModel_area)

# get associations for manual plotting - replicate
assoc_mean <- reshape2::melt(OmegaCor_area[[1]]$mean)
assoc_support <- reshape2::melt(OmegaCor_area[[1]]$support)

associations_area <- cbind.data.frame(assoc_mean, support = assoc_support$value)

colnames(associations_area) <- c("species1", "species2", "mean", "support")

## Select association that are significant in the 95% CI -->> NOTHING

associations_area %>%
    filter(support < 0.025 | support > 0.975) %>%
    filter(species1 != species2)

## use the 75% confidence interval for visualisation
associations_area %>%
    filter(support < 0.25 | support > 0.75)%>%
    filter(species1 != species2) %>%
    dplyr::rename(from = species1,
                  to = species2) ->
    associations_area75

## ## plot as hierarchical bundle
edges <- data.frame(from = "origin", to = associations_area$species1)
vertices <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))

## calculate the ANGLE of the labels
vertices$id <- NA
myleaves <- which(is.na(match(vertices$name, edges$from))) # Select the rows with the final species, not the groups
nleaves <- length(myleaves) # This should be the number of species
vertices$id[myleaves] <- seq(1:nleaves)
vertices$angle <- 110 - 360*vertices$id/nleaves
# text adjustment
vertices$hjust<-ifelse( vertices$angle < -90, 0, 1)
vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

## Create a graph object
mygraph <- graph_from_data_frame(d = edges, vertices = vertices)

# The connection object must refer to the ids of the leaves:
from <- match(associations_area75$from, vertices$name)
to <-  match(associations_area75$to, vertices$name)

graph_area <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to, values = associations_area75$mean), 
                   alpha=0.5, width=1.8, tension = 0.8, aes(colour=values), 
                   show.legend = TRUE )+
  scale_edge_colour_distiller(palette = "Blues",
                              direction = +1,
                              guide = "edge_colourbar",
                              limits = c(0.65, 0.999),
                              name = "Correlation\n",
                              ) +
    geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label = name,
                       angle = angle, hjust=hjust), size=4, alpha=1) +
    geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, size=1, alpha=0.5),
                    colour = "grey30") +
  expand_limits(x = c(-2, 2), y = c(-2, 2)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        plot.background=element_rect(fill = "white"),
        legend.key = element_rect(color = "gray", fill = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black")) +
  guides(alpha = "none", 
         colour = "none",
         size = "none")

## NOT clearly assigned yet! -> NOT SIGNIFICANT ANYMORE!!!
ggsave(plot = graph_area, "./figures/suppl/PAModel_area_sp_assoc.png", 
       width = 6, height = 6, dpi = 600)


## Gradient model
OmegaCor_grad <- computeAssociations(PAModel_grad)

# get associations for manual plotting - replicate
assoc_mean <- reshape2::melt(OmegaCor_grad[[1]]$mean)
assoc_support <- reshape2::melt(OmegaCor_grad[[1]]$support)

associations_grad <- cbind.data.frame(assoc_mean, support = assoc_support$value)

colnames(associations_grad) <- c("species1", "species2", "mean", "support")

# Select association that are significant in the 95% CI
associations_grad %>%
  filter(support < 0.025 | support > 0.975) %>%
  filter(species1 != species2)

## ## empty! Nothing to select from anymore!!
##   select(from = species1, to = species2, mean = mean, support = support) 

head(associations_grad)
## head(associations_grad95)

## plot as hierarchical bundle
edges <- data.frame(from = "origin", to = associations_grad$species1)
vertices <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))

## calculate the ANGLE of the labels
vertices$id <- NA
myleaves <- which(is.na(match(vertices$name, edges$from))) # Select the rows with the final species, not the groups
nleaves <- length(myleaves) # This should be the number of species
vertices$id[myleaves] <- seq(1:nleaves)
vertices$angle <- 110 - 360*vertices$id/nleaves
# text adjustment
vertices$hjust<-ifelse( vertices$angle < -90, 0, 1)
vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

## Create a graph object
mygraph <- graph_from_data_frame(d = edges, vertices = vertices)

## ## This does not exit anymore1
## ##The connection object must refer to the ids of the leaves:
## from <-  match(associations_grad95$from, vertices$name)
## to <-  match(associations_grad95$to, vertices$name)

## graph_grad <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
##   geom_conn_bundle(data = get_con(from = from, to = to, values = associations_grad95$mean), 
##                    alpha=0.5, width=1.8, tension = 0.8, aes(colour=values), 
##                    show.legend = TRUE )+
##   scale_edge_colour_distiller(palette = "Blues",
##                               direction = +1,
##                               guide = "edge_colourbar",
##                               limits = c(0.85, 0.999),
##                               name = "Correlation\n",
##   ) +
##   geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label = name, angle = angle, hjust=hjust), size=4, alpha=1) +
##   geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, size=1, alpha=0.5), colour = "grey30") +
##   expand_limits(x = c(-2, 2), y = c(-2, 2)) +
##   theme_minimal() +
##   theme(panel.grid.major = element_blank(),
##         panel.grid.minor = element_blank(), 
##         axis.title = element_blank(),
##         axis.text = element_blank(), 
##         plot.background=element_rect(fill = "white"),
##         legend.key = element_rect(color = "gray", fill = "black"),
##         legend.title = element_text(color = "black"),
##         legend.text = element_text(color = "black")) +
##   guides(alpha = "none", 
##          colour = "none",
##          size = "none")

## ## not clearly assigned yet  -> NOT SIGNIFICANT anymore!!!
## ggsave(plot = graph_grad, "./JSDM_models/figures_PA/PAModel_grad_sp_assoc.png", 
##        width = 6, height = 6, dpi = 600)


############################
## Plot Variance partitioning

## area model
head(PAModel_area$X)

## No grouping necessary anymore!
VP_area <- computeVariancePartitioning(PAModel_area)

## ## first view?!
## plotVariancePartitioning(PAModel_area, VP_area)

# Extract the values for the manual plot
VP_vals_area <- as.data.frame(VP_area$vals) 
VP_vals_area

# get mean variance explained
mean_vp_area <- as.data.frame(rowSums(VP_vals_area)/ncol(VP_vals_area))
colnames(mean_vp_area) <- "mean"
mean_vp_area <- mean_vp_area %>% 
    mutate(percent = round(mean * 100, 2), 
           Variable = factor(rownames(mean_vp_area),
                             levels = c("area", "weight_kg", "season",
                                        "sex", "Random: site")))

# set species names
my_species <- colnames(PAModel_area$Y) 

colnames(VP_vals_area) <- my_species

# Melt the data for plotting and add a column with the variable group 
VP_toplot_area <- VP_vals_area %>% 
  pivot_longer(everything(), names_to = "Species") %>% 
  mutate(Variable = rep(rownames(VP_vals_area), each = length(my_species))) %>% 
    mutate(Variable = factor(Variable,
                             levels = c("area", "weight_kg", "season",
                                        "sex", "Random: site")))

## get the species in descending order of host-intrinsic values
species_order <- VP_toplot_area %>% 
  filter(Variable == "area") %>% 
  arrange(desc(value))
  
## order and add colours to plot
VP_toplot_area2 <- VP_toplot_area %>% 
  mutate(Species_ord = factor(Species, levels = species_order$Species)) 

## add colour to the species based on one-host or multi-host trait
multihost_sp <- species_order %>% 
  left_join(PAModel_area$Tr %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "Species"), by = "Species") %>% 
    as.data.frame() %>% 
    mutate(multihost = case_when(
               lifecyclethree.host == 1 ~ "Three",
               lifecycletwo.host == 1 ~ "Two",
               TRUE ~ "One"
           )) %>% 
    mutate(colour_text = case_when(
               multihost == "Three" ~ "orange",
               multihost == "Two" ~ "darkgreen",
               multihost == "One" ~ "grey40"
           )) %>% 
    dplyr::select(Species, colour_text)

## load theme for the plot
source("./R/plot_setup.R")

## plot
vp_plot_area <- ggplot(VP_toplot_area2, aes(x = Species_ord, y = value, fill = Variable)) +
    geom_bar(stat = 'identity', colour = "grey40", alpha = 0.3) +
    scale_fill_manual(values = alpha(c("lightyellow", "darkgreen", "darkorange3",
                                       "firebrick4", "darkblue"), 0.7),
                      name = "Variable", 
                      labels=c(paste0("Admin. area\n(mean = ", 
                                      mean_vp_area$percent[mean_vp_area$Variable ==
                                                           "area"], ")"),
                               paste0("Weight\n(mean = ", 
                                      mean_vp_area$percent[mean_vp_area$Variable ==
                                                           "weight_kg"], ")"),
                               paste0("Season\n(mean = ", 
                                      mean_vp_area$percent[mean_vp_area$Variable ==
                                                           "season"], ")"),
                               paste0("sex\n(mean = ", 
                                      mean_vp_area$percent[mean_vp_area$Variable ==
                                                           "sex"], ")"),
                               paste0("Random: site\n(mean = ", 
                                      mean_vp_area$percent[mean_vp_area$Variable ==
                                                           "Random: site"], ")")
                               )) +
    labs(x = "\nHelminth genus", 
         y = "Variance partitioning (%)\n", col = "black") +
    scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)) +
    theme_custom() +
    theme(panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, 
                                     size=12, colour = multihost_sp$colour_text,
                                     face = "italic"),
          plot.margin = ggplot2::margin(4,1, 1, 1),
          legend.text = element_text(margin = margin(t = 5, b = 5), size = 12),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
          axis.title.y = ggplot2::element_text(angle = 90, margin = ggplot2::margin(r = 6)), 
          axis.text.y = element_text(colour = "black", size = 12)
          )

ggsave(plot = vp_plot_area, "./figures/PAModel_area_varpart.png",  
       dpi = 600, width = 6, height = 5,
       bg = "white")





## gradient model
head(PAModel_grad$X)

VP_grad <- computeVariancePartitioning(PAModel_grad)

# Extract the values for the manual plot
VP_vals_grad <- as.data.frame(VP_grad$vals) 
VP_vals_grad

# get mean variance explained
mean_vp_grad <- as.data.frame(rowSums(VP_vals_grad)/ncol(VP_vals_grad))
colnames(mean_vp_grad) <- "mean"
mean_vp_grad <- mean_vp_grad %>% 
  mutate(percent = round(mean * 100, 2), 
         Variable = factor(rownames(mean_vp_grad),
                           levels = c("human_fpi_1000m", "tree_cover_1000m",
                                      "weight_kg", "season",
                                      "sex", "Random: site")))

mean_vp_grad

# set species names
my_species <- colnames(PAModel_grad$Y) 

colnames(VP_vals_grad) <- my_species

# Melt the data for plotting and add a column with the variable group 
VP_toplot_grad <- VP_vals_grad %>% 
  pivot_longer(everything(), names_to = "Species") %>% 
  mutate(Variable = rep(rownames(VP_vals_grad), each = length(my_species))) %>% 
    mutate(Variable = factor(Variable,
                             levels = c("human_fpi_1000m", "tree_cover_1000m",
                                        "weight_kg", "season",
                                        "sex", "Random: site")))
                             

## get the species in descending order of host-intrinsic values
species_order_grad <- VP_toplot_grad %>% 
  filter(Variable == "human_fpi_1000m") %>% 
  arrange(desc(value)) %>% 
  ungroup() 


## order and add colours to plot
VP_toplot_grad2 <- VP_toplot_grad %>% 
  mutate(Species_ord = factor(Species, levels = species_order_grad$Species)) 

## add colour to the species based on one-host or multi-host trait
multihost_sp_grad <- species_order_grad %>% 
  left_join(PAModel_grad$Tr %>% 
              as.data.frame() %>% 
              rownames_to_column(var = "Species"), by = "Species") %>% 
  as.data.frame() %>% 
  mutate(multihost = case_when(
    lifecyclethree.host == 1 ~ "Three",
    lifecycletwo.host == 1 ~ "Two",
    TRUE ~ "One"
  )) %>% 
  mutate(colour_text = case_when(
    multihost == "Three" ~ "orange",
    multihost == "Two" ~ "darkgreen",
    multihost == "One" ~ "grey40"
  )) %>% 
  dplyr::select(Species, colour_text)


vp_plot_grad <- ggplot(VP_toplot_grad2, aes(x = Species_ord, y = value, fill = Variable)) +
    geom_bar(stat = 'identity', colour = "grey40", alpha = 0.3) +
    scale_fill_manual(values = alpha(c("lightyellow", "lightblue", "darkgreen",
                                       "darkorange3","firebrick4", "darkblue"), 0.7),
                    name = "Variable", 
                    labels=c(paste0("Human FPI\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable ==
                                                         "human_fpi_1000m"], ")"),
                             paste0("Tree cover\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable ==
                                                         "tree_cover_1000m"], ")"),
                             paste0("Weight\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable ==
                                                         "weight_kg"], ")"),
                             paste0("Season\n(mean = ", 
                                    mean_vp_area$percent[mean_vp_grad$Variable ==
                                                         "season"], ")"),
                             paste0("Sex\n(mean = ", 
                                    mean_vp_area$percent[mean_vp_grad$Variable ==
                                                         "sex"], ")"),
                             paste0("Random: site\n(mean = ", 
                                    mean_vp_grad$percent[mean_vp_grad$Variable ==
                                                         "Random: site"], ")"))) +
  labs(x = "\nHelminth taxa",
       y = "Variance partitioning (%)\n", col = "black") +
  scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)) +
  theme_custom() +
  theme(panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, 
                                   size=12, colour = multihost_sp$colour_text, face = "italic"),
        plot.margin = ggplot2::margin(4,1, 1, 1),
        legend.text = element_text(margin = margin(t = 5, b = 5), size = 12),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(angle = 90, margin = ggplot2::margin(r = 6)), 
        axis.text.y = element_text(colour = "black", size = 12)
  )

ggsave(plot = vp_plot_grad, "./figures/suppl/VarPart_PAModel_grad.png",  
       dpi = 600, width = 6, height = 5,
       bg = "white")

