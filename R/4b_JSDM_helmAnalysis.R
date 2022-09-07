library(ggplot2)
library(patchwork)

recomputejSDMModels <- FALSE

if(!exists("COModel") | !exists("PAModel")){
    if(recomputejSDMModels & recomputejSDMModels){
        source("R/3_JSDM_helminths.R")
    } else {
        PAModel <- readRDS(file="/SAN/Metabarcoding/AA_Fox/PAModel_jSDM.rds")
##      COModel <- readRDS(file="intermediate_data/COModel_jSDM.rds")
    }
}

## set theme for plots
theme_set(theme_minimal(base_family = "Roboto", base_size = 12))
theme_update(
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    strip.text = element_text(face = "bold", color = "black", size = 12, margin = margin(b = 10)),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.spacing.x = unit(2, "lines"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(rep(12, 4))
)

## font for numeric label
font_num <- "Roboto Condensed"



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

PAMpost <- convertToCodaObject(PAModel)
##COMpost <- convertToCodaObject(COModel)

PAConv <- getConvergenceStats(PAMpost)
## COConv <- getConvergenceStats(COMpost)

PAConv[, c("ESS", "GELMAN.est", "GELMAN.CI")] <-
    apply(PAConv[, c("ESS", "GELMAN.est", "GELMAN.CI")], 2, as.numeric)

## COConv[, c("ESS", "GELMAN.est", "GELMAN.CI")] <-
##     apply(COConv[, c("ESS", "GELMAN.est", "GELMAN.CI")], 2, as.numeric)


ESShist <- ggplot(PAConv, aes(ESS)) +
    geom_histogram() +
    facet_wrap(~variable, scales = "free_y")


GELMhist <- ggplot(PAConv, aes(GELMAN.est)) +
    geom_histogram() +
    facet_wrap(~variable, scales = "free_y")


## Graphical output (which I assume is a supplementary file)
png("figures/suppl/jSDM_PA_convergence.png", width = 800, height = 1000,
    pointsize = 20)
ESShist / GELMhist
dev.off()


## ESSCOhist <- ggplot(COConv, aes(ESS)) +
##     geom_histogram() +
##     facet_wrap(~variable, scales = "free_y")


## GELMCOhist <- ggplot(COConv, aes(GELMAN.est)) +
##     geom_histogram() +
##     facet_wrap(~variable, scales = "free_y")



## ## Graphical output (which I assume is a supplementary file)
## png("figures/suppl/jSDM_CO_convergence.png", width = 800, height = 1000,
##     pointsize = 20)
## ESSCOhist / GELMCOhist
## dev.off()


MCMCtrace(PAMpost$Beta, 
          pdf = FALSE,
          plot = TRUE,
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_beta",
          wd= "figures/suppl/"
          )

MCMCtrace(PAMpost$Gamma, 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_gamma",
          wd= "figures/suppl/")


MCMCtrace(PAMpost$Omega[[1]], 
          pdf = TRUE, o
          open_pdf = FALSE,
          filename = "jSDM_PA_MCMCtrace_omega",
          wd = "figures/suppl/")

## MCMCtrace(COMpost$Beta, 
##           pdf = TRUE, 
##           open_pdf = FALSE,
##           filename = "jSDM_CO_MCMCtrace_beta",
##           wd= "figures/suppl/")

## MCMCtrace(COMpost$Gamma, 
##           pdf = TRUE, 
##           open_pdf = FALSE,
##           filename = "jSDM_CO_MCMCtrace_gamma",
##           wd = "figures/suppl/")

## MCMCtrace(COMpost$Omega[[1]], 
##           pdf = TRUE, 
##           open_pdf = FALSE,
##           filename = "jSDM_CO_MCMCtrace_omega",
##           wd = "figures/suppl/")

### -> No proper convergence for the Count (CO) models will continue
### -> only with the presence/absence (PA models) for now!!!


PApreds <- computePredictedValues(PAModel, expected = TRUE)

## Median of the predictions
PApreds.values <- apply(abind(PApreds, along=3), c(1,2), median)
# Mean of the predictions
PApreds.values.mean <- apply(abind(PApreds, along = 3), c (1,2), mean)


# R2 with the built in function
modelr2.explanatory <- evaluateModelFit(hM = PAModel, predY = PApreds)
modelr2.explanatory

# AUC of the model
mean(modelr2.explanatory$AUC, na.rm=TRUE)
## 0.8826433 [was 0.8524357 in previous versions)



postGamma <- getPostEstimate(PAModel, parName = "Gamma")


pdf("figures/suppl/jSDM_PA_GammaEffects.pdf", width=10, height=20)
## ## don't understand this strange "heatmap"
##plotGamma(hM = PAModel, post = postGamma, param = "Support", supportLevel = 0.95)
MCMCplot(PAMpost$Gamma, ref_ovl = TRUE)
dev.off()





pdf("figures/suppl/jSDM_PA_BetaEffects.pdf", width=10, height=20)
## ## don't understand this strange "heatmap"
MCMCplot(PAMpost$Beta, ref_ovl = TRUE)
dev.off()


## MCMCplot(PAMpost$Omega, ref_ovl = TRUE)

## MCMCplot(PAMpost$Beta, ref_ovl = TRUE)






