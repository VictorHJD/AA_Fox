library(phyloseq)
library(tidyverse)
library(MASS)

recomputeBioinfo <- FALSE

recomputeDiversity <- FALSE

if(!exists("PSG")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSG <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
    }
}


PSGHelm <- subset_taxa(PSG, category%in%"Helminth")

PSGHelm  %>% otu_table() %>%
    as.data.frame -> HelmCounts

colnames(HelmCounts) <- tax_table(PSGHelm)[, "genus"]

Sdat <- as.data.frame(sample_data(PSGHelm))
class(Sdat) <- "data.frame"

Sdat$seq.depth <- rowSums(otu_table(PSG)[rownames(Sdat), ])

## Helminth traits
traits <- as.data.frame(tax_table(PSGHelm))
class(traits) <- "data.frame"
rownames(traits) <- NULL

traits %>%
    column_to_rownames("genus")  -> traits 

HelmCounts.t <- t(HelmCounts)[rownames(traits), ]

by(HelmCounts.t, traits$pet.infecting, sum)

by(HelmCounts.t, traits$fox.parasite, sum)

by(HelmCounts.t, traits$zoonotic, sum)




### Multiple testing betwwen Areas: Prevalence!
Fishing <- apply(HelmCounts.t, 1, function (x) {fisher.test(x>0, Sdat$area)})
Fishing.pval <- lapply(Fishing, "[[", "p.value")

Fishing[p.adjust(Fishing.pval)<0.05]


## Mesocesstoides more prevalent in Brandenburg
table(HelmCounts.t["Mesocestoides", ] >0, Sdat$area)


## Angiostrongylus more prevalent in Berlin
table(HelmCounts.t["Angiostrongylus", ] >0, Sdat$area)

## Alaria more prevalent in Brandenburg
table(HelmCounts.t["Alaria", ] >0, Sdat$area)


### Multiple testing betwwen Areas: Count data (non-parametric)!
FishingC <- apply(HelmCounts.t, 1, function (x) {
    wilcox.test(x[Sdat$area%in%"Brandenburg"],
                x[Sdat$area%in%"Berlin"])
})

FishingC.pval <- lapply(FishingC, "[[", "p.value")

table(p.adjust(FishingC.pval)<0.05)

FishingC[p.adjust(FishingC.pval)<0.05]

## Eucolus more prevalent AND ABUNDANT in Brandenburg
tapply(HelmCounts.t["Eucoleus", ], Sdat$area, median)
tapply(HelmCounts.t["Eucoleus", ], Sdat$area, quantile, 0.75)
tapply(HelmCounts.t["Eucoleus", ], Sdat$area, mean)

## Clonorchis more prevalent AND ABUNDANT in Brandenburg
tapply(HelmCounts.t["Clonorchis", ], Sdat$area, median)
tapply(HelmCounts.t["Clonorchis", ], Sdat$area, quantile, 0.75)
tapply(HelmCounts.t["Clonorchis", ], Sdat$area, mean)

## Angiostrongylus more prevalent AND ABUNDANT in Berlin
tapply(HelmCounts.t["Angiostrongylus", ], Sdat$area, median)
tapply(HelmCounts.t["Angiostrongylus", ], Sdat$area, mean)

## Alaria more prevalent AND ABUNDANT in Brandenburg
tapply(HelmCounts.t["Alaria", ], Sdat$area, median)
tapply(HelmCounts.t["Alaria", ], Sdat$area, quantile, 0.75)
tapply(HelmCounts.t["Alaria", ], Sdat$area, mean)

## Strongyloides more prevalent AND ABUNDANT in Brandenburg
tapply(HelmCounts.t["Strongyloides", ], Sdat$area, median)
tapply(HelmCounts.t["Strongyloides", ], Sdat$area, quantile, 0.75)
tapply(HelmCounts.t["Strongyloides", ], Sdat$area, mean)

## COUNT models with all the fox-covariates
FishingM <- apply(HelmCounts.t, 1, function (x) {
    d <- cbind(count=x, Sdat)
    glm(count~ area + condition + I(as.numeric(weight_kg)) + sex + age
         + season + year,
        offset = log(seq.depth),
        data=d, family="poisson")
})

## only the converged models
FishingM <- FishingM[unlist(lapply(FishingM, "[[", "converged"))]

MBrand.pvals <- lapply(FishingM, function(x) {
    c <- coefficients(summary(x))
    c["areaBrandenburg","Pr(>|z|)"]
})

MBrand.effects <- lapply(FishingM, function(x) {
    c <- coefficients(summary(x))
    c["areaBrandenburg","Estimate"]
})


MBrand.pvals.adj <- p.adjust(unlist(MBrand.pvals))


AreaModells <- data.frame(cbind(effect.size=
                                    exp(unlist(MBrand.effects[MBrand.pvals.adj<0.1])),
                                adj.p.value=
                                    unlist(MBrand.pvals.adj[MBrand.pvals.adj<0.1])))


AreaModells <- merge(AreaModells, traits, by=0)
  
AMtab <- AreaModells[order(AreaModells["effect.size"]),
                     c("Row.names", "effect.size", "adj.p.value",
                       "exact.tax", "zoonotic", "transmission.fox", "lifecycle", "host.range")]

write.csv(AMtab, "tables/Pre_models.csv", row.names=FALSE)

## presence/absence models with all the fox-covariates
FishingPA <- apply(HelmCounts.t, 1, function (x) {
    d <- cbind(PA=x>0, Sdat)
    glm(PA~ area + condition + I(as.numeric(weight_kg)) + sex + age  +
            season + year,
        data=d, family="binomial")
})


## only the converged models
FishingPA <- FishingPA[unlist(lapply(FishingPA, "[[", "converged"))]

PABrand.pvals <- lapply(FishingPA, function(x) {
    c <- coefficients(summary(x))
    c["areaBrandenburg","Pr(>|z|)"]
})

PABrand.effects <- lapply(FishingPA, function(x) {
    c <- coefficients(summary(x))
    c["areaBrandenburg","Estimate"]
})


PABrand.pvals.adj <- p.adjust(unlist(PABrand.pvals))


AreaModellsPA <- data.frame(cbind(effect.size=
                                      unlist(MBrand.effects[PABrand.pvals.adj<0.1]),
                                  adj.p.value=
                                      unlist(MBrand.pvals.adj[PABrand.pvals.adj<0.1])))


AreaModellsPA <- merge(AreaModellsPA, traits, by=0)

AreaModellsPA[order(AreaModellsPA["effect.size"]),
              c("Row.names", "effect.size", "adj.p.value",
                "exact.tax", "zoonotic", "transmission.fox", "lifecycle")]

## Prevalences:
apply(HelmCounts.t, 1, function (x){
    cbind(length(x[x>0]),
    (length(x[x>0])/length(x))*100)
})
