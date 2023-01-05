library(phyloseq)
library(tidyverse)
library(MASS)
library(stargazer)

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

foo <- cbind(HelmCounts, Sdat)

foo$weight_kg <- as.numeric(foo$weight_kg)

helminths <- c("Ancylostoma", "Aelurostrongylus", "Opisthorchis", "Taenia", 
               "Eucoleus", "Clonorchis", "Mesocestoides", "Angiostrongylus", 
               "Uncinaria", "Crenosoma", "Alaria", "Strongyloides", 
               "Toxocara", "Toxascaris", "Brachylaima", "Pearsonema")


FishingM <- lapply(helminths, function (x) {    
    fo <- formula(paste0(x, "/seq.depth ~ area + condition + weight_kg +
                       sex + age +  season + DNAng.ul + DNA260.230 + DNA260.280"))
    print(fo)
    tryCatch({
        glm(fo, data=foo, family="quasipoisson", weights=seq.depth)
    }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

names(FishingM) <- helminths

helminths.converged <- unlist(lapply(FishingM, "[[", "converged"))

FishingM <- FishingM[names(helminths.converged)[helminths.converged]]

lapply(FishingM, summary)


## only the converged models
FishingM <- FishingM[unlist(lapply(FishingM, "[[", "converged"))]


## to figure a p-value correction out for stargrazer

lapply(FishingM, function(M) {
    summary(M)$coefficients[, "Pr(>|t|)"]
}) %>% do.call(rbind, .)  %>% apply(., 2, p.adjust, method="fdr")



pCor <- function (x) p.adjust(x, method="fdr", n=length(FishingM))

stargazer(FishingM, out="tables/IndHelmAbu.html", type="html",
          column.labels=names(FishingM), dep.var.caption="",
          dep.var.labels="", dep.var.labels.include=FALSE,
          apply.p = pCor, notes.append=TRUE,
          notes= "fals discoverey rate (FDR) corrected for multiple testing")

## Prevalences:
apply(HelmCounts.t, 1, function (x){
    cbind(Numb = round(length(x[x>0])),
          Prev = round(length(x[x>0])/length(x)*100, 2))
})
