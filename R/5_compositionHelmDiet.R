library(vegan)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ggnewscale)

recomputeBioinfo <- FALSE

if(!exists("PSG")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSG <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
    }
}

recomputeDiversity <- FALSE

if(!"FunM_Species_richness"%in%colnames(sample_data(PSG))|
   recomputeDiversity){
    source("R/2_iNEXT_fox.R")
}

## Helminth traits
traits <- read.csv("input_data/helminth_traits.csv")
traits %>%
    column_to_rownames("t.genus")  -> traits 

PSGHelm <- phyloseq::subset_taxa(PSG, category%in%"Helminth")

PSGHelm <- subset_samples(PSGHelm, !is.na(sample_data(PSGHelm)[, "condition"]) &
                                   !is.na(sample_data(PSGHelm)[, "weight_kg"]) &
                                   !is.na(sample_data(PSGHelm)[, "season"]) &
                                   !is.na(sample_data(PSGHelm)[, "tree_cover_1000m"]))

PSGHelm ### 114 foxes without NA anywhere

## only taxa with reads
PSGHelm <- prune_taxa(taxa_sums(PSGHelm) > 0, PSGHelm)
PSGHelm <- prune_samples(sample_sums(PSGHelm) > 0, PSGHelm)

PSGHelm ## 102 foxes with all samples any taxa

HelmData <- otu_table(PSGHelm)
colnames(HelmData) <- tax_table(PSGHelm)[, "genus"]

EnvData <- sample_data(PSGHelm)
class(EnvData) <- "data.frame"
EnvData$weight_kg <- as.numeric(EnvData$weight_kg)
EnvData$tree_cover_1000m <- as.numeric(EnvData$tree_cover_1000m)
EnvData$imperv_1000m <- as.numeric(EnvData$imperv_1000m)
EnvData$human_fpi_1000m <- as.numeric(EnvData$human_fpi_1000m)

### This shoud be the same now after removing all the NAs already
### above
EnvDataNA <- na.omit(EnvData)
HelmDataNA <- HelmData[rownames(EnvDataNA), ]

### NO OTHER environmental variables are better explaining composition!

PERMA <- vegan::adonis2(HelmDataNA ~ area + weight_kg + age +
                     sex  + season + year +
                     Diet_Species_richness + BacM_Species_richness +
                     FunM_Species_richness,
                     data=EnvDataNA, 
                     na.action = na.fail, by="margin")

PERMA

### still area is the best model, everything below (the other
### environmental predictors) is not as good!

adonis2(HelmDataNA ~ tree_cover_1000m + weight_kg + age +
            sex  + season + year +
            Diet_Species_richness + BacM_Species_richness + FunM_Species_richness,
        data=EnvDataNA, 
        na.action = na.fail, by="margin")

adonis2(HelmDataNA ~ area + tree_cover_1000m + weight_kg + age +
            sex  + season + year +
            Diet_Species_richness + BacM_Species_richness + FunM_Species_richness,
        data=EnvDataNA, 
        na.action = na.fail, by="margin")

adonis2(HelmDataNA ~ human_fpi_1000m + weight_kg + age +
            sex  + season + year +
            Diet_Species_richness + BacM_Species_richness + FunM_Species_richness,
        data=EnvDataNA, 
        na.action = na.fail, by="margin")


## ## Phyloseq basic ordination
## bray_dist <- phyloseq::distance(PSGHelm,
##                                 method="bray", weighted=F)

## ordination <- ordinate(microbiome::transform(PSGHelm, "log"),
##                        method="PCoA", distance="bray")

## HelmOrdination <- merge(ordination$vectors, sample_data(PSGHelm), by=0)

## ggplot(HelmOrdination, aes(Axis.1, Axis.2, shape=area, color=season)) +
##     geom_point(size=4)

## ## now nMDS and envfit! For the real deal

nMDSHelm <- metaMDS(HelmData, distance = "jaccard", weakties = FALSE,
                    try=1500, trymax=1500, k=3,
                    center = TRUE)

HelmEnvFit <- envfit(nMDSHelm, EnvData[ , c("area", "age", "weight_kg",
                                            "sex", "condition", "season", "year")])



### AMAZING! ThIs MAKeS SeNSe!!!!
HelmEnvFit

## to see which helminths drive this
HelmHelmFit <- envfit(nMDSHelm, HelmData, na.rm=TRUE)

### Haha now seeing that envfit has a "tidy" argument to produce
### ggplot compatible output... all the below would likely not have
### been necessary

faclevels <- sapply(names(HelmEnvFit$factors$pvals),
                    function (x) {
                        thatoften <- sum(grepl(x, rownames(scores(HelmEnvFit, "factors"))))
                        rep(x, thatoften)
                    })

HelmEnvFitDf <-
    as.data.frame(
        rbind(cbind(scores(HelmEnvFit, "vectors") * ordiArrowMul(HelmEnvFit),
                    pvals=HelmEnvFit$vectors$pvals),
              cbind(scores(HelmEnvFit, "factors") * ordiArrowMul(HelmEnvFit),
                    ## repeat times the numbers of factor levels
                    pvals=HelmEnvFit$factors$pvals[unlist(faclevels)])))
def.off() ## this had somehow opened a graphics device?!

HelmHelmDf <-
    as.data.frame(
        cbind(scores(HelmHelmFit, "vectors") * ordiArrowMul(HelmHelmFit),
              pvals=HelmHelmFit$vectors$pvals)
    )
def.off() ## this had somehow opened a graphics device?!


ScoresHelm <-  as.data.frame(scores(nMDSHelm)$sites)
ScoresHelm <- cbind(ScoresHelm, EnvData)

## theme_set(theme_minimal(base_family = "Roboto", base_size = 12))
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


ggplot(data = ScoresHelm, aes(x = NMDS1, y = NMDS2)) +
    geom_point(data = ScoresHelm, aes(colour = area, shape = season), size = 3) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    new_scale_color()+
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, color = log(pvals)),
                 data = subset(HelmEnvFitDf, pvals<0.1),
                 arrow = arrow(length = unit(0.04, "npc"), angle=23),
                 size=1.5) +
    scale_color_viridis_c(option = "cividis") + 
    geom_text(data = subset(HelmEnvFitDf, pvals<0.1), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(HelmEnvFitDf, pvals<0.1)), size=5.5,
              color="darkgrey") +
    theme_bw() ->
    ggHelmEnv

ggHelmEnv +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, color = pvals),
                 data = subset(HelmHelmDf, pvals<=0.1),
                 arrow = arrow(length = unit(0.04, "npc"), angle=23), size=1.5) +
    geom_text(data = subset(HelmHelmDf, pvals<0.1), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(HelmHelmDf, pvals<0.1)), size=4.5,
              color="blue")+
    theme_bw() -> ggHelmEnvHelm

pdf("figures/composition_Env_Helm.pdf", width=14, height=7)
ggHelmEnvHelm
dev.off()



