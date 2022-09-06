library(vegan)
library(tidyverse)
library(phyloseq)
library(patchwork)

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
                                   !is.na(sample_data(PSGHelm)[, "weight_kg"]))


## only taxa with reads
PSGHelm <- prune_taxa(taxa_sums(PSGHelm) > 0, PSGHelm)
PSGHelm <- prune_samples(sample_sums(PSGHelm) > 0, PSGHelm)

HelmData <- otu_table(PSGHelm)
colnames(HelmData) <- tax_table(PSGHelm)[, "genus"]

EnvData <- sample_data(PSGHelm)
class(EnvData) <- "data.frame"
EnvData$weight_kg <- as.numeric(EnvData$weight_kg)


### NO OTHER environmental variables are better explaining composition!

adonis2(HelmData[!is.na(EnvData$tree_cover_1000m), ] ~ age + weight_kg +
        + sex + condition + tree_cover_1000m,
        data=EnvData[!is.na(EnvData$tree_cover_1000m), ],
        na.action = na.omit, by="margin")


adonis2(HelmData ~ area + FunM_Species_richness  +
            Diet_Species_richness + BacM_Species_richness,
        data=EnvData,
        na.action = na.omit, by="margin")


adonis2(HelmData[!is.na(EnvData$imperv_1000m), ] ~ imperv_1000m + age + weight_kg
        + sex + condition + DDiv + BacDiv,
        data=EnvData[!is.na(EnvData$imperv_1000m), ],
        na.action = na.omit, by="margin")

### Even human fpi is barely as good as rural/urban!
adonis2(HelmData[!is.na(EnvData$human_fpi_1000m), ] ~ human_fpi_1000m + age + weight_kg
        + sex + condition + DDiv + BacDiv,
        data=EnvData[!is.na(EnvData$human_fpi_1000m), ],
        na.action = na.omit, by="margin")


## how does diet influence Helminth occurence?
PSDiet <- subset_taxa(PSF, phylum %in% "Chordata" &
                           !genus%in%c("Vulpes", "Homo", "Globicephala", "Hylobates", "Procyon", "Canis"))


## Collapse Chordata by family
PSDietFam <- tax_glom(PSDiet, "family")
PSDietFam <- prune_taxa(taxa_sums(PSDietFam) > 10, PSDietFam)

## Collapse invertebrates by class
PSDietCollapsed <- merge_phyloseq(PSDietFam, 
                                  tax_glom(subset_taxa(PSF, phylum %in%
                                                            c("Annelida", "Arthropoda","Mollusca")), "class"))

PSDietCollapsed <- prune_taxa(taxa_sums(PSDietCollapsed) > 10, PSDietCollapsed)

PSDietCollapsed <- prune_samples(sample_names(PSDietCollapsed)%in%sample_names(PSHelm), PSDietCollapsed)

DietCollapsed <- otu_table(PSDietCollapsed)

## names from faimily in Chordata from class in invertebrates
colnames(DietCollapsed) <- ifelse(tax_table(PSDietCollapsed)[, "phylum"] %in%"Chordata",
                                  as.vector(tax_table(PSDietCollapsed)[, "family", drop=T]),
                                  as.vector(tax_table(PSDietCollapsed)[, "class", drop=T]))

set.seed(123)
nMDSHelm <- metaMDS(HelmData, distance = "jaccard", weakties = FALSE, try=1500, trymax=1500, k=3,
                    center = TRUE)

HelmDietEnvFit <- envfit(nMDSHelm, cbind(DietCollapsed,
                                         EnvData[ , c("area", "age", "weight_kg", "sex", "condition")]),
                         na.rm=TRUE)
### AMAZING! ThIs MAKeS SeNSe!!!!

## to see which helminths drive this
HelmHelmFit <- envfit(nMDSHelm, HelmData, na.rm=TRUE)


HelmDietEnvFitDf <-
    as.data.frame(
        rbind(cbind(scores(HelmDietEnvFit, "vectors") * ordiArrowMul(HelmDietEnvFit),
                    pvals=HelmDietEnvFit$vectors$pvals),
              cbind(scores(HelmDietEnvFit, "factors") * ordiArrowMul(HelmDietEnvFit),
                    ### works only because each factor has two levels!!!
                    pvals=rep(HelmDietEnvFit$factors$pvals, each=2)))
    )

HelmHelmDf <-
    as.data.frame(
        cbind(scores(HelmHelmFit, "vectors") * ordiArrowMul(HelmHelmFit),
              pvals=HelmHelmFit$vectors$pvals)
    )



HelmDietEnvFitDf$Cat <- ifelse(rownames(HelmDietEnvFitDf)%in%colnames(DietCollapsed), "Diet",
                               "Environmental")

ScoresHelm <-  as.data.frame(scores(nMDSHelm))
ScoresHelm <- cbind(ScoresHelm, EnvData[, c("area", "age", "weight_kg", "sex", "condition")])

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
## font_num <- "Roboto Condensed"


ggHelmEnv <-  ggplot(data = ScoresHelm, aes(x = NMDS1, y = NMDS2)) +
    geom_point(data = ScoresHelm, aes(colour = area), size = 3) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") 

ggHelmEnv +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, size = 0.1-pvals),
                 data = subset(HelmDietEnvFitDf, pvals<0.1),
                 arrow = arrow(length = unit(0.04, "npc"), angle=23)) +
    geom_text(data = subset(HelmDietEnvFitDf, pvals<0.1), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(HelmDietEnvFitDf, pvals<0.1)), size=4.5,
              color=ifelse(subset(HelmDietEnvFitDf, pvals<0.1)$Cat%in%"Diet", "red", "green"))+
              scale_size_continuous(range=c(1, 3))  -> ggHelmEnvDiet


ggHelmEnv +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, size = pvals),
                 data = subset(HelmHelmDf, pvals<=0.1),
                 arrow = arrow(length = unit(0.04, "npc"), angle=23)) +
    scale_size_binned(trans = 'reverse', range=c(0,3)) +
    geom_text(data = subset(HelmHelmDf, pvals<0.1), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(HelmHelmDf, pvals<0.1)), size=4.5,
              color="blue")+
    theme_bw() -> ggHelmEnvHelm

pdf("figures/composition_Diet_Helm.pdf", width=21, height=7)
ggHelmEnvDiet + ggHelmEnvHelm
dev.off()
