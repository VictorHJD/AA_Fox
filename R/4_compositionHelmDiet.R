library(vegan)


if(!exists("PS")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PS <- readRDS(file="intermediate_data/PhyloSeqCombi.Rds")
    }
}


## Helminth traits
traits <- read.csv("input_data/helminth_traits.csv")
traits %>%
    column_to_rownames("t.genus")  -> traits 

BADtaxa <- rownames(traits)[!traits$BlastEvaluation%in%"Okay"]
NOPara <- rownames(traits)[traits$fox.parasite%in%"No"]

## collapse to genus level
PSG <- phyloseq::tax_glom(PS, "genus")

PSF <- subset_taxa(PSG, !genus%in%NOPara &
                        !genus%in%BADtaxa)

PSF <- subset_samples(PSG, !is.na(sample_data(PSF)[, "condition"]) &
                           !is.na(sample_data(PSF)[, "weight_kg"]))

PSHelm <- subset_taxa(PSF, phylum %in% c("Nematoda", "Platyhelminthes"))

## only taxa with at lest 10 reads
PSHelm <- prune_taxa(taxa_sums(PSHelm) > 0, PSHelm)
PSHelm <- prune_samples(sample_sums(PSHelm) > 0, PSHelm)

HelmData <- otu_table(PSHelm)

EnvData <- sample_data(PSHelm)
class(EnvData) <- "data.frame"
EnvData <- EnvData[, c("area", "weight_kg", "age", "sex", "condition")]
EnvData$weight_kg <- as.numeric(EnvData$weight_kg)

adonis2(HelmData ~ age + weight_kg + sex + condition,
        data=EnvData, na.action = na.omit, by="margin")


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
nMDSHelm <- metaMDS(HelmData, distance = "jaccard", weakties = FALSE, try=250, trymax=250, k=4,
                    center = TRUE)

HelmDietFit <- envfit(nMDSHelm, cbind(DietCollapsed, EnvData), na.rm=TRUE)
HelmEnvFit <- envfit(nMDSHelm, EnvData, na.rm=TRUE)

### AMAZING! ThIs MAKeS SeNSe!!!!


