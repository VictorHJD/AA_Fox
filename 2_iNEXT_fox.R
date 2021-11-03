library(phyloseq)
library(iNEXT)
library(vegan)
library(tidyverse)


recompute <- FALSE

if(!exists("PS")){
    if(recompute){
        source("1_Fox_general_MA.R")
    } else {
        PS <- readRDS(file="data/PhyloSeqCombi.Rds")
    }
}


## Subsetting for Helminths:

## WHY do we only look a helminths not acutally at protozoa
## (e.g. Coccidia) too?

PSHelm <- phyloseq::subset_taxa(PS, phylum%in%c("Nematoda", "Platyhelminthes"))

PSHelmS <- phyloseq::tax_glom(PSHelm, "species")

PSHelmG <- phyloseq::tax_glom(PSHelm, "genus")

HelmCounts <- as.data.frame(otu_table(PSHelmG))

colnames(HelmCounts) <- tax_table(PSHelmG)[, "genus"]
rownames(HelmCounts) <- paste("Fox", rownames(HelmCounts))

## Per species sequence numbers
colSums(HelmCounts)
## Per Fox sequence numbers
rowSums(HelmCounts)
## Per Fox species numbers
rowSums(HelmCounts>0)

### First a rarefaction analysis per sample!!!
## We need to keep only samples with more than one species
indHelmCounts <- HelmCounts[rowSums(HelmCounts>0)>1, ]

OTU_inext_imp <- iNEXT(t(indHelmCounts), q =0, datatype = "abundance")

### Now the plot gets a bit messy redo by hand
zet <- fortify(OTU_inext_imp)

## merge it with the sample data
Sdat <- as.data.frame(sample_data(PSHelmG))
class(Sdat) <- "data.frame"

Sdat$IZW_ID <- paste("Fox", Sdat$IZW_ID)

zet <- merge(zet, Sdat, by.x="site", by.y="IZW_ID")

divFoxPlot2000 <- ggplot(zet, aes(x = x, y = y, colour = area, group=site)) +
    geom_point(data=subset(zet, method%in%"observed"), size=4) + 
    geom_line(data=subset(zet, method%in%"interpolated"),
              lwd = 1.5, alpha=0.3) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    ## guides(linetype = guide_legend(title = "Method"),
    ##        colour = guide_legend(title = "Guides"), 
    ##        fill = guide_legend(title = "Guides"),
    ##        shape = guide_legend(title = "Guides")) + 
    theme(legend.position = "bottom", legend.title = element_blank(), 
          text = element_text(size = 18), legend.key.width = unit(1.2, 
                                                                  "cm")) +
    theme(legend.position = "bottom", legend.title = element_blank(), 
          text = element_text(size = 18), legend.key.width = unit(1.2, 
                                                                  "cm")) +
    ylab("helminth diversity\n") +
    xlab("number of sequence reads\n") +
##    scale_x_continuous(limits=c(0, 2000)) +
    theme_bw()

divFoxPlot2000
### 52 observed values above the sequencing count shown onx!

## Now: the way Caro designed the analysis it distinguishes between
## Berlin and Brandenburg as a whole (and between male and female,
## maybe?)

## Metadata annotated presence species number per fox

HelmCountsAnn <- merge(HelmCounts, Sdat, by.x=0, by.y="IZW_ID")

HelmPresAbsAnn1 <- by(HelmCountsAnn, HelmCountsAnn$area,
                     function (x) {colSums(x[, 2:ncol(HelmCounts)]>0)})

HelmPresAbsAnn1 <- do.call(rbind, HelmPresAbsAnn1)
    
OTU_inext_area1 <- iNEXT(t(HelmPresAbsAnn1), q =0, datatype = "abundance")

HelmPresAbsAnn100 <- by(HelmCountsAnn, HelmCountsAnn$area,
                     function (x) {colSums(x[, 2:ncol(HelmCounts)]>100)})

HelmPresAbsAnn100 <- do.call(rbind, HelmPresAbsAnn100)
    
OTU_inext_area100 <- iNEXT(t(HelmPresAbsAnn100), q =0, datatype = "abundance")

plot_area1 <- ggiNEXT(OTU_inext_area1) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    theme_minimal() +
    xlab("number of sampled foxes") +
    ylab("helminth diversity\n") +
    theme(legend.position="bottom",
          legend.box = "vertical", 
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16),
          legend.title = element_text(color = "#5c564b", size = 12, face = "bold"),
          legend.text = element_text(color = "#5c564b", size = 12)) +
    guides(size = guide_legend(override.aes = list(shape = 1)))

plot_area1


plot_area100 <- ggiNEXT(OTU_inext_area100) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    theme_minimal() +
    xlab("number of sampled foxes") +
    ylab("helminth diversity\n") +
    theme(legend.position="bottom",
          legend.box = "vertical", 
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16),
          legend.title = element_text(color = "#5c564b", size = 12, face = "bold"),
          legend.text = element_text(color = "#5c564b", size = 12)) +
    guides(size = guide_legend(override.aes = list(shape = 1)))

plot_area100

HelmReadAnn <- by(HelmCountsAnn, HelmCountsAnn$area,
                  function (x) {colSums(x[, 2:ncol(HelmCounts)])})

HelmReadAnn <- do.call(rbind, HelmReadAnn)
    
Read_inext_area <- iNEXT(t(HelmReadAnn), q =0, datatype = "abundance")

plot_areaSeq <- ggiNEXT(Read_inext_area) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    theme_minimal() +
    xlab("number of sequence reads") +
    ylab("helminth diversity\n") +
    theme(legend.position="bottom",
          legend.box = "vertical", 
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16),
          legend.title = element_text(color = "#5c564b", size = 12, face = "bold"),
          legend.text = element_text(color = "#5c564b", size = 12)) +
    guides(size = guide_legend(override.aes = list(shape = 1)))

plot_areaSeq

shannonFox <- vegan::diversity(HelmCounts, index="shannon")

shannonPresAbs1 <- vegan::diversity(HelmPresAbsAnn1, index="shannon")

shannonPresAbs100 <- vegan::diversity(HelmPresAbsAnn100, index="shannon")

shannonRead <- vegan::diversity(HelmReadAnn, index="shannon")
