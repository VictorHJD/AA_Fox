library(phyloseq)
library(iNEXT)
library(vegan)
library(betapart)
library(tidyverse)
library(patchwork)

library(sjPlot)
library(sjmisc)
library(sjlabelled)

recomputeBioinfo <- FALSE

if(!exists("PS")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PS <- readRDS(file="intermediate_data/PhyloSeqCombi.Rds")
    }
}


## Subsetting for Helminths:

## WHY do we only look a helminths not acutally at protozoa
## (e.g. Coccidia) too?

PSHelm <- phyloseq::subset_taxa(PS, phylum%in%c("Nematoda", "Platyhelminthes"))
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


### For inext diversity analysis we !
## We need to keep only samples with at least two species
indHelmCounts <- HelmCounts[rowSums(HelmCounts>0)>1, ]


## Sample data in data frame
Sdat <- as.data.frame(sample_data(PSHelmG))
class(Sdat) <- "data.frame"

## sample data for samples with helminth present
SdatHPres <- Sdat[rowSums(HelmCounts>0)>1, ]

## how many lost for each area
table(Sdat$area, MoreOne=rowSums(HelmCounts>0)>1)

## 21 samples with less than 2 species. No significant differences in
## the areas
fisher.test(table(Sdat$area, MoreOne=rowSums(HelmCounts>0)>1))

OTU_inext_imp <- iNEXT(t(indHelmCounts), q =0, datatype = "abundance")

### Now the plot gets a bit messy redo by hand
zet <- fortify(OTU_inext_imp)

Sdat$IZW_ID <- paste("Fox", Sdat$IZW_ID)

zet <- merge(zet, Sdat, by.x="site", by.y="IZW_ID")

alphaDivFox <- ggplot(zet, aes(x = x, y = y, colour = area, group=site)) +
    geom_point(data=subset(zet, method%in%"observed"), size=4) + 
    geom_line(data=subset(zet, method%in%"interpolated"),
              lwd = 1.5, alpha=0.3) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    ylab("helminth diversity\n") +
    xlab("number of sequence reads\n") +
    theme_minimal() +
    theme(legend.position="none")


## compare the asymptotic estimates between Berlin and Brandenburg
## But iNext loses the rownames and all extra information, add back:
EstimatesAsy <- cbind(OTU_inext_imp$AsyEst, SdatHPres[rep(1:nrow(SdatHPres), each=3), ])

EstimatesAsy <- filter(EstimatesAsy,
                       !is.na(area) &
                       !is.na(tree_cover_1000m) &
                       !is.na(imperv_1000m)&
                       !is.na(human_fpi_1000m))
                       

alphaCompared <- ggplot(EstimatesAsy, aes(area, Estimator, color=area)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    scale_y_continuous("") +
    scale_x_discrete("") +
    facet_wrap(~Diversity, scales="free_y")+
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    theme_minimal() +
    theme(legend.position="none")

### Models for Shannon
EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelShannonArea

summary(DivModelShannonArea)

EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelShannonTree

summary(DivModelShannonTree)

EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelShannonHum

summary(DivModelShannonHum)


EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelShannonImp

summary(DivModelShannonImp)

AIC(DivModelShannonArea, DivModelShannonTree, DivModelShannonImp, DivModelShannonHum)

### Models for Simpson

EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelSimpsonArea

summary(DivModelSimpsonArea)

EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelSimpsonTree

summary(DivModelSimpsonTree)

EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelSimpsonHum

summary(DivModelSimpsonHum)


EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelSimpsonImp

summary(DivModelSimpsonImp)

AIC(DivModelSimpsonArea, DivModelSimpsonTree, DivModelSimpsonImp, DivModelSimpsonHum)


### Models for Species richness

EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelHillArea

summary(DivModelHillArea)

EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelHillTree

summary(DivModelHillTree)

EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelHillHum

summary(DivModelHillHum)


EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    DivModelHillImp

summary(DivModelHillImp)

AIC(DivModelHillArea, DivModelHillTree, DivModelHillImp, DivModelHillHum)


tab_model(DivModelHillArea, DivModelHillTree, DivModelHillImp, DivModelHillHum,
          file="tables/HillDiv.html", show.aic = TRUE)

tab_model(DivModelSimpsonArea, DivModelSimpsonTree, DivModelSimpsonImp,
          DivModelSimpsonHum, file="tables/SimpsonDiv.html", show.aic = TRUE)

tab_model(DivModelShannonArea, DivModelShannonTree, DivModelShannonImp,
          DivModelShannonHum, file="tables/ShannonDiv.html", show.aic = TRUE)


## Now: the way Caro designed the analysis it distinguishes between
## Berlin and Brandenburg as a whole (and between male and female,
## maybe?)


HelmPres1 <- as.data.frame(HelmCounts>0)
HelmPres1 <- apply(HelmPres1, 1, as.numeric)

HelmPresArea <- by(t(HelmPres1), Sdat$area, function (x) x)


## for iNext all sites have to have more than one species!
Fox_inext_area1 <- iNEXT(list(Berlin=t(HelmPresArea$Berlin),
                              Brandenburg=t(HelmPresArea$Brandenburg)), 
                         q=0, datatype = "incidence_raw")

gammaDivFox <- ggiNEXT(Fox_inext_area1) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    theme_minimal() +
    xlab("number of sampled foxes") +
    ylab("helminth diversity\n") +
    theme(legend.position="none")



## beta diversity
JaccPairsDist <- beta.pair(t(apply(indHelmCounts>0, 1, as.numeric)),
                           index.family="jaccard")

JaccGrups <- betadisper(JaccPairsDist[[3]], SdatHPres$area)

## plot(JaccGrups, col=c("#e7b800", "#2e6c61"))

### this is the proper way to do it!
anova(JaccGrups)
## nonparametric test shows distances in groups don't differe significantly
## wilcox.test(JaccGrups$distances ~ JaccGrups$group)

data.frame(distances=JaccGrups$distances,
             area=JaccGrups$group) %>%
    ggplot(aes(area, distances, color=area)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    scale_y_continuous("distance to area centroid") +
    scale_x_discrete("") +
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    theme_minimal() + 
    theme(legend.position="top")->
    betaDivJac


betaDivJacMulti <- 
    wrap_elements(full =
                      ~(plot(JaccGrups, col=c("#e7b800", "#2e6c61"), main="",
                             label = FALSE, sub="")))


pdf("figures/Diversity.pdf", width=12, height=8)
(alphaDivFox + alphaCompared) /
    (betaDivJacMulti + betaDivJac  + gammaDivFox) +
    plot_annotation(tag_levels = 'a')
dev.off()    

