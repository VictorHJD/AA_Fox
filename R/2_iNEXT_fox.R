library(phyloseq)
library(iNEXT)
library(vegan)
library(betapart)
library(tidyverse)
library(patchwork)
library(colorspace)
library(MASS)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

## You can leave this at true here, as the "bioinformatics" script has
## finer controls within itself for now. Re-executing this for now
## adds only the (potentially new) sample data (created in
## 0_Extract_Einvir_Covariates.R) to the PS object.
recomputeBioinfo <- TRUE

if(!exists("PS")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PS <- readRDS(file="intermediate_data/PhyloSeqCombi.Rds")
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


## Subsetting for Helminths:


## "parasites" plus other nematoda
helminths <- c("Nematoda", "Platyhelminthes")

## WHY do we only look a helminths not acutally at protozoa
## (e.g. Coccidia) too? --> NOW see below they are also more diverse
## in Brandenburg

## unicellular "parasites" (incl. Gregarina)
protist.phyla <- c("Microsporidia", "Apicomplexa")
## and more by genus (see
## http://tolweb.org/accessory/Parasitic_Protists?acc_id=53) searched
## the list of "weird (no identity developed" protist genera there
## against the taxon table")
protist.genera <- c("Acanthamoeba",
                    "Dermocystidium",
                    "Ichthyophonus")

## or here an idea for diet:
diet <- c("Annelida", "Arthropoda", "Chordata", "Mollusca")


PSHelm <- phyloseq::subset_taxa(PS, phylum%in%helminths)

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
table(Sdat$area, MoreOne=rowSums(HelmCounts>0)>1)
fisher.test(table(Sdat$area, MoreOne=rowSums(HelmCounts>0)>1))

OTU_inext_imp <- iNEXT(t(indHelmCounts), q =0, datatype = "abundance")

### Now the plot gets a bit messy redo by hand
zet <- fortify(OTU_inext_imp)

Sdat$IZW_ID <- paste("Fox", Sdat$IZW_ID)

zet <- merge(zet, Sdat, by.x="site", by.y="IZW_ID")

alphaDivFox <- ggplot(zet, aes(x = x, y = y, colour = area, group=site)) +
    geom_line(data=subset(zet, method%in%"interpolated"),
              lwd = 1.5, alpha=0.3) +
    geom_point(data=subset(zet, method%in%"observed"), shape = 21, fill = "white", size=2.7, stroke = .8) + 
    geom_point(data=subset(zet, method%in%"observed"), shape = 21, fill = "transparent", size=2.7, stroke = .8) + 
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    ylab("Helminth diversity") +
    xlab("Number of sequence reads") +
    theme(legend.position="none", axis.text = element_text(family = font_num))


## compare the asymptotic estimates between Berlin and Brandenburg
## But iNext loses the rownames and all extra information, add back:
EstimatesAsy <- cbind(OTU_inext_imp$AsyEst, SdatHPres[rep(1:nrow(SdatHPres), each=3), ])

EstimatesAsy <- filter(EstimatesAsy,
                       !is.na(area) &
                       !is.na(tree_cover_1000m) &
                       !is.na(imperv_1000m)&
                       !is.na(human_fpi_1000m))
                       

alphaCompared <- ggplot(EstimatesAsy, aes(area, Estimator, color=area, fill = after_scale(lighten(color, .7)))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021), fill = "white", size = 1.3, stroke = .7) +
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL) +
    facet_wrap(~Diversity, scales="free_y")+
    scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
    scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
    theme(legend.position="none", axis.text.y = element_text(family = font_num))

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
    glm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
       data=., family="poisson") ->
    DivModelHillArea

summary(DivModelHillArea)

EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    glm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
           data=., family="poisson") ->
    DivModelHillTree

summary(DivModelHillTree)

EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    glm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=., family="poisson") ->
    DivModelHillHum

summary(DivModelHillHum)


EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    glm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
       data=., family="poisson") ->
    DivModelHillImp

summary(DivModelHillImp)

AIC(DivModelHillArea, DivModelHillTree, DivModelHillImp, DivModelHillHum)


tab_model(DivModelHillArea, DivModelHillTree, DivModelHillImp, DivModelHillHum,
          file="tables/HillDiv.html", show.aic = TRUE)

tab_model(DivModelSimpsonArea, DivModelSimpsonTree, DivModelSimpsonImp,
          DivModelSimpsonHum, file="tables/SimpsonDiv.html", show.aic = TRUE)

tab_model(DivModelShannonArea, DivModelShannonTree, DivModelShannonImp,
          DivModelShannonHum, file="tables/ShannonDiv.html", show.aic = TRUE)


### Let's see whether there's a non-linear "ecotone" effect of any of
### the continuous environmental variables!

EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    ggplot(aes(human_fpi_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    ggplot(aes(human_fpi_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    ggplot(aes(human_fpi_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    ggplot(aes(imperv_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    ggplot(aes(imperv_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    ggplot(aes(imperv_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()


EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
    ggplot(aes(tree_cover_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    ggplot(aes(tree_cover_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    ggplot(aes(tree_cover_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()



## Now: the way Caro designed the analysis it distinguishes between
## Berlin and Brandenburg as a whole (and between male and female,
## maybe?). I'd call this level gamma-diversity!


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
    xlab("Number of sampled foxes") +
    ylab("Helminth diversity") +
    ## need to repeat theme because it is overwritten by the wrapper
    theme_minimal(base_family = "Roboto", base_size = 12) +
    theme(legend.position="none", 
          axis.text = element_text(family = font_num),
          axis.title.x = element_text(margin = margin(t = 12)),
          axis.title.y = element_text(margin = margin(r = 12)),
          panel.grid.minor = element_blank(),
          plot.margin = margin(rep(12, 4)))



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
    ggplot(aes(area, distances, color=area, fill = after_scale(lighten(color, .7)))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021), fill = "white", size = 2, stroke = .7) +
    #geom_point(position = position_jitter(width = .3, seed = 2021)) +
    scale_y_continuous(name = "Distance to area centroid") +
    scale_x_discrete(name = NULL) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    guides(fill = guide_legend(title.position = "top", title.hjust = .5),
           color = guide_legend(title.position = "top", title.hjust = .5)) +
    theme(legend.position="top", axis.text.y = element_text(family = font_num))->
    betaDivJac


betaDivJacMulti <- 
    wrap_elements(full =
                      ~(plot(JaccGrups, col=c("#e7b800", "#2e6c61"), main="",
                             label = FALSE, sub="")))


wrap_plots(
    ## place first two plots
    alphaDivFox, alphaCompared, ## -> A + B
    ## place the legend in the middle
    guide_area(), ## -> C
    ## ... then the other three plots
    betaDivJacMulti, betaDivJac, gammaDivFox, ## -> D, E + F
    ## you can build more complex layouzts by providing simple letters that are
    ## then filled accordingly by the plots you defined in the previous step;
    ## the plots are "named" as the appear here: the first one is A, the next B and so on...
    design = "AAABBB\n##CC##\nDDEEFF",
    ## by default all rows and columns have similar widths and heights but we
    ## don't want our legend to fill up 1/3 of the plot height
    heights = c(21/45, 3/45, 21/45),
    ## this tells ggplot to just draw the legend once, not six times
    guides = "collect"
) +
plot_annotation(tag_levels = 'a')

ggsave("figures/Diversity.pdf", width=13, height=9, device=cairo_pdf)


## checking Apicomplexa
PSApico <- phyloseq::subset_taxa(PS, phylum%in%c("Apicomplexa"))
PSApicoG <- phyloseq::tax_glom(PSApico, "genus")
ApicoCounts <- as.data.frame(otu_table(PSApicoG))

Sdat <- as.data.frame(sample_data(PSApicoG))

## 33 samples with less than 2 species
table(Sdat$area, MoreOne=rowSums(ApicoCounts>0)>1)
## Many more of those in Berlin (significantly so)
fisher.test(table(Sdat$area, MoreOne=rowSums(ApicoCounts>0)>1))

## 13 samples with 0 species
table(Sdat$area, MoreOne=rowSums(ApicoCounts)>0)
## All in Berlin, that's significant
fisher.test(table(Sdat$area, MoreOne=rowSums(ApicoCounts)>0))


PSArthro <- phyloseq::subset_taxa(PS, phylum%in%c("Arthropoda"))
PSArthroG <- phyloseq::tax_glom(PSArthro, "genus")
ArthroCounts <- as.data.frame(otu_table(PSArthroG))

fisher.test(table(Sdat$area, MoreOne=rowSums(ArthroCounts>0)>1))
fisher.test(table(Sdat$area, MoreOne=rowSums(ArthroCounts)>0))


remove_geom <- function(ggplot2_object, geom_type) {
    ## Delete layers that match the requested type.
    layers <- lapply(ggplot2_object$layers, function(x) {
        if (class(x$geom)[1] == geom_type) {
            NULL
        } else {
            x
        }
    })
    ## Delete the unwanted layers.
    layers <- layers[!sapply(layers, is.null)]
    ggplot2_object$layers <- layers
    ggplot2_object
}


apicoRichness <- plot_richness(PSApicoG, measures=c("Observed", "Shannon"),
                               color="area", x="area") + geom_boxplot(outlier.shape = NA)

remove_geom(apicoRichness, "GeomPoint")  +
    geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021),
               fill = "white", size = 2, stroke = .7) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_y_continuous("Study area") + 
    ggtitle("Apicomplexa (incl. gregarina) diversity")
ggsave("figures/suppl/ApicoDiversity.pdf", width=13, height=9, device=cairo_pdf)

cbind(estimate_richness(PSApico, measures=c("Observed", "Shannon")),
      sample_data(PSApicoG))  %>% 
    lm(Shannon~ area + condition + I(as.numeric(weight_kg)) + sex + age,
        data=.) ->
    ModDivApicoShannon

summary(ModDivApicoShannon)

library(MASS)
cbind(estimate_richness(PSApico, measures=c("Observed", "Shannon")),
      sample_data(PSApicoG))  %>% 
    glm.nb(Observed~ area + condition + I(as.numeric(weight_kg)) + sex + age,
        data=.) ->
    ModDivApicoObserved

summary(ModDivApicoObserved)

arthroRichness <- plot_richness(PSArthroG, measures=c("Observed", "Shannon"),
                               color="area", x="area") + geom_boxplot(outlier.shape = NA)

remove_geom(arthroRichness, "GeomPoint")  +
    geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021),
               fill = "white", size = 2, stroke = .7) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    ggtitle("Arthropoda diversity")


cbind(estimate_richness(PSArthro, measures=c("Observed", "Shannon")),
      sample_data(PSApicoG))  %>% 
    lm(Shannon~ area + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    ModDivArthroShannon

summary(ModDivArthroShannon)


cbind(estimate_richness(PSArthro, measures=c("Observed", "Shannon")),
      sample_data(PSApicoG))  %>% 
    glm.nb(Observed~ area + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    ModDivArthroObserved

summary(ModDivArthroObserved)

