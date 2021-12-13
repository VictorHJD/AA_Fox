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
recomputeBioinfo <- FALSE

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


### As we now want diversity for different taxonomic (phyla) subsets
### I've put all off this in on giant function 
getAllDiversity <- function (PS, output_string) {
    PSG <- phyloseq::tax_glom(PS, "genus")
    Counts <- as.data.frame(otu_table(PSG))

    colnames(Counts) <- tax_table(PSG)[, "genus"]
    rownames(Counts) <- paste("Fox", rownames(Counts))

    ## For inext diversity analysis we need to keep only samples with
    ## at least two species
    indCounts <- Counts[rowSums(Counts>0)>1, ]    
    
    ## Sample data in data frame
    Sdat <- as.data.frame(sample_data(PSG))
    class(Sdat) <- "data.frame"

    ## same for the data
    SdatHPres <- Sdat[rowSums(Counts>0)>1, ]
    
    ## produce ouptut for interactive review
    message("\n Significance of removed data:") 
    print(table(Sdat$area, MoreOne=rowSums(Counts>0)>1))
    print(fisher.test(table(Sdat$area, MoreOne=rowSums(Counts>0)>1)))
    message("\n")
    
    OTU_inext_imp <- iNEXT(t(indCounts), q =0, datatype = "abundance")

### Now the plot gets a bit messy redo by hand
    zet <- fortify(OTU_inext_imp)

    Sdat$IZW_ID <- paste("Fox", Sdat$IZW_ID)

    zet <- merge(zet, Sdat, by.x="site", by.y="IZW_ID")

    alphaDivFox <- ggplot(zet, aes(x = x, y = y, colour = area, group=site)) +
        geom_line(data=subset(zet, method%in%"interpolated"),
                  lwd = 1.5, alpha=0.3) +
        geom_point(data=subset(zet, method%in%"observed"), shape = 21,
                   fill = "white", size=2.7, stroke = .8) + 
        geom_point(data=subset(zet, method%in%"observed"), shape = 21,
                   fill = "transparent", size=2.7, stroke = .8) + 
        scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
        scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
        ylab(paste0(deparse(substitute(output_string)), " diversity")) +
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
    

    alphaCompared <- ggplot(EstimatesAsy,
                            aes(area, Estimator, color=area,
                                fill = after_scale(lighten(color, .7)))) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021), fill = "white", size = 1.3, stroke = .7) +
        scale_y_continuous(name = NULL) +
        scale_x_discrete(name = NULL) +
        facet_wrap(~Diversity, scales="free_y")+
        scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
        scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
        theme(legend.position="none", axis.text.y = element_text(family = font_num))

    ## Now: the way Caro designed the analysis it distinguishes between
    ## Berlin and Brandenburg as a whole (and between male and female,
    ## maybe?). I'd call this level gamma-diversity!

    Pres1 <- as.data.frame(Counts>0)
    Pres1 <- apply(Pres1, 1, as.numeric)

    PresArea <- by(t(Pres1), Sdat$area, function (x) x)


    ## for iNext all sites have to have more than one species!
    Fox_inext_area1 <- iNEXT(list(Berlin=t(PresArea$Berlin),
                                  Brandenburg=t(PresArea$Brandenburg)), 
                             q=0, datatype = "incidence_raw")

    gammaDivFox <- ggiNEXT(Fox_inext_area1) +
        scale_colour_manual(values = c("#e7b800", "#2e6c61")) +
        scale_fill_manual(values = c("#e7b800", "#2e6c61")) +
        xlab("Number of sampled foxes") +
        ylab(paste0(deparse(substitute(output_string)), " diversity")) +
        ## need to repeat theme because it is overwritten by the wrapper
        theme_minimal(base_family = "Roboto", base_size = 12) +
        theme(legend.position="none", 
              axis.text = element_text(family = font_num),
              axis.title.x = element_text(margin = margin(t = 12)),
              axis.title.y = element_text(margin = margin(r = 12)),
              panel.grid.minor = element_blank(),
              plot.margin = margin(rep(12, 4)))

    ## beta diversity
    JaccPairsDist <- beta.pair(t(apply(indCounts>0, 1, as.numeric)),
                               index.family="jaccard")

    JaccGrups <- betadisper(JaccPairsDist[[3]], SdatHPres$area)

    message("\n Significance of the beta-diversity differences\n")
    print(anova(JaccGrups))
    message("\n")

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
    
    f4 <- paste0("figures/Diversity", output_string, ".pdf")
    ggsave(f4, width=13, height=9, device=cairo_pdf)

    ## we plot an return: 1. AsymptoticAlpaEstimates, 2. the beta Diversity Anova
    return(EstimatesAsy)
}


## HelmEstimateAsy <- getAllDiversity(
##     subset_taxa(PS, phylum %in% c("Nematoda", "Platyhelminthes")),
##     "Helminth")

## DietEstimateAsy <- getAllDiversity(
##     subset_taxa(PS, phylum %in% c("Annelida", "Arthropoda", "Chordata", "Mollusca")),
##     "Diet")


## ### Models for Shannon
##     EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
##         lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelShannonArea

##     summary(DivModelShannonArea)

##     EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
##         lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelShannonTree

##     summary(DivModelShannonTree)

##     EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
##         lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelShannonHum

##     summary(DivModelShannonHum)


##     EstimatesAsy %>% filter(Diversity %in% "Shannon diversity") %>%
##         lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelShannonImp

##     summary(DivModelShannonImp)

##     AIC(DivModelShannonArea, DivModelShannonTree, DivModelShannonImp, DivModelShannonHum)

## ### Models for Simpson

##     EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
##         lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelSimpsonArea

##     summary(DivModelSimpsonArea)

##     EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
##         lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelSimpsonTree

##     summary(DivModelSimpsonTree)

##     EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
##         lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelSimpsonHum

##     summary(DivModelSimpsonHum)


##     EstimatesAsy %>% filter(Diversity %in% "Simpson diversity") %>%
##         lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##            data=.) ->
##         DivModelSimpsonImp

##     summary(DivModelSimpsonImp)

##     AIC(DivModelSimpsonArea, DivModelSimpsonTree, DivModelSimpsonImp, DivModelSimpsonHum)


## ### Models for Species richness

##     EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
##         glm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
##             data=., family="poisson") ->
##         DivModelHillArea

##     summary(DivModelHillArea)

##     EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
##         glm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##             data=., family="poisson") ->
##         DivModelHillTree

##     summary(DivModelHillTree)

##     EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
##         glm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##             data=., family="poisson") ->
##         DivModelHillHum

##     summary(DivModelHillHum)


##     EstimatesAsy %>% filter(Diversity %in% "Species richness") %>%
##         glm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) + sex + age,
##             data=., family="poisson") ->
##         DivModelHillImp

##     summary(DivModelHillImp)

##     AIC(DivModelHillArea, DivModelHillTree, DivModelHillImp, DivModelHillHum)

##     f1 <- paste0("tables/HillDiv", output_string, ".html")

##     tab_model(DivModelHillArea, DivModelHillTree, DivModelHillImp, DivModelHillHum,
##               file=deparse(substitute(f1)), show.aic = TRUE)

##     f2 <- paste0("tables/SimpsonDiv", output_string, ".html")
    
##     tab_model(DivModelSimpsonArea, DivModelSimpsonTree, DivModelSimpsonImp,
##               DivModelSimpsonHum, file=deparse(substitute(f2)), show.aic = TRUE)

##     f3 <- paste0("tables/ShannonDiv", output_string, ".html")
    
##     tab_model(DivModelShannonArea, DivModelShannonTree, DivModelShannonImp,
##               DivModelShannonHum, file=deparse(substitute(f3)), show.aic = TRUE)





## WHY do we only look a helminths not acutally at protozoa
## (e.g. Coccidia) too? --> NOW see below they are also more diverse
## in Brandenburg

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




### Let's see whether there's a non-linear "ecotone" effect of any of
### the continuous environmental variables!

HelmEstimateAsy %>% filter(Diversity %in% "Species richness") %>%
    ggplot(aes(human_fpi_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

HelmEstimateAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    ggplot(aes(human_fpi_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

HelmEstimateAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    ggplot(aes(human_fpi_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

HelmEstimateAsy %>% filter(Diversity %in% "Species richness") %>%
    ggplot(aes(imperv_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

HelmEstimateAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    ggplot(aes(imperv_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

HelmEstimateAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    ggplot(aes(imperv_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()


HelmEstimateAsy %>% filter(Diversity %in% "Species richness") %>%
    ggplot(aes(tree_cover_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

HelmEstimateAsy %>% filter(Diversity %in% "Shannon diversity") %>%
    ggplot(aes(tree_cover_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

HelmEstimateAsy %>% filter(Diversity %in% "Simpson diversity") %>%
    ggplot(aes(tree_cover_1000m, Estimator, color=area)) +
    geom_point() +
    stat_smooth()

