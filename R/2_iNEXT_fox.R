library(phyloseq)
library(iNEXT)
library(vegan)
library(betapart)
library(tidyverse)
library(patchwork)
library(colorspace)
library(MASS)
library(broom)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

## You can leave this at true here, as the "bioinformatics" script has
## finer controls within itself for now. Re-executing this for now
## adds only the (potentially new) sample data (created in
## 0_Extract_Einvir_Covariates.R) to the PS object.
recomputeBioinfo <- FALSE

if(!exists("PSG")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSG <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
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
getAllDiversity <- function (ps, output_string) {
    Counts <- as(otu_table(ps), "matrix")
    Counts <- as.data.frame(t(Counts))
    rownames(Counts) <- tax_table(ps)[, "genus"]
    
    ## For inext diversity analysis we need to keep only samples with
    ## at least two species
    indCounts <- Counts[, colSums(Counts>0)>1]    
    
    ## Sample data in data frame
    Sdat <- as.data.frame(as(sample_data(ps), "matrix"))

    ## same for the data
    SdatHPres <- Sdat[rowSums(Counts>0)>1, ]


    
    OTU_inext_imp <- iNEXT(indCounts, datatype = "abundance", q = 0)

    ## UGLY set theme again within function
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
    
   ### Now the plot gets a bit messy redo by hand
    zet <- fortify(OTU_inext_imp)

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

    ##  get the the asymptotic diversity estimates
    EstimatesAsy <- OTU_inext_imp$AsyEst

    ## But iNext loses the rownames and all extra information, add back the IDs
    ## only the sites "present" with more than 2 taxa 
    EstimatesAsy$IZW_ID <- SdatHPres[rep(1:nrow(SdatHPres), each=3), "IZW_ID"]
    
    ## now also add back the pure observed diversity for all the excluded samples
    SdatObs <- merge(data.frame(rawObs=rowSums(Counts>0), IZW_ID=rownames(Counts)), Sdat)

    EstimatesAsy <- merge(EstimatesAsy, SdatObs, by="IZW_ID", all.y=TRUE)

    alphaCompared <- ggplot(EstimatesAsy,
                            aes(area, Estimator, color=area,
                                fill = after_scale(lighten(color, .7)))) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021),
                   fill = "white", size = 1.3, stroke = .7) +
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
        geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021),
                   fill = "white", size = 2, stroke = .7) +
        ##geom_point(position = position_jitter(width = .3, seed = 2021)) +
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

    ## we plot an return AsymptoticAlphaEstimates
    return(EstimatesAsy)
}


HelmEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"Helminth"),
                                   "Helminth")

BacterialEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"Microbiome"),
                                        "Microbiome")

DietEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"Diet"),
                                   "Diet")

## ### Models for all diversity indices HELMINTHS
HelmEstimateAsy %>% filter(!Diversity %in% "Species richness") %>% group_by(Diversity) %>%
    drop_na() %>%
    do(modelArea = lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
                      data=.),
       modelImperv = lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) +
                            sex + age,
                        data=.),
       modelTree = lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) +
                          sex + age,
                      data=.),
       modelHFPI = lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                          sex + age,
                      data=.)
       ) -> lmHelm

HelmEstimateAsy %>% filter(Diversity %in% "Species richness") %>% group_by(Diversity) %>%
    drop_na() %>%
    do(modelArea = glm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
                      data=., family="poisson"),
       modelImperv = glm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) +
                             sex + age,
                         data=., family="poisson"),
       modelTree = glm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) +
                           sex + age,
                       data=., family="poisson"),
       modelHFPI = glm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                           sex + age,
                       data=., family="poisson")
       ) -> glmHelm

HelmModels<- rbind(lmHelm, glmHelm)


HelmModels %>%
    pivot_longer(!Diversity, names_to = "predictor", values_to="model") %>%
    mutate(tidied = map(model, tidy),
           glanced = map(model, glance)) %>%
    mutate(envPvals = unlist(map(tidied, ~ dplyr::select(.x[2,], p.value)))) ->
    HelmModels

HelmModels

rawGLM <- glm.nb(rawObs ~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                     sex + age,
              data=HelmEstimateAsy[!duplicated(HelmEstimateAsy$IZW_ID), ])

summary(rawGLM)

rawGLM <- glm.nb(rawObs ~ area + condition + I(as.numeric(weight_kg)) +
                     sex + age,
              data=HelmEstimateAsy[!duplicated(HelmEstimateAsy$IZW_ID), ])

summary(rawGLM)

## ### Models for all diversity indices DIET
DietEstimateAsy %>% filter(!Diversity %in% "Species richness") %>% group_by(Diversity) %>%
    drop_na() %>%
    do(modelArea = lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
                      data=.),
       modelImperv = lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) +
                            sex + age,
                        data=.),
       modelTree = lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) +
                          sex + age,
                      data=.),
       modelHFPI = lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                          sex + age,
                      data=.)
       ) -> lmDiet

DietEstimateAsy %>% filter(Diversity %in% "Species richness") %>% group_by(Diversity) %>%
    drop_na() %>%
    do(modelArea = glm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
                      data=., family="poisson"),
       modelImperv = glm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) +
                             sex + age,
                         data=., family="poisson"),
       modelTree = glm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) +
                           sex + age,
                       data=., family="poisson"),
       modelHFPI = glm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                           sex + age,
                       data=., family="poisson")
       ) -> glmDiet

DietModels<- rbind(lmDiet, glmDiet)

DietModels %>%
    pivot_longer(!Diversity, names_to = "predictor", values_to="model") %>%
    mutate(tidied = map(model, tidy),
           glanced = map(model, glance)) %>%
    mutate(envPvals = unlist(map(tidied, ~ dplyr::select(.x[2,], p.value)))) ->
    DietModels

DietModels

rawGLM <- glm.nb(rawObs ~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                     sex + age,
              data=DietEstimateAsy[!duplicated(DietEstimateAsy$IZW_ID), ])

summary(rawGLM)

rawGLM <- glm.nb(rawObs ~ area + condition + I(as.numeric(weight_kg)) +
                     sex + age,
              data=DietEstimateAsy[!duplicated(DietEstimateAsy$IZW_ID), ])

summary(rawGLM)


BacterialEstimateAsy %>% filter(!Diversity %in% "Species richness") %>% group_by(Diversity) %>%
    drop_na() %>%
    do(modelArea = lm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
                      data=.),
       modelImperv = lm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) +
                            sex + age,
                        data=.),
       modelTree = lm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) +
                          sex + age,
                      data=.),
       modelHFPI = lm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                          sex + age,
                      data=.)
       ) -> lmBacterial

BacterialEstimateAsy %>% filter(Diversity %in% "Species richness") %>% group_by(Diversity) %>%
    drop_na() %>%
    do(modelArea = glm(Estimator~ area + condition + I(as.numeric(weight_kg)) + sex + age,
                      data=., family="poisson"),
       modelImperv = glm(Estimator~ imperv_1000m + condition + I(as.numeric(weight_kg)) +
                             sex + age,
                         data=., family="poisson"),
       modelTree = glm(Estimator~ tree_cover_1000m + condition + I(as.numeric(weight_kg)) +
                           sex + age,
                       data=., family="poisson"),
       modelHFPI = glm(Estimator~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                           sex + age,
                       data=., family="poisson")
       ) -> glmBacterial

BacterialModels<- rbind(lmBacterial, glmBacterial)


BacterialModels %>%
    pivot_longer(!Diversity, names_to = "predictor", values_to="model") %>%
    mutate(tidied = map(model, tidy),
           glanced = map(model, glance)) %>%
    mutate(envPvals = unlist(map(tidied, ~ dplyr::select(.x[2,], p.value)))) ->
    BacterialModels

BacterialModels


rawGLM <- glm.nb(rawObs ~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                     sex + age,
                 data=BacterialEstimateAsy[!duplicated(BacterialEstimateAsy$IZW_ID), ])

summary(rawGLM)

rawGLM <- glm.nb(rawObs ~ area + condition + I(as.numeric(weight_kg)) +
                     sex + age,
                 data=BacterialEstimateAsy[!duplicated(BacterialEstimateAsy$IZW_ID), ])

summary(rawGLM)






## CHECKING APICOMPLEXA 


ApicoPEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"ApicoParasites"),
                                     "ApicoParasites")

rawGLM <- glm.nb(rawObs ~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                     sex + age,
                 data=ApicoPEstimateAsy[!duplicated(ApicoPEstimateAsy$IZW_ID), ])

summary(rawGLM)


rawGLM <- glm.nb(rawObs ~ area + condition + I(as.numeric(weight_kg)) +
                     sex + age,
                 data=ApicoPEstimateAsy[!duplicated(ApicoPEstimateAsy$IZW_ID), ])

summary(rawGLM)


ApicoXEstimateAsyy <- getAllDiversity(subset_taxa(PSG, category%in%"noClue" & 
                                                      phylum%in%"Apicomplexa"),
                                     "ApicoX")

rawGLM <- glm.nb(rawObs ~ human_fpi_1000m + condition + I(as.numeric(weight_kg)) +
                     sex + age,
                 data=ApicoXEstimateAsy[!duplicated(ApicoXEstimateAsy$IZW_ID), ])

summary(rawGLM)


rawGLM <- glm.nb(rawObs ~ area + condition + I(as.numeric(weight_kg)) +
                     sex + age,
                 data=ApicoXEstimateAsy[!duplicated(ApicoXEstimateAsy$IZW_ID), ])

summary(rawGLM)






## ApicoEstimateAsy <- getAllDiversity(
##     subset_taxa(PS, phylum %in% c("Apicomplexa")),
##     "Apicomplexa")

## problem is that there are more with ZERO AND ONE in BRANDENBURG

##  Significance of removed data:
##              MoreOne
##               FALSE TRUE
##   Berlin         30   79
##   Brandenburg     3   39

## 	Fisher's Exact Test for Count Data

## data:  table(Sdat$area, MoreOne = rowSums(Counts > 0) > 1)
## p-value = 0.00749
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##   1.388368 26.606115
## sample estimates:
## odds ratio 
##   4.895872 

### Therefore we here compute diversity estimates without iNEXT hill
### numbers

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


PSApicoG <- phyloseq::tax_glom(subset_taxa(PS, phylum%in%"Apicomplexa"), "genus")

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

cbind(estimate_richness(PSApicoG, measures=c("Observed", "Shannon")),
      sample_data(PSApicoG))  %>% 
    lm(Shannon~ area + condition + I(as.numeric(weight_kg)) + sex + age,
       data=.) ->
    ModDivApicoShannon

summary(ModDivApicoShannon)

library(MASS)
cbind(estimate_richness(PSApicoG, measures=c("Observed", "Shannon")),
      sample_data(PSApicoG))  %>% 
    glm.nb(Observed~ area + condition + I(as.numeric(weight_kg)) + sex + age,
        data=.) ->
    ModDivApicoObserved

summary(ModDivApicoObserved)

### Let's see whether there's a non-linear "ecotone" effect of any of
### the continuous environmental variables!

HelmEstimateAsy %>% dplyr::select(Diversity, Estimator, area,
                                  tree_cover_1000m, imperv_1000m, human_fpi_1000m) %>%
    pivot_longer(cols=contains("1000m")) %>%
    ggplot(aes(value, Estimator, color=area)) +
    geom_point() +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    facet_wrap(name~Diversity, scales = "free") + 
    stat_smooth() +
    geom_smooth(aes(value, Estimator), color="black") +
    ggtitle("Helminth diversity")
ggsave("figures/suppl/HelmDiv_Conti_Env.pdf", width=25, height=15, device=cairo_pdf)

DietEstimateAsy %>% dplyr::select(Diversity, Estimator, area,
                                  tree_cover_1000m, imperv_1000m, human_fpi_1000m) %>%
    pivot_longer(cols=contains("1000m")) %>%
    ggplot(aes(value, Estimator, color=area)) +
    geom_point() +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    facet_wrap(name~Diversity, scales = "free") + 
    stat_smooth() +
    geom_smooth(aes(value, Estimator), color="black") +
    ggtitle("Diet diversity")
ggsave("figures/suppl/DietDiv_Conti_Env.pdf", width=25, height=15, device=cairo_pdf)


#### Should do the same for bacteria!!!!




## ## Obtaining diet diversity as a predictor for later models
## ## Will add this to sample data!!!
DietEstimateAsy %>% dplyr::filter(Diversity%in%"Species richness") %>%
    dplyr::transmute(EstDietD = Estimator) %>% data.frame() %>%
    merge(., rowSums(otu_table(PSG.diet) > 0), by=0, all=TRUE) -> Div

colnames(Div)[colnames(Div)%in%"y"] <- "ObsDietD"

HelmEstimateAsy %>% dplyr::filter(Diversity%in%"Species richness") %>%
    dplyr::transmute(EstHelmD = Estimator) %>% data.frame()  %>%
        merge(., Div , by.x= 0, by.y="Row.names", all=TRUE) %>%
        merge(., rowSums(otu_table(PSG.helm) > 0), by.x="Row.names",
              by=0, all=TRUE) -> Div

colnames(Div)[colnames(Div)%in%"y"] <- "ObsHelmD"

BacterialEstimateAsy %>% dplyr::filter(Diversity%in%"Species richness") %>%
    dplyr::transmute(EstBacD = Estimator) %>% data.frame()  %>%
        merge(., Div , by.x= 0, by.y="Row.names", all=TRUE) %>%
        merge(., rowSums(otu_table(PSG.bac) > 0), by.x="Row.names",
              by=0, all=TRUE) -> Div
colnames(Div)[colnames(Div)%in%"y"] <- "ObsBacD"

cor(Div[, -1], use="pairwise.complete.obs")

### observed and Estimated diversities are so close that we can use
### observed to include low numers... now: add all this to sample-data!


#### TABLE 1 : WHATCHOUT!!! ###

## ## ## this is removed from the PS/PSG
## traits %>% dplyr::filter(BlastEvaluation%in%"Okay") %>%

## ## ## this is done on the talbe from tax_table data.
##     dplyr::select(zoonotic, fox.parasite, transmission.fox) ->
##         outHelm

### this is not needed anymore
## helmData <- as.data.frame(otu_table(PSG))
## colnames(helmData) <- unname(tax_table(PSG)[, "genus"])

## helmData%>%dplyr::select(any_of(rownames(traits))) %>% t() -> helmData

### thi will be done in 2_ from the tax_table data
## outHelm <- cbind(outHelm,
##                  foxesN=rowSums(helmData>0),
##                  prevalence=round(rowSums(helmData>0)/ncol(helmData) *100, 2),
##                  meanAbundance=round(rowMeans(helmData), 2))

## outHelm$meanIntensity <- round(apply(helmData, 1, function (x) mean(x[x>0])), 2)

## head(outHelm[order(outHelm$prevalence, decreasing=TRUE), ], n=50)

## table(colSums(helmData)==0) / ncol(helmData)
