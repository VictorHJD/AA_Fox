library(phyloseq)
library(iNEXT)
library(vegan)
library(betapart)
library(tidyverse)
library(patchwork)
library(colorspace)
library(MASS)
library(broom)
library(stargazer)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

## The "bioinformatics" script "1_Fox_general_MA.R" has finer controls
## within itself for now. But we can still read it's single final
## output object.
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
    rownames(Counts) <- make.unique(tax_table(ps)[, "genus"])
    
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
    
    zet <- fortify(OTU_inext_imp)

    ## now add the sample data and also the zero and 1 "sites" (foxes
    ## back into this, by "all.y")

    zet <- merge(zet, Sdat, by.x="site", by.y=0, all.y=TRUE)

    ## Now the plot gets a bit messy redo by hand
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

    ## ## ## observed diverstiy of 1 and 0 div samples
    NullOne <- as.data.frame(cbind(Observed=colSums(Counts[, colSums(Counts>0)<2]>0),
                                   Estimator=colSums(Counts[, colSums(Counts>0)<2]>0)))
    NullOne$Site <- rownames(NullOne)
    NullOne  <- NullOne[rep(1:nrow(NullOne), each=3),]
    NullOne$Diversity <-  c("Species richness",
                            "Shannon diversity",
                            "Simpson diversity")
    NullOne$s.e. <- NullOne$LCL <- NullOne$UCL <- NA

    EstimatesAsy <- rbind(EstimatesAsy, NullOne[, colnames(EstimatesAsy)])
    
    ## now  add back the pure observed diversity for all the excluded samples
    EstimatesAsy <- merge(EstimatesAsy, Sdat, by.x="Site", by.y=0)
    
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

    Pres1 <- as.data.frame(t(Counts)>0)
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
    JaccPairsDist <- beta.pair(t(apply(Counts>0, 2, as.numeric)),
                               index.family="jaccard")

    JaccGrups <- betadisper(JaccPairsDist[[3]], Sdat$area)

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



gimmeModels <- function (EA){
    ## all models on the Asymptotic estimate tables from the function
    ## above (wrapping iNEXT)
    EA %>% mutate_at(c("weight_kg","tree_cover_1000m" ,"imperv_1000m",
                       "human_fpi_1000m"), as.numeric)%>%
    mutate(REstimator = round(Estimator)) -> EA

    EA %>%
        filter(!Diversity %in% "Species richness") %>% group_by(Diversity) %>%
        do(modelArea = lm(Estimator~ area + condition + weight_kg + sex + age +
                              season + year, data=.),
           modelImperv = lm(Estimator~ imperv_1000m + condition + weight_kg +
                                sex + age  + season + year, 
                            data=.),
           modelTree = lm(Estimator~ tree_cover_1000m + condition + weight_kg +
                              sex + age + season + year,
                          data=.),
           modelHFPI = lm(Estimator~ human_fpi_1000m + condition + weight_kg +
                              sex + age + season + year,
                          data=.)
           ) -> lmEA

    EA %>% filter(Diversity %in% "Species richness") %>% group_by(Diversity) %>%
        do(modelArea = glm(REstimator ~ area + condition + weight_kg +
                               sex + age  + season + year,
                           data=., family="poisson"),
           modelImperv = glm(REstimator ~ imperv_1000m + condition + weight_kg +
                                 sex + age + season + year,
                             data=., family="poisson"),
           modelTree = glm(REstimator ~ tree_cover_1000m + condition + weight_kg +
                               sex + age + season + year,
                           data=., family="poisson"),
           modelHFPI = glm(REstimator ~ human_fpi_1000m + condition + weight_kg +
                               sex + age + season + year,
                           data=., family="poisson")
           ) -> glmEA

    EAModels<- rbind(lmEA, glmEA)

    EAModels %>%
        pivot_longer(!Diversity, names_to = "predictor", values_to="model") %>%
        mutate(tidied = map(model, tidy),
               glanced = map(model, glance)) %>%
        mutate(envPvals = unlist(map(tidied, ~ dplyr::select(.x[2,], p.value))))
}




HelmEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"Helminth"),
                                   "Helminth")

gimmeModels(HelmEstimateAsy)

## Parasitic Helminths not more diverse in brandenburg

DietEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"Diet"),
                                   "Diet")

gimmeModels(DietEstimateAsy)

###  -> Diet clearly more diverse in Brandenburg

WormDietEA <- getAllDiversity(subset_taxa(PSG, category%in%"Diet"&
                                          phylum%in%c("Platyhelminthes", "Nematoda")),
                                   "WormDiet")

gimmeModels(WormDietEA)

### -> "Diet-worms" more diverse in Brandenburg


BacterialEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"Microbiome"),
                                        "Microbiome")

gimmeModels(BacterialEstimateAsy)

### -> Bacterial microbiome clearly more diverse in Brandenburg

FungalEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"FungalMicrobiome"),
                                   "FungalMicrobiome")

gimmeModels(FungalEstimateAsy)

### -> Fungal microbiome more diverse in Brandenburg


ApicoPEstimateAsy <- getAllDiversity(subset_taxa(PSG, category%in%"ApicoParasites"),
                                     "ApicoParasites")

gimmeModels(ApicoPEstimateAsy)


ApicoEnvirEstimateAsy <- getAllDiversity(subset_taxa(PSG, !category%in%"ApicoParasites"&
                                                          phylum%in%"Apicomplexa"), 
                                         "EnvironmApicomplexa")

gimmeModels(ApicoEnvirEstimateAsy)

### -> Environmental Apicomplexans are more diverse in Brandenburg,
### -> parasitic Apicomplexans not really. 


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


## ## Obtaining diet diversity as a predictor for later models (could
## ## have done this before the modelling above for cleaner code, but
## ## now move on ;-))
DietEstimateAsy %>%
    dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
    pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix="Diet ") ->
    DietDiversity

HelmEstimateAsy %>%
    dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
    pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix="Helm ") ->
    HelminthDiversity

ApicoPEstimateAsy %>%
    dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
    pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix="ApicoP ") ->
    ApicoPDiversity

BacterialEstimateAsy %>%
    dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
    pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix="BacM ") ->
    BacterialDiversity

FungalEstimateAsy %>%
    dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
    pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix="FunM ") ->
    FungalDiversity

## only append the diversity statustics to the sample data if it's not
## alredy there
if(!"Helm_Species_richness"%in%colnames(sample_data(PSG))){
    Reduce(merge, list(DietDiversity, HelminthDiversity,
                       ApicoPDiversity, BacterialDiversity,
                       FungalDiversity)) %>%
        rename_with(~ gsub(" ", "_", .x, fixed = TRUE)) ->
        AllDiv
    ## all(AllDiv$Site == rownames(sample_data(PSG)))
    rownames(AllDiv) <- AllDiv$Site
    sample_data(PSG) <- AllDiv
    ### write new phylseq object only if it's now in sample data 
    if("Helm_Species_richness"%in%colnames(sample_data(PSG))){
        saveRDS(PSG, file="intermediate_data/PhyloSeqGenus.Rds")
    }
}
    
gimmeModels(HelmEstimateAsy) %>%
    ## filter(predictor%in%"modelArea")%>%
    dplyr::select(model) %>% .[["model"]]%>%
    stargazer(type="html", out="./tables/HelmDiversity.html")

