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
library(ggeffects)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

source("./R/plot_setup.R")

## The "bioinformatics" script "1_Fox_general_MA.R" has finer controls
## within itself for now. But we can still read it's single final
## output object.

iNV <- packageVersion("iNEXT")
if (iNV < 3){
    stop(paste("We need iNEXT version < 3, your version is",
               iNV))
}
    


recomputeBioinfo <- FALSE
## recomputeBioinfo <- TRUE

if(!exists("PSG")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSG <- readRDS(file = "intermediate_data/PhyloSeqGenus.Rds")
    }
}

### WARNING: most of the functions in this script are more complicatd
### than necessary for the present (Dec 2022) manuscript. The
### functions allow a general assessment of diversity (q=0, 1 and 2)
### for different taxonomic (phyla) subsets. I've put all off this in
### two big functions. One runs iNEXT, compiles the data and generates
### plots (getAllDiversity) the other one (gimmeModels) runs all the
### statistical models on them


getAllDiversity <- function (ps, output_string, plot=FALSE) {
    Counts <- as(otu_table(ps), "matrix")
    Counts <- as.data.frame(t(Counts))
    rownames(Counts) <- make.unique(tax_table(ps)[, "genus"])
    
    ## For inext diversity analysis we need to keep only samples with
    ## at least two species
    indCounts <- Counts[, colSums(Counts > 0) > 1]    
    
    ## Sample data in data frame
    Sdat <- as.data.frame(as(sample_data(ps), "matrix"))

    ## same for the data
    SdatHPres <- Sdat[rowSums(Counts > 0) > 1, ]

    OTU_inext_imp <- iNEXT(indCounts, datatype = "abundance", q = 0)
    
    zet <- fortify(OTU_inext_imp)

    ## now add the sample data and also the zero and 1 "Assemblages" (foxes
    ## back into this, by "all.y")
    
    zet <- merge(zet, Sdat, by.x = "Assemblage", by.y = 0, all.y = TRUE) ## TODO: doesn't work for me, no col assemblage

    
    ## Now the plot gets a bit messy redo by hand
    alphaDivFox <- 
      ggplot(zet, aes(x = x, y = y, colour = area, group=Assemblage)) +
      geom_line(data = subset(zet, Method %in% "Rarefaction"),
                lwd = 1.5, alpha = 0.3) +
      geom_point(data = subset(zet, Method %in% "Observed"), shape = 21,
                 fill = "white", size = 2.7, stroke = .8) + 
      geom_point(data = subset(zet, Method %in% "Observed"), shape = 21,
                 fill = "transparent", size = 2.7, stroke = .8) + 
      coord_cartesian(clip = "off") +
      scale_colour_manual(values = colors_regions, guide = "none") +
      scale_x_continuous(labels = scales::label_comma(), expand = c(0, 0), limits = c(0, NA)) +
      scale_y_continuous(breaks = 2:10, expand = c(0, 0)) +
      labs(x = "Number of sequence reads",
           #y = paste0(deparse(substitute(output_string)), " diversity"))
           y = paste0('"', output_string, '" diversity')) +
      theme(plot.margin = margin(rep(20, 4)))

    ##  get the the asymptotic diversity estimates
    EstimatesAsy <- OTU_inext_imp$AsyEst

    ## ## ## observed diversity of 1 and 0 div samples
    NullOne <- as.data.frame(cbind(Observed =  colSums(Counts[, colSums(Counts > 0) < 2] > 0),
                                   Estimator = colSums(Counts[, colSums(Counts > 0) < 2] > 0)))
    NullOne$Assemblage <- rownames(NullOne)
    NullOne <- NullOne[rep(1:nrow(NullOne), each = 3),]
    NullOne$Diversity <-  c("Species richness",
                            "Shannon diversity",
                            "Simpson diversity")
    NullOne$s.e. <- NullOne$LCL <- NullOne$UCL <- NA

    EstimatesAsy <- rbind(EstimatesAsy, NullOne[, colnames(EstimatesAsy)])
    
    ## now  add back the pure observed diversity for all the excluded samples
    EstimatesAsy <- merge(EstimatesAsy, Sdat, by.x = "Assemblage", by.y = 0)
    
    alphaCompared <- 
      ggplot(EstimatesAsy,
             aes(area, Estimator, color = area,
                 fill = after_scale(lighten(color, .7)))) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021),
                 fill = "white", size = 1.3, stroke = .7) +
      facet_wrap(~Diversity, scales = "free_y") +
      coord_cartesian(clip = "off") +
      scale_x_discrete(expand = c(.2, .3), name = NULL) +
      scale_y_continuous(expand = c(.005, .005), name = NULL) +
      scale_colour_manual(values = colors_regions, guide = "none") +
      theme(panel.grid.major.x = element_blank(),
            plot.margin = margin(rep(20, 4)))

    ## Now: the way Caro designed the analysis it distinguishes between
    ## Berlin and Brandenburg as a whole (and between male and female,
    ## maybe?). I'd call this level gamma-diversity!

    Pres1 <- as.data.frame(t(Counts) > 0)
    Pres1 <- apply(Pres1, 1, as.numeric)

    PresArea <- by(t(Pres1), Sdat$area, function (x) x)


    ## for iNext all assemblages have to have more than one species!
    Fox_inext_area1 <- iNEXT(list(Berlin = t(PresArea$Berlin),
                                  Brandenburg = t(PresArea$Brandenburg)), 
                             q = 0, datatype = "incidence_raw")

    gammaDivFox <- 
      ggiNEXT(Fox_inext_area1) +
      coord_cartesian(expand = FALSE, clip = "off") +
      scale_x_continuous(breaks = c(1, seq(50, 200, by = 50))) +
      scale_y_continuous(breaks = seq(0, 20, by = 5)) + 
      scale_colour_manual(values = colors_regions) +
      scale_fill_manual(values = colors_regions) +
      labs(x = "Number of sampled foxes",
           #y = paste0(deparse(substitute(output_string)), " diversity")) +
           y = paste0('"', output_string, '" diversity')) +
      ## need to repeat theme because it is overwritten by the wrapper
      theme_custom() +
      theme(legend.position = "none",
            plot.margin = margin(rep(20, 4)))

    ## beta diversity
    JaccPairsDist <- beta.pair(t(apply(Counts > 0, 2, as.numeric)),
                               index.family = "jaccard")

    JaccGrups <- betadisper(JaccPairsDist[[3]], Sdat$area)

    message("\n Significance of the beta-diversity differences\n")
    print(anova(JaccGrups))
    message("\n")

    betaDivJac <- 
      data.frame(distances = JaccGrups$distances,
                 area = JaccGrups$group) %>%
      ggplot(aes(area, distances, color = area, fill = after_scale(lighten(color, .7)))) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(shape = 21, position = position_jitter(width = .25, seed = 2021),
                 fill = "white", size = 2, stroke = .7) +
      ##geom_point(position = position_jitter(width = .3, seed = 2021)) +
      coord_cartesian(clip = "off") +
      scale_x_discrete(expand = c(.2, .3), name = NULL) +
      scale_y_continuous(expand = c(0, 0), name = "Distance to area centroid") +
      scale_colour_manual(values = colors_regions, name = "Study area:") +
      guides(fill = guide_legend(title.position = "top", title.hjust = .5),
             color = guide_legend(title.position = "top", title.hjust = .5)) +
      theme(legend.position = "top",  ## TODO: why do we keep the legend here? Why boxplots as legend keys?
            panel.grid.major.x = element_blank(),
            plot.margin = margin(rep(20, 4)))

    # betaDivJacMulti <- 
    #   wrap_elements(full =
    #                   ~(plot(JaccGrups, col = colors_regions, main = "",
    #                          label = FALSE, sub = "")))

    dat <- data.frame(vec_x = JaccGrups$vectors[,1],
                      vec_y = JaccGrups$vectors[,2],
                      region = JaccGrups$group)
    
    dat$centr_x <- ifelse(dat$region == "Berlin", JaccGrups$centroids[1,1], JaccGrups$centroids[2,1])
    dat$centr_y <- ifelse(dat$region == "Berlin", JaccGrups$centroids[1,2], JaccGrups$centroids[2,2])
    
    hull_dat <- 
      dat %>%
      group_by(region) %>% 
      slice(chull(vec_x, vec_y))
    
    betaDivJacMulti <- 
      ggplot(dat, aes(vec_x, vec_y, color = region, fill = region, shape = region)) +
      #ggConvexHull::geom_convexhull(alpha = .2, lwd = 1.5) +
      geom_polygon(data = hull_dat, alpha = .2, lwd = 1.5) +
      #geom_segment(aes(xend = centr_x, yend = centr_y), color = "grey65", lwd = .4) +
      geom_segment(aes(xend = centr_x, yend = centr_y), color = "white", lwd = .4) +
      geom_segment(aes(xend = centr_x, yend = centr_y), alpha = .4, lwd = .4) +
      geom_point(aes(centr_x, centr_y), shape = 21, color = "white", stroke = .7, size = 4, shape = 16) +
      geom_point(fill = "white", size = 2.5, stroke = .8) +
      coord_cartesian(clip = "off") +
      scale_x_continuous(expand = c(.005, .005), name = "PCoA 1", 
                         breaks = -3:3 * .2) +
      scale_y_continuous(expand = c(.005, .005), name = "PCoA 2",
                         breaks = -3:3 * .2) +
      scale_shape_manual(values = c(21, 23), guide = "none") +
      scale_colour_manual(values = colors_regions, guide = "none") +
      scale_fill_manual(values = colors_regions, guide = "none") +
      theme(plot.margin = margin(rep(20, 4)))
      
    
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
    plot_annotation(tag_levels = 'a', theme = theme(legend.title = element_text(hjust = .5)))
    
    if(plot %in% c("pdf", "png")){
        if(plot %in% "pdf"){
            f4 <- paste0("figures/suppl/Diversity", output_string, ".pdf")
            ggsave(f4, width = 16.5, height = 10.5, device = cairo_pdf)
        }
        if(plot %in% "png"){
            f4 <- paste0("figures/suppl/Diversity", output_string, ".png")
            ggsave(f4, width = 16.5, height = 10.5, bg = "white", dpi = 600)
        }
    } else{
        message("to produce plots give \"pdf\" or \"png\" as argument")
    }
    ## we return AsymptoticAlphaEstimates
    return(EstimatesAsy)
}



gimmeModels <- function (EA){
    ## all models on the Asymptotic estimate tables from the function
    ## above (wrapping iNEXT)
    EA %>% 
      mutate_at(c("weight_kg","tree_cover_1000m" ,"imperv_1000m", "human_fpi_1000m", 
                  "DNAng.ul", "DNA260.230", "DNA260.280"), as.numeric) %>%
      mutate(REstimator = round(Estimator)) -> EA

    EA %>%
        filter(!Diversity %in% "Species richness") %>% group_by(Diversity) %>%
        do(modelArea = lm(Estimator~ area + condition + weight_kg + sex + age +
                              season +
                              DNAng.ul + DNA260.230 + DNA260.280, 
                          data = .),
           modelImperv = lm(Estimator~ imperv_1000m + condition + weight_kg +
                                sex + age  + season +
                                DNAng.ul + DNA260.230 + DNA260.280, 
                            data = .),
           modelTree = lm(Estimator~ tree_cover_1000m + condition + weight_kg +
                              sex + age + season +
                              DNAng.ul + DNA260.230 + DNA260.280,
                          data = .),
           modelHFPI = lm(Estimator~ human_fpi_1000m + condition + weight_kg +
                              sex + age + season +
                              DNAng.ul + DNA260.230 + DNA260.280,
                          data = .)
           ) -> lmEA

    EA %>% filter(Diversity %in% "Species richness") %>% group_by(Diversity) %>%
        do(modelArea = glm(REstimator ~ area + condition + weight_kg +
                               sex + age  + season
                           + DNAng.ul + DNA260.230 + DNA260.280,
                           data = ., family = "poisson"),
           modelImperv = glm(REstimator ~ imperv_1000m + condition + weight_kg +
                                 sex + age + season +
                                 DNAng.ul + DNA260.230 + DNA260.280,
                             data = ., family = "poisson"),
           modelTree = glm(REstimator ~ tree_cover_1000m + condition + weight_kg +
                               sex + age + season + 
                               DNAng.ul + DNA260.230 + DNA260.280,
                           data = ., family = "poisson"),
           modelHFPI = glm(REstimator ~ human_fpi_1000m + condition + weight_kg +
                               sex + age + season +
                               DNAng.ul + DNA260.230 + DNA260.280,
                           data = ., family = "poisson")
           ) -> glmEA

    EAModels <- rbind(lmEA, glmEA)

    EAModels %>%
        pivot_longer(!Diversity, names_to = "predictor", values_to = "model") %>%
        mutate(tidied = map(model, tidy),
               glanced = map(model, glance)) %>%
        mutate(envPvals = unlist(map(tidied, ~ dplyr::select(.x[2,], p.value)))) %>%
        arrange(factor(Diversity,
                       levels = c("Species richness",
                                  "Shannon diversity",
                                  "Simpson diversity")))
}




HelmEstimateAsy <- getAllDiversity(subset_taxa(PSG, category %in% c("Helminth")),
                                   "Helminth", plot="png")

HelmModels <- gimmeModels(HelmEstimateAsy)

AreaRich <- HelmModels %>%
  filter(Diversity %in% "Species richness")%>%
  dplyr::select(model) %>% .[["model"]] %>% .[[1]]

plot(ggeffect(AreaRich, terms = c("weight_kg", "season", "area")), rawdata = TRUE) +
  labs(x = "Weight (kg)", y = "Species richness (Hill number q=0)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  scale_color_manual(values = colors_seasons, 
                     labels = c("Spring", "Summer + Autumn", "Winter"),
                     name = "Season:") +
  scale_fill_manual(values = colors_seasons, guide = "none") +
  ## need to repeat theme because it is overwritten by the wrapper
  theme_custom() +
  theme(plot.margin = margin(.5, .5, 2, 0))

ggsave("figures/Div_Model.png", width = 11, height = 7, bg = "white", dpi = 600)



HelmEstimateAsy %>% 
  mutate_at(c("weight_kg","tree_cover_1000m" ,"imperv_1000m", "human_fpi_1000m", 
              "DNAng.ul", "DNA260.230", "DNA260.280"), as.numeric) %>%
    filter(Diversity %in% "Species richness")%>%
    mutate(REstimator = round(Estimator)) -> EA


ContiRich <- glm(REstimator ~ tree_cover_1000m + 
                     imperv_1000m +  human_fpi_1000m +
                     weight_kg + sex + age  + season +
                     + condition + 
                     DNAng.ul + DNA260.230 + DNA260.280,
                 data = EA, family = "poisson")

### The problem is the Continous model has a better AIC but it is not
### any better at explaining anything. 

## Parasitic Helminths not more diverse in brandenburg

DietEstimateAsy <- getAllDiversity(subset_taxa(PSG, category %in% "Diet"),
                                   "Diet", plot=FALSE)

gimmeModels(DietEstimateAsy)

###  -> Diet clearly more diverse in Brandenburg

WormDietEA <- getAllDiversity(subset_taxa(PSG, category %in% "Diet"&
                                          phylum %in% c("Platyhelminthes", "Nematoda")),
                                   "WormDiet", plot=FALSE)

gimmeModels(WormDietEA)

### -> "Diet-worms" more diverse in Brandenburg


BacterialEstimateAsy <- getAllDiversity(subset_taxa(PSG, category %in% "Microbiome"),
                                        "Microbiome")

gimmeModels(BacterialEstimateAsy)

### -> Bacterial microbiome clearly more diverse in Brandenburg

FungalEstimateAsy <- getAllDiversity(subset_taxa(PSG, category %in% "FungalMicrobiome"),
                                   "FungalMicrobiome", plot=FALSE)

gimmeModels(FungalEstimateAsy)

### -> Fungal microbiome more diverse in Brandenburg


ApicoPEstimateAsy <- getAllDiversity(subset_taxa(PSG, category %in% "ApicoParasites"),
                                     "ApicoParasites", plot=FALSE)

gimmeModels(ApicoPEstimateAsy)


ApicoEnvirEstimateAsy <- getAllDiversity(subset_taxa(PSG, !category %in% "ApicoParasites" &
                                                          phylum %in% "Apicomplexa"), 
                                         "EnvironmApicomplexa", plot=FALSE)

gimmeModels(ApicoEnvirEstimateAsy)

### -> Environmental Apicomplexans are more diverse in Brandenburg,
### -> parasitic Apicomplexans not really. 


### Let's see whether there's a non-linear "ecotone" effect of any of
### the continuous environmental variables!

## ## Diet diversity along a gradient 
## DietEstimateAsy %>% dplyr::select(Diversity, Estimator, area,
##                                   tree_cover_1000m, imperv_1000m, human_fpi_1000m) %>%
##     pivot_longer(cols = contains("1000m")) %>%
##     ggplot(aes(value, Estimator, color = area)) +
##     geom_point() +
##     scale_colour_manual(values = colors_regions, name = "Study area:") +
##     scale_fill_manual(values = colors_regions, name = "Study area:") +
##     facet_wrap(name~Diversity, scales = "free") + 
##     stat_smooth() +
##     geom_smooth(aes(value, Estimator), color = "black") +
##     ggtitle("Diet diversity")

## ggsave("figures/suppl/DietDiv_Conti_Env.pdf", width = 25, height = 15, device = cairo_pdf)

### As a reviewer figure only the q0

HelmEstimateAsy %>% dplyr::select(Diversity, Estimator, area,
                                  tree_cover_1000m, imperv_1000m, human_fpi_1000m) %>%
    dplyr::filter(Diversity %in%"Species richness") %>%
    pivot_longer(cols = contains("1000m")) %>%
    ggplot(aes(value, Estimator, color = area)) +
    geom_point() +
    scale_colour_manual(values = colors_regions, name = "Study area:") +
    scale_fill_manual(values = colors_regions, name = "Study area:") +
    facet_wrap(~name, nrow=1, scales = "free") + 
    stat_smooth() +
    geom_smooth(aes(value, Estimator), color = "black") +
    scale_y_continuous("Helminth species richness (Hill number q=0)")

ggsave("figures/suppl/HelmRich_Conti_Env.png", 
       width = 21, height = 7, bg = "white", dpi = 600)



## ## Obtaining diet diversity as a predictor for later models
## ## Decided to remove this for now from the manuscript (December 2022). 

## DietEstimateAsy %>%
##     dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
##     pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix = "Diet ") ->
##     DietDiversity

## HelmEstimateAsy %>%
##     dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
##     pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix = "Helm ") ->
##     HelminthDiversity

## ApicoPEstimateAsy %>%
##     dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
##     pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix = "ApicoP ") ->
##     ApicoPDiversity

## BacterialEstimateAsy %>%
##     dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
##     pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix = "BacM ") ->
##     BacterialDiversity

## FungalEstimateAsy %>%
##     dplyr::select(-c(Observed, LCL, s.e.,UCL)) %>%
##     pivot_wider(names_from = Diversity, values_from = Estimator, names_prefix = "FunM ") ->
##     FungalDiversity

## ## only append the diversity statustics to the sample data if it's not
## ## alredy there
## if(!"Helm_Species_richness" %in% colnames(sample_data(PSG))){
##     Reduce(merge, list(DietDiversity, HelminthDiversity,
##                        ApicoPDiversity, BacterialDiversity,
##                        FungalDiversity)) %>%
##         rename_with(~ gsub(" ", "_", .x, fixed = TRUE)) ->
##         AllDiv
##     ## all(AllDiv$Assemblage == rownames(sample_data(PSG)))
##     rownames(AllDiv) <- AllDiv$Assemblage
##     sample_data(PSG) <- AllDiv
##     ### write new phylseq object only if it's now in sample data 
##     if("Helm_Species_richness" %in% colnames(sample_data(PSG))){
##         saveRDS(PSG, file = "intermediate_data/PhyloSeqGenus.Rds")
##     }
## }
    
HelmModels %>%
    filter(Diversity %in% "Species richness")%>%
    dplyr::select(model) %>% .[["model"]] %>%
    append(. , list(ContiRich)) %>%
    stargazer(type = "html", out="./tables/HelmDiversityq0.html",
              dep.var.caption = "Species Richness (q=0)",
              dep.var.labels  = c("", ""),
              intercept.top=TRUE, intercept.bottom=FALSE)

HelmModels %>%
    filter(Diversity %in% c("Shannon diversity", "Simpson diversity"))%>%
    dplyr::select(model) %>% .[["model"]] %>%
    stargazer(type = "html", out="./tables/suppl/HelmDiversityq1aq2.html",
              column.labels = rep(c("Shannon diversity (q=1)",
                                    "Simpson diversity (q=2)"),
                                  each=4),
              dep.var.caption = "",
              dep.var.labels  = "",
              intercept.top=TRUE, intercept.bottom=FALSE
              )

## No pattern here (in the non-infected)
table(sample_data(PSG)[
    rowSums(otu_table(subset_taxa(PSG, category %in% "Helminth")))==0,
    c("sex", "age", "season")])


table(sample_data(PSG)[
    rowSums(otu_table(subset_taxa(PSG, category %in% "Helminth")))==0,
    c("sex", "area", "season")])

