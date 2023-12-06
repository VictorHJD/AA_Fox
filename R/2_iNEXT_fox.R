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
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

## extrafont::font_import() ## only run once
extrafont::loadfonts(device = "all") ## run every time

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

if(!exists("PSGHelm")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSGHelm <- readRDS(file = "intermediate_data/PhyloSeqGenus.Rds")
    }
}


getAllDiversity <- function (ps) {
    Counts <- as.data.frame(t(unclass(otu_table(ps))))
    rownames(Counts) <- make.unique(tax_table(ps)[, "genus"])
    
    ## For inext diversity analysis we need to keep only samples with
    ## at least two species
    indCounts <- Counts[, colSums(Counts > 0) > 1]    
        
    ## Sample data in data frame
    Sdat <- as.data.frame(unclass(sample_data(ps)))

    ## same for the data
    SdatHPres <- Sdat[colSums(Counts > 0) > 1, ]
    rownames(SdatHPres) <- SdatHPres$IDb
    
    ## a Hack including othr reads as an additional animal
    ## (assemblage), which I then remove

    if(!all(rownames(SdatHPres)== colnames(indCounts))){
        stop("Sample data and counts do not match")
    }
    
    indCounts <- rbind(indCounts, ALLSEQ = SdatHPres$nSeq - 
                                      colSums(indCounts))
    
    OTU_inext_imp <- iNEXT(indCounts, datatype = "abundance", q = 0)

    
    zet <- fortify(OTU_inext_imp)
    ## remove estimates for the "ALLSEQ" pseudo sample
    zet <- zet[!zet$Assemblage%in%"ALLSEQ", ]

    ## now add the sample data and also the zero and 1 "Assemblages" (foxes
    ## back into this, by "all.y")
    zet <- merge(zet, Sdat, by.x = "Assemblage", by.y = "IDb", all.y = TRUE) 
    
    ## Now the plot gets a bit messy redo by hand
    ##  get the the asymptotic diversity estimates
    EstimatesAsy <- OTU_inext_imp$AsyEst
    EstimatesAsy <- EstimatesAsy[!EstimatesAsy$Assemblage%in%"ALLSEQ", ]
    
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
    EstimatesAsy <- merge(EstimatesAsy, Sdat, by.x = "Assemblage", by.y = "IDb")
    return(list(EstimatesAsy, zet))
}


gimmeModels <- function (EA){
    ## all models on the Asymptotic estimate tables from the function
    ## above (wrapping iNEXT)
    EA %>% 
        mutate_at(c("weight_kg","tree_cover_1000m" ,"imperv_1000m", "human_fpi_1000m"),
                  as.numeric) %>%
        mutate(REstimator = round(Estimator)) -> EA

    EA %>%
        filter(Diversity %in% "Species richness") %>% 
        do(modelArea = glm(Estimator~ area + weight_kg + sex +
                              season, 
                          data = ., family = "poisson"),
           modelImperv = lm(Estimator~ imperv_1000m + weight_kg +
                                sex  + season, 
                            data = .),
           modelTree = lm(Estimator~ tree_cover_1000m + weight_kg +
                              sex + season,
                          data = .),
           modelHFPI = lm(Estimator~ human_fpi_1000m + weight_kg +
                              sex + season,
                          data = .)
           ) -> lmEA
    lmEA
}
  


HelmEstimateAsy <- getAllDiversity(PSGHelm)
AsyEst <- HelmEstimateAsy[[1]]
zet <- HelmEstimateAsy[[2]]

models <- gimmeModels(AsyEst)

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
    scale_x_continuous(labels = scales::label_comma(),
                       expand = c(0, 0), limits = c(0, NA)) +
    scale_y_continuous(breaks = 2:10, expand = c(0, 0)) +
    labs(x = "Number of sequence reads",
         y = "Helminth diversity",
         title = "") +
    theme(plot.margin = margin(rep(20, 4)))


AreaRich <- models %>%
    .[["modelArea"]] %>% .[[1]]

modelFig <- plot(ggeffect(model = AreaRich, terms = c("weight_kg", "season", "area")),
                 rawdata = TRUE) +
    labs(x = "Weight (kg)", y = "Species richness (Hill number q=0)",
         title = "") +
    coord_cartesian(expand = FALSE, clip = "off") +
    scale_color_manual(values = colors_seasons, 
                       labels = c("Spring", "Summer + Autumn", "Winter"),
                       name = "Season:") +
    scale_fill_manual(values = colors_seasons, guide = "none") +
    ## need to repeat theme because it is overwritten by the wrapper
    theme_custom() +
    theme(plot.margin = margin(.5, .5, 2, 0))

wrap_plots(
  ## place first two plots
      alphaDivFox,  ## -> A
      ## ... then the other plot
      modelFig, ## -> B
      ## you can build more complex layouts by providing simple letters that are
      ## then filled accordingly by the plots you defined in the previous step;
      ## the plots are "named" as the appear here: the first one is A, the next B and so on...
      design = "A\nB"
      ## by default all rows and columns have similar widths and heights but we
      ## don't want our legend to fill up 1/3 of the plot height
    ## heights = c(21/45, 3/45, 21/45),
    ## this tells ggplot to just draw the legend once, not six times
    ##  guides = "collect"
) +
    plot_annotation(tag_levels = 'a', theme = theme(legend.title = element_text(hjust = .5)))

ggsave("figures/Fig2DivModel.png", width = 12.5, height = 10.5, units = "in")


stargazer(lapply(models, "[[", 1), type="html",
          out="tables/HelmDiversityq0.html")
