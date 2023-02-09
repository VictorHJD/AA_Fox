library(vegan)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(ggnewscale)

source("./R/plot_setup.R")

recomputeBioinfo <- FALSE

if(!exists("PSG")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSG <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
    }
}

## ## Dec 2022 we don't need this anymore -> Excluded!
## recomputeDiversity <- FALSE

## if(!"FunM_Species_richness"%in%colnames(sample_data(PSG))|
##    recomputeDiversity){
##     source("R/2_iNEXT_fox.R")
## }

## ## We also don't us this here!!
## ## Helminth traits
## traits <- read.csv("input_data/helminth_traits.csv")
## traits %>%
##     column_to_rownames("t.genus")  -> traits 

PSGHelm <- phyloseq::subset_taxa(PSG, category%in%"Helminth")

PSGHelm <- subset_samples(PSGHelm, !is.na(sample_data(PSGHelm)[, "condition"]) &
                                   !is.na(sample_data(PSGHelm)[, "weight_kg"]) &
                                   !is.na(sample_data(PSGHelm)[, "season"]) &
                                   !is.na(sample_data(PSGHelm)[, "tree_cover_1000m"]) &
                                   !is.na(sample_data(PSGHelm)[, "DNAng.ul"]) &
                                   !is.na(sample_data(PSGHelm)[, "DNA260.230"]) &
                                   !is.na(sample_data(PSGHelm)[, "DNA260.280"]) 
                          )

PSGHelm ### 150 foxes without NA anywhere

## only taxa with reads
PSGHelm <- prune_taxa(taxa_sums(PSGHelm) > 0, PSGHelm)
PSGHelm <- prune_samples(sample_sums(PSGHelm) > 0, PSGHelm)

PSGHelm ## 139 foxes with all samples any taxa

HelmData <- otu_table(PSGHelm)
colnames(HelmData) <- tax_table(PSGHelm)[, "genus"]

EnvData <- sample_data(PSGHelm)[, !colnames(sample_data(PSGHelm))%in%c("date",
                                                                       "date_found")]
class(EnvData) <- "data.frame"
EnvData$weight_kg <- as.numeric(EnvData$weight_kg)
EnvData$tree_cover_1000m <- as.numeric(EnvData$tree_cover_1000m)
EnvData$imperv_1000m <- as.numeric(EnvData$imperv_1000m)
EnvData$human_fpi_1000m <- as.numeric(EnvData$human_fpi_1000m)
EnvData$DNAng.ul <- as.numeric(EnvData$DNAng.ul)
EnvData$DNA260.280 <- as.numeric(EnvData$DNA260.280)
EnvData$DNA260.230 <- as.numeric(EnvData$DNA260.230)


### This shoud be the same now after removing all the NAs already
### above
EnvDataNA <- na.omit(EnvData)
HelmDataNA <- HelmData[rownames(EnvDataNA), ]

### NO OTHER environmental variables are better explaining composition!

PERMA <- vegan::adonis2(HelmDataNA ~ area + weight_kg + age +
                            sex  + season + condition +
                            DNAng.ul + DNA260.230 + DNA260.280,
                     data=EnvDataNA, 
                     na.action = na.fail, by="margin",
                     method="jaccard")

PERMA

### still area is the best model, everything below (the other
### environmental predictors) is not as good!

PERMAimp <- adonis2(HelmDataNA ~  imperv_1000m + weight_kg + age +
                        sex  + season + condition +
                        DNAng.ul + DNA260.230 + DNA260.280,
                    data=EnvDataNA, 
                    na.action = na.fail, by="margin",
                    method="jaccard")

PERMAtree <- adonis2(HelmDataNA ~ tree_cover_1000m + weight_kg + age +
                         sex  + season + condition +
                         DNAng.ul + DNA260.230 + DNA260.280,
                     data=EnvDataNA, 
                     na.action = na.fail, by="margin",
                     method="jaccard")

PERMAhuman <- adonis2(HelmDataNA ~ human_fpi_1000m + weight_kg + age +
                          sex  + season + condition +
                          DNAng.ul + DNA260.230 + DNA260.280, 
                      data=EnvDataNA, 
                      na.action = na.fail, by="margin",
                      method="jaccard")

PERMAall <- adonis2(HelmDataNA ~ human_fpi_1000m + tree_cover_1000m +
                        imperv_1000m + weight_kg + age + condition +
                        sex  + season +
                        DNAng.ul + DNA260.230 + DNA260.280, 
                      data=EnvDataNA, 
                      na.action = na.fail, by="margin",
                      method="jaccard")


PERMAallX <- adonis2(HelmDataNA ~ area + tree_cover_1000m +
                         imperv_1000m + weight_kg + age + condition +
                         sex  + season +
                         DNAng.ul + DNA260.230 + DNA260.280, 
                     data=EnvDataNA, 
                     na.action = na.fail, by="margin",
                     method="jaccard")



## this does not work, let's stick with the csv
## stargazer(list(PERMA, PERMAtree, PERMAimp, PERMAhuman), type="html",
   ##       out="tables/Permanova.html")

write.csv(round(PERMA, 3), "tables/Permanova.csv")

write.csv(round(rbind(PERMAimp, PERMAtree, PERMAhuman, PERMAallX), 2),
          "tables/suppl/PermanovaConti.csv")

nMDSHelm <- metaMDS(HelmData, distance = "jaccard", weakties = FALSE,
                    try=1500, trymax=1500, k=3,
                    center = TRUE)

HelmEnvFit <- envfit(nMDSHelm, EnvData[ , c("area", "human_fpi_1000m", "tree_cover_1000m",
                                            "imperv_1000m",
                                            "age", "weight_kg", "sex",
                                            "condition", "season",
                                            "DNAng.ul", "DNA260.230", "DNA260.280")])


### AMAZING! ThIs MAKeS SeNSe!!!!
HelmEnvFit

## to see which helminths drive this
HelmHelmFit <- envfit(nMDSHelm, HelmData, na.rm=TRUE)

### Haha now seeing that envfit has a "tidy" argument to produce
### ggplot compatible output... all the below would likely not have
### been necessary

faclevels <- sapply(names(HelmEnvFit$factors$pvals),
                    function (x) {
                        thatoften <- sum(grepl(x, rownames(scores(HelmEnvFit, "factors"))))
                        rep(x, thatoften)
                    })

HelmEnvFitDf <-
    as.data.frame(
        rbind(cbind(scores(HelmEnvFit, "vectors") * ordiArrowMul(HelmEnvFit),
                    pvals=HelmEnvFit$vectors$pvals),
              cbind(scores(HelmEnvFit, "factors") * ordiArrowMul(HelmEnvFit),
                    ## repeat times the numbers of factor levels
                    pvals=HelmEnvFit$factors$pvals[unlist(faclevels)])))
dev.off() ## this had somehow opened a graphics device?!

HelmHelmDf <-
    as.data.frame(
        cbind(scores(HelmHelmFit, "vectors") * ordiArrowMul(HelmHelmFit),
              pvals=HelmHelmFit$vectors$pvals)
    )
dev.off() ## this had somehow opened a graphics device?!


ScoresHelm <-  as.data.frame(scores(nMDSHelm))
ScoresHelm <- cbind(ScoresHelm, EnvData)
ScoresHelm$season <- factor(ScoresHelm$season, levels = c("spring", "S_autumn", "winter"))

## arrow head size and offsetting for the graphic
arrowhead <- 0.04
offset <- 0.04

ggHelmEnv <- 
  ggplot(data = ScoresHelm, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = ScoresHelm, aes(fill = area, color = after_scale(fill), shape = season), 
             size = 3, alpha = .5, stroke = .8) +
  #scale_colour_manual(values = colors_regions, name = "Study area:") +
  scale_fill_manual(values = colors_regions, name = "Study area:") +
  ### ordiArrowMul didn't work somehow have to scale manually
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2,
                   color = log(pvals)),
               data = subset(HelmEnvFitDf, pvals < 0.1),
               arrow = arrow(length = unit(arrowhead, "npc"), angle = 23, type = "closed"),
               size = 1.5) +
  coord_cartesian(clip = "off") +
  scale_color_viridis_c(option = "cividis", name = "log(pvals)") +
  scale_shape_manual(values = c(21, 22, 23), name = "Season:", 
                     labels = c("Spring", "Summer + Autumn", "Winter"),
                     guide = guide_legend(override.aes = list(alpha = 1, color = "black", fill = "black")))
    
pal <- scico::scico(palette = "batlowK", n = 100, begin = .1, end = .8)

ggHelmEnvHelm <-
  ggHelmEnv +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, color = pvals),
               data = subset(HelmHelmDf, pvals < 0.1),
               arrow = arrow(length = unit(arrowhead, "npc"), angle = 23, type = "closed"), size = 1.5) +
  shadowtext::geom_shadowtext( ## TODO: moved down here to avoid text overplotting
    data = subset(HelmEnvFitDf, pvals < 0.1),
    aes(x = NMDS1, y = NMDS2 + offset),
    label = row.names(subset(HelmEnvFitDf, pvals < 0.1)), size = 5.5,
    color = "darkgrey", family = "Open Sans", fontface = "bold"
  ) +
  shadowtext::geom_shadowtext( ## TODO: why are the text labels not colored with the same gradient?
    data = subset(HelmHelmDf, pvals < 0.1), aes(x = NMDS1, y = NMDS2 + offset),
    label = row.names(subset(HelmHelmDf, pvals < 0.1)), size = 4.5,
    color = pal[1], bg.colour = "white", family = "Open Sans", fontface = "bold"
  ) +
  scale_color_gradientn(colors = pal, name = "log(pvals)") +
  theme(legend.margin = margin(l = 20),
        legend.title = element_text(margin = margin(b = 4)))

ggsave("figures/CompositionEnvHelm.png", ggHelmEnvHelm, 
       width = 18, height = 7, bg = "white", dpi = 600)


## now for the table

HelmTabEF <- cbind(R2=HelmHelmFit$vectors$r, Pval=HelmHelmFit$vectors$pvals)

EnvTabEF <- rbind(
    cbind(R2=HelmEnvFit$factors$r, Pval=HelmEnvFit$factors$pvals),
    cbind(R2=HelmEnvFit$vectors$r, Pval=HelmEnvFit$vectors$pvals)
    )

write.csv(round(
    rbind(EnvTabEF[order(EnvTabEF[, "Pval"]), ], 
          cbind(R2=0, Pval=0),
          HelmTabEF[order(HelmTabEF[, "Pval"]), ]),
    3), "tables/suppl/EnvFitnMDS.csv", 
    )
      
