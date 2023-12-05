library(vegan)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(ggnewscale)

source("./R/plot_setup.R")

recomputeBioinfo <- FALSE

if(!exists("PSGHelm")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSGHelm <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
    }
}


PSGHelm <- prune_taxa(taxa_sums(PSGHelm) > 0, PSGHelm)
PSGHelm <- prune_samples(sample_sums(PSGHelm) > 0, PSGHelm)

PSGHelm ### 131 foxes without NA anywhere

HelmData <- as.data.frame(otu_table(PSGHelm))
colnames(HelmData) <- tax_table(PSGHelm)[, "genus"]
class(HelmData) <- "data.frame"

EnvData <- as.data.frame(
    sample_data(PSGHelm)[,
                         !colnames(sample_data(PSGHelm))%in%c("date",
                                                              "date_found")])
class(EnvData) <- "data.frame"

EnvData <- EnvData[rownames(HelmData), ]

EnvData$weight_kg <- as.numeric(EnvData$weight_kg)
EnvData$tree_cover_1000m <- as.numeric(EnvData$tree_cover_1000m)
EnvData$imperv_1000m <- as.numeric(EnvData$imperv_1000m)
EnvData$human_fpi_1000m <- as.numeric(EnvData$human_fpi_1000m)
EnvData$DNAng.ul <- NULL
EnvData$DNA260.280 <- NULL
EnvData$DNA260.230 <- NULL
EnvData$age <- NULL

HelmData$Aelurostrongylus <- NULL

### NO OTHER environmental variables are better explaining composition!
PERMA <- vegan::adonis2(HelmData>0 ~ area + weight_kg +
                            sex  + season, 
                        data=EnvData, 
                        na.action = na.fail, by="margin",
                        method="jaccard")

PERMA

### still area is the best model, everything below (the other
### environmental predictors) is not as good!
PERMAimp <- adonis2(HelmData>0 ~  imperv_1000m + weight_kg + 
                        sex  + season,
                    data=EnvData, 
                    na.action = na.fail, by="margin",
                    method="jaccard")

PERMAtree <- adonis2(HelmData>0 ~ tree_cover_1000m + weight_kg +
                         sex  + season,
                     data=EnvData, 
                     na.action = na.fail, by="margin",
                     method="jaccard")

PERMAhuman <- adonis2(HelmData>0 ~ human_fpi_1000m + weight_kg + 
                          sex  + season,
                      data=EnvData, 
                      na.action = na.fail, by="margin",
                      method="jaccard")

PERMAall <- adonis2(HelmData>0 ~ human_fpi_1000m + tree_cover_1000m +
                        imperv_1000m + weight_kg +
                        sex  + season, 
                      data=EnvData, 
                      na.action = na.fail, by="margin",
                      method="jaccard")


PERMAallX <- adonis2(HelmData>0 ~ area + tree_cover_1000m +
                         imperv_1000m + weight_kg + 
                         sex  + season, 
                     data=EnvData, 
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

HelmEnvFit <- envfit(nMDSHelm, EnvData[ , c("area", "human_fpi_1000m",
                                            "tree_cover_1000m",
                                            "imperv_1000m",
                                            "weight_kg", "sex",
                                            "season"
                                            )])


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


ScoresHelm <-  as.data.frame(scores(nMDSHelm)$sites)
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
               linewidth = 1.5) +
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
               arrow = arrow(length = unit(arrowhead, "npc"), angle = 23, type = "closed"),
               linewidth = 1.5) +
  shadowtext::geom_shadowtext( ## TODO: moved down here to avoid text overplotting
    data = subset(HelmEnvFitDf, pvals < 0.1),
    aes(x = NMDS1, y = NMDS2 + offset),
    label = row.names(subset(HelmEnvFitDf, pvals < 0.1)), size = 5.5,
    color = "#757575", bg.colour = "white", family = "Open Sans", fontface = "bold"
  ) +
  shadowtext::geom_shadowtext( ## TODO: why are the text labels not colored with the same gradient?
    data = subset(HelmHelmDf, pvals < 0.1), aes(x = NMDS1, y = NMDS2 + offset),
    label = row.names(subset(HelmHelmDf, pvals < 0.1)), size = 4.5,
    color = pal[1], bg.colour = "white", family = "Open Sans", fontface = "bold"
  ) +
  scale_x_continuous(limits = c(-2.9, NA)) +
  scale_color_gradientn(colors = pal, name = "log(pvals)") +
  theme(legend.box.margin = margin(l = 30),
        legend.title = element_text(margin = margin(t = 1, b = 3)))

ggsave("figures/CompositionEnvHelm.png", ggHelmEnvHelm, 
       width = 14, height = 6, bg = "white", dpi = 600)


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
      
