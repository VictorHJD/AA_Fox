library(dplyr)
library(raster)
library(sf)
library(readr)
library(patchwork)


## the required raw data
rawdata_dir <- "input_data/tifs/"
foxsample_data <- read.csv("input_data/Fox_data.csv")

## remove foxes with geocoordinate information 
filter(foxsample_data, !is.na(location_long) & !is.na(location_lat)) %>%
    st_as_sf(coords = c("location_long", "location_lat"), crs = 4236) %>%
    st_transform(crs = 3035) ->
fox_sp 


## read the raw raster data
tree_cover <- raster(paste0(rawdata_dir, "/tcd_bb_mv_b_20m_3035.tif")) 
tree_cover <- clamp(tree_cover, upper = 200, useValues = FALSE)

imperv <- raster(paste0(rawdata_dir, "/imp_bb_mv_b_20m_3035.tif"))
imperv <- clamp(imperv, upper = 200, useValues = FALSE)

human_fpi <- raster(paste0(rawdata_dir, "/HFP2009_int_3035.tif"))
#human_fpi <- clamp(human_fpi, upper = 200, useValues = FALSE)


## extract environmental summary variables
BE_BB_envcov1 <- extract(tree_cover, as_Spatial(fox_sp), buffer = 1000,
                         fun = mean, na.rm = T, sp = T)

BE_BB_envcov2 <- extract(tree_cover, BE_BB_envcov1, buffer = 100,
                         fun = mean, na.rm = T, sp = T)

BE_BB_envcov3 <- extract(imperv, BE_BB_envcov2, buffer = 1000,
                         fun = mean, na.rm = T, sp = T)

BE_BB_envcov4 <- extract(imperv, BE_BB_envcov3, buffer = 100,
                         fun = mean, na.rm = T, sp = T)

BE_BB_envcov5 <- extract(human_fpi, BE_BB_envcov4, buffer = 1000,
                         fun = mean, na.rm = T, sp = T)

BE_BB_alldata <- extract(human_fpi, BE_BB_envcov5, buffer = 100,
                         fun = mean, na.rm = T, sp = T)



as.data.frame(BE_BB_alldata) %>% 
    rename(tree_cover_1000m = tcd_bb_mv_b_20m_3035,
           tree_cover_100m = tcd_bb_mv_b_20m_3035.1, 
           imperv_1000m = imp_bb_mv_b_20m_3035,
           imperv_100m = imp_bb_mv_b_20m_3035.1,
           human_fpi_1000m = HFP2009_int_3035,
           human_fpi_100m = HFP2009_int_3035.1) ->
    fox_variables

### now that's something I don't understand! Why is this multipied (differently!)?
fox_variables$tree_cover_1000m <- fox_variables$tree_cover_1000m * 100
fox_variables$tree_cover_100m <- fox_variables$tree_cover_100m * 100
fox_variables$human_fpi_1000m <- fox_variables$human_fpi_1000m * 2
fox_variables$human_fpi_100m <- fox_variables$human_fpi_100m * 2

readr::write_rds(fox_variables,
                 "intermediate_data/Fox_data_envir.RDS")

fox_variables %>% group_by(area) %>%
    summarize(treeCoverCor=cor(tree_cover_1000m, tree_cover_100m,
                               use="pairwise.complete.obs"),
              impervCor=cor(imperv_100m, imperv_1000m,
                            use="pairwise.complete.obs"),
              hfpiCor=cor(human_fpi_100m, human_fpi_1000m, 
                          use="pairwise.complete.obs")) 

tree_cover <- ggplot(fox_variables, aes(tree_cover_1000m, tree_cover_100m, color=area)) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    theme(legend.position="none", axis.text.y = element_text(family = font_num))+
    geom_point() 

imperv <- ggplot(fox_variables, aes(imperv_1000m, imperv_100m, color=area)) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    theme(legend.position="none", axis.text.y = element_text(family = font_num))+
    geom_point() 

hfpi <- ggplot(fox_variables, aes(human_fpi_1000m, human_fpi_100m, color=area)) +
    scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
    geom_point()

pdf("figures/suppl/Env100_1000Cors.pdf", width=15, height=5)
tree_cover + imperv + hfpi
dev.off()

