library(dplyr)
library(raster)
library(sf)
library(stars)
library(readr)
library(patchwork)
library(ggplot2)
library(systemfonts)
library(scico)
library(ggspatial)
#library(tmap)

## devtools::install_github("EcoDynIZW/d6berlin")

source("./R/plot_setup.R")
extrafont::loadfonts(device = "all") ## run every time


## the required raw data
rawdata_dir <- "input_data/tifs/"
foxsample_data <- read.csv("input_data/Fox_data.csv")

## remove foxes with geocoordinate information 
filter(foxsample_data, !is.na(location_long) & !is.na(location_lat)) %>%
  st_as_sf(coords = c("location_long", "location_lat"), crs = 4236) %>%
  st_transform(crs = 3035) ->
  fox_sp 


######################################################################################################
## THIS IS THE CODE TO CREATE THE TREE COVER RASTER FOR THE STUDY AREA - NOT INCLUDED IN REPOSITORY ##
######################################################################################################

## Original raster not included in github because it is ~600Mb
## to download it got to:


### Load Copernicus raster for the area
## tree_cover_cop <- raster(paste0(rawdata_dir, "TCD_2015_020m_eu_03035_d05_E40N30.tif"))
## tree_cover_cop
## plot(tree_cover_cop)
## summary(tree_cover_cop)
## unique(values(tree_cover_cop))


## ### To make the raster smaller, there are multiple options, here use the extent of another raster
## crop_ext <- extent(imperv)

## # cut to selected extent
## tree_cover_bb <- crop(tree_cover_cop, crop_ext)
## tree_cover_bb <- clamp(tree_cover_bb, upper = 200, useValues = FALSE) # set NA values
## summary(tree_cover_bb)
## plot(tree_cover_bb)
## plot(st_geometry(fox_sp), add = TRUE, col = "red", pch = 16)


## # INTERACTIVE PLOTTING TO SEE BETTER THE POINTS
## tmap_mode("view")
## tm_shape(tree_cover_bb) +
##   tm_raster(palette = "Greens") +
##   tm_shape(fox_sp) +
##   tm_dots("red")

## # save raster for further use
## writeRaster(tree_cover_bb, filename = paste0(rawdata_dir, "NEW_TCD_2015_bb_020m_03035.tif"), overwrite = TRUE)


#####################################
## END PREPARING TREE COVER RASTER ##
#####################################

#######################################
### Load environmental rasters

## from here (smaller raster files this is included directly in the repository). 

tree_cover_bb <- raster(paste0(rawdata_dir, "NEW_TCD_2015_bb_020m_03035.tif"))

imperv <- raster(paste0(rawdata_dir, "/imp_bb_mv_b_20m_3035.tif"))

human_fpi <- raster(paste0(rawdata_dir, "/HFP2009_int_3035.tif"))

#######################################
## extract environmental variables

filepath <- "intermediate_data/Fox_data_envir.RDS"

if(!file.exists(filepath)) {
  
  # 1000m buffer 
  envcov_1 <- raster::extract(tree_cover_bb, fox_sp, buffer = 1000,
                              fun = mean, na.rm = TRUE, sp = TRUE)
  envcov_2 <- raster::extract(imperv, envcov_1, buffer = 1000,
                              fun = mean, na.rm = TRUE, sp = TRUE)
  envcov_1000m_df <- raster::extract(human_fpi, envcov_2, buffer = 1000,
                                     fun = mean, na.rm = TRUE, sp = TRUE) %>% 
    as.data.frame() %>% 
    rename(tree_cover_1000m = NEW_TCD_2015_bb_020m_03035,
           imperv_1000m = imp_bb_mv_b_20m_3035, 
           human_fpi_1000m = HFP2009_int_3035)
  
  ## 100m buffer 
  envcov_3 <- raster::extract(tree_cover_bb, fox_sp, buffer = 100,
                              fun = mean, na.rm = TRUE, sp = TRUE)
  envcov_4 <- raster::extract(imperv, envcov_3, buffer = 100,
                              fun = mean, na.rm = TRUE, sp = TRUE)
  envcov_100m_df <- raster::extract(human_fpi, envcov_4, buffer = 100,
                                    fun = mean, na.rm = TRUE, sp = TRUE) %>% 
    as.data.frame() %>% 
    rename(tree_cover_100m = NEW_TCD_2015_bb_020m_03035,
           imperv_100m = imp_bb_mv_b_20m_3035, 
           human_fpi_100m = HFP2009_int_3035) %>% 
    dplyr::select(tree_cover_100m, imperv_100m, human_fpi_100m, IZW_ID)
  
  ## Put together in one table
  fox_envcov <-  left_join(envcov_1000m_df, envcov_100m_df, 
                           by = "IZW_ID")
  
  ## Save environmental values 
  readr::write_rds(fox_envcov, filepath)
} else {
  fox_envcov <- readr::read_rds(filepath)
}


###################################################
### Plots comparing values at the two buffers 


fox_envcov %>% 
  group_by(area) %>%
  summarize(
    treeCoverCor = cor(tree_cover_1000m, tree_cover_100m,
                       use = "pairwise.complete.obs"),
    impervCor = cor(imperv_100m, imperv_1000m,
                    use = "pairwise.complete.obs"),
    hfpiCor = cor(human_fpi_100m, human_fpi_1000m, 
                  use = "pairwise.complete.obs")) 

plot_corr_buffer <- function(x, y, label) {
  g <- 
    ggplot(fox_envcov, aes_string(x, y, fill = "area")) +
    geom_point(aes(color = area), shape = 21, size = 2.5, stroke = .8, fill = "white") +
    geom_point(shape = 21, size = 2.5, alpha = .3, stroke = .8, color = "NA") +
    scale_color_manual(values = colors_regions, name = "Study area:") +
    scale_fill_manual(values = colors_regions, name = "Study area:") +
    labs(x = paste(label, "(1000 m buffer)"), y = paste(label, "(100 m buffer)"))
  return(g)
}

tree_cover <- plot_corr_buffer(
  x = "tree_cover_1000m", y = "tree_cover_100m", label = "Tree cover"
)

imperv <- plot_corr_buffer(
  x = "imperv_1000m", y = "imperv_100m", label = "Imperviousness"
)

hfpi <- plot_corr_buffer(
  x = "human_fpi_1000m", y = "human_fpi_100m", label = "Human FPI"
)

panel <- tree_cover + imperv + hfpi + plot_layout(guides = "collect")

## Not explored in the currently submitted version (December 2023)
## ggsave("figures/suppl/Env100_1000Cors.png", width = 18, height = 7, bg = "white", dpi = 600)
