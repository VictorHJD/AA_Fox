library(dplyr)
library(raster)
library(sf)
library(readr)

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
BE_BB_envcov <- extract(tree_cover, as_Spatial(fox_sp), buffer = 1000,
                        fun = mean, na.rm = T, sp = T)

BE_BB_envcov2 <- extract(imperv, BE_BB_envcov, buffer = 1000,
                         fun = mean, na.rm = T, sp = T)

BE_BB_alldata <- extract(human_fpi, BE_BB_envcov2, buffer = 1000,
                         fun = mean, na.rm = T, sp = T)

as.data.frame(BE_BB_alldata) %>% 
rename(tree_cover_1000m = tcd_bb_mv_b_20m_3035, 
       imperv_1000m = imp_bb_mv_b_20m_3035,
       human_fpi_1000m = HFP2009_int_3035) ->
    fox_variables

fox_variables$tree_cover_1000m <- fox_variables$tree_cover_1000m * 100
fox_variables$human_fpi_1000m <- fox_variables$human_fpi_1000m * 2


readr::write_rds(fox_variables,
                 "intermediate_data/Fox_data_envir.RDS")
