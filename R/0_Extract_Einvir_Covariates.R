library(dplyr)
library(raster)
library(sf)
library(stars)
library(readr)
library(patchwork)
library(ggplot2)
library(systemfonts)
library(scico)
#library(tmap)

source("./R/plot_setup.R")

#devtools::install_github("EcoDynIZW/d6berlin")


## the required raw data
rawdata_dir <- "input_data/tifs/"
foxsample_data <- read.csv("input_data/Fox_data.csv")

## remove foxes with geocoordinate information 
filter(foxsample_data, !is.na(location_long) & !is.na(location_lat)) %>%
  st_as_sf(coords = c("location_long", "location_lat"), crs = 4236) %>%
  st_transform(crs = 3035) ->
  fox_sp 


#########################################################################
## THIS IS THE CODE TO CREATE THE TREE COVER RASTER FOR THE STUDY AREA ##
#########################################################################

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

## double check values make sense in a map
## make env cov spatial
fox_envcov_sf <- st_as_sf(fox_envcov, coords = c("coords.x1", "coords.x2"), crs = 3035)


## static map with ggplot2 + sf 
## External circle represents the values at the 1000m buffer, the inner circle 
## represents the values at the 100m buffer
b <- as(extent(4400000, 4700000, 3100000, 3400000), 'SpatialPolygons')
crs(b) <- crs(imperv)
human_fpi_crop <-  st_as_stars(crop(human_fpi, b))
human_fpi_agg_1000m <- st_as_stars(terra::aggregate(crop(human_fpi, b), fact = 10, fun = "mean"))

## shape of federal states
boundaries <- 
  st_read("input_data/VG250_Bundeslaender_esri.geojson") %>% 
  st_transform(crs = st_crs(fox_envcov_sf)) %>% 
  filter(GEN %in% c("Berlin", "Brandenburg"))


## map study area
map_study_base <- 
  ggplot(fox_envcov_sf) +
  geom_stars(data = human_fpi_crop) +
  ## state boundaries
  geom_sf(data = boundaries, fill = NA, color = "black") +
  ## 1000m buffer
  geom_sf(size = 5, shape = 21, stroke = 1.2, fill = "white", color = "white") +
  geom_sf(size = 5, shape = 21, stroke = 0, fill = "white") +
  geom_sf(aes(color = human_fpi_1000m), shape = 16, size = 5, alpha = .7) +
  ## 100m buffer
  geom_sf(size = 1.2, shape = 21, stroke = .8, fill = "white", color = "black") +
  geom_sf(size = 1.2, shape = 21, stroke = 0, fill = "white") +
  geom_sf(aes(color = human_fpi_100m), shape = 16, size = 1.2, alpha = .7) +
  coord_sf(xlim = c(4410000, 4650000), ylim = c(3150000, 3387000)) +
  scale_fill_gradient(low = "grey30", high = "grey96", guide = "none") +
  scale_color_scico(
    palette = "batlow", begin = .1,
    name = "Human Footprint Index (2009)", limits = c(0, 50), breaks = 1:9*5, 
    guide = guide_colorsteps(barwidth = unit(18, "lines"), barheight = unit(.6, "lines"),
                             title.position = "top", title.hjust = 0, show.limits = TRUE)
  ) +
  labs(x = NULL, y = NULL) +
  theme_map()

map_study <- map_study_base +
  ggspatial::annotation_scale(
    location = "bl", text_family = "Open Sans", text_cex = 1.2
  ) +
  ggspatial::annotation_north_arrow(location = "tr")

ggsave("figures/raw/map_study_area.png", width = 6.6, height = 7, bg = "white", dpi = 600)


## map Berlin
map_berlin <- map_study_base +
  coord_sf(xlim = c(4531042, 4576603), ylim = c(3253866, 3290780)) +
  theme_void() + 
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = .8))

ggsave("figures/raw/map_berlin.png", width = 4, height = 3.7, dpi = 600)


## overview map
sf_world <- 
  st_as_sf(rworldmap::getMap(resolution = "low")) %>% 
  st_transform(crs = st_crs(fox_envcov_sf)) %>% 
  st_buffer(dist = 0) %>% 
  dplyr::select(ISO_A2, SOVEREIGNT, LON, continent) %>% 
  mutate(area = st_area(.))

map_europe <- 
  ggplot(sf_world) +
  geom_sf(fill = "grey80", color = "grey96", lwd = .1) +
  geom_rect(
    xmin = 4430000, xmax = 4640000, ymin = 3160000, ymax = 3385000,
    color = "#212121", fill = "#a4cbb6", size = .7
  ) +
  geom_sf_text(
    data = filter(sf_world, ISO_A2 %in% c(
      "DE", "SE", "FR", "PL", "CZ", "IT", "ES", "AT", "CH", "GB", "PT", "NL", "BE", "IR", "IS"
    )),
    aes(label = ISO_A2),
    family = "Open Sans", color = "grey40", fontface = "bold", size = 4.5,
    nudge_x = 20000, nudge_y = -10000
  ) +
  ggspatial::annotation_scale(
    location = 'tr', text_family = "Open Sans", text_cex = 1.2
  ) +
  coord_sf(xlim = c(2650000, 5150000), ylim = c(1650000, 5100000)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-10, 30, by = 10)) +
  labs(x = NULL, y = NULL) +
  theme_map() +
  theme(panel.ontop = FALSE,
        panel.grid.major = element_line(color = "grey75", linetype = "15", size = .3))

#map_europe
ggsave("figures/raw/map_europe.png", width = 5, height = 7, bg = "white", dpi = 600)


map_globe <- d6berlin::globe(col_earth = "grey80", col_water = "grey96", bg = TRUE)

#map_globe
ggsave("figures/raw/map_globe.png", width = 2.2, height = 2.2, dpi = 600)


### combined map
# map_overview <- map_europe + labs(tag = "A.") + inset_element(map_globe, .14, .75, .59, 1.2, align_to = "plot")
# m <- map_overview + (map_study + labs(tags = "B."))
# 
# ggsave(""figures/raw/study_overview.png", width = 11.5, height = 7, bg = "white", dpi = 600)




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

## Not explored in the currently submitted version (Februrary 2023)
## ggsave("figures/suppl/Env100_1000Cors.png", width = 18, height = 7, bg = "white", dpi = 600)
