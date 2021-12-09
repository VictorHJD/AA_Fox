library(dplyr)
library(raster)
library(sf)
library(readr)
library(patchwork)
library(ggplot2)
library(tmap)


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

### Load Copernicus raster for the area
tree_cover_cop <- raster(paste0(rawdata_dir, "TCD_2015_020m_eu_03035_d05_E40N30.tif"))
tree_cover_cop
plot(tree_cover_cop)
summary(tree_cover_cop)
unique(values(tree_cover_cop))


### To make the raster smaller, there are multiple options, here use the extent of another raster
crop_ext <- extent(imperv)

# cut to selected extent
tree_cover_bb <- crop(tree_cover_cop, crop_ext)
tree_cover_bb <- clamp(tree_cover_bb, upper = 200, useValues = FALSE) # set NA values
summary(tree_cover_bb)
plot(tree_cover_bb)
plot(st_geometry(fox_sp), add = TRUE, col = "red", pch = 16)


# INTERACTIVE PLOTTING TO SEE BETTER THE POINTS
tmap_mode("view")
tm_shape(tree_cover_bb) +
  tm_raster(palette = "Greens") +
  tm_shape(fox_sp) +
  tm_dots("red")

# save raster for further use
writeRaster(tree_cover_bb, filename = paste0(rawdata_dir, "NEW_TCD_2015_bb_020m_03035.tif"), overwrite = TRUE)


#####################################
## END PREPARING TREE COVER RASTER ##
#####################################

#######################################
### Load environmental rasters

tree_cover_bb <- raster(paste0(rawdata_dir, "NEW_TCD_2015_bb_020m_03035.tif"))

imperv <- raster(paste0(rawdata_dir, "/imp_bb_mv_b_20m_3035.tif"))
imperv
summary(imperv)
unique(values(imperv)) # this layer is correct, with values from 0 to 100.
# imperv <- clamp(imperv, upper = 200, useValues = FALSE) # I removed this because it doesn't do anything.

human_fpi <- raster(paste0(rawdata_dir, "/HFP2009_int_3035.tif"))
#human_fpi <- clamp(human_fpi, upper = 200, useValues = FALSE)



#######################################
## extract environmental variables

# 1000m buffer 
envcov_1 <- raster::extract(tree_cover_bb, fox_sp, buffer = 1000,
                            fun = mean, na.rm = T, sp = TRUE)
envcov_2 <- raster::extract(imperv, envcov_1, buffer = 1000,
                            fun = mean, na.rm = T, sp = TRUE)
envcov_1000m_df <- raster::extract(human_fpi, envcov_2, buffer = 1000,
                                   fun = mean, na.rm = T, sp = TRUE) %>% 
  as.data.frame() %>% 
  rename(tree_cover_1000m = NEW_TCD_2015_bb_020m_03035,
         imperv_1000m = imp_bb_mv_b_20m_3035, 
         human_fpi_1000m = HFP2009_int_3035)

# 100m buffer 
envcov_3 <- raster::extract(tree_cover_bb, fox_sp, buffer = 100,
                            fun = mean, na.rm = T, sp = TRUE)
envcov_4 <- raster::extract(imperv, envcov_3, buffer = 100,
                            fun = mean, na.rm = T, sp = TRUE)
envcov_100m_df <- raster::extract(human_fpi, envcov_4, buffer = 100,
                                  fun = mean, na.rm = T, sp = TRUE) %>% 
  as.data.frame() %>% 
  rename(tree_cover_100m = NEW_TCD_2015_bb_020m_03035,
         imperv_100m = imp_bb_mv_b_20m_3035, 
         human_fpi_100m = HFP2009_int_3035) %>% 
  dplyr::select(tree_cover_100m, imperv_100m, human_fpi_100m, IZW_ID)

# Put together in one table
fox_envcov <-  left_join(envcov_1000m_df, envcov_100m_df, 
                         by = "IZW_ID")

# checking everything looks fine
nrow(envcov_1000m_df)
nrow(envcov_100m_df)
nrow(fox_envcov)

## Save environmental values 
readr::write_rds(fox_envcov,
                 "intermediate_data/Fox_data_envir.RDS")

## double check values make sense in a map
# make env cov spatial
fox_envcov_sf <- st_as_sf(fox_envcov, coords = c("coords.x1", "coords.x2"), crs = 3035)

# interactive map. External circle represents the values at the 1000m buffer, the inner circle represents the values at the 100m buffer
tm_shape(tree_cover_bb) +
  tm_raster(palette = "Greens", alpha = 0.5) +
  tm_shape(fox_envcov_sf) +
  tm_dots("tree_cover_1000m", palette = "Greens", size = 0.1) +
  tm_shape(fox_envcov_sf) +
  tm_dots("tree_cover_100m", palette = "Greens", size = 0.04) 


# Map to save to figure
tmap_mode("plot")
map_tc <- tm_shape(tree_cover_bb, bbox = ) +
  tm_raster(palette = "Greens", alpha = 0.5, legend.show = FALSE) + #remove the legend because it is the same for all layers
  tm_shape(fox_envcov_sf) +
  tm_symbols(col = "tree_cover_1000m", palette = "Greens", border.col = "black", 
             size = 1.5, legend.col.show = TRUE, title.col =  "Tree cover \npercentage") + 
  tm_shape(fox_envcov_sf) +
  tm_symbols(col = "tree_cover_100m", palette = "Greens", border.col = "black", 
             size = 0.5, legend.col.show = FALSE) + #remove the legend because it is the same 
  tm_scale_bar(position = c("left", "BOTTOM"), text.size = 1) +
  tm_compass(position = c("right", "TOP")) +
  tm_layout(legend.frame = TRUE)

map_tc

tmap_save(map_tc, filename = "figures/suppl/Map_treecover.pdf")



###################################################
### Plots comparing values at the two buffers 

fox_envcov %>% group_by(area) %>%
  summarize(treeCoverCor=cor(tree_cover_1000m, tree_cover_100m,
                             use="pairwise.complete.obs"),
            impervCor=cor(imperv_100m, imperv_1000m,
                          use="pairwise.complete.obs"),
            hfpiCor=cor(human_fpi_100m, human_fpi_1000m, 
                        use="pairwise.complete.obs")) 

tree_cover <- ggplot(fox_envcov, aes(tree_cover_1000m, tree_cover_100m, color=area)) +
  scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
  scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
  theme(legend.position="none", #axis.text.y = element_text(family = font_num) # missing the font_num object
  )+
  geom_point() 

imperv <- ggplot(fox_variables, aes(imperv_1000m, imperv_100m, color=area)) +
  scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
  scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
  theme(legend.position="none", #axis.text.y = element_text(family = font_num)
  )+
  geom_point() 

hfpi <- ggplot(fox_variables, aes(human_fpi_1000m, human_fpi_100m, color=area)) +
  scale_colour_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
  scale_fill_manual(values = c("#e7b800", "#2e6c61"), name = "Study area:") +
  geom_point()

pdf("figures/suppl/Env100_1000Cors_2.pdf", width=15, height=5)
tree_cover + imperv + hfpi
dev.off()


