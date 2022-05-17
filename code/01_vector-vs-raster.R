# Script that tests differences in area based metrics when we use a vector or a raster
# TL:DR; differences emerge if we are running patch based metrics (especially for very
# small patch sizes)

# We will use three objects for testing
# Use the 1870 nilgiris vector file
# Rasterize the above 1870 nilgiris vector file 
# Rasterize and resample the above 1870 nilgiris vector file

# To load these files and libraries, please use the land-cover-analysis script 
# which is located in the main Git Repo

# Use the nil1870 shapefile as an example
areaVect <- nil1870 %>%
  mutate(areaSqKm = as.numeric(st_area(nil1870)/1000000)) # in sq. km.

sumAreaVect <- areaVect %>%
  group_by(class) %>%
  mutate(sumArea = sum(areaSqKm)) %>%
  distinct(class, sumArea)

# use the nil1870 as a raster
vectnil1870 <- terra::vect(nil1870)
emptyRast <- terra::rast(res = 12,xmin=661012.738348484, xmax= 718286.685853398, ymin=1240744.4111945, ymax=1275965.51753469, crs="+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs")
rastnil1870 <- terra::rasterize(vectnil1870, emptyRast,"class")
sz <- cellSize(rastnil1870, unit="m")
areaRast <- zonal(sz, rastnil1870, sum)
areaRastSqKm <- (areaRast$area/1000000)
areaRast <- cbind(areaRast,areaRastSqKm)

# Since the scale has not changed and the rasterization was done at 12m resolution, you will notice that we did not lose any data in terms of overall area calculations. 

# convert to 30m resolution and calculate areas
# To do so, let's use the nil1973 file to resample it to that resolution:
  
vectnil1973 <- vect(nil1973)
emptyRast1973 <- terra::rast(res = 30,xmin=652338.464550138, xmax= 718740, ymin=1226658.15395737, ymax=1275965.51753469, crs="+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs")
rastnil1973 <- terra::rasterize(vectnil1973, emptyRast1973,"class")
  
# resampling the 1870 raster to the same resolution as the 1973 raster
resamNil1870 <- terra::resample(rastnil1870,rastnil1973, method="near")
szAgg1870 <- cellSize(resamNil1870, unit="m")
areaResam1870 <- zonal(szAgg1870, resamNil1870, sum)
areaResamSqKm <- (areaResam1870$area/1000000)
areaResam1870 <- cbind(areaResam1870,areaResamSqKm)
  
# We observed no significant differences in the overall area calculated 
# between a vector file and a vector file that was rasterized and resampled.

  
# Testing if there are differences in area based statistics at the patch level

# for the sake of this analysis we create size classes to 
# test for differences in the number of patches by area and by land cover type
areaVect <- areaVect %>%
    mutate(size_class = case_when(areaSqKm>=0 & areaSqKm<0.5 ~ "0-0.5",
                                  areaSqKm>=0.5 & areaSqKm<5 ~"0.5-5",
                                  areaSqKm>=5 & areaSqKm<50 ~ "5-50",
                                  areaSqKm>=50 & areaSqKm<100 ~ "50-100",
                                  areaSqKm>=100 & areaSqKm <1000 ~ "100-1000"))
  
size_class_vect <- areaVect %>%
    group_by(class, size_class) %>%
    count() %>%
    filter(!is.na(class))
  
# Plot the above result
fig_sizeClass_vect <- ggplot(size_class_vect, aes(x=size_class, y=n, fill=class)) + geom_bar(stat="identity", position=position_dodge(), fill="#883107", alpha=0.9) +
    geom_text(aes(label=n, hjust="middle", vjust=0.2), 
              position = position_dodge(), angle=0, size=5) +
    facet_wrap(~class) +
    theme_bw() +
    labs(x="\nLand cover type and size class", 
         y="nPatches\n") +
    theme(axis.title = element_text(family = "Century Gothic",
                                    size = 14, face = "bold"),
          axis.text = element_text(family="Century Gothic",size = 14),
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          legend.position = "none")
  
ggsave(fig_sizeClass_vect, filename = "figs/fig_size_vect.png", width=17, height=7,
         device = png(), units="in", dpi = 300)
  
# Here, we notice a >1000 forest patches (0-0.5 sqkm in size) a
# and two very large grassland patches (due to autocompletion)
  
# Using a raster of the vector file above at the same resolution:
rastpatch12 <- lsm_p_area(rastnil1870) # 12 refers to the 12m resolution
rastpatch12$value <- (rastpatch12$value)/100
  
# mutate on the tibble
areaRast12 <- rastpatch12 %>%
    mutate(size_class = case_when(value>=0 & value<0.5 ~ "0-0.5",
                                  value>=0.5 & value<5 ~"0.5-5",
                                  value>=5 & value<50 ~ "5-50",
                                  value>=50 & value<100 ~ "50-100",
                                  value>=100 & value <1000 ~ "100-1000")) %>%
    mutate(class = case_when(class == 0 ~ "NA",
                             class == 1 ~ "Shola Forest",
                             class == 2 ~ "Shola Grassland",
                             class == 3 ~ "cultivated tracts",
                             class == 4 ~ "plantations",
                             class == 5 ~ "settlements",
                             class == 6 ~ "waterbodies"))
  
size_class_rast12 <- areaRast12 %>%
    group_by(class, size_class) %>%
    count() %>%
    filter(class != "NA")
  
# Plot the above result
fig_sizeClass_rast12 <- ggplot(size_class_rast12, aes(x=size_class, y=n, fill=class)) + geom_bar(stat="identity", position=position_dodge(), fill="#883107", alpha=0.9) +
    geom_text(aes(label=n, hjust="middle", vjust=0.2), 
              position = position_dodge(), angle=0, size=5) +
    facet_wrap(~class) +
    theme_bw() +
    labs(x="\nLand cover type and size class", 
         y="nPatches\n") +
    theme(axis.title = element_text(family = "Century Gothic",
                                    size = 14, face = "bold"),
          axis.text = element_text(family="Century Gothic",size = 14),
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          legend.position = "none")
  
ggsave(fig_sizeClass_rast12, filename = "figs/fig_size_rast12.png", width=17, height=7, device = png(), units="in", dpi = 300)
  
# repeating the above for the resampled raster
rastPatch <- lsm_p_area(resamNil1870)
rastPatch$value <- (rastPatch$value)/100
  
# mutate on the tibble
areaRast <- rastPatch %>%
    mutate(size_class = case_when(value>=0 & value<0.5 ~ "0-0.5",
                                  value>=0.5 & value<5 ~"0.5-5",
                                  value>=5 & value<50 ~ "5-50",
                                  value>=50 & value<100 ~ "50-100",
                                  value>=100 & value <1000 ~ "100-1000")) %>%
    mutate(class = case_when(class == 0 ~ "NA",
                             class == 1 ~ "Shola Forest",
                             class == 2 ~ "Shola Grassland",
                             class == 3 ~ "cultivated tracts",
                             class == 4 ~ "plantations",
                             class == 5 ~ "settlements",
                             class == 6 ~ "waterbodies"))
  
size_class_rast <- areaRast %>%
    group_by(class, size_class) %>%
    count() %>%
    filter(class != "NA")
  
# Plot the above result
fig_sizeClass_rast <- ggplot(size_class_rast, aes(x=size_class, y=n, fill=class)) + geom_bar(stat="identity", position=position_dodge(), fill="#883107", alpha=0.9) +
    geom_text(aes(label=n, hjust="middle", vjust=0.2), 
              position = position_dodge(), angle=0, size=5) +
    facet_wrap(~class) +
    theme_bw() +
    labs(x="\nLand cover type and size class", 
         y="nPatches\n") +
    theme(axis.title = element_text(family = "Century Gothic",
                                    size = 14, face = "bold"),
          axis.text = element_text(family="Century Gothic",size = 14),
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          legend.position = "none")
  
ggsave(fig_sizeClass_rast, filename = "figs/fig_size_rast.png", width=17, height=7,
         device = png(), units="in", dpi = 300)

# Following completion of this analysis, we notice that the number of small
# patches increased dramatically for Shola Grassland and waterbodies when we
# rasterized as well as rasterized+resampled.

# However, marginal increases in the number of patches were seen for other land cover
# classes.

  