#MARBEFES
#IRISH SEA
#Vulnerability Analysis
#18/06/2025

pkgs = c("tidyverse", "vegan", "patchwork", "lubridate", "sf", 
         "rnaturalearth", "rnaturalearthdata", "raster")
for(p in pkgs){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
} 

## Plot parameters**************************************************************************************************
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")

## Raster from Mitchell et al. 2019
EUNIS <- raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/substrate/EUNIS_class.tif")# 6 missing
BATHY<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/Mitchell_2018/Input data 1/Bathymetry.tif")# 6 missing
CURRENT<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/Mitchell_2018/Input data 1/DistanceCo.tif")# 6 missing
DISTCO<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/Mitchell_2018/Input data 1/Current_Sp.tif")# 6 missing
GRAVEL<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/Mitchelletal_2019/Sediment predictions 1/Predicted_Gravel_Fraction.tif")# 6 missing
MUD<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/Mitchelletal_2019/Sediment predictions 1/Predicted_Mud_Fraction.tif")# 6 missing
SAND<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/Mitchelletal_2019/Sediment predictions 1/Predicted_Sand_Fraction.tif")# 6 missing

## Raster from mNCEA
DISTCO<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/B_distcoast.tif")# 8 missing
avCURRENT<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/L_av_current.tif")# 44 missing
pkCURRENT<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/M_peak_current.tif")# 44 missing
CDM<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/CDM_mean.tif")# No missing
CHL<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/CHL2_mean.tif")# No missing
pkORBVEL<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/N_peak_wave_orb.tif")# 44 missing
pkWVCURR<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/O_peak_wave_current_stress.tif")# 44 missing
ampSAL<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/Q_ann_amp_bottom_sal.tif")# 36 missing
ampTEMP<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/S_ann_amp_bottom_temp.tif")# 36 missing
ROUGHNESS<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/roughness.tif")# 17 missing
GRAVEL<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/G_gravel.tif")# 8 missing
MUD<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/E_mud.tif")# 8 missing
SAND<-raster("C:/Users/cg05/OneDrive - CEFAS/Science/Operation/Data/Environment/mNCEA_GIS_rasters/F_sand.tif")# 8 missing

****************************************************************************************************************************
  
  
setwd("C:/Users/cg05/OneDrive - CEFAS/Science/Project - EU/MARBEFES/WP3/CVA")
irish.dat<-read.csv("./Irish Sea/input/Nov24_benthic_data.csv")
coordIrish<-irish.dat %>% 
  distinct(SiteNumber, Latitude, Longitude)

#write.csv(coordIrish, "C:/Users/cg05/OneDrive - CEFAS/Science/Project - EU/MARBEFES/WP3/CVA/coordIrish.csv")

# Extracting environment rasters from coordinate
response<-coordIrish
coordinates(response) <- ~Longitude+Latitude
mypoints <- SpatialPoints(response, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

rasStack <- stack(BATHY, CURRENT, DISTCO, GRAVEL, MUD, SAND)
rasStack1 <- stack(EUNIS)
rasStack2 <- stack(CDM)
rasStack3 <- stack(CHL)

Plat.data <- extract(rasStack, mypoints)
Plat.data1 <- extract(rasStack1, mypoints)
Plat.data2 <- extract(rasStack2, mypoints)
Plat.data3 <- extract(rasStack3, mypoints)

resp <- cbind(response, Plat.data, Plat.data1, Plat.data2, Plat.data3)
resp<-as.data.frame(resp)

irish.meta<-resp
#*************************************************************

#Adding back some information on survey and sieve mesh size
site<-irish.dat %>% 
  distinct(SiteNumber, id, Latitude, Longitude, Coordinates, EventDate, Year)

#write.csv(site, "C:/Users/cg05/OneDrive - CEFAS/Science/Project - EU/MARBEFES/WP3/CVA/Irish.site.csv")
irish.site.add<-read.csv("./Irish Sea/input/Irish.site.add.csv")
colnames(irish.site.add)[2]<-"id"

site.add<-left_join(site, irish.site.add, by = "id")

#write.csv(site.add, "C:/Users/cg05/OneDrive - CEFAS/Science/Project - EU/MARBEFES/WP3/CVA/Irish.site.added.csv")

# This Irish.site data has some info on survey and other things added to it
# The potential rep are the sited where "id" and "siteNumber" mismatch, sometimes they arguably seem to be genuine rep (i.e. same coordinates and date),
# sometimes they have different date but the same "SiteNumber" (not many are like that, I might just pick one id to align with the rest either way)
irish.site<-read.csv("./Irish Sea/input/Irish.site.csv")

save(irish.meta, irish.site, file="./IrishSea_BBT.RData")

#Checking the biology

