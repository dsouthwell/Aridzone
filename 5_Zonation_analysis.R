################################################################
#Preparing spatial analysis for Zonation analysis  
################################################################

library(raster)
library(rgdal)
library(raster)
library(lubridate)
library(rgeos)
library(oz)
library(spatstat)
library(viridis)
library(stringr)
library(sp)
library(rgdal)
options(stringsAsFactors = FALSE)

###############################################################
#Load in SDMs as crop to study area/mask out lakes etc
###############################################################

setwd("......")

cat <- raster("./input_data/unclipped_SDMs/Cat.tif", sep = "")
fox <- raster("./input_data/unclipped_SDMs/Fox.tif", sep = "")
dingo <- raster("./input_data/unclipped_SDMs/Dingo.tif", sep = "")
cattle <- raster("./input_data/unclipped_SDMs/Cattle.tif", sep = "")
rabbit <- stack("./input_data/clipped_SDMs/species_old.tif")[[5]]
ampurta <- stack("./input_data/clipped_SDMs/species_old.tif")[[6]]
roo <- raster("./input_data/unclipped_SDMs/macropus.sp.tif", sep = "")
dusky.mouse <- stack("./input_data/clipped_SDMs/species_old.tif")[[8]]
spinifex.mouse <- raster("./input_data/unclipped_SDMs/Hopping.mouse.Spinifex.tif", sep = "")
camel <- raster("./input_data/unclipped_SDMs/Camel.tif", sep = "")
emu <- raster("./input_data/unclipped_SDMs/Emu.tif", sep = "")
goanna <- raster("./input_data/unclipped_SDMs/Goanna.tif", sep = "")

species <- stack(cat, fox, dingo, cattle, rabbit, ampurta, roo, dusky.mouse, spinifex.mouse, camel, emu, goanna)

data <- read.csv("./input_data/SA_data_combined_covariates.csv")

index<-48

dat <- as.data.frame(cbind(data$x, data$y, data[,index], data$year,
                           #Bioclim climatic variables
                           data$Annual_mean_temp, data$Temp_annual_range, data$Temperature_seasonality, data$Isothermality, data$Mean_diurnal_range,            
                           data$Precipitation_seasonality, data$Max_temp_warmest_month, data$Min_temp_coldest_month, data$Precipitation_driest_month, data$Precipitation_wettest_month,
                           data$Mean_temp_coldest_quarter, data$Mean_temp_driest_quarter, data$Mean_temp_wettest_quarter, data$Mean_temp_warmest_quarter, data$Precipitation_coldest_quarter,
                           data$Precipitation_driest_quarter, data$Precipitation_warmest_quarter, data$Precipitation_wettest_quarter,
                           #Time dependent rainfall and NDVI values
                           data$rain.yr, data$rain.1yrlag, data$NDVI.yr,
                           #Topographic variables
                           data$elevation, data$roughness, data$slope,
                           #Environmental variables
                           data$clay, data$distanywat, data$grazing_dist, data$bdensity, data$calcrete, data$NVIS, data$Dogfence))
#Other info for analysis
#data$visits, data$idsite))


#colnames(dat)<-c("x", "y","presence", "elevation", "roughness","slope","clay","distwater","Rain","Rainlag1","Rainlag2","MeanTemp","NDVI","Min_temp","Rain_wettest","GrazingDist","bdensity","calcrete","NVIS","visits","idsite")

colnames(dat)<-c("lon", "lat","presence","Year", 
                 #Bioclim climatic variables
                 "Annual_mean_temp", "Temp_annual_range", "Temperature_seasonality", "Isothermality", "Mean_diurnal_range",            
                 "Precipitation_seasonality", "Max_temp_warmest_month", "Min_temp_coldest_month", "Precipitation_driest_month", "Precipitation_wettest_month",
                 "Mean_temp_coldest_quarter", "Mean_temp_driest_quarter", "Mean_temp_wettest_quarter", "Mean_temp_warmest_quarter", "Precipitation_coldest_quarter",
                 "Precipitation_driest_quarter", "Precipitation_warmest_quarter", "Precipitation_wettest_quarter",
                 #Time dependent rainfall and NDVI values
                 "rain.yr", "rain.1yrlag", "NDVI.yr",
                 #Topographic variables
                 "elevation", "roughness", "slope",
                 #Environmental variables
                 "clay", "distanywater", "grazing_dist", "bdensity", "calcrete", "NVIS", "fence"
                 #Other info for analysis
                  )#"visits", "idsite")

################################################################################################################
#Make polygon around sites to clip SDM predictions 
###############################################################################################################

sites <- SpatialPoints(cbind(dat$lon,dat$lat), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
sites <- spTransform(sites, CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
coords <- coordinates(sites)
ch <- chull(coords)
coords <- coords[c(ch, ch[1]),]
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)),proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#sp_poly <- spTransform(sp_poly, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
# e.g. CRS("+proj=longlat +datum=WGS84")
sp_poly_df <- SpatialPolygonsDataFrame(sp_poly, data=data.frame(ID=1))

#study.area <- gBuffer(sp_poly, width = 100000)
#plot(study.area, add=TRUE, col="red")

study_raster <- raster("./input_data/clipped_SDMs/Cat_clipped.tif", sep = "")

species.cr <- species
species.cr[is.na(study_raster)] <- NA

plot(species.cr[[2]])
plot(aus, add=TRUE)
points(sites, pch=19)

#Now mask out lakes and waterbodies
lakes <- raster("./input_data/lakes.tif")
species.cr[lakes==110] <- NA


writeRaster(species.cr[[1]], "./input_data/clipped_SDMs/Cat_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[2]], "./input_data/clipped_SDMs/Fox_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[3]], "./input_data/clipped_SDMs/Dingo_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[4]], "./input_data/clipped_SDMs/Cattle_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[5]], "./input_data/clipped_SDMs/Rabbit_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[6]], "./input_data/clipped_SDMs/Ampurta_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[7]], "./input_data/clipped_SDMs/Kangaroo_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[8]], "./input_data/clipped_SDMs/Dusky_hopping_mouse.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[9]], "./input_data/clipped_SDMs/Spinifex_hopping_mouse.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[10]], "./input_data/clipped_SDMs/Camel_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[11]], "./input_data/clipped_SDMs/Emu_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(species.cr[[12]], "./input_data/clipped_SDMs/Goanna_clipped.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

writeRaster(species.cr, "./input_data/clipped_SDMs/species_new.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

###############################################################
#Prepare Zonation files
###############################################################

#Load raster layer of each state
rr=raster(list.files("./input_data/clipped_SDMs",pattern=".tif$",full.names = T)[1]) 
mask <- rr
mask[!is.na(mask)] <- 1

fireBuffer <- condition <- mu <- mask

#Set new working directory and save raster layers for Zonation analysis
dir.create("./input_data/masks")

#condition layer
writeRaster(condition,file="masks/condition.tif",format="GTiff",options=c("TFW=YES"),overwrite=T)

#Admin units mask
writeRaster(fireBuffer,file="masks/analysis_mask_5kBuffer.tif",format="GTiff",options=c("TFW=YES"),overwrite=T)
writeRaster(mu,file="masks/ADMU_fireRegions_5kBuffer.tif",format="GTiff",options=c("TFW=YES"),overwrite=T)


#NOW SET UP INPUTS FOR SPECIFIC ZONATION RUNS: START HERE IF MAST AND CONDITION LAYERS MADE
##########################################################################################

library(raster)
library(lubridate)
library(rgeos)
options(stringsAsFactors = FALSE)

#species list

method<-"Revisions"

tifnames <- list.files(".input_files/clipped_SDMs",pattern=".tif$",full.names = T)
tifnames <- tifnames[c(1:11)] #Get rid of stack

#tifnames <- tifnames[c(1,6)]
wgt <- rep(1,11)

# writing Zonation species list file
nr=length(tifnames)
splist.table=data.frame(weight=wgt,alpha=rep(0,nr),c3=rep(0,nr),c4=rep(0,nr),exp=0.25,file=tifnames)
write.table(splist.table,file=paste0("specieslist_",method,".spp"),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

#Admin units weights..


zonvals=sort(na.omit(unique(values(mu))))

#may want to explore different admin unit weights..(how important are regional occurrences versus state/nation-wide distributions)
admu.table=data.frame(ID=zonvals[zonvals>0],G_A=round(1/max(zonvals),3),beta_A=0.5,name=c("AUS"))
write.table(admu.table,file="AMDU_descriptions.txt",quote=FALSE,sep=" ",row.names=FALSE)

admuWeights=matrix(1,nrow=nr,ncol=max(zonvals))
write.table(admuWeights,file="ADMU_weights.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

#groups file
groups=matrix(-1,nrow=nr,ncol=5)
groups[,2]=1
write.table(groups,file=paste0("groupsfile_",method,".txt"),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
#condition layer file

sink(file="condition_layers.txt")
cat("1  masks/condition.tif")
sink(NULL)

##create Z run settings file

sink(file=paste0("Zrun_settings_",method,".dat"))



cat(fill=TRUE,"
    [Settings]
    
    removal rule = 1
    warp factor = 1000
    edge removal = 1
    add edge points = 10000
    
    mask missing areas = 1
    area mask file = masks/analysis_mask_5kBuffer.tif
    
    use groups = 1")
cat("groups file =",paste0("groupsfile_",method,".txt"),"\n")
cat(fill=TRUE,
    "use condition layer = 1
    condition file = condition_layers.txt
    
    use ADMUs = 1
    ADMU descriptions file = ADMU_descriptions.txt
    ADMU layer file = masks/ADMU_fireRegions_5kBuffer.tif
    ADMU weight matrix = Inputs/ADMU_weights.txt
    ADMU mode = 2
    Mode 2 global weight = 0.5
    row count for per ADMU output curves = 0
    output weighted range size corrected richness = 0
    ")
sink(NULL)

#write Zonation batch file (assumes Zonation exe (zig4.exe) is in R working directory (where all input files have been written to)

outname=paste0("Zout/Z_",method,"_5k") #base name of Z outputs

dir.create("Zout")
sink(file=paste0("Z_",method,"_5K.bat"))  # name the .bat file
cat("call zig4 -r ",paste0("Zrun_settings_",method,".dat"),paste0("specieslist_",method,".spp"), outname ,"0 0 1 0")
sink()


# Now run the batch file in Zonation gui, or by double clicking .bat file created above
# can also call the .bat directly from R with 'system' function (but then your R session is busy while Z runs, which may be a long time!!)

####

#setwd("C:/Darren/12_Postdoc/7_Desert_monitoring_project/Paper2_SDMs_power/SDMs_clipped/")  # set working directory
outname="Z_Revisions_5k"  # the base name of the Zonation outputs, see line 220 of Zonation inputs script
Zoutputs=list.files("Zout",pattern=outname,full.names=T) #get the names of all ouputs for this analysis (assuming in a ...'/Zout' directory)
r1=raster(Zoutputs[grep("rank.compressed.tif",Zoutputs)])
writeRaster(r1, "./input_layers/Zout/Zonation_equal_weights_raw.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
r1[r1<0.9] <- NA
writeRaster(r1, "./input_layers/Zout//Zonation_equal_weights.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


