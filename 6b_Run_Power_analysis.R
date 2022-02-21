
#SCENARIO 1 - #VARY NUMBER OF SITES, 1 REPEAT VISIT

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/park.tif") #Load in remote layer
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- remote <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
remote[] <- 1
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

SITES <- seq(50,700,50)
#SITES <- c(50,75)

for (zz in 1:length(SITES)) {
  
  prespecify.sites <- TRUE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- SITES[zz]} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- FALSE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  #n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(1,1,1,1) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "increasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_increasing_sites",n.sites,"_Freq1_repeats1",".RData", sep = ""))
}






#SCENARIO 2 - #VARY NUMBER OF SITES, 2 REPEAT VISITS

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/park.tif") #Load in remote layer
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- remote <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
remote[] <- 1
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

SITES <- seq(50,700,50)
#SITES <- c(50,75)

for (zz in 1:length(SITES)) {
  
  prespecify.sites <- TRUE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- SITES[zz]} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- FALSE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "increasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_increasing_sites",n.sites,"_Freq1_repeats2",".RData", sep = ""))
}





#SCENARIO 3 - #VARY NUMBER OF SITES, 3 REPEAT VISITS

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/park.tif") #Load in remote layer
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- remote <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
remote[] <- 1
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

SITES <- seq(50,700,50)
#SITES <- c(50,75)

for (zz in 1:length(SITES)) {
  
  prespecify.sites <- TRUE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- SITES[zz]} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- FALSE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  #n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(3,3,3,3) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "increasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_increasing_sites",n.sites,"_Freq1_repeats3",".RData", sep = ""))
}




#SCENARIO 1 - #VARY NUMBER OF SITES, 1 REPEAT VISIT, DECREASING TREND

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/park.tif") #Load in remote layer
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- remote <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
remote[] <- 1
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

SITES <- seq(50,700,50)
#SITES <- c(50,75)

for (zz in 1:length(SITES)) {
  
  prespecify.sites <- TRUE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- SITES[zz]} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- FALSE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  #n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(1,1,1,1) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "decreasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_decreasing_sites",n.sites,"_Freq1_repeats1",".RData", sep = ""))
}






#SCENARIO 2 - #VARY NUMBER OF SITES, 2 REPEAT VISITS

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/park.tif") #Load in remote layer
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- remote <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
remote[] <- 1
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

SITES <- seq(50,700,50)
#SITES <- c(50,75)

for (zz in 1:length(SITES)) {
  
  prespecify.sites <- TRUE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- SITES[zz]} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- FALSE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "decreasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_decreasing_sites",n.sites,"_Freq1_repeats2",".RData", sep = ""))
}





#SCENARIO 3 - #VARY NUMBER OF SITES, 3 REPEAT VISITS

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/park.tif") #Load in remote layer
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- remote <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
remote[] <- 1
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

SITES <- seq(50,700,50)
#SITES <- c(50,75)

for (zz in 1:length(SITES)) {
  
  prespecify.sites <- TRUE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- SITES[zz]} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- FALSE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  #n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(3,3,3,3) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "decreasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_decreasing_sites",n.sites,"_Freq1_repeats3",".RData", sep = ""))
}



#SCENARIO 4 - #OPTIMISE PLACEMENT OF SITES USING ZONATION UNEQUAL WEIGHTS

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/Zonation_two_sp_weights.tif") #Load in remote layer
remote[remote>0] <- 1
remote[is.na(remote)] <- 0
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

#SITES <- seq(150)
REPEATS <- c(1,2,3)

for (zz in 2:2) {
  
  prespecify.sites <- FALSE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  #sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- 200} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- TRUE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- rep(REPEATS[zz],4) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "decreasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  occ.time[is.na(occ.time)] <- 0
  det.method1.time[is.na(det.method1.time)] <- 0
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_decreasing_optimal_unequalweight_Freq1_repeats",REPEATS[zz],".RData", sep = ""))
  
}


#SCENARIO 4 - #OPTIMISE PLACEMENT OF SITES USING ZONATION EQUAL WEIGHTS

#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/Zonation_equal_weights.tif") #Load in remote layer
remote[remote>0] <- 1
remote[is.na(remote)] <- 0
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

#SITES <- seq(150)
REPEATS <- c(1,2,3)

for (zz in 2:2) {
  
  prespecify.sites <- FALSE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  #sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- 200} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- TRUE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- rep(REPEATS[zz],4) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "decreasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_NATURAL_VARIATION.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  occ.time[is.na(occ.time)] <- 0
  det.method1.time[is.na(det.method1.time)] <- 0
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_decreasing_optimal_equalweight_Freq1_repeats",REPEATS[zz],".RData", sep = ""))
  
}



###############################################################################
#Test the effect of only doing repeats on a fraction of sites each year (note save a new functions file)
###############################################################################


#####################################################################################################################
#SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES
###################################################################################################################
# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au
# This package calculates the statistical power to detect simulated trends in occupancy for multiple species in dynamic landscapes
# Simulations require an occupancy raster layer for each species and estimates of detectability for up to three different types of detection methods
# Users can specify the length of the monitoring program, the years in which monitoring occurs and the number of repeat visits to a site within a survey year
# Pre-specified monitoring sites can be loaded into the package or selected randomly across space at the start of each simulation
# Monitoring sites can also be startified across environmental strata, be positioned on cells with the highest species richness, or selected manually from the raster layers
# Monitoring sites can be divided between remote and non-remote areas across the landscape
# Occupancy and detectability layers can remain static over time or change dynamically in response to simulated fire events
# Fire events are determined by a hazard function which depends on the time since a cell last burned
# The relationship between fire events and species occupancy and detectability is determined by statistical model
# Users can calculate power at a landscape level as well as power at smaller regional level management units

#######################################################################################################################################
#STEP 1 - Load packages

#Loading this package will automatically load other packages required for simulations
rm(list=ls()) #Remove any pre-existing files in the global environment
library(unmarked)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
#profvis({

##########################################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species

#Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
#should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
#For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
#Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
#All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch 
#Cells not considered for monitoriong should be assigned an NA value
#The package contains a sample dataset from Kakadu and Nitmiluk National Parks in the Northern Territory, Australia. 
#To get an idea about how your raster layers should look, type out either occ, det.method1, det.method2 or det.method3 before loading in your own raster stacks

setwd("~/Desert_monitoring/Data3")
occ <- stack("species_new.tif") #Load in species occupancy layers combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method1 <- stack("species_new.tif") #Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method2 <- stack("species_new.tif") #Load in species detectability layers for method 2 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method3 <- stack("species_new.tif") #Load in species detectability layers for method 3 combined into a raster stack. Note, if you only assess one species this should still be a stack
det.method4 <- stack("species_new.tif")

is.na(det.method1) <- is.na(occ) <- 0

det.method1[[1]][det.method1[[1]] > -1000] <- 0.36 #cat
det.method1[[2]][det.method1[[2]] > -1000] <- 0.85 #fox
det.method1[[3]][det.method1[[3]] > -1000] <- 0.28 #dingo
det.method1[[4]][det.method1[[4]] > -1000] <- 0.44 #cattle
det.method1[[5]][det.method1[[5]] > -1000] <- 0.79 #rabbit
det.method1[[6]][det.method1[[6]] > -1000] <- 0.92#ampurta
det.method1[[7]][det.method1[[7]] > -1000] <- 0.41 #roo
det.method1[[8]][det.method1[[8]] > -1000] <- 0.92 #dusky mouse
det.method1[[9]][det.method1[[9]] > -1000] <- 0.92 #spinifex mouse
det.method1[[10]][det.method1[[10]] > -1000] <- 0.44 #camel
det.method1[[11]][det.method1[[11]] > -1000] <- 0.44 #emu
det.method1[[12]][det.method1[[12]] > -1000] <- 0.44 #goanna
det.method2[det.method2 > -1000] <- 0.7
det.method3[det.method3 > -1000] <- 0.7
det.method4[det.method4 > -1000] <- 0.7

plot(occ[[1]]) #Plot occupancy for the first species. Change the number to plot a different species

##########################################################################################################################################
#STEP 3 - Load files specifying which methods apply to each species and a statistical model relating occupancy and detectability to covariates

#Simulations require a list that specifies which detection methods are relevant to each species. 
#A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
#If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
#The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites
#The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates
#Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's

species.list <- read.csv("~/Desert_monitoring/Data3/Species_list.csv",stringsAsFactors=FALSE) #Load in which methods are relevant to each species
species.cov <- read.csv("~/Desert_monitoring/Data3/Species_covariates.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

##########################################################################################################################################
#STEP 4 - Load in additional raster layers for simulations

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas
#Type 'remote to have a look at the example remote layer for Kakadu. Here, remote areas are all cells greater than 1 km from roads
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks <- raster("~/Desert_monitoring/Data3/park.tif") #Load in parks layer
crs(parks) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
remote <- raster("~/Desert_monitoring/Data3/park.tif") #Load in remote layer
veg <- raster("~/Desert_monitoring/Data3/park.tif") #Load in veg layer

parks <- remote <- occ[[1]]
#parks[parks<0.5] <- 1
#parks[parks<1] <- 2
remote[] <- 1
parks[] <- 1
####################################################################################################################################
#Step 4 - Decide whether to model fire and load in fire history layer 

#if you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
#The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
#In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
#The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
#The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

model.fire <- FALSE #Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise
if (model.fire == TRUE) {
  covar <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack
  fire <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in fire history raster stack
}

############################################################################################################################################
#Step 5: Define where to monitor

#A wide variety of monitoring scenarios can be evaulated. Power can be assessed at:
#1) Pre-selected sites - that is, monitoring can be repeated at the same locations each simulation. 
#2) A subset of pre-selected sites - a subset of pre-selected sites can be randomly selected for monitoring each simulation
#3) randomly selected sites - sites are selected randomly throughout the landscape each simulation
#4) randomly stratified sites - sites are randomly selected across the landscape, but evenly across environmental strata
#5) regions with highest species richness - sites can be located on cells with the highest combined occupancy
#6) maunally selcted sites - users can click on a raster layer to identify monitoring sites
#Note, when sites are selected randomly, they can be divided between remote and non-remote areas

SITES <- 200
#SITES <- c(50,75)

for (zz in 1:length(SITES)) {
  
  prespecify.sites <- TRUE #If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE
  sites<-read.csv("~/Desert_monitoring/Data3/Sites_EqualArea.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data
  sites<-sites[-c(200,604,605),]
  all.prespecify.sites <- FALSE #If you want to monitor all of the pre-selected sites set to TRUE, otherwise set to FALSE. If FALSE a subset will be selected at random each simulation 
  if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list
    n.sites <- SITES[zz]} #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape
  R <- 1 #This defines the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. If you don't have remote or non-remote areas, make all cells in your remote layer equal to 1 and make R=1 
  new.site.selection <- "random" #If you don't pre-select sites, specify how sites will be selected. Choose from 'random', 'stratified','maxocc' or ''manual'
  #random - sites randomly selected given ratio of remote to non-remote 
  #stratified - sites stratified across layers of a chosen raster layer
  #maxocc - sites positioned on cells with highest species richness (i.e. highest total occupancy for all species)
  #manual - select sites by clicking on the map
  plotter <- FALSE #Set to TRUE to plot monitoring sites at the start of each simuation
  
  ##################################################################################################################################################
  #Step 6: Define when to monitor
  
  #Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each sie for each detection method
  
  Tmax <- 15 #Define length of monitoring program in years
  s.years <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  #n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #Define the years when monitoring occurs. Note, the last value must be equal to Tmax.
  n.method <- c(2,2,2,2) #Number of repeat visits to sites using method 1
  #n.method2 <- 4 #Number of repeat visits to sites using method 2
  #n.method3 <- 6 #Number of repeat visits to sites using method 3
  
  ################################################################################################################################################
  #Step 7: Define power analysis
  
  #Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
  #whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations
  
  park.level <- FALSE #Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level 
  decline <- "random" #Random is the only option - it means we simulate a constant decline in occupancy across space
  trend <- "decreasing" #Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend
  model.variation <- FALSE
  variation <- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376)
  alpha <- 0.05 #Set significance level. Choose from 0.01, 0.05 or 0.1
  two.tailed <- TRUE #Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test
  effect.size <- c(0.1,0.3,0.5,0.7,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes
  nsims <- 100 #Define the number of simulations used to calculate power
  
  ###############################################################################################################################################
  #Step 8: Create empty arrays to store results of simualtions
  
  n.species <- nlayers(occ) #Calculates the number of species
  n.park <- unique(parks) #Calculates the number of parks
  det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays
  
  source("~/Desert_monitoring/Rcode/2_Desert_monitoring_functions_vary_repeats.R") #Make sure the most recent source file is loaded
  
  ###############################################################################################################################################
  #Step 10: Set up monitoring sites
  
  #This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
  #Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed, or if 'maxocc' is selected this function is not called again
  
  xy.sites <- select.sites(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) #Define monitoring sites
  xy.sites2 <- SpatialPoints(xy.sites, proj4string=CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  park.ID <- parks[cellFromXY(parks, xy.sites2)] #Extract parks values at monitoring sites
  park.ID[is.na(park.ID)] <- 3
  xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
  occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
  for (ss in 1:n.species){
    if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
    if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
    if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
    if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
  }
  
  #################################################################################################################################################
  #Step 11: Check raster layers and run power analysis
  
  #The last step us to check all raster layers have the same dimensions and then run the power analysis. 
  #A warning message will be displayed if there is a mismatch between rasters
  #This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
  #Simulations can be sped up by running them in parallel. Excessively large datasets should be run on a computer with at least 10 cores   
  
  
  check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 
  
  n.cores <- 10 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
  cl <- makeCluster(n.cores, outfile="log.txt") #initiate clusters
  registerDoParallel(cl) #initiate clusters
  #clusterEvalQ(cl, source("~/Desert_monitoring/Rcode/1_Desert_monitoring_functions_NATURAL_VARIATION.R"))
  start.time<-proc.time() #Record start time
  pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("raster","unmarked")) %dopar%{ #Runs the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation
    sapply(effect.size, run.power, 
           nsims=nsims, 
           alpha=alpha, 
           decline=decline, 
           Tmax=Tmax, 
           s.years=s.years, 
           trend=trend,
           sites=sites, 
           model.fire=model.fire, 
           species.list=species.list, 
           park.level=park.level, 
           xy.sites=xy.sites, 
           R=R,
           two.tailed=two.tailed, 
           plotter=plotter, 
           n.sites=n.sites, 
           n.species=n.species, 
           n.park=n.park, 
           all.prespecify.sites=all.prespecify.sites,
           prespecify.sites=prespecify.sites, 
           new.site.selection=new.site.selection,
           occ.time=occ.time, 
           det.method1.time=det.method1.time, 
           det.method2.time=det.method2.time, 
           det.method3.time=det.method3.time, 
           det.method4.time=det.method4.time,
           n.method=n.method,
           occ=occ,
           det.method1=det.method1,
           det.method2=det.method2,
           det.method3=det.method3,
           det.method4=det.method4,
           veg=veg,
           remote=remote,
           parks=parks,
           variation=variation,
           model.variation=model.variation)
  }
  stopCluster(cl)
  time.elapsed<-proc.time()-start.time #Record end time
  time.elapsed #Print elapsed time
  
  ################################################################################################################################
  #Step 12: PLot results
  
  Results <- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 
  
  #################################################################################################################################
  #Step 13: Save results
  
  setwd("~/Desert_monitoring/Results")
  save(pwr, file = paste("New_AZM_Natvar_sims1000_decreasing_sites",n.sites,"_Freq1_repeats2_1sites",".RData", sep = ""))
}






