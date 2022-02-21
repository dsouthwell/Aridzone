#########################################################################################################
# Assessing the univariate variable importance of each covariate using GAMs in 5-fold spatial block cross-validation design.
##########################################################################################################

#----------------------------------------------------------------------------------------------------------

# load required packages

library(sperrorest)
library(dismo)
library(mgcv)
library(dplyr)
library(gamm4)

#----------------------------------------------------------------------------------------------------------

#Load in data
setwd("......")
data <- read.csv("./input_data/SA_data_combined_covariates.csv") 

#Identify target species
species <- c("Camel",
               "Cat",
               "Cattle",
               "Dingo",
               "Emu",
               "Fox",
               "Goanna",
               "Hopping.mouse.Dusky",
               "Hopping.mouse.Spinifex",
               "Mulgara.Crest.tailed",
               "Rabbit",
               "Macropus.sp")

#Create empty results matrix
results <- matrix(NA,ncol=length(species), nrow=25)
colnames(results) <- species

#Loop trough species
for (i in 1:length(species)){
  print(i)
  #Find column for species i
  index <- which(colnames(data)==species[i]) 
  
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
                             data$clay, data$distanywat, data$grazing_dist, data$bdensity, data$calcrete, data$Dogfence, data$NVIS))#,
                             #Other info for analysis
                             #data$visits, data$idsite))
  
  colnames(dat)<-c("x", "y","presence","Year", 
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
                   "clay", "distanywater", "grazing_dist", "bdensity", "calcrete", "fence", "NVIS"
                   #Other info for analysis
                   #"visits", "idsite")
  )
  
  #Identify how many sites were surveyed more than once in the same year
  dups<-duplicated(dat[,c('x','y','Year')])
  sum(dups) 
  
  #Collapse within year repeat surveys  
  dat <- distinct(dat) # remove one of the sites with 0,0 / 1,1
  dat <- dat[order(-dat$presence),] # reorder so 1s come first
  dat<-dat[!duplicated(dat[,c('x','y','Year')]),]
  
  #Remove any rows with NA values
  row.has.na <- apply(dat, 1, function(x){any(is.na(x))})
  sum(row.has.na) #Record number of rows with NA values
  dat <- dat[!row.has.na,] #Remove rows with NA values
  
  #Remove Year column
  dat <- dat[,-4]
  dat[,4:32] <- apply(dat[,4:32],2,function(x){scale(x)})
  
  #Split the data to make training and testing datasets
  sample <- sample.int(n = nrow(dat), size = floor(.8*nrow(dat)), replace = F)

  avinames.50presences <- colnames(dat)[3] 	
  #Identify columns with predictors
  avi.pred.names <- colnames(dat)[4:(ncol(dat))]  	
  avi.pred.names <-  c("elevation","roughness","clay","calcrete","distanywater","Annual_mean_temp","Temp_annual_range","Temperature_seasonality",      
                       "Isothermality","Mean_diurnal_range","Precipitation_seasonality","Max_temp_warmest_month","Min_temp_coldest_month",
                       "Precipitation_wettest_month","Mean_temp_coldest_quarter","Mean_temp_driest_quarter",
                       "Mean_temp_wettest_quarter","Mean_temp_warmest_quarter","Precipitation_coldest_quarter",
                       "Precipitation_warmest_quarter","Precipitation_wettest_quarter","grazing_dist","NDVI.yr","rain.yr","rain.1yrlag","fence","NVIS")
  
  avi.pred.names <-  c("elevation","roughness","clay","calcrete","distanywater","Annual_mean_temp","Temp_annual_range","Temperature_seasonality",      
                       "Isothermality","Mean_diurnal_range","Max_temp_warmest_month","Min_temp_coldest_month",
                       "Precipitation_wettest_month",
                       "grazing_dist","NDVI.yr","rain.yr","rain.1yrlag","fence","NVIS")
  
  # Compute univariate variable importance
  ## Univariate GAMs for loop for creating R squareds and AIC values or each variable and response to pres-abs data
  #names(data)
  aicoutput.allpredictors <- data.frame()
  rsqoutput.allpredictors <- data.frame()
  mypredictors <- c("s(TimeSinceFire)", "s(DistanceToRoad)", "s(DistanceToWater)", "s(bio1.var)", "s(bio5.var)", "s(bio6.var)", "s(bio10.var)", "s(bio11.var)", "s(bio12.var)", "s(Slope)", "s(EVC)", "s(Elevation)")
  
  mypredictors <- c("s(elevation)","s(roughness)","s(clay)","s(calcrete)","s(distanywater)","s(Annual_mean_temp)","s(Temp_annual_range)","s(Temperature_seasonality)",      
                    "s(Isothermality)","s(Mean_diurnal_range)","s(Precipitation_seasonality)","s(Max_temp_warmest_month)","s(Min_temp_coldest_month)",
                    "s(Precipitation_wettest_month)","s(Mean_temp_coldest_quarter)","s(Mean_temp_driest_quarter)",
                    "s(Mean_temp_wettest_quarter)","s(Mean_temp_warmest_quarter)","s(Precipitation_coldest_quarter)",
                    "s(Precipitation_warmest_quarter)","s(Precipitation_wettest_quarter)","s(grazing_dist)","s(NDVI.yr)","s(rain.yr)","s(rain.1yrlag)")
  
  mypredictors <- c("s(elevation)","s(roughness)","s(clay)","s(calcrete)","s(distanywater)","s(Annual_mean_temp)","s(Temp_annual_range)","s(Temperature_seasonality)",      
                    "s(Isothermality)","s(Mean_diurnal_range)","s(Max_temp_warmest_month)","s(Min_temp_coldest_month)",
                    "s(Precipitation_wettest_month)","s(grazing_dist)","s(NDVI.yr)","s(rain.yr)","s(rain.1yrlag)")
  
  
  for (m in mypredictors)
  {
    speciesmod <- paste0("presence ~ ", m)
    PA.gamm.m <- gamm4(as.formula(speciesmod), family = binomial, data = dat)
    aic.output <- AIC(PA.gamm.m$mer)
    rsq.output <- summary(PA.gamm.m$gam)$r.sq
    aic.output.table<- as.data.frame(aic.output)
    rsq.output.table <- as.data.frame(rsq.output)
    aic.output.table$predictor <- m
    rsq.output.table$predictor <- m
    aicoutput.allpredictors <- rbind(aicoutput.allpredictors, aic.output.table)
    rsqoutput.allpredictors <- rbind(rsqoutput.allpredictors, rsq.output.table)
  }
  
  modeloutput.allpredictors <- full_join(aicoutput.allpredictors, rsqoutput.allpredictors, by = "predictor")
  
  select07 <- function(imp, X, threshold=0.7, method="spearman")
  {
    cm <- cor(X, method=method)
    sort.imp <- colnames(X)[order(imp,decreasing=T)]
    
    pairs <- which(abs(cm)>= threshold, arr.ind=T) # identifies correlated variable pairs
    index <- which(pairs[,1]==pairs[,2])           # removes entry on diagonal
    pairs <- pairs[-index,]                        # -"-
    
    exclude <- NULL
    for (i in seq_len(length(sort.imp)))
    {
      if ((sort.imp[i] %in% row.names(pairs))&
          ((sort.imp[i] %in% exclude)==F)) {
        cv<-cm[setdiff(row.names(cm),exclude),sort.imp[i]]
        cv<-cv[setdiff(names(cv),sort.imp[1:i])]
        exclude<-c(exclude,names(which((abs(cv)>=threshold)))) }
    }
    
    sort.imp[!(sort.imp %in% unique(exclude)),drop=F]
  }
  
  avi.univar.cv.d2 <- data.frame(modeloutput.allpredictors$rsq.output)
  names(avi.univar.cv.d2) <- avinames.50presences
  rownames(avi.univar.cv.d2) <- avi.pred.names[1:17]
  
  avi.select07.cv <- apply(avi.univar.cv.d2, 2,FUN=function(x){select07(x, dat[, avi.pred.names])[1:25]})
  results[,i] <-avi.select07.cv 
  
} #End loop

#Save ranked univariate responses
write.csv(results, file = "./input_data/Predator_predictors.csv", row.names=FALSE)


