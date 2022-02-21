#########################################################################################
#FIT ENSEMBLE SDMS TO SPECIES
#########################################################################################

library(RColorBrewer)
library(raster)
library(sp)
library(rgdal)
library(plyr)
library(dplyr)
library(dismo)

####################################################################
#Load in top ranked predictors
setwd(".......")
predictors  <- read.csv("./input_data/Predator_predictors.csv", stringsAsFactors=FALSE)

#Load in map of Australia
aus <- readOGR(dsn="......a",layer="cstauscd_r")
aus <- spTransform(aus, CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

envdat <- stack("./input_data/environmental_covariates.tif")

names(envdat) <- c(#Bioclim climatic variables
  "Annual_mean_temp", "Temp_annual_range", "Temperature_seasonality", "Isothermality", "Mean_diurnal_range",            
  "Precipitation_seasonality", "Max_temp_warmest_month", "Min_temp_coldest_month", "Precipitation_driest_month", "Precipitation_wettest_month",
  "Mean_temp_coldest_quarter", "Mean_temp_driest_quarter", "Mean_temp_wettest_quarter", "Mean_temp_warmest_quarter", "Precipitation_coldest_quarter",
  "Precipitation_driest_quarter", "Precipitation_warmest_quarter", "Precipitation_wettest_quarter",
  #Time dependent rainfall and NDVI values
  "rain.yr", "rain.1yrlag", "NDVI.yr",
  #Topographic variables
  "elevation", "roughness", "slope",
  #Environmental variables
  "clay", "distanywater", "grazing_dist", "bdensity", "calcrete", "NVIS", "fence",
  #Other info for analysis
  "visits", "lon","lat")            


###############################################################################################################
#Load in data to fit model
###############################################################################################################

data <- read.csv("./input_data/SA_data_combined_covariates.csv")

predators <- colnames(predictors)
predators[12] <- "Macropus.sp"


train.performance <- list()
crossval.performance <- list()
ensemble.performance <- list() 

for (i in 1:length(predators)){ 
  
  index <- which(colnames(data)==predators[i]) 
  
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
  
  #Remove duplicates
  dups<-duplicated(dat[,c('lon','lat','Year')])
  sum(dups) #175 dups, multiple recordings for same location (all have the same PA value)
  
  #dat <- add_count(dat, x, y)
  
  #Collapse within year repeat surveys  
  dat <- distinct(dat) # remove one of the sites with 0,0 / 1,1
  dat <- dat[order(-dat$presence),] # reorder so 1s come first
  dat<-dat[!duplicated(dat[,c('lon','lat','Year')]),]
  
  #Remove rows with NA values
  row.has.na <- apply(dat, 1, function(x){any(is.na(x))})
  sum(row.has.na) #Record number of rows with NA values
  dat <- dat[!row.has.na,] #Remove rows with NA values
  
  dat <- dat[,-4]
  
  sum(dat$presence)
  
  #Scale continuous covariates
  dat[,4:32] <- apply(dat[,4:32],2,function(x){scale(x)})
  
  #Split the data to make training and testing datasets
  sample <- sample.int(n = nrow(dat), size = floor(.80*nrow(dat)), replace = F)
  train <- dat[sample, ]
  test  <- dat[-sample, ]
  
  # Load data. Our datasets were called the following:
  aviCH.forest.train <- train		# birds - training data for SDM and JSDM calibration
  aviCH.forest.test <- test		# birds - test data for evaluating community predictions
  
  ####################################################################
  #Fit models for species i
  ####################################################################
  
  ttetrix.env <- dat
  colnames(ttetrix.env)[c(1:3)] <- c("lon","lat","occ")
  
  sum(ttetrix.env$occ)
  
  n.pred <- floor(sum(ttetrix.env$occ)/40)
  
  
  if (n.pred < 20) {pred <- as.character(predictors[,i])[1:n.pred]} else {pred <- as.character(predictors[,i])}
  
  pred <- pred[!is.na(pred)]
  
  if (i==10) {pred <- pred[-2]}
  
  if (i %in% c(1,2,3,4,5,6,12)) {pred <- c(pred, "fence")} #for cattle, dingo, emus, kangaroos, foxes and cats
  
  
  #Fit models
  
  # Fit generalised linear model (GLM)
  m.glm <- step(glm( as.formula(
    paste('occ ~',paste(pred,paste0('+ I(',pred,'^2)'),collapse=' + '))),
    family='binomial', data=ttetrix.env))
  
  library(rpart)
  m.cart <- rpart( as.formula(
    paste('occ ~',paste(pred,collapse=' + '))),
    data=ttetrix.env, control=rpart.control(minsplit=20,xval=10))
  
  # Fit Random Forest (RF)
  library(randomForest)
  m.rf <- randomForest( x=ttetrix.env[,pred], y=as.factor(ttetrix.env$occ),
                        ntree=1000, importance =T, nodesize=20)
  
  # Fit boosted regression tree (BRT)
  library(gbm)
  m.brt <- gbm.step(data = ttetrix.env,
                    gbm.x = pred,
                    gbm.y = 'occ',
                    family = 'bernoulli',
                    tree.complexity = 2,
                    bag.fraction = 0.75,
                    learning.rate = 0.01)
  
  }
  
  ####################################
  #Functions
  
  make.preds <- function(model, newdata) {
    require(dismo)
    require(gam)
    require(rpart)
    require(randomForest)
    require(gbm)
    # require(maxnet)
    switch(class(model)[1],
           Bioclim = predict(model, newdata),
           Domain = predict(model, newdata),
           glm = predict(model, newdata, type='response'),
           Gam = predict(model, newdata, type='response'),
           rpart = predict(model, newdata),
           randomForest = predict(model, newdata, type='prob')[,2],
           gbm = predict.gbm(model, newdata,
                             n.trees=model$gbm.call$best.trees, type="response")
           #maxnet = predict(model, newdata, type="logistic")
    )
  }
  
  # Second, a function for deriving cross-validated predictions.
  # The function partitions the data into k folds, determines the model
  # algorithm, updates the model for the new training data and makes
  # predictions to the hold-out data.
  crossval.preds <- function(model, traindat, colname.species, colname.pred,
                             env.r, colname.coord, kfold) {
    require(dismo)
    require(gam)
    require(rpart)
    require(randomForest)
    require(gbm)
    #require(maxnet)
    # Make k-fold data partitions
    ks <- kfold(traindat, k = kfold, by = traindat[,colname.species])
    cross.val.preds <- data.frame(row = row.names(traindat),
                                  cross.val.preds = numeric(length = nrow(traindat)))
    for(i in seq_len(kfold)){
      cv.train <- traindat[ks!=i,]
      cv.test <- traindat[ks==i,]
      # Because we used the gbm.step() for BRTs, we need a small work-around:
      if (class(model)[1]=='gbm') {
        cv.train.gbm <- cv.train;
        names(cv.train.gbm)[names(cv.train.gbm)==colname.species] <-
          model$response.name
      }
      
      # We update the model for the new training data
      modtmp <- switch(class(model)[1],
                       Bioclim = bioclim(env.r[[colname.pred]],
                                         cv.train[cv.train[, colname.species]==1, colname.coord]),
                       Domain = domain(env.r[[colname.pred]],
                                       cv.train[cv.train[, colname.species]==1, colname.coord]),
                       glm = update(model, data=cv.train),
                       Gam = update(model, data=cv.train),
                       rpart = update(model, data=cv.train),
                       randomForest = update(model, data=cv.train),
                       gbm = gbm(model$call, 'bernoulli',
                                 data=cv.train.gbm[,c(colname.pred,model$response.name)],
                                 n.trees=model$gbm.call$best.trees)
                       #maxnet = maxnet(p= cv.train[,colname.species],
                       #data= cv.train[,colname.pred])
      )
      # We make predictions for k-fold:
      if (class(model)[1]=='gbm') {
        cross.val.preds[which(ks==i),2] <-
          predict.gbm(modtmp, cv.test[, colname.pred],
                      n.trees=model$gbm.call$best.trees, type="response")
      } else {
        cross.val.preds[which(ks==i),2] <-
          make.preds(modtmp, cv.test[, colname.pred])
      }
    }
    cross.val.preds[,2]
  }
  
  
  # Third, a function for calculating performance measures:
  calc.eval <- function(dat, colname.species, preds, thresh.method='MaxSens+Spec'){
    require(PresenceAbsence)
    require(dismo)
    # Helper functions - we had defined them before in "SDM_functions.r". However,
    # if we want this code to be independent, we have to insert the functions here again.
    # True Skill Statistic:
    TSS = function(cmx){
      PresenceAbsence::sensitivity(cmx, st.dev=F) +
        PresenceAbsence::specificity(cmx, st.dev=F) - 1
    }
    
    # Explained deviance:
    d.square <- function(obs, pred, family='binomial'){
      if (family=='binomial') {pred <-
        ifelse(pred<.00001,.00001,ifelse(pred>.9999,.9999,pred))}
      null.pred <- rep(mean(obs), length(obs))
      1 - (calc.deviance(obs, pred, family=family) /
             calc.deviance(obs, null.pred, family=family))
    }
    
    thresh.dat <- data.frame(ID=seq_len(nrow(dat)),
                             obs = dat[, colname.species],
                             pred = preds)
    thresh <- optimal.thresholds(DATA= thresh.dat)
    cmx.maxSSS <- cmx(DATA= thresh.dat, threshold=thresh[thresh$Method==thresh.method,2])
    data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
               TSS = TSS(cmx.maxSSS),
               Sens = PresenceAbsence::sensitivity(cmx.maxSSS, st.dev=F),
               Spec = PresenceAbsence::specificity(cmx.maxSSS, st.dev=F),
               PCC = PresenceAbsence::pcc(cmx.maxSSS, st.dev=F),
               D2 = d.square(thresh.dat$obs, thresh.dat$pred),
               thresh = thresh[thresh$Method==thresh.method,2])
  }
  #------------------
  
  # Now we make predictions on training data and assess internal model performance
  our.models <- c('m.glm', 'm.rf', 'm.brt')
  
  train.preds <- sapply(our.models, FUN=function(m){make.preds(eval(parse(text=m)), ttetrix.env[, pred])})
  
  train.perf <- sapply(our.models, FUN=function(x){calc.eval(ttetrix.env,'occ',train.preds[,x])})
  
  # Now we make cross-validated predictions and assess cross-validated model performance:
  crossval.preds <- sapply(our.models, FUN=function(m){crossval.preds(eval(parse(text=m)), ttetrix.env,
                                                                      colname.species='occ', colname.pred=pred, env.r=envdat,
                                                                      colname.coord=c('lon','lat'), kfold=10)})
  
  crossval.perf <- sapply(our.models, FUN=function(x){calc.eval(ttetrix.env,'occ',crossval.preds[,x])})
  
  ###############################################################################################################3
  #Making ensembles
  
  # Mean of probabilities
  mean.prob <- rowMeans(crossval.preds)
  
  # Median of probabilities
  median.prob <- apply(crossval.preds,1,median)
  
  # Weighted mean of probabilities, weighted by TSS
  #wmean.prob <- apply(crossval.preds,1,weighted.mean, w=crossval.perf['TSS',])
  wmean.prob <- apply(crossval.preds,1,weighted.mean, w=as.numeric(crossval.perf['TSS',]))
  
  # Committee averaging of binary predictions: calculates the proportion of models that predict the species to be present.
  committee.av <- rowSums(sapply(our.models, FUN=function(x){ ifelse(crossval.preds[,x]>=crossval.perf['thresh',x],1,0)}))/length(our.models)
  
  # We can also calculate uncertainty measures, e.g. the standard deviation when
  # making ensembles of mean probabilities.
  sd.prob <- apply(crossval.preds,1,sd)
  
  # We can also put these into a function:
  make.ensemble <- function(preds, eval.metric, thresh){
    # "preds" is a data.frame containing predictions for different algorithms.
    # "eval.metric" is a vector with same length as number of columns in preds. It
    # provides the evaluation metric used for weighting probabilities.
    # "thresh" is a vector with same length as number of columns in preds. It provides
    # the algorithm-specific threshold for making binary predictions.
    data.frame(mean.prob = rowMeans(preds),
               median.prob = apply(preds,1,median),
               #wmean.prob = apply(preds,1,weighted.mean, w=eval.metric),
               wmean.prob = apply(preds,1,weighted.mean, w=as.numeric(eval.metric)),
               committee.av = rowSums(sapply(seq_len(ncol(preds)), FUN=function(x){
                 ifelse(preds[,x]>=thresh[x],1,0)}))/ncol(preds),
               sd.prob = apply(preds,1,sd))
  }
  
  # Make ensemble predictions:
  ensemble.preds <- make.ensemble(crossval.preds,
                                  crossval.perf['TSS',], crossval.perf['thresh',])
  # Evaluate ensemble predictions:
  ensemble.perf <- sapply(names(ensemble.preds)[1:4], FUN=function(x){
    calc.eval(ttetrix.env, 'occ', ensemble.preds[,x])})
  
  ##############################################################################################################
  #Visualising response surfaces
  
  #library(RColorBrewer)
  #library(lattice)
  #cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)
  # We prepare our grid of environmental predictors:
  #xyz <- expand.grid(
  #seq(min(ttetrix.env[,pred[1]]),max(ttetrix.env[,pred[1]]),length=50),
  #seq(min(ttetrix.env[,pred[2]]),max(ttetrix.env[,pred[2]]),length=50))
  #names(xyz) <- pred[1:2]
  
  xyz <- data.frame(1:50)
  for (j in 1:length(pred)) {
    xyz <- data.frame(xyz, rep(mean(ttetrix.env[,pred[j]]),length=50))
  }
  xyz <- xyz[-1]
  names(xyz) <- pred
  
  setwd("./input_data/unclipped_SDMs")
  tiff(paste(predators[i],"_response_curves.tiff",sep=""), width = 18, height = 22, units = 'cm', res = 100, compression="lzw") 
  par(mfrow=c(5,3), mar=c(4,2,1.5,2), oma=c(4,4,2,1), xpd=FALSE)
  
  for (j in 1:length(pred)){
    test <- xyz
    test[,j] <- seq(min(ttetrix.env[,pred[j]]),max(ttetrix.env[,pred[j]]),length=50)
    
    # We make predictions of all models and make ensembles:
    xyz.preds <- sapply(our.models, FUN=function(m){make.preds(eval(parse(text=m)), test)})
    
    xyz.ensemble <- make.ensemble(xyz.preds, crossval.perf['TSS',], crossval.perf['thresh',])
    
    plot(seq(min(ttetrix.env[,pred[j]]),max(ttetrix.env[,pred[j]]),length=50), xyz.ensemble$wmean.prob, ylim=c(0,1), xlab=pred[j], cex=0.01, las=1)
    polygon(c(seq(min(ttetrix.env[,pred[j]]),max(ttetrix.env[,pred[j]]),length=50),rev(seq(min(ttetrix.env[,pred[j]]),max(ttetrix.env[,pred[j]]),length=50))),
            c(xyz.ensemble$wmean.prob + xyz.ensemble$sd.prob,rev(xyz.ensemble$wmean.prob - xyz.ensemble$sd.prob)), col="lightgrey", border=NA)
    
    lines(seq(min(ttetrix.env[,pred[j]]),max(ttetrix.env[,pred[j]]),length=50), xyz.ensemble$wmean.prob, lwd=2.5)
  }
  
  dev.off()
  
  ###################################################################################
  #Mapping ensemble predictions
  
  #ext = c(129, 141, -39, -25)
  ext = c(-304362.2, 902637.8, -4297576, -2697576)
  r.mean.prob <- r.median.prob <- r.wmean.prob <- r.committee.av <-
    r.sd.prob <- raster(resolution=res(envdat),ext=extent(ext))
  
  crs(r.mean.prob) <- crs(r.median.prob) <- crs(r.wmean.prob) <- crs(r.committee.av) <-
    crs(r.sd.prob) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "
  
  envdat[is.na(envdat)] <- 0
  
  # We make predictions of all models:
  env.preds <- sapply(our.models, FUN=function(m){make.preds(eval(parse(text=m)), na.exclude(data.frame(values(crop(envdat[[pred]],ext)))) )})
  
  #ext = c(129, 141, -39, -25)
  ext = c(-304362.2, 902637.8, -4297576, -2697576)
  glm.occ <- gam.occ <- rf.occ <- brt.occ <-
    brt.occ <- raster(resolution=res(envdat),ext=extent(ext))
  
  crs(glm.occ) <- crs(gam.occ) <- crs(rf.occ) <- 
    crs(brt.occ) <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "
  
  values(glm.occ)[!is.na(values(crop(envdat[[1]],ext)))] <- env.preds[,'m.glm']
  values(gam.occ)[!is.na(values(crop(envdat[[1]],ext)))] <- env.preds[,'m.gam']
  values(rf.occ)[!is.na(values(crop(envdat[[1]],ext)))] <- env.preds[,'m.rf']
  values(brt.occ)[!is.na(values(crop(envdat[[1]],ext)))] <- env.preds[,'m.brt']
  
  #env.preds <- sapply(our.models, FUN=function(m){make.preds(eval(parse(text=m)), data.frame(values(envdat))  )})
  
  #test <- na.exclude(data.frame(values(crop(envdat[[pred]],ext))))
  #test <- data.frame(values(envdat))
  
  # We make ensembles:
  env.ensemble <- make.ensemble(env.preds, crossval.perf['TSS',], crossval.perf['thresh',])
  
  # We assign the ensemble predictions to the raster:
  values(r.mean.prob)[!is.na(values(crop(envdat[[1]],ext)))] <- env.ensemble[,'mean.prob']
  values(r.median.prob)[!is.na(values(crop(envdat[[1]],ext)))] <- env.ensemble[,'median.prob']
  values(r.wmean.prob)[!is.na(values(crop(envdat[[1]],ext)))] <- env.ensemble[,'wmean.prob']
  values(r.committee.av)[!is.na(values(crop(envdat[[1]],ext)))] <- env.ensemble[,'committee.av']
  values(r.sd.prob)[!is.na(values(crop(envdat[[1]],ext)))] <- env.ensemble[,'sd.prob']
  
  r.mean.prob[is.na(envdat[[1]])] <- NA
  r.median.prob[is.na(envdat[[1]])] <- NA
  r.wmean.prob[is.na(envdat[[1]])] <- NA
  r.committee.av[is.na(envdat[[1]])] <- NA
  r.sd.prob[is.na(envdat[[1]])] <- NA
  
  r.mean.prob[is.na(elevation)] <- NA
  r.median.prob[is.na(elevation)] <- NA
  r.wmean.prob[is.na(elevation)] <- NA
  r.committee.av[is.na(elevation)] <- NA
  r.sd.prob[is.na(elevation)] <- NA
  
  glm.occ[is.na(elevation)] <- NA
  gam.occ[is.na(elevation)] <- NA
  rf.occ[is.na(elevation)] <- NA
  brt.occ[is.na(elevation)] <- NA
  
  plot(stack(r.mean.prob, r.median.prob, r.wmean.prob, r.committee.av, r.sd.prob), main=names(env.ensemble))
  
  #Idenitfy coordinates with presences
  Sites <- SpatialPoints(cbind(ttetrix.env[,1],ttetrix.env[,2]),proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  Sites <- spTransform(Sites, CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  Presence <- subset(ttetrix.env, ttetrix.env$occ==1)
  Presence <- SpatialPoints(cbind(Presence[,1],Presence[,2]),proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  Presence <- spTransform(Presence, CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  ###############################################################################
  #Plot ensemble model
  
  setwd("./input_data/unclipped_SDMs")
  
  par(mfrow=c(1,1), oma = c(3, 4, 2, 0), mar=c(0.5,2,1.5,1)) 
  tiff(paste(predators[i],".tiff",sep=""), width = 18, height = 22, units = 'cm', res = 100, compression="lzw")
  
  cuts=seq(0,1,0.1) #set break
  
  plot(r.wmean.prob, xaxt="n", yaxt="n", main=predators[i], breaks = cuts, col=colorRampPalette(c("lightgrey", "blue"), space="rgb")(16), axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5))#, legend.shrink=1.1) #, axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5), bty="n", box=FALSE)
  plot(aus, add=TRUE)
  #Plot incidentals
  
  #Plot presence-absence from 2 ha plots
  points(Sites, pch=19, col="khaki", cex=0.9)
  points(Presence, pch=19, col="black", cex=0.9)
  
  dev.off()
  
  #################################################################################
  #Plot 4 different model algorithms
  ext <- c(-304362.2, 902637.8, -4297576, -2697576)
  glm.occ <- crop(glm.occ, ext)
  gam.occ <- crop(gam.occ, ext) 
  rf.occ <- crop(rf.occ, ext) 
  brt.occ <- crop(brt.occ, ext)
  r.wmean.prob2 <- crop(r.wmean.prob, ext) 
  
  tiff(paste(predators[i],"allmodels.tiff",sep=""), width = 18, height = 22, units = 'cm', res = 100, compression="lzw") 
  par(mfrow=c(3,2), oma = c(3, 4, 2, 0), mar=c(0.5,2,1.5,1)) 
  plot(glm.occ, xaxt="n", yaxt="n", main="GLM", breaks = cuts, col=colorRampPalette(c("lightgrey", "blue"), space="rgb")(16), axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5))#, legend.shrink=1.1) #, axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5), bty="n", box=FALSE)
  plot(aus, add=TRUE)
  plot(gam.occ, xaxt="n", yaxt="n", main="GAM", breaks = cuts, col=colorRampPalette(c("lightgrey", "blue"), space="rgb")(16), axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5))#, legend.shrink=1.1) #, axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5), bty="n", box=FALSE)
  plot(aus, add=TRUE)
  plot(rf.occ, xaxt="n", yaxt="n", main="RF", breaks = cuts, col=colorRampPalette(c("lightgrey", "blue"), space="rgb")(16), axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5))#, legend.shrink=1.1) #, axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5), bty="n", box=FALSE)
  plot(aus, add=TRUE)
  plot(brt.occ, xaxt="n", yaxt="n", main="BRT", breaks = cuts, col=colorRampPalette(c("lightgrey", "blue"), space="rgb")(16), axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5))#, legend.shrink=1.1) #, axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5), bty="n", box=FALSE)
  plot(aus, add=TRUE)
  plot(r.wmean.prob2, xaxt="n", yaxt="n", main="Ensemble model", breaks = cuts, col=colorRampPalette(c("lightgrey", "blue"), space="rgb")(16), axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5))#, legend.shrink=1.1) #, axis.args=list(at=seq(0,1,0.1), tick=TRUE, labels=seq(0,1,0.1), lwd=0.5), bty="n", box=FALSE)
  plot(aus, add=TRUE)
  dev.off()
  
  #Save ensemble model
  writeRaster(r.wmean.prob, file = paste(predators[i],".tif", sep = ""), overwrite=TRUE)
  #writeRaster(r.wmean.prob, file = paste(predators[i],".asc", sep = ""), format ='ascii',overwrite=TRUE)
  save(Sites, Presence, file = paste(predators[i],"_sites.RData", sep = ""))
  
  train.performance[[i]] <- train.perf
  crossval.performance[[i]] <- crossval.perf
  ensemble.performance[[i]] <- ensemble.perf
  
}

#write.csv(unlist(train.performance), file = "C:/Darren/12_Postdoc/7_Desert_monitoring_project/Maps/SDM_paper/train.performance.csv")
setwd("./input_data/unclipped_SDMs")
train <- as.matrix(do.call("rbind", lapply(train.performance, as.data.frame)))
crossval <- as.matrix(do.call("rbind", lapply(crossval.performance, as.data.frame)))
ensemble <- as.matrix(do.call("rbind", lapply(ensemble.performance, as.data.frame)))
write.csv(train,file="train.performance_predators.csv",row.names=F)
write.csv(crossval,file="crossval.performance_predators.csv",row.names=F)
write.csv(ensemble,file="ensemble.performance_predators.csv",row.names=F)


