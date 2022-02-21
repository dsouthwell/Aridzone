#######################################################################################
# ANALYSIS OF THE OCCURRENCE AND DETECTABILITY OF SPECIES FROM TRACK BASED SURVEY PLOTS
#######################################################################################

# set pathway
rm(list=ls())
getwd()
setwd(".......")

library(unmarked)
library(MuMIn)
library(parallel) #Runs the dredge function across multiple cores to speed it up
library(snow) #Runs the dredge function across multiple cores to speed it up
library(stats4) #Runs the dredge function across multiple cores to speed it up

library(AICcmodavg)
library(data.table)
library(ggplot2)

# get data
detection.data <- read.csv("AWDetectionData2019.csv", header=T, sep=",") #this works3.csv file has the substrate variable updated
head(detection.data)
dim(detection.data)
summary(detection.data)


#Scale/standardise the continuous vars to avoid non-convergence
dd.scaled <- as.data.frame(scale(detection.data[,c(9:13,28)], center=TRUE, scale=TRUE)) #just scales the continuous variables, columns indicated as new dataframe
dd.scaled.full <- cbind(detection.data[,c(1:8,14:27,29:52)], dd.scaled) # combines the new dataframe with the unchanged columns = 45 variables

#Make categorical vars factors again
dd.scaled.full$rain_wind <- as.factor(dd.scaled.full$rain_wind)
dd.scaled.full$spinifex <- as.factor(dd.scaled.full$spinifex_dom) #this adds another variable of spinifex
dd.scaled.full$visit<- as.factor(dd.scaled.full$visit)
dd.scaled.full$plot_observer<- as.factor(dd.scaled.full$plot_observer)
dd.scaled.full$year<- as.factor(dd.scaled.full$year)

#make vectors of site and visit 
idsite <- as.vector(dd.scaled.full[,"idsite"])
vis2 <- as.vector(dd.scaled.full[,"vis2"])
year<- as.vector(dd.scaled.full[,"year"])

# make the observer into a couple of values
levels(dd.scaled.full$plot_observer)[levels(dd.scaled.full$plot_observer)=="rc,nw"] <- "rc"
levels(dd.scaled.full$plot_observer)[levels(dd.scaled.full$plot_observer)=="rcnw"] <- "rc"
levels(dd.scaled.full$plot_observer)[levels(dd.scaled.full$plot_observer)=="rs,cm"] <- "rs"
levels(dd.scaled.full$plot_observer)[levels(dd.scaled.full$plot_observer)=="rs,cm,cq,cy"] <- "rs"
levels(dd.scaled.full$plot_observer)[levels(dd.scaled.full$plot_observer)=="rs,cm,cy"] <- "rs"
levels(dd.scaled.full$plot_observer)[levels(dd.scaled.full$plot_observer)=="rs,nathan,con,plus"] <- "rs"

indexYear1<- c(1:334)  #not including 2007
indexYear2<- c(335:715)
indexYear3<- c(716:979)

indexAllYears <- c(1-979)

Sites <- 150
Repeats <-3   # secondary sample periods
Years <- 3   # primary sample periods - not including 2007

######################
# Detection covariates
#####################

# Create list matrix for each observaton covariate, it loops through and creates seperate matrices for each of the listed variables
# this will then be fed into the dataframe


det_data <- with(dd.scaled.full, list(plot_observer, shadow_intensity, continuity, start_time, rain_wind, visit))
obs.covs <- vector('list', length = length(det_data))    # create empty vector
names(obs.covs) <- c('observer.cov', 'shadow.cov','continuity.cov','start.cov','rainwind.cov','visit.cov' )

for (i in 1:length(det_data))  {

  det.mat1 <- matrix(NA, nrow = Sites, ncol = Repeats) 
  det.mat2 <- matrix(NA, nrow = Sites, ncol = Repeats)
  det.mat3 <- matrix(NA, nrow = Sites, ncol = Repeats)  
  
  for(j in indexYear1){
    row.ind <- idsite[j] #jth occurrence
    col.ind<- vis2[j]
    det.mat1[row.ind, col.ind]<- det_data[[i]][j] # vector i, jth occurrence
  }
  
  for(j in indexYear2){
    row.ind <- idsite[j]
    col.ind<- vis2[j]
    det.mat2[row.ind, col.ind]<- det_data[[i]][j]
  }
  
  for(j in indexYear3){
    row.ind <- idsite[j]
    col.ind<- vis2[j]
    det.mat3[row.ind, col.ind]<- det_data[[i]][j]
  }
  
  obs.covs[[i]] <- cbind(det.mat1, det.mat2, det.mat3) #combines the matrices
  
}

str(obs.covs)

###############################################
## Year covariates (the year of the record)
###############################################

yearly_data <- with(dd.scaled.full, list(rainfall_preceeding, year))
yearly.site.covs <- vector('list', length = length(yearly_data))    # create empty vector
names(yearly.site.covs) <- c('rain','year' )


for (i in 1:length(yearly_data))  {
  
  yearly.mat <- matrix(NA, nrow = Sites, ncol = 3) 

  for(j in 1:334){
    row.ind <- idsite[j] #jth occurrence
    yearly.mat[row.ind, 1]<- yearly_data[[i]][j] # vector i, jth occurrence
  }
  
  for(j in indexYear2){
    row.ind <- idsite[j]
    yearly.mat[row.ind, 2]<- yearly_data[[i]][j]
  }
  
  for(j in indexYear3){
    row.ind <- idsite[j]
    yearly.mat[row.ind, 3]<- yearly_data[[i]][j]
  }
  
  yearly.site.covs[[i]] <- yearly.mat
  
}  

sum(is.na(dd.scaled.full$rainfall_preceeding))
sum(is.na(dd.scaled.full$year))
#No NAs in yearly data


#there are some sites that didn't have surveys in some years, to run colext we need to impute some year values below to remove NAs
m<- mean(yearly.site.covs$rain[,1],na.rm=TRUE)
# year three rainfall all the same value
yearly.site.covs$rain[,3]<- -1.186267
yearly.site.covs$rain[91,2]<-0.8588763
yearly.site.covs$rain[146:150,1]<--0.345518  # mean value for the year

yearly.site.covs$year[,1]<-"1"
yearly.site.covs$year[,2]<-"2"
yearly.site.covs$year[,3]<-"3"

#####################
#### Site covariates 
#####################

test <- dd.scaled.full[order(idsite),] #order the dataset by site

landform.vec<-rep(NA,150) #create an empty vector for 'landform' var
landform <- as.vector(test[,'map_landform_cat4']) #vector of all 'rock_ha' ordered by site ID no
id.vec <- as.vector(test[,'idsite']) #vector of site IDs (in order)
for(i in indexAllYears){
  landform.vec[id.vec[i]] <- landform[i] #fill rock.vec with 'rocks' value for each individual 'id.vec'
}
landform.vec[1:10]

#map_veg_name
veg.vec<-rep(NA,150) 
veg <- as.vector(test[,'map_veg_name']) 
id.vec <- as.vector(test[,'idsite']) 
for(i in indexAllYears){
  veg.vec[id.vec[i]] <- veg[i] 
}

#spinifex (this had been already changed from spinifex_dom)
spini.vec<-rep(NA,150) 
spini <- as.vector(test[,'spinifex']) 
id.vec <- as.vector(test[,'idsite']) 
for(i in indexAllYears){
  spini.vec[id.vec[i]] <- spini[i] 
}


#obs_shrub_dom
shrubs.vec<-rep(NA,150) 
shrubs <- as.vector(test[,'obs_shrub_dom']) 
id.vec <- as.vector(test[,'idsite']) 
for(i in indexAllYears){
  shrubs.vec[id.vec[i]] <- shrubs[i] 
}

plotname.vec<-rep(NA,150) 
plotname <- as.vector(dd.scaled.full[,'plot_name']) 
id.vec <- as.vector(dd.scaled.full[,'idsite']) 
for(i in indexAllYears){      # for first year only
  plotname.vec[id.vec[i]] <- plotname[i] 
}

site.covs <- data.frame(landform=landform.vec, veg=veg.vec, spini=spini.vec, shrubs=shrubs.vec, plot=plotname.vec)

#####################################################################
## Section 2: Detection matrices for species
####################################################################

# add in red_kangaroo ANY
dd.scaled.full$redkangaroo_recmid<-(dd.scaled.full$redkangaroo_mid+dd.scaled.full$redkangaroo_rec)
dd.scaled.full$redkangaroo_recmid<-replace(dd.scaled.full$redkangaroo_recmid, dd.scaled.full$redkangaroo_recmid==2, 1)


# add in grey_kangaroo ANY
dd.scaled.full$greykangaroo_recmid<-(dd.scaled.full$greykangaroo_mid+dd.scaled.full$greykangaroo_rec)
dd.scaled.full$greykangaroo_recmid<-replace(dd.scaled.full$greykangaroo_recmid, dd.scaled.full$greykangaroo_recmid==2, 1)

#####################
## CREATE DATAFRAMES
#####################

#setup data to predict occupancy
newData.occ <- data.frame(landform=factor("Colluvium", levels = c("Colluvium", "Residual","Sand dunes", "Sand plain")),
                          veg= factor("Maralinga", levels = c("Maralinga","Nurrairi", "Purdu", "Victoria Desert", "Yellabina")),
                          spini = factor("0", levels = c("0","1")),
                          shrubs = factor("mallee", levels = c("mallee", "mix", "mulga")),
                          rainwind.cov = factor("1", levels = c("1","2")),
                          continuity.cov = 0,
                          shadow.cov=0,
                          start.cov=0,
                          #visit.cov=factor("1", levels= c("1","2","3")),
                          observer.cov= factor("1", levels = c("1","2","3","4")))

#Set up data to predict detectability
newData.det <- data.frame(landform=factor("Colluvium", levels = c("Colluvium", "Residual","Sand dunes", "Sand plain")),
                          veg= factor("Maralinga", levels = c("Maralinga","Nurrairi", "Purdu", "Victoria Desert", "Yellabina")),
                          spini = factor("0", levels = c("0","1")),
                          shrubs = factor("mallee", levels = c("mallee", "mix", "mulga")),
                          rainwind.cov = 0,
                          continuity.cov = 0,
                          shadow.cov=0,
                          start.cov=0,
                          #visit.cov=0,
                          observer.cov= factor("1", levels = c("1","2","3","4"))) 

newData.det <- data.frame(landform=factor("Colluvium", levels = c("Colluvium", "Residual","Sand dunes", "Sand plain")),
                          veg= factor("Maralinga", levels = c("Maralinga","Nurrairi", "Purdu", "Victoria Desert", "Yellabina")),
                          spini = factor("0", levels = c("0","1")),
                          shrubs = factor("mallee", levels = c("mallee", "mix", "mulga")),
                          rainwind.cov = 0,
                          continuity.cov = c(1,2,3),
                          shadow.cov=0,
                          start.cov=0,
                          #visit.cov=0,
                          observer.cov= factor("1", levels = c("1","2","3","4")))  

# Load function needed for the bootstrapping
fitstats <- function(model.name, method = "nonparboot") {
  observed <- getY(model.name@data)
  expected <- fitted(model.name)
  resids <- residuals(model.name, method = "nonparboot")
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

save.parboot <- function(object){
  t.star <- object@t.star
  t0 <- object@t0
  nsim <- nrow(t.star)
  biasMat <- pMat <- matrix(NA, nsim, length(t0))
  for(i in 1:nsim) {
    biasMat[i,] <- t0 - t.star[i,]
    pMat[i,] <- abs(t.star[i,] - 1) > abs(t0 - 1)
  }
  mean.tB <- mean(t.star)
  tB.se <- apply(t.star, 2, sd)
  bias <- colMeans(biasMat)
  bias.se <- apply(biasMat, 2, sd)
  p.val <- colSums(pMat) / (1 + nsim)
  stats <- data.frame("Pr(t_B > t0)" = p.val,	
                      check.names = FALSE)     
  return(stats)
}

modAvgCoefs <- function(ms) {
  
  theta.hat <- coef(ms)
  vars <- SE(ms)^2
  wts <- ms@Full$AICwt
  
  theta.hat[is.na(theta.hat)] <- vars[is.na(vars)] <- 0
  theta.hatbar <- colSums(theta.hat*wts)
  
  theta.hatbar.rk <- matrix(theta.hatbar, nrow(theta.hat), ncol(theta.hat), 
                            byrow=TRUE)
  theta.hatbar.se <- colSums(wts*sqrt(vars + (theta.hat-theta.hatbar.rk)^2))
  return(data.frame(Ests=theta.hatbar, SE=theta.hatbar.se))
}

#Set up working directory where the results should be saved
setwd("......") 

  
######################################
# Running the models

#Err <- rep(0, length(species))

species<- as.vector(c('hopm_rec','dog_rec','fox_rec','cat_rec','camel_rec','rabbit_rec',
 'rabbit_recmid','dog_recmid', 'fox_recmid', 'cat_recmid', 'camel_recmid',
'greykangaroo_rec', 'greykangaroo_recmid','malleefowl_recmid')) 

Output <- c("Species","No. occ sites","Naive occ","Naive det", "Mean occ","SE Occ","Mean det","SE det","Model averaging","No. models averaged")

Species.coeffs<-list()

OUT <- list()

k <- 1
i<- "fox_rec"
#ratio<-10  #what does this do??

par(mfrow=c(4,4), mar=c(5,2,2,2), oma=c(4,4,2,2), xpd=FALSE)

for (i in species) {
  
  output <- c()
  mdp <- dd.scaled.full[[i]] #Extract data for species i
  
# set up the detection matrix  
 species.y1 <- matrix(NA, nrow = Sites, ncol = Repeats) 
  species.y2 <- matrix(NA, nrow = Sites, ncol = Repeats)
  species.y3 <- matrix(NA, nrow = Sites, ncol = Repeats)  
  
    for(j in indexYear1){
    row.ind <- idsite[j] #jth occurrence
    col.ind<- vis2[j]
    species.y1[row.ind, col.ind]<- mdp[j] #ith occurrence
    }
  
    for(j in indexYear2){
    row.ind <- idsite[j]
    col.ind<- vis2[j]
    species.y2[row.ind, col.ind]<- mdp[j]
    }
  
    for(j in indexYear3){
     row.ind <- idsite[j]
    col.ind<- vis2[j]
    species.y3[row.ind, col.ind]<- mdp[j]
    }
  
  y <- cbind(species.y1, species.y2, species.y3) #combines the matrices
  
  #calculate number of occupied sites
  Obs <- ncol(y)
  sites <- nrow(y)
  occ.sites <- y
  occ.sites[is.na(occ.sites)] <- 0
  occ.sites<-apply(occ.sites,1,function(x)sum(x))
  occ.sites[occ.sites>0] <- 1
  n.sites <- sum(occ.sites)
  #n.cov <- floor(n.sites/ratio)
  #if (n.cov == 0) {n.cov <- 1}
  #output[6] <- n.cov
  output[2]<- n.sites
  
  # calculate naiive occupancy
  output[3]<-n.sites/150
  
  # calculate naiive detectability
  #output [4] <- 
  
  

    #Define special unmarked dataframe
    umf_out <- unmarkedMultFrame(y=y, siteCovs=site.covs,
                                 yearlySiteCovs=yearly.site.covs, obsCovs=obs.covs,                
                                 numPrimary=3) # number of seasons/years organize data

    print(i) #Print the name of the species to track the loop
    output[1] <- paste(i)
    
    #Run models   
    ##adding "-1" to the factors so that the first value is not the intercept
    
    null <- colext(~1, ~ 1, ~1, ~ 1, data=umf_out) 
    
    m1 <- nonparboot(null, B = 100)
    occ <- cbind(smoothed=smoothed(m1)[2,], SE=m1@smoothed.mean.bsse[2,])
    
    years <- c("2009", "2011", "2013")
    
    nd <- data.frame(year=c("2009", "2011", "2013"), date=0)
    E.det <- predict(null, type='det',newdata=nd) #Detection estimate for year 1 and 2
    print(E.det$Predicted[1])
    
    plot(years,occ[,1], ylim=c(0,1), pch=19,main=paste(species[k]), xlab="Year", ylab="Occupancy probability",cex=0.8)
    arrows(as.numeric(years), occ[,1]-occ[,2]*1.96, as.numeric(years), occ[,1]+occ[,2]*1.96, code=3, angle=90, length=0.03, col=4, lwd=2.5)
    
    OUT[[k]] <- occ
    k <- k+1
    
}
    
   model.list<-list(null)
   
   names(model.list)<- c("null")
   
   conv=NULL
   se.nan=NULL
   
   for (p in 1:length(model.list))
   {
     conv[p]= model.list[[p]]@opt$convergence==0 # these are just the models that converged
     se.nan[p]=any(is.nan(SE(model.list[[p]])))
   }
   
   print(se.nan)
   
   model.list.2 = model.list[conv]  
   Cand.mods=model.list.2[se.nan != "TRUE"]   # models not including those that have TRUE NaNs
   
  #list of model groups that don't need to be inspected
   names(Cand.mods)
   
   ms<- modSel(fitList(fits = Cand.mods))
   ms
   
# export to excel docs with models in order using the coeffs
   toExport <- as(ms, "data.frame")
   write.csv(toExport, file = paste(paste(i),".csv", sep=""))
   
   # extract models that have delta <2
   n.delta<- sum(ms@Full$delta <= 2)
   cand.set<- list(ms@Full$model[1:n.delta])
   cand.set<-unlist(cand.set)
   output[10]<- n.delta
   
   ##### loops depending on number of models in <2 delta AIC
   
   if(n.delta == 1){  # if one top model then
     
     output[9]<- "no"
     
     best.model<- null
     
     best.model <- get(cand.set)
     
     occu_model<-predict(best.model, type='psi',newData.occ) # what is type here?
     #print(summary(occu_model, newData))
     mean.occ <- summary(occu_model, umf)
     test <- gsub(":","", mean.occ[4,1], perl=TRUE)
     test2 <- gsub("[A-z]","", test)
     output[5] <- paste(as.numeric(test2))
     test <- gsub(":","", mean.occ[4,2], perl=TRUE) #how do I deal when its a -e number?
     test2 <- gsub("[A-z]","", test)
     output[6] <- paste(as.numeric(test2))
     
     det_model<-predict(best.model, type='det',newData.det)   # with this you'll get predictions for detectability for every site
     #print(summary(det_model, newData))
     mean.det <- summary(det_model, newData)
     test <- gsub(":","", mean.det[4,1], perl=TRUE)
     test2 <- gsub("[A-z]","", test)
     output[7] <- paste(as.numeric(test2))
     test <- gsub(":","", mean.det[4,2], perl=TRUE)
     test2 <- gsub("[A-z]","", test)
     output[8] <- as.numeric(test2)
     
     
     #Estimate occupancy in year 1 2 and 3
     m1 <- nonparboot(model1, B = 100)
     occ <- cbind(smoothed=smoothed(m1)[2,], SE=m1@smoothed.mean.bsse[2,])
     
     #estimate the coeffs from the top model
     coeff.best<-coef(best.model)
     coeff.best<-data.frame(as.list(coeff.best))
     coeff.best$Species<-paste(i)
     coeff.best$Coeff<-"Ests"
     
     #second the SEs
     SE.best<-SE(best.model)
     SE.best<-data.frame(as.list(SE.best))
     SE.best$Species<-paste(i)
     SE.best$Coeff<-"SE"
     
     #reformat and rearrange to long format for graphing later
     melt.Est<-melt(coeff.best)
     names(melt.Est)[4]<-"Ests"
     
     melt.SE<-melt(SE.best)
     names(melt.SE)[4]<-"SE"
     
     Est.coeffs<-merge(melt.SE,melt.Est, by=c("Species","variable"))
     
     Est.coeffs$Coeff.x<-NULL
     Est.coeffs$Coeff.y<-NULL
     names(Est.coeffs)[2]<-"rn"
     Est.coeffs$Averaged<-"no"
     
     Species.coeffs[[i]]<- Est.coeffs  #add to the list
    
      }
   

   if (n.delta > 1) { # if more than 1 model in the delta <2 group
     
      output[9]<- "yes"
      
      # get the cand.set of models to be listed together for model averaging
      mods.to.av<-list()
    
      
       for(j in (1:n.delta)){
         mods.to.av[[j]]<- get(cand.set[j])
        
       }
      
      names(mods.to.av)<-cand.set
        
       modav_spp<- modavgPred(mods.to.av, modnames = NULL,
                              newdata=newData.occ, second.ord = TRUE, nobs = 979, uncond.se = "revised",
                              conf.level = 0.95, type = "response", c.hat = 1,
                              parm.type = "psi")  
       
       output[5]<-paste(as.numeric(modav_spp$matrix.output[1])) #occ estimate
       output[6]<-paste(as.numeric(modav_spp$matrix.output[2])) # SE 

       
       #capture the mean det plus se for each type of sign where there are multiple models in <2AIC
       modav_spp<- modavgPred(mods.to.av, modnames = NULL,
                                  newdata=newData.det, second.ord = TRUE, nobs = 979, uncond.se = "revised",
                                  conf.level = 0.95, type = "response", c.hat = 1,
                                  parm.type = "detect")  # type says if occ or det side
       output[7]<-paste(as.numeric(modav_spp$matrix.output[1]))
       output[8]<-paste(as.numeric(modav_spp$matrix.output[2]))
  
       # average the coefficients where there is more than one top model
          ms.s<- modSel(fitList(fits = mods.to.av))
          
          Est.coeffs<-modAvgCoefs(ms.s)      #apply the function to average the coeffs in the top models
       
          #print the coeffs to a large table
          setDT(Est.coeffs,keep.rownames = TRUE)[]
          
          # add some extra info
          Est.coeffs$Species<- paste(i)
          Est.coeffs$Averaged <- "yes"
          Species.coeffs[[i]]<- Est.coeffs  # add the coeffs to the list
         
          }
   


  
   #Save results to file
  #save(best.model, file = paste(i,"_Best_model",".RData", sep = ""))
  save(ms, file = paste(i,"_Model_list",".RData", sep = ""))
    #save(pb, file = paste(i,"_GOF_test",".RData", sep = "")) 
  
   #Output has the mean estimates of occ and det for the sign type 
  Output <- rbind(Output, output)

  big_data = do.call(rbind, Species.coeffs)


write.csv(Output, file = "Allsign_occ_det.csv")
write.csv(big_data, file="Model_coefficients.csv")

#####################################

# code for the occ and det predictions
Output.graph<-as.data.frame(Output)
colnames(Output.graph)<- c("Sign","Occ.sites","Naive.Occ", "Naive.Det", "Occ", "OccSE","Det","DetSE", "Mod.av", "n.models")
Output.graph=Output.graph[-1,] #remove first row with names
Output.graph$Occ<- as.numeric(as.character(Output.graph$Occ)) # R automatically makes it factor levels, need to tell it to use the chr as numbers as.character
Output.graph$OccSE<- as.numeric(as.character(Output.graph$OccSE))
Output.graph$Det<- as.numeric(as.character(Output.graph$Det))
Output.graph$DetSE<- as.numeric(as.character(Output.graph$DetSE))
Output.graph$ID<-(1:length(species))

#make extra columns for the Occ and Det CIs
Output.graph$CIs <- (Output.graph$OccSE * 1.96)
Output.graph$DetCIs <- (Output.graph$DetSE * 1.96)


## Rearrange data to plot
Occ_det_long<- reshape(Output.graph,
                       varying=list(c("Occ","Det"), c("CIs","DetCIs")),
                       timevar = "Type", 
                       times = c("Occupancy","Detectability"),
                       v.names= c("Value","CI"),
                       direction ="long",
                       idvar="Species")

Occ_det_long$Type<-as.factor(as.character(Occ_det_long$Type))

