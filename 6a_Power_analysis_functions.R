#############################################################################################################################
#START BY RUNNING CHECKS ON INPUT RASTER LAYERS
#############################################################################################################################
#' Checks input raster layers
#'
#' This function checks the occupancy and detectability raster layers before running the power analysis
#' @param occ  The occ raster stack
#' @keywords 
#' @export
#' @examples
check.inputs <- function(occ){
  
  #Check dimensions of raster layers and print warning if they don't match
  if (ncol(occ) != ncol(det.method1) & ncol(det.method2) & ncol(det.method3) & ncol(parks) & ncol(remote)) {
    cat('\n',"Warning: Raster layers have unequal number of columns")
  }
  
  if (model.fire == TRUE){
    #if ((nrow(fire) != nrow(veg)) | (ncol(fire) != ncol(veg))){
    #cat('\n',"Warning: Raster layers have unequal number of columns")}
  }
  
  #Check dimensions of raster layers and print warning if they don't match
  if (nrow(occ) != nrow(det.method1) & nrow(det.method2) & nrow(det.method3) & nrow(parks) & nrow(remote)) {
    cat('\n',"Warning: Raster layers have unequal number of columns")
  }
  
  #Check number of species match
  if (nrow(species.list) != nlayers(occ) & nrow(species.cov) & nlayers(det.method1) & nlayers(det.method2) & nlayers(det.method3)) {
    cat('\n',"Warning: Number of species in species/list and/or species.cov not equal to number of raster layers in occ stack")
  }
  
  if (Tmax != max(s.years)) {
    cat('\n',"Warning: Monitoring must occur in final year (last year of monitoring must be equal to Tmax")}
  
  if (park.level == TRUE & length(n.park) == 1) {
    cat('\n',"Warning: Park level analysis selected but only one park identified in raster layer")
  }
  
  if (n.sites == 0) {
    cat('\n',"Warning: No monitoring sites identified")
  }
  
  if (prespecify.sites == TRUE & all.prespecify.sites == FALSE & new.site.selection == "maxocc") {
    cat('\n',"Warning: If new.site.selection == maxocc, set prespecify.sites == FALSE")
  }
  
  if (prespecify.sites == TRUE & all.prespecify.sites == FALSE & new.site.selection == "stratified") {
    cat('\n',"Warning: If new.site.selection == stratified, set prespecify.sites == FALSE")
  }
  
  if (prespecify.sites == TRUE & all.prespecify.sites == FALSE & new.site.selection == "manual") {
    cat('\n',"Warning: If new.site.selection == manual, set prespecify.sites == FALSE")
  }
}


#' Site selection function
#'
#' This function selects sites in the landscape for monitoring. Sites can be pre-selected by loading in file of the XY coordinates, randomly selected, randomly stratified across environmental layers, positioned on cells
#' with the highest expected species richness, or manually selected by clicking on the raster map. Sites can also be divided between remote and non-remote areas 
#' @param sites  Contains the XY coordinates of sites loaded into the package for analysis
#' @param n.sites  The number of sites to be monitored. This is equal to the number of rows in 'sites' if all.prespecified.sites == TRUE. Otherwise, it is specified by the user      
#' @param R The ratio of remote to non-remote sites to be monitored. Setting R=1 means all sites will be in remote areas. Setting R=0 means no sites will be in remote areas
#' @keywords 
#' @export
#' @examples
#' prespecify.sites <- TRUE
#' all.prespecify.sites <- TRUE
#' randomise.sites <- FALSE
#' new.site.selection <- 'random'
#' R <- 0.5
#' if (all.prespecify.sites == TRUE) {n.sites <- nrow(sites)} else {n.sites <- 50}
#' plotter <- TRUE
#' xy.sites <- select.sites(sites, n.sites, R)

select.sites <-  function(sites, n.sites, R, all.prespecify.sites, prespecify.sites, new.site.selection, plotter) { #I should really change the name of function and add new.site.selection + plotter as arguments in main function. Note, doesn't apply when sites are randomly selected - need to add, plots maxocc if selected
  
  #Identify and then separate remote and non-remote sites
  rem.ras <- acc.ras <- remote
  rem.ras[rem.ras==0] <- NA
  acc.ras[acc.ras==1] <- NA
  Rem <- remote[cellFromXY(remote, cbind(sites$POINT_X,sites$POINT_Y))] #Extract remote ID from remote layer
  n.rem <- floor(n.sites*R) #Calculate number of remote sites
  n.acc <- n.sites - n.rem #Calculate number of non-remote sites
  XY.rem <- cbind(sites$POINT_X,sites$POINT_Y)[which(Rem==1),] #Separate remote and non-remote sites
  XY.acc <- cbind(sites$POINT_X,sites$POINT_Y)[which(Rem==0),]
  
  #Select all pre-selected sites
  if  (all.prespecify.sites == TRUE) {
    XY.sites <- cbind(sites$POINT_X,sites$POINT_Y)
  }
  
  #Select a random subset of pre-specified sites given the number of remote and non-remote sites is less then what is available
  if (prespecify.sites == TRUE & all.prespecify.sites == FALSE & n.rem <= nrow(XY.rem)  &  n.acc <= nrow(XY.acc)) {
    XY.rem <- XY.rem[sample(nrow(XY.rem), n.rem, replace=FALSE), ] #Randomly select remote and non-remote sites
    XY.acc <- XY.acc[sample(nrow(XY.acc), n.acc, replace=FALSE), ]
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  #The number of remote sites to survey is less than what is available. But the number of non-remote sites for survey is greater than what is available
  if (prespecify.sites == TRUE & all.prespecify.sites == FALSE & n.rem <= nrow(XY.rem)  &  n.acc > nrow(XY.acc)) {
    XY.rem <- XY.rem[sample(nrow(XY.rem), n.rem, replace=FALSE), ] #Randomly select remote and non-remote sites
    XY.acc <- rbind(XY.acc, sampleRandom(acc.ras,(n.acc-nrow(XY.acc)),na.rm=TRUE, xy=TRUE)[,1:2])
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  #The number of remote sites to survey is greater than what is available. But the number of 4WD sites for survey is less than what is available
  if (prespecify.sites == TRUE & all.prespecify.sites == FALSE & n.rem > nrow(XY.rem)  &  n.acc <= nrow(XY.acc)) {
    XY.acc <- XY.acc[sample(nrow(XY.acc), n.acc, replace=FALSE), ] #Randomly select remote and non-remote sites
    XY.rem <- rbind(XY.rem, sampleRandom(rem.ras,(n.rem-nrow(XY.rem)),na.rm=TRUE, xy=TRUE)[,1:2])
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  #The number of remote and 4WD sites to be surveyed is greater than the number of pre-existing sites
  if (prespecify.sites == TRUE & all.prespecify.sites == FALSE & n.rem > nrow(XY.rem)  &  n.acc > nrow(XY.acc)) {
    XY.acc <- rbind(XY.acc, sampleRandom(acc.ras, (n.acc-nrow(XY.acc)), na.rm=TRUE, xy=TRUE)[,1:2]) #Randomly select remote and non-remote sites
    XY.rem <- rbind(XY.rem, sampleRandom(rem.ras, (n.rem-nrow(XY.rem)), na.rm=TRUE, xy=TRUE)[,1:2])
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  if (prespecify.sites == FALSE & new.site.selection == "random") {
    XY.acc <- XY.rem <- NULL
    if (n.acc > 0) {XY.acc <- sampleRandom(acc.ras, n.acc, na.rm=TRUE, xy=TRUE)[,1:2]} #Randomly select remote and non-remote sites
    if (n.rem > 0){XY.rem <- sampleRandom(rem.ras, n.rem, na.rm=TRUE, xy=TRUE)[,1:2]}
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  if (prespecify.sites == FALSE & new.site.selection == "stratified") {
    classes <- unique(values(veg))
    classes <- na.omit(classes)
    strata <- length(classes)
    XY.sites <- sampleStratified(veg, size=floor(n.sites/strata), na.rm=TRUE, xy=TRUE)[,2:3] 
    if (nrow(XY.sites) < n.sites) {
      extra <- sampleRandom(veg, (n.sites-nrow(XY.sites)), na.rm=TRUE, xy=TRUE)[,1:2]
      XY.sites <- rbind(XY.sites, extra)
    }
  }
  
  if (prespecify.sites == FALSE & new.site.selection == "manual") {
    X <- sum(occ)
    plot(X)
    cat('\n',"Click on the map to select sites then press Esc") 
    XY.sites <- "click"(X, n=Inf, id=FALSE, xy=TRUE, type="p", col="red")
    XY.sites <- XY.sites[,1:2]
    rm(X)
  }
  
  if (prespecify.sites == FALSE & new.site.selection == "maxocc") {
    X <- sum(occ)
    test2 <- unlist(sort(as.vector(X), decreasing = TRUE, index.return=TRUE)$x)[n.sites]
    XY.sites <- xyFromCell(X, Which(X >= test2, cells=TRUE,na.rm=TRUE))
    rm(X)
  }
  
  if (plotter == TRUE) {
    plot(remote)
    points(XY.sites, col="red", pch=19)
  }
  
  rm(acc.ras, rem.ras)
  return(XY.sites)
}


#' Significance level function
#'
#' This function returns a critical value needed to calculate confidence intervals around the trend parameter given the Type I error rate and whether it is a one or two-tailed significance test.
#' @param two.tailed Set to TRUE for a two-tailed test, FALSE for a one-tailed test. 
#' @param alpha Set the type one error rate. Choose between 0.1, 0.05 and 0.01
#' @keywords 
#' @export
#' @examples
#' two.tailed <- TRUE
#' alpha <- 0.05
#' sig.test(two.tailed, alpha)

sig.test <- function(two.tailed, alpha) {
  values <- matrix(c(1.28, 1.645, 1.65, 1.96, 2.33, 2.58), ncol=3, nrow=2)
  if (two.tailed == TRUE) {ind <- 2} else {ind <- 1}
  if (alpha == 0.1) {value <- values[ind,1]}
  if (alpha == 0.05) {value <- values[ind,2]}
  if (alpha == 0.01) {value <- values[ind,3]}
  return(value)
}


#' Model fire function
#'
#' This function models the incidence of fire at monitoring sites based on a hazard function that depends on time since fire in each cell.
#' At present, separate hazard functions are defined for three different vegetation classes. It returns a vector specifying whether each site burns (1) or not (0) in a given year jj
#' @param time.hist An array specifying the time since fire at each monitoring site prior to simulations
#' @param time.fire An array specifying the time since fire at each monitoring site during simulations
#' @param veg.ID Specifies the different layers with different hazard functions
#' @param jj The year of simulation. Can be a value between 1 and Tmax
#' @keywords 
#' @export
#' @examples
#' fire.point.model()

#Model the incidence of fire at selected sites according using a hazard function
fire.point.model <- function(time.hist, time.fire, veg.ID, jj){
  burn <- veg.ID
  OW <- which(veg.ID==2) #Record which cells are open woodland
  CF <- which(veg.ID==3) #Record which cells are open forest
  SW <- which(veg.ID==1) #Record which cells are sandstone and woodland
  if (jj==1) {
    burn[OW] <- 1-((exp(-(time.hist[OW]/2.24)^1.22))/(exp(-((time.hist[OW]-1)/2.24)^1.22))) #Discrete Weibull hazard function for open woodland
    burn[CF] <- 1-((exp(-(time.hist[CF]/4.34)^0.96))/(exp(-((time.hist[CF]-1)/4.34)^0.96))) #Discrete Weibull hazard function for closed forest
    burn[SW] <- 1-((exp(-(time.hist[SW]/3.56)^1.26))/(exp(-((time.hist[SW]-1)/3.56)^1.26))) #Discrete Weibull hazard function for sandstone and woodland
  } else {
    burn[OW] <- 1-((exp(-(time.fire[OW,jj-1]/2.24)^1.22))/(exp(-((time.fire[OW,jj-1]-1)/2.24)^1.22))) #Discrete Weibull hazard function for open woodland
    burn[CF] <- 1-((exp(-(time.fire[CF,jj-1]/4.34)^0.96))/(exp(-((time.fire[CF,jj-1]-1)/4.34)^0.96))) #Discrete Weibull hazard function for closed forest
    burn[SW] <- 1-((exp(-(time.fire[SW,jj-1]/3.56)^1.26))/(exp(-((time.fire[SW,jj-1]-1)/3.56)^1.26))) #Discrete Weibull hazard function for sandstone and woodland
  }
  burn <- ifelse(runif(n.sites)<burn,1,0)
  return(burn)
}

#Model the incidence of fire at selected sites using a hazard function
fire.point.model2 <- function(time.hist, jj){
  PP <- c(0.48,0.53,0.42,0.40,0.39,0.26,0.32,0.25,0.28,0.30,0.26,0.21,0.22,0.26,0.05)
  #PP[] <- 1
  burn <- PP[time.hist]
  burn <- ifelse(runif(n.sites)<burn,1,0)
  return(burn)
}

#' Refit occupancy raster layers function
#'
#' This function is called only if fire is modelled at sites. It re-fits occupancy raster layers for fire sensitive species given simulated time since fire and number of fire layers
#' It returns a new raster stack with the updated occupancy raster layers for each species
#' @param occ.new A raster stack of occupancy maps. Raster layers for fire sensitive species are updated depending on the simulated fire history at time jj
#' @param layers An array of covararite values at each site. These are used to update the occ.new raster stack   
#' @param time.fire A vector specifying the time since fire at each site in year jj given the simulated fire history
#' @param fire.freq A vector specifying the number of fires at sites during a 15 year moving window
#' @param jj The year of the monitoring program
#' @keywords 
#' @export
#' @examples
#' refit.occ()

#Function to refit the occupancy models for fire sensitive species should the user choose to model fire at sites
refit.occ <- function(occ.new, layers, time.fire, fire.freq, jj, veg.ID) {
  for (i in 1:n.species) {
    covr <- as.numeric(species.cov[i,-1])
    if (covr[14] | covr[16] !=0) { #Only remap occupancy if fire influences occupancy
      freq.scale <- array(scale(fire.freq[,jj])[,1])
      time.scale <- array(scale(time.fire[,jj])[,1])
      #occ.new[,i,jj] <- rep(covr[1],n.sites) + covr[3]*layers[,1] + covr[4]*layers[,1]^2 + covr[5]*layers[,2] + 
      #covr[7]*freq.scale + covr[8]*freq.scale^2 + covr[9]*layers[,7] + covr[10]*layers[,7]^2 + 
      #covr[11]*layers[,3] + covr[12]*layers[,3]^2 + covr[14]*layers[,6] + covr[15]*time.scale + covr[15]*time.scale^2
      
      occ.new[,i,jj] <- rep(covr[1],n.sites) + 
        covr[2]*layers[,8] + #Clay
        covr[3]*layers[,8]^2 + 
        covr[4]*layers[,11] + #Veg
        covr[5]*layers[,11]^2 + 
        covr[6]*layers[,1] + #DEM
        covr[7]*layers[,1]^2 + 
        covr[8]*layers[,2] + #Creek
        covr[9]*layers[,9] + #Ruggedness
        covr[10]*layers[,5] + #Temp
        covr[11]*layers[,5]^2 + 
        covr[12]*layers[,10] + #Rain
        covr[13]*layers[,10]^2 +
        covr[14]*freq.scale + #Fire freq
        covr[15]*freq.scale^2 +
        covr[16]*time.scale + #Time since fire
        covr[17]*time.scale^2
      #covr[18]*layers[,1] + #Extent
      #covr[19]*layers[,1] #Patch
      
      occ.new[,i,jj] <- exp(occ.new[,i,jj])/(1+exp(occ.new[,i,jj]))
      
      if (i == 5 & nrow(species.list == 20)) {occ.new[veg.ID==0,i,jj] <- 0}
      
    }
  }
  return(occ.new)
}

#' Refit detectability raster layers function
#'
#' This function is called only if fire is modelled at sites. It re-fits the detectability raster layers for fire sensitive species given simulated time since fire and number of fire layers
#' A new raster stack is returned for each detection method
#' @param det.method1.new A raster stack of detectability under method 1. Raster layers for fire sensitive species are updated depending on the simulated fire history at time jj
#' @param det.method2.new A raster stack of detectability under method 2. Raster layers for fire sensitive species are updated depending on the simulated fire history at time jj
#' @param det.method3.new A raster stack of detectability under method 3. Raster layers for fire sensitive species are updated depending on the simulated fire history at time jj
#' @param layers An array of covararites at each site. These are used to update the detectability raster stacks 
#' @param time.fire A vector specifying the time since fire at each site in year jj given the simulated fire history
#' @param fire.freq A vector specifying the number of fires at sites during a 15 year moving window
#' @param jj The year of the monitoring program
#' @keywords 
#' @export
#' @examples
#' refit.det()

#Function to refit the detection models should the user choose to model fire in the landscape 
refit.det <- function(det.method1.new, det.method2.new, det.method3.new, det.method4.new, layers, time.fire, fire.freq, jj) {
  det.method <- rep(1,n.sites)
  freq.scale <- array(scale(fire.freq)[,jj])
  time.scale <- array(scale(time.fire)[,jj])
  obj <- list(det.method1.new, det.method2.new, det.method3.new, det.method4.new)
  
  for (i in 1:n.species) {
    
    covr <- as.numeric(species.cov[i,-1])
    det.covr <- covr[20:27]
    
    #Species detected using 1 method
    if ((sum(species.list[i,-1]) == 1) & (covr[22] | covr[23] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      logit <- det.covr[1] + det.covr[2]*layers[,8] + det.covr[3]*freq.scale + det.covr[4]*time.scale 
      obj[[pos]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
    #Species detected using 2 methods
    if ((sum(species.list[i,-1]) == 2) & (covr[22] | covr[23] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      
      det.method[] <- 0  #Method 1
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[2]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale  
      obj[[pos[1]]][,i,jj] <- exp(logit)/(1+exp(logit))
      
      det.method[] <- 1  #Method 2
      logit <- covr[2] + covr[17]*freq.scale + det.covr[pos[2]+3]*det.method + covr[22]*layers[,6] + covr[23]*time.scale 
      obj[[pos[2]]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
    #Species detected using 3 methods
    if ((sum(species.list[i,-1]) == 3) & (covr[22] | covr[23] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      
      det.method[] <- 0  #Method 1
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[2]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale  
      obj[[pos[1]]][,i,jj] <- exp(logit)/(1+exp(logit))
      
      det.method[] <- 1  #Method 2
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[2]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale 
      obj[[pos[2]]][,i,jj] <- exp(logit)/(1+exp(logit))
      
      det.method[] <- 1  #Method 3
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[3]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale
      obj[[pos[3]]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
    #Species detected using 4 methods
    if ((sum(species.list[i,-1]) == 4) & (covr[22] | covr[23] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      
      det.method[] <- 0  #Method 1
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[2]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale 
      obj[[pos[1]]][,i,jj] <- exp(logit)/(1+exp(logit))
      
      det.method[] <- 1  #Method 2
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[2]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale 
      obj[[pos[2]]][,i,jj] <- exp(logit)/(1+exp(logit))
      
      det.method[] <- 1  #Method 3
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[3]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale
      obj[[pos[3]]][,i,jj] <- exp(logit)/(1+exp(logit))
      
      det.method[] <- 1  #Method 3
      logit <- det.covr[1] + det.covr[3]*freq.scale + det.covr[pos[4]+3]*det.method + det.covr[2]*layers[,8] + det.covr[4]*time.scale
      obj[[pos[4]]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
    
    #if (sp.lst$Method1 == 1 & sp.lst$Method2 == 0 & sp.lst$Method3 == 1 & covr[17] | covr[23] != 0){ #Species detected using method 1 and method 3
    
    #det.method[] <- 1  #Method 1 first
    #det.method1.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[18]*det.method + covr[22]*layers[,6] + covr[23]*time.scale 
    #det.method1.new[,i,jj] <- exp(det.method1.new[,i,jj])/(1+exp(det.method1.new[,i,jj]))
    
    #det.method[] <- 0  #Method 3 second
    # det.method3.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[18]*det.method + covr[22]*layers[,6] + covr[23]*time.scale  
    #det.method3.new[,i,jj] <- exp(det.method3.new[,i,jj])/(1+exp(det.method3.new[,i,jj]))
    # }
    
    # if (sp.lst$Method1 == 0 & sp.lst$Method2 == 1 & sp.lst$Method3 == 1 & covr[17] | covr[23] != 0){ #Species detected using method 2 and method 3
    
    #det.method[] <- 1  #Method 2
    # det.method2.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[18]*det.method + covr[22]*layers[,6] + covr[23]*time.scale 
    # det.method2.new[,i,jj] <- exp(det.method2.new[,i,jj])/(1+exp(det.method2.new[,i,jj]))
    
    # det.method[] <- 0  #Method 3
    # det.method3.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[18]*det.method + covr[22]*layers[,6] + covr[23]*time.scale  
    # det.method3.new[,i,jj] <- exp(det.method3.new[,i,jj])/(1+exp(det.method3.new[,i,jj]))
    # }
    
    # if (sp.lst$Method1 == 1 & sp.lst$Method2 == 1 & sp.lst$Method3 == 1 & covr[17] | covr[23] != 0){ #Species detected using methods 1 and 2 and 3
    
    # det.method[] <- 1  #Method 1
    # det.method1.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[20]*det.method + covr[22]*layers[,6] + covr[23]*time.scale 
    #  det.method1.new[,i,jj] <- exp(det.method1.new[,i,jj])/(1+exp(det.method1.new[,i,jj]))
    
    #  det.method[] <- 1  #Method 2
    #  det.method2.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[19]*det.method + covr[22]*layers[,6] + covr[23]*time.scale 
    #  det.method2.new[,i,jj] <- exp(det.method2.new[,i,jj])/(1+exp(det.method2.new[,i,jj]))
    
    #  det.method[] <- 0  #Repeat 3
    #  det.method3.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[19]*det.method + covr[22]*layers[,6] + covr[23]*time.scale  
    #  det.method3.new[,i,jj] <- exp(det.method3.new[,i,jj])/(1+exp(det.method3.new[,i,jj]))
    #}
  }
  return(obj)
}


#Species detected using 1 method
#if (sp.lst$Method1 == 1 & sp.lst$Method2 == 0 & sp.lst$Method3 == 0 & covr[17] | covr[23] != 0){ 
#det.method1.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[22]*layers[,6] + covr[23]*time.scale 
#det.method1.new[,i,jj] <- exp(det.method1.new[,i,jj])/(1+exp(det.method1.new[,i,jj]))
#  }
#Species detected using method 2 only
# if (sp.lst$Method1 == 0 & sp.lst$Method2 == 1 & sp.lst$Method3 == 0 & covr[17] | covr[23] != 0){ 
#det.method2.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[22]*layers[,6] + covr[23]*time.scale 
#det.method2.new[,i,jj] <- exp(det.method2.new[,i,jj])/(1+exp(det.method2.new[,i,jj]))
# }
#Species detected using method 3 only
# if (sp.lst$Method1 == 0 & sp.lst$Method2 == 0 & sp.lst$Method3 == 1 & covr[17] | covr[23] != 0){ 
#det.method3.new[,i,jj] <- covr[2] + covr[17]*freq.scale + covr[22]*layers[,6] + covr[23]*time.scale 
#det.method3.new[,i,jj] <- exp(det.method3.new[,i,jj])/(1+exp(det.method3.new[,i,jj]))
# }


#' Fit occupancy model with one detection method
#'
#' This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected using one detection method
#' A trend in occupancy is estimated and confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
#' A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
#' intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from each regional level management unit. 
#' @param method The detection method relevant to species ss
#' @param repeats The number of repeat visits for the specified detection method
#' @param s.years A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)
#' @param n.sites The number of sites monitored
#' @param xy.sites The XY coordinates of monitored sites
#' @param park.ID A vector specifying that location of each site with sub-level parks
#' @param park.level Set to TRUE is power is estimated within regional level management unit, FALSE otherwise
#' @param powcnt A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape
#' @param fail A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a landscape level in unmarked
#' @param pow.park A vector that keeps track of how many times a significant trend in occupancy is detected within each park
#' @param fail.park A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a park level in unmarked
#' @param value The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test
#' @param ss An index to loop through each species 
#' @param two.tailed Set to TRUE if conducting a two-tailed test, FALSE otherwise
#' @keywords 
#' @export
#' @examples
#' fit.occ.1method()

#Function to fit occupancy model when only 1 method is used to detect species s
fit.occ.1method <- function(method, repeats, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats == 1) {
    y <- aperm(method[,ss,,s.years],c(1,2))
    dim(y)<- c(n.sites*length(s.years), repeats)
  } else {
    y <- aperm(method[,ss,,s.years],c(1,3,2))
    dim(y)<- c(n.sites*length(s.years), repeats)
  }
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  
  #if 1 repeat use logistic regression  
  if (repeats == 1) {
    Time <- siteCovs[,1]
    y <- as.vector(y)
    mod <- glm(formula = y ~ Time, family=binomial)  
    if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
      Time.mean <- summary(mod)$coefficients[2,1]
      Time.CI <- summary(mod)$coefficients[2,2] * value
      Time.coeff<-summary(mod)$coefficients[2,1]
      Time.SE<-summary(mod)$coefficients[2,2]
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { 
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { 
        powcnt[ss] <- powcnt[ss] + 1}
    }  
  }
  
  
  if (repeats > 1) {
    umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = NULL)
    inits<- rbind(c(0,0,0),c(-1,-1,-1),c(1,1,1))
    for (vv in 1:nrow(inits)) {
      mod<- try(occu(~1 ~Time, umf, control=list(maxit=1000000), starts=inits[vv,]), silent=TRUE)
      if (!is(mod,"try-error")){break}
    }
  #}
  
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  }
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        
        inits<- rbind(c(0,0,0),c(-1,-1,-1),c(1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~1 ~Time, y.park, control=list(maxit=1000000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}

#' Fit occupancy model with two detection methods
#'
#' This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected using two detection methods
#' A trend in occupancy is estimated and confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
#' A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
#' intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from each regional level management unit. 
#' @param method1 The first detection method relevant to species ss
#' @param method2 The second detection method relevant to species ss
#' @param repeats1 The number of repeat visits for the first detection method
#' @param repeats2 The number of repeat visits for the second detection method
#' @param s.years A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)
#' @param n.sites The number of sites monitored
#' @param xy.sites The XY coordinates of monitored sites
#' @param park.ID A vector specifying that location of each site with sub-level parks
#' @param park.level Set to TRUE is power is estimated within regional level management unit, FALSE otherwise
#' @param powcnt A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape
#' @param fail A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a landscape level in unmarked
#' @param pow.park A vector that keeps track of how many times a significant trend in occupancy is detected within each park
#' @param fail.park A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a park level in unmarked
#' @param value The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test
#' @param ss An index to loop through each species 
#' @param two.tailed Set to TRUE if conducting a two-tailed test, FALSE otherwise
#' @keywords 
#' @export
#' @examples
#' fit.occ.2method()

#Function to fit occupancy model when only 2 methods are used to detect species s
fit.occ.2method <- function(method1, method2, repeats1, repeats2, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats1 == 1) {
    y1 <- aperm(method1[,ss,,s.years],c(1,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  } else {
    y1 <- aperm(method1[,ss,,s.years],c(1,3,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  }
  if (repeats2 == 1) {
    y2 <- aperm(method2[,ss,,s.years],c(1,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  } else {
    y2 <- aperm(method2[,ss,,s.years],c(1,3,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  }
  
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  
  m1.code <- matrix(factor(1),nrow=nrow(xy.sites)*length(s.years), ncol=repeats1)
  m2.code <- matrix(factor(2),nrow=nrow(xy.sites)*length(s.years), ncol=repeats2)
  obsCovs <- data.frame(as.matrix(cbind(m1.code,m2.code)))
  colnames(obsCovs) <- c(rep("M1",repeats1),rep("M2",repeats2))
  
  umf <- unmarkedFrameOccu(y = cbind(y1, y2), siteCovs = siteCovs, obsCovs = list(METHOD = obsCovs))
  
  inits<- rbind(c(0,0,0,0),c(-1,-1,-1,-1),c(1,1,1,1))
  for (vv in 1:nrow(inits)) {
    mod<- try(occu(~METHOD ~Time, umf, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
    if (!is(mod,"try-error")){break}
  }
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        
        inits<- rbind(c(0,0,0,0),c(-1,-1,-1,-1),c(1,1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~METHOD ~Time, y.park, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}

#' Fit occupancy model with three detection methods
#'
#' This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected using 3 detection methods
#' A trend in occupancy is estimated and confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
#' A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
#' intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from each regional level management unit. 
#' @param method1 The first detection method relevant to species ss
#' @param method2 The second detection method relevant to species ss
#' @param method3 The third detection method relevant to species ss
#' @param repeats1 The number of repeat visits for the first detection method
#' @param repeats2 The number of repeat visits for the second detection method
#' @param repeats3 The number of repeat visits for the third detection method
#' @param s.years A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)
#' @param n.sites The number of sites monitored
#' @param xy.sites The XY coordinates of monitored sites
#' @param park.ID A vector specifying that location of each site with sub-level parks
#' @param park.level Set to TRUE is power is estimated within regional level management unit, FALSE otherwise
#' @param powcnt A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape
#' @param fail A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a landscape level in unmarked
#' @param pow.park A vector that keeps track of how many times a significant trend in occupancy is detected within each park
#' @param fail.park A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a park level in unmarked
#' @param value The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test
#' @param ss An index to loop through each species 
#' @param two.tailed Set to TRUE if conducting a two-tailed test, FALSE otherwise
#' @keywords 
#' @export
#' @examples
#' fit.occ.3method()

#Function to fit occupancy model when 3 methods are used to detect species s
fit.occ.3method <- function(method1, method2, method3, repeats1, repeats2, repeats3, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats1 == 1) {
    y1 <- aperm(method1[,ss,,s.years],c(1,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  } else {
    y1 <- aperm(method1[,ss,,s.years],c(1,3,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  }
  if (repeats2 == 1) {
    y2 <- aperm(method2[,ss,,s.years],c(1,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  } else {
    y2 <- aperm(method2[,ss,,s.years],c(1,3,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  }
  if (repeats3 == 1) {
    y3 <- aperm(method3[,ss,,s.years],c(1,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  } else {
    y3 <- aperm(method3[,ss,,s.years],c(1,3,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  }
  
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  
  m1.code <- matrix(factor(1),nrow=nrow(xy.sites)*length(s.years), ncol=repeats1)
  m2.code <- matrix(factor(2),nrow=nrow(xy.sites)*length(s.years), ncol=repeats2)
  m3.code <- matrix(factor(3),nrow=nrow(xy.sites)*length(s.years), ncol=repeats3)
  obsCovs <- data.frame(as.matrix(cbind(m1.code, m2.code, m3.code)))
  colnames(obsCovs) <- c(rep("M1",repeats1), rep("M2",repeats2), rep("M3",repeats3))
  
  umf <- unmarkedFrameOccu(y = cbind(y1, y2, y3), siteCovs = siteCovs, obsCovs = list(METHOD = obsCovs))
  
  inits<- rbind(c(0,0,0,0,0),c(-1,-1,-1,-1,-1),c(1,1,1,1,1))
  for (vv in 1:nrow(inits)) {
    mod<- try(occu(~METHOD ~Time, umf, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
    if (!is(mod,"try-error")){break}
  }
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        inits<- rbind(c(0,0,0,0,0),c(-1,-1,-1,-1,-1),c(1,1,1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~METHOD ~Time, y.park, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}


#Function to fit occupancy model when 4 methods are used to detect species s
fit.occ.4method <- function(method1, method2, method3, method4, repeats1, repeats2, repeats3, repeats4, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats1 == 1) {
    y1 <- aperm(method1[,ss,,s.years],c(1,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  } else {
    y1 <- aperm(method1[,ss,,s.years],c(1,3,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  }
  if (repeats2 == 1) {
    y2 <- aperm(method2[,ss,,s.years],c(1,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  } else {
    y2 <- aperm(method2[,ss,,s.years],c(1,3,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  }
  if (repeats3 == 1) {
    y3 <- aperm(method3[,ss,,s.years],c(1,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  } else {
    y3 <- aperm(method3[,ss,,s.years],c(1,3,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  }
  if (repeats4 == 1) {
    y4 <- aperm(method4[,ss,,s.years],c(1,2))
    dim(y4)<- c(n.sites*length(s.years), repeats4)
  } else {
    y4 <- aperm(method4[,ss,,s.years],c(1,3,2))
    dim(y4)<- c(n.sites*length(s.years), repeats4)
  }
  
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  
  m1.code <- matrix(factor(1),nrow=nrow(xy.sites)*length(s.years), ncol=repeats1)
  m2.code <- matrix(factor(2),nrow=nrow(xy.sites)*length(s.years), ncol=repeats2)
  m3.code <- matrix(factor(3),nrow=nrow(xy.sites)*length(s.years), ncol=repeats3)
  m4.code <- matrix(factor(4),nrow=nrow(xy.sites)*length(s.years), ncol=repeats4)
  obsCovs <- data.frame(as.matrix(cbind(m1.code, m2.code, m3.code, m4.code)))
  
  colnames(obsCovs) <- c(rep("M1",repeats1), rep("M2",repeats2), rep("M3",repeats3), rep("M4",repeats4))
  
  umf <- unmarkedFrameOccu(y = cbind(y1, y2, y3, y4), siteCovs = siteCovs, obsCovs = list(METHOD = obsCovs))
  
  inits<- rbind(c(0,0,0,0,0,0),c(-1,-1,-1,-1,-1,-1),c(1,1,1,1,1,1))
  for (vv in 1:nrow(inits)) {
    mod<- try(occu(~METHOD ~Time, umf, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
    if (!is(mod,"try-error")){break}
  }
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        inits<- rbind(c(0,0,0,0,0,0),c(-1,-1,-1,-1,-1,-1),c(1,1,1,1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~METHOD ~Time, y.park, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}

#' Run power analysis
#'
#' This function runs the main power analysis. It starts by creating empty vectors to record the number of times a significant trend is detected in the simulated detection histories
#' It then simulates fire at sites if model.fire = TRUE, and models a decline in occupancy over time depending on the specified effect size. 
#' Detection histories are simulated at each site given the occupancy status at each site, the number of survey methods used to detect a species and the detection probability of each relevent method.
#' Simulated detection histories are loaded into the package unmarked, and a trend in occupancy is modelled over time. The proportion of times a significant trend is detected is recorded.
#' This process is repeated for each effect size.
#' @param effect.size An integer specifying the proportional reduction in occupancy between the start and end of the monitoring program
#' @param nsims An integer specifying the number of simulations 
#' @param alpha The Type I error rate. Must be either 0.1, 0.05 or 0.01
#' @param decline Specifies the pattern of decline across space. Set to 'random' to model a constant 'blanket' decline across space
#' @param Tmax An integer specifying the length of the monitoring program
#' @param s.years A vector specifying the years in which monitoring occurs. Note, the final year or monitoring must be equal to Tmax
#' @param sites An array containing the XY coordinates of pre-specified monitoring sites
#' @param fire.hist Delete
#' @param model.fire Set to TRUE if the incidence of fire is modelled at sites, FALSE to model a static landscape over time
#' @param species.list An array containing the name of each species in the first column, and a zero or one describing how many methods are used to detect the species
#' @param park.level Set to TRUE if power is to be estimated within sub-level management units
#' @param fire.fr Delete
#' @param combined.effect An array that records the reduction in occupancy across sites due to the effect size and fire
#' @param xy.sites The XY coordinates of monitoring sites
#' @param two.tailed Set to TRUE is conducting a two-tailed significance test, FALSE if conducting a one-tailed test
#' @param plotter Set to TRUE to plot the location of sites at the start of each simulation, FALSE otherwise
#' @param n.sites The number of sites to be monitored.  
#' @param n.species The number of species for which monitoring is simulated. Is equal to the number of layers in the occ raster stack
#' @param n.park The number of sub-level management units in which power is estimated
#' @param randomise.sites Set to TRUE if sites are to randomly selected for monitoring, FALSE if pre-specified sites are surveyed 
#' @param all.prespecify.sites Set to TRUE if all of the prespecified sites are to be monitored, FALSE if fewer or more sites are to be monitored
#' @param occ.new An array containing the occupancy value for each species at monitoring sites over time 
#' @param occ.time An array containing the occupancy value for each species at monitoring sites over time 
#' @param det.method1.new An array containing the detectability values of method 1 for each species at monitoring sites over time
#' @param det.method2.new An array containing the detectability values of method 2 for each species at monitoring sites over time
#' @param det.method3.new An array containing the detectability values of method 3 for each species at monitoring sites over time
#' @param det.method1.time An array containing the detectability values of method 1 for each species at monitoring sites over time
#' @param det.method2.time An array containing the detectability values of method 2 for each species at monitoring sites over time
#' @param det.method3.time An array containing the detectability values of method 3 for each species at monitoring sites over time
#' @param det1.method1 An array that records simulated detection histories at sites using method 1
#' @param det1.method2 An array that records simulated detection histories at sites using method 2
#' @param det1.method3 An array that records simulated detection histories at sites using method 3
#' @param fire.freq Delete!
#' @param n.method1 The number of repeat visits to sites using method 1
#' @param n.method2 The number of repeat visits to sites using method 2
#' @param n.method3 The number of repeat visits to sites using method 3
#' @keywords 
#' @export
#' @examples
#' run.power()

#Run pwr analysis at a regional level
#run.power <- function(effect.size, nsims, alpha, det1, working.raster, decline, plot.decline, plot.fire, Tmax, s.years, n.remote, n.4WD, sites, fire.hist, burn.area, study.area, model.fire, species.list, park.level, fire.fr, combined.effect, xy.sites) {
#run.power <- function(effect.size, nsims, alpha, working.raster, decline, Tmax, s.years, sites, fire.hist, model.fire, species.list, park.level, fire.fr, combined.effect, xy.sites, two.tailed, plotter, n.sites, n.species, n.park, randomise.sites, all.prespecify.sites) {
run.power <- function(effect.size, nsims, alpha, decline, Tmax, s.years, sites, model.fire, species.list, park.level, xy.sites, R, two.tailed, plotter, n.sites, n.species, n.park, all.prespecify.sites, prespecify.sites, new.site.selection, occ.time, det.method1.time, det.method2.time, det.method3.time, det.method4.time, n.method, occ, det.method1, det.method2, det.method3, det.method4, veg, remote, parks, trend, variation, model.variation) {
  
  det1.method1 <- array(NA, dim=c(n.sites,n.species,n.method[1],Tmax))
  det1.method2 <- array(NA, dim=c(n.sites,n.species,n.method[2],Tmax))
  det1.method3 <- array(NA, dim=c(n.sites,n.species,n.method[3],Tmax)) 
  det1.method4 <- array(NA, dim=c(n.sites,n.species,n.method[4],Tmax))
  det.method1.new <- det.method2.new <- det.method3.new <- det.method4.new <- occ.new <- array(0, dim=c(n.sites,n.species,Tmax))
  combined.effect<- matrix(NA,nsims,n.species)
  fr.hist <- time.hist <-veg.ID <- park.ID <- rep(NA,n.sites)
  time.fire <- array(NA, dim=c(n.sites,Tmax))
  fire.freq <- array(NA, dim=c(n.sites,Tmax))
  
  #Set up vectors and matrices to record the number of instances when we detect a change in occupancy
  powcnt <-rep(0,n.species) #REMOVE FROM LOOP
  pow.park <- matrix(0,ncol=n.species,nrow=length(n.park)) #REMOVE FROM LOOP
  fail <- rep(0,n.species) #REMOVE FROM LOOP
  fail.park <- matrix(0,ncol=n.species,nrow=length(n.park)) #REMOVE FROM LOOP
  
  for (ii in 1:nsims){ #Loop through number of simulations from 1 to nsims
    cat('\n',"Effect size = ", effect.size, " Simulation = ", ii) 
    
    #Select survey sites
    if (all.prespecify.sites == FALSE & new.site.selection != "maxocc" & new.site.selection != "manual") {
      cat('\n',"Selecting survey sites.....") 
      xy.sites <- select.sites(sites=sites, n.sites=n.sites, R=R, all.prespecify.sites=all.prespecify.sites, prespecify.sites=prespecify.sites, new.site.selection=new.site.selection, plotter=plotter) #Returns matrix of XY coordinates plus description for which park each point belongs in 	
      
      #Extract occ and det values from selected sites
      cat('\n',"Extracting values from selected sites.....") 
      occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)]
      
      for (ss in 1:n.species){
        if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]}
        if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]}
        if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]}
        if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
      }
    }
    
    #Extract park ID, veg type and fire histories from selected sites
    park.ID <- parks[cellFromXY(parks, xy.sites)]
    xy.sites <- cbind(xy.sites,park.ID)
    
    #Make a copy of the occ and det values at each site so they can be manipulated by fire and the effect size over time
    det.method1.new[] <- det.method1.time
    det.method2.new[] <- det.method2.time
    det.method3.new[] <- det.method3.time
    det.method4.new[] <- det.method4.time
    occ.new[] <- occ.time
    
    #Model fire at monitoring sites if the model.fire == TRUE
    if (model.fire == TRUE) { #If we specified that fire be modelled then simulate fire spread during the year
      cat('\n',"Modelling fire and updating occ and det layers.....")
      fire.ID <- matrix(NA, nrow = n.sites, ncol=nlayers(fire) + Tmax)
      veg.ID <- veg[cellFromXY(veg, xy.sites)]
      layers <- covar[cellFromXY(covar, xy.sites)]
      fire.ID[,1:nlayers(fire)] <- fire[cellFromXY(fire, xy.sites)]
      
      time.hist <- apply(fire.ID[,1:nlayers(fire)], 1, function(x) (nlayers(fire)+1)-max(which(x==1))) #time.hist <- round(fire.hist[cellFromXY(fire.hist, xy.sites)]) #time.fire <- 
      fr.hist <- rowSums(fire.ID[,1:nlayers(fire)]) #fire.freq <- rowSums(fire.ID[,1:nlayers(fire)])
      
      for (jj in 1:Tmax) { #For each simulation loop ii, loop through time from 1 to Tmax
        #burn <- fire.point.model(time.hist, time.fire, veg.ID, jj) #Simulate whether monitoring sites burn
        burn <- fire.point.model2(time.hist, jj)
        fire.ID[,nlayers(fire) + jj] <- burn #Record fire history at sites for year jj
        if (jj==1) {
          time.fire[,jj] <- ifelse(burn == 1, 1, time.hist+1)
          fire.freq[,jj] <- rowSums(fire.ID[, (jj+1):(jj+nlayers(fire))]) #Sum number of fires in 15 year moving window
        } 
        if (jj>1) {
          time.fire[,jj] <- ifelse(burn == 1, 1, time.fire[,jj-1]+1)
          fire.freq[,jj] <- rowSums(fire.ID[, (jj+1):(jj+nlayers(fire))]) #Sum number of fires in 15 year moving window
        }	
        if (jj %in% s.years) { 
          occ.new <- refit.occ(occ.new, layers, time.fire, fire.freq, jj, veg.ID) 
          det.out <- refit.det(det.method1.new=det.method1.new, det.method2.new=det.method2.new, det.method3.new=det.method3.new, det.method4.new=det.method4.new, layers, time.fire, fire.freq, jj)
        }
      }
      det.method1.new[] <- det.out[[1]]
      det.method2.new[] <- det.out[[2]]
      det.method3.new[] <- det.out[[3]]
      det.method4.new[] <- det.out[[4]]
    }
    
    if (trend == "decreasing") { #Model decreasing trend
      effect.time <- 1-(effect.size/Tmax*c(1:Tmax))
      if (model.variation == TRUE) {
        effect.time <- runif(Tmax,effect.time-variation[ss]/2,effect.time+variation[ss]/2)
        effect.time[effect.time>1] <- 1
        effect.time[effect.time<0] <- 0
      }
      for (i in 1:Tmax) {occ.new[,,i] <- occ.new[,,i]*effect.time[i]}
    } 
    
    if (trend == "increasing") { #Model increasing trend
      effect.time <- effect.size/Tmax*c(1:Tmax)
      if (model.variation == TRUE) {
        effect.time <- runif(Tmax,effect.time-variation[ss]/2,effect.time+variation[ss]/2)
        effect.time[effect.time>1] <- 1
        effect.time[effect.time<0] <- 0
      }
      for (i in 1:Tmax) {occ.new[,,i] <- occ.new[,,i] + ((1-occ.new[,,i])*effect.time[i])}
    } 
    
    
    combined.effect[ii,] <- (occ.time[1,,1]-occ.new[1,,Tmax])/occ.time[1,,1]
    
    if (decline== "random") {
      for (i in 1:Tmax) {
        occ.new[,,i] <- ifelse(runif(ncell(occ.new[,,i])) < occ.new[,,i],1,0)
      }
    }	
    
    cat('\n',"Simulating monitoring at sites.....")
    for (jj in 1:n.method[1]){
      for (tt in 1:Tmax) { #MIGHT BE ABLE TO REMOVE TT LOOP 
        det1.method1[,,jj,tt] <- ifelse(runif(ncell(det.method1.new[,,tt])) < det.method1.new[,,tt] & occ.new[,,tt] == 1, 1, 0)
      }}
    
    for (jj in 1:n.method[2]){
      for (tt in 1:Tmax) { #MIGHT BE ABLE TO REMOVE TT LOOP
        det1.method2[,,jj,tt] <- ifelse(runif(ncell(det.method2.new[,,tt])) < det.method2.new[,,tt] & occ.new[,,tt] == 1,1,0)
      }
    }
    
    for (jj in 1:n.method[3]){
      for (tt in 1:Tmax) { #MIGHT BE ABLE TO REMOVE TT LOOP
        det1.method3[,,jj,tt] <- ifelse(runif(ncell(det.method3.new[,,tt])) < det.method3.new[,,tt] & occ.new[,,tt] == 1 ,1,0)
      }
    }
    
    for (jj in 1:n.method[4]){
      for (tt in 1:Tmax) { #MIGHT BE ABLE TO REMOVE TT LOOP
        det1.method4[,,jj,tt] <- ifelse(runif(ncell(det.method4.new[,,tt])) < det.method4.new[,,tt] & occ.new[,,tt] == 1 ,1,0)
      }
    }
    
    #Combine the results of each detection method for each species, create detection histories and fit occupancy model 
    
    cat('\n',"Fitting occ model.....")
    value <- sig.test(two.tailed, alpha)
    methods <- list(det1.method1,det1.method2,det1.method3,det1.method4)
    
    for (ss in 1:n.species) {
      
      #Fit occ model for species detected using method 1
      if (sum(species.list[ss,-1]) == 1) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.1method(method=methods[[pos]], repeats=n.method[pos], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      #Fit occ model for species detected using 2 methods
      if (sum(species.list[ss,-1]) == 2) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.2method(method1=methods[[pos[1]]], method2=methods[[pos[2]]], repeats1=n.method[pos[1]], repeats2=n.method[pos[2]], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      #Fit occ model for species detected using 3 methods
      if (sum(species.list[ss,-1]) == 3) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.3method(method1=methods[[pos[1]]], method2=methods[[pos[2]]], method3=methods[[pos[3]]], repeats1=n.method[pos[1]], repeats2=n.method[pos[2]], repeats3=n.method[pos[3]], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      #Fit occ model for species detected using 4 methods
      if (sum(species.list[ss,-1]) == 4) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.4method(method1=methods[[pos[1]]], method2=methods[[pos[2]]], method3=methods[[pos[3]]], method4=methods[[pos[4]]], repeats1=n.method[pos[1]], repeats2=n.method[pos[2]], repeats3=n.method[pos[3]], repeats4=n.method[pos[4]], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      #Fit occ model for species only detected using method 1
      #if (species.list[ss,2] == 1 & species.list[ss,3] == 0 & species.list[ss,4] == 0) {
      #occ.fit <- fit.occ.1method(method=det1.method1, repeats=n.method[1], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      # }
      
      #Sit occ model for species only detected using method 2
      #if (species.list[ss,2] == 0 & species.list[ss,3] == 1 & species.list[ss,4] == 0) {
      #occ.fit <- fit.occ.1method(method=det1.method2, repeats=n.method[2], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      #}
      
      #Fit occ model for species only detected using method 3
      #if (species.list[ss,2] == 0 & species.list[ss,3] == 0 & species.list[ss,4] == 1) {
      # occ.fit <- fit.occ.1method(method=det1.method3, repeats=n.method[3], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      #}
      
      #Fit occ model for species detected using method 1 and 2
      #if (species.list[ss,2] == 1 & species.list[ss,3] == 1 & species.list[ss,4] == 0) {
      #occ.fit <- fit.occ.2method(method1=det1.method1, method2=det1.method2, repeats1=n.method[1], repeats2=n.method[2], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      # }
      
      #Fit occ model for species detected using method 1 and 3
      # if (species.list[ss,2] == 1 & species.list[ss,3] == 0 & species.list[ss,4] == 1) {
      # occ.fit <- fit.occ.2method(method1=det1.method1, method2=det1.method3, repeats1=n.method[1], repeats2=n.method[3], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      #}
      
      #Fit occ model for species detected using method 2 and 3
      # if (species.list[ss,2] == 0 & species.list[ss,3] == 1 & species.list[ss,4] == 1) {
      #occ.fit <- fit.occ.2method(method1=det1.method2, method2=det1.method3, repeats1=n.method[2], repeats2=n.method[3], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      # }
      
      #Fit occ model for species detected using methods 1, 2 and 3
      # if (species.list[ss,2] == 1 & species.list[ss,3] == 1 & species.list[ss,4] == 1) {
      #  occ.fit <- fit.occ.3method(method1=det1.method1, method2=det1.method2, method3=det1.method3, repeats1=n.method[1], repeats2=n.method[2], repeats3=n.method[3], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      # }
      
      powcnt <- occ.fit[[1]]
      fail <- occ.fit[[2]]
      pow.park <- occ.fit[[3]]
      fail.park <- occ.fit[[4]]
    } #End species loop
    
  } #End sims loop (ii)
  
  #Calculate pwr from simulations for each species
  combined.effect <- colMeans(combined.effect)
  output <- rbind(powcnt,fail,pow.park,fail.park,combined.effect) #If running in parallel
  #pwr <- powcnt/(nsims-fail)
  #pwr.park <- pow.park/(nsims-fail.park)
  #print(fail)
  #output <- rbind(matrix(pwr,nrow=1,ncol=n.species),pwr.park,combined.effect)
  #cat("power = ", pwr , " seconds\n\n\n") 
  return(output)
}

#' Plot results of power analysis
#'
#' This function plots results of the power analysis. It first unpacks the list generated by the foreach function and estimates power given the number of simulations 
#' and the number of cores used. Power is plotted for each species on a single figure. Users can then scroll through separate plots for power at a park level if it park.power = TRUE
#' @param pwr An array containing the results of the power analysis
#' @param n.species The number of species
#' @param nsims The number of simulations
#' @param effect.size The proportional decline in occupancy
#' @param n.park The number of parks in which to restimate power
#' @param species.list A vector specifying the name of each species analysed
#' @param save.wd A directory to save results and figures
#' @keywords cats
#' @export
#' @examples
#' plot.results()

plot.results <- function(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) {
  Results <- array(dim=c(length(n.park)+1,n.species,length(effect.size)))
  new.effect <- matrix(NA, length(effect.size), n.species)
  pwr_all <- pwr
  for (i in 1:length(effect.size)) {
    pwr1 <- pwr_all[,i]
    dim(pwr1) <- c(length(pwr1)/n.species, n.species)
    results <- matrix(NA,length(n.park)+1,n.species)
    results[1,] <- pwr1[1,]/((nsims*n.cores)-pwr1[2,])
    results[2:(length(n.park)+1),] <- pwr1[3:(length(n.park)+2),]/((nsims*n.cores)-pwr1[(length(n.park)+3):(length(n.park)*2+2),])
    new.effect[i,] <- pwr1[nrow(pwr1),]/n.cores
    Results[,,i] <- results
  }
  
  if (park.level==TRUE) { #Change back to TRUE
    for (v in 1:(length(n.park)+1)) {
      cl <- rainbow(n.species)
      par(mfcol=c(1,1), mar=c(0.1,0.1,0.1,0.1), oma=c(4,4,4,4), mai=c(0.1,0.1,0.1,2),xpd=NA)
      new.effect[,1] <- c(0.1,0.3,0.5,0.7,0.9)
      #plot(Results[v,1,]~new.effect[,1], type="l",ylim=c(0,1), xlim=c(0,1), lwd=1, main="", col="grey", ylab="Statistical power", xlab="Effect size")
      for (i in 1:n.species) {
        new.effect[,i] <- c(0.1,0.3,0.5,0.7,0.9)
        #lines(Results[v,i,]~new.effect[,i], type="l",ylim=c(0,1), xlim=c(0,1), lwd=2, col=cl[i], lty=1)
      }
      #legend("topright", c(as.character(species.list[,1])), inset=c(-0.3,0), lwd=c(2,2,2,2,2,2), cex=0.5, col=cl[1:n.species])
      if (v==1) {mtext("Power: landscape level", side = 3, outer=TRUE, adj=0.5, line = 1, cex=1.1)} 
      else {mtext(paste("Power: Park ", v-1), side = 3, outer=TRUE, adj=0.5, line = 1, cex=1.1)}
      locator(1)
    }
  } else {
    cl <- rainbow(n.species)
    #par(mfcol=c(1,1), mar=c(0.1,0.1,0.1,0.1), oma=c(4,4,4,4), mai=c(0.1,0.1,0.1,2),xpd=NA)
    new.effect[,1] <- c(0.1,0.3,0.5,0.7,0.9)
    #plot(Results[1,1,]~new.effect[,1], type="l",ylim=c(0,1), xlim=c(0,1), lwd=1, main="Power: landscape level", col="grey", ylab="Statistical power", xlab="Effect size")
    for (i in 1:n.species) {
      #new.effect[,i] <- c(0.1,0.3,0.5,0.7,0.9)
      #lines(Results[1,i,]~new.effect[,i], type="l",ylim=c(0,1), xlim=c(0,1), lwd=2, col=cl[i], lty=1)
    }
    #legend("topright", c(as.character(species.list[,1])), inset=c(-0.3,0), lwd=c(2,2,2,2,2,2), cex=0.5, col=cl[1:n.species])
  }
  
  
  cat('\n',"##########################################################") 
  cat('\n',"OUTPUT SUMMARY.....") 
  cat('\n',"##########################################################")  
  cat('\n',"Number of species = ", n.species)
  cat('\n',"Number of parks = ", max(n.park))
  cat('\n',"Fire modelled at sites = ", model.fire)
  cat('\n',"Analysis based on prespecified sites = ", prespecify.sites)
  cat('\n',"All prespecified sites surveyed = ", all.prespecify.sites)
  cat('\n',"Number of sites surveyed = ", n.sites)
  cat('\n',"Ratio of remote to non-remote sites = ", R)
  cat('\n',"Time horizon = ", Tmax, " years")
  cat('\n',"Survey years = ", s.years)
  cat('\n',"Number of repeat visits (method 1) = ", n.method[1])
  cat('\n',"Number of repeat visits (method 2) = ", n.method[2])
  cat('\n',"Number of repeat visits (method 3) = ", n.method[3])
  cat('\n',"Number of repeat visits (method 4) = ", n.method[4])
  cat('\n',"Effect size(s) = ", effect.size) 
  cat('\n',"Two tailed test = ", two.tailed)
  cat('\n',"Type I error rate = ", alpha)
  cat('\n',"Number of simulations = ", nsims*n.cores)
  
  return(Results) 
}  




