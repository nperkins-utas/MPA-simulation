#Please see the readme file for more details
#Set the working directory to the correct location

#Clear the workspace
rm(list = ls())

#Note that the below code is given as an example only, with one MPA site and reference pair (Bodega Head) used as an example
#Code will need to be modified for specific applications elsewhere.

#Load in required packages
require(spatstat)
require(raster)
require(maptools)
require(fields)
require(rgdal)
require(rgeos)
require(RandomFields)

###################################################################################################
#READ IN DATA AND MODIFY TO FORMAT REQUIRED

#Read in the example fish location data from Bodega Head
#Note that site "BB1" is the MPA site and site "BB4" is the reference site
bbf <- read.csv("BB_Fish.csv", header = TRUE, stringsAsFactors = FALSE)

#Read in the bathymetry covariate data for Bodega Head
bbc <- read.csv("BB_Covs.csv", header = TRUE, stringsAsFactors = FALSE)

#Read in output of the IPM.
#Note - Please see IPM model for details
load("brf_bb_ipm_output.Rdata")
#Extract the data for the MPA
brf.mpa <- brf.bb.theoretical[[1]]
#Extract the data for the reference site
brf.ref <- brf.bb.theoretical[[2]]
#Extract the abundance ratio change for the MPA site. This represents tha abundance change relative to the starting abundance for each year
#in the simulation
brf.mpa.N <- brf.mpa$N
#Extract the MPA size structure. This is a matrix where rows represent the size class (in cm)
#of fish, columns represent years, and cell entries are the proportion of fish that fall into 
#the relevant size class for that year
brf.mpa.size <- brf.mpa$str
#Extract the reference site size structure
brf.ref.size <- brf.ref$str

#Recalibrate size structure to fish lengths >= the observable size classes (20 cm for Brown rf)
colnames(brf.mpa.size) <- paste("year", 1:25, sep = "_")
rownames(brf.mpa.size) <- 1:61
brf.mpa.size <- brf.mpa.size[-c(1:19), ] 
brf.mpa.size <- sweep(brf.mpa.size, 2, colSums(brf.mpa.size), FUN = "/")
rownames(brf.ref.size) <- 1:61
brf.ref.size <- brf.ref.size[-c(1:19), ] 
brf.ref.size <- brf.ref.size/sum(brf.ref.size)


#Define sample locations (ROV survey sites) from shapefiles
rov.shp <- readOGR(dsn = ".", layer = "BBsites") #read in shapefile of sites as spatial points polygon
rov.shp@data$SiteID #site names - not order (maintained below)
#Extract shapefiles of site boundaries and convert to owin objects
rp <- as(rov.shp, "SpatialPolygons")
rov.sites <- slot(rp, "polygons")
rov.sites <- lapply(rov.sites,  function(x) { SpatialPolygons(list(x)) })
rov.wins <- lapply(rov.sites, as.owin)
bb.4 <- rov.wins[[1]] #REF
bb.1 <- rov.wins[[3]] #SMR

#Read in parameters from spatial model
#This reads in a list that contains beta covariate estimates, spatial range parameter estimates and spatial variance estimates
all.pars <- readRDS("spatial_pars")


#Add a depth^2 column to the covariate file
bbc$depth <- bbc$depth*-1
bbc$Depth2 <- (bbc$depth)^2


#Scale covariates by centering on the mean and dividing by the standard deviation using the "scale" function
bbc[, c(4:12, 14:19, 21)] <- scale(bbc[, c(4:12, 14:19, 21)])
bbc$SF_char_Indicator_of_mixed <- rep(0, nrow(bbc))
bbc$SF_char_Indicator_of_mixed[which(bbc$sf_character == "mixed")] <- 1
bbc[, "SF_char_Indicator_of_mixed"] <- scale(bbc[, "SF_char_Indicator_of_mixed"])
bbc$SF_char_Indicator_of_soft <- rep(0, nrow(bbc))
bbc$SF_char_Indicator_of_soft[which(bbc$sf_character == "soft")] <- 1
bbc[, "SF_char_Indicator_of_soft"] <- scale(bbc[, "SF_char_Indicator_of_soft"])

#Create raster stack of covariates for each site for the simulation
#Note: these only contain "significant" covariates from the spatial modelling
###MPA SITE - BB1###
spg.bb1 <- bbc[which(bbc$site == "BB1"), c("x", "y", "Depth2", "bpi_2_10", "bpi_3_15", "bpi_5_25", "bpi_10_50", "curvature", "slope", "vrm_5",
                                           "vrm_15", "sf_character", "eastness", "northness", "dist_hard", "rdmv")]
spg.bb1$sf_character <- as.factor(spg.bb1$sf_character)
coordinates(spg.bb1) <- ~ x + y
gridded(spg.bb1) <- TRUE
spg.bb1.ras <- stack(spg.bb1)

###REFERENCE SITE - BB4###
spg.bb4 <- bbc[which(bbc$site == "BB4"), c("x", "y", "Depth2", "bpi_2_10", "bpi_3_15", "bpi_5_25", "bpi_10_50", "curvature", "slope", "vrm_5",
                                           "vrm_15", "sf_character", "eastness", "northness", "dist_hard", "rdmv")]
spg.bb4$sf_character <- as.factor(spg.bb4$sf_character)
coordinates(spg.bb4) <- ~ x + y
gridded(spg.bb4) <- TRUE
spg.bb4.ras <- stack(spg.bb4)

###FUNCTIONS USED IN THE SIMULATIONS####

##########################################################################################################
img.func <- function(pars){
  #This function takes the site and covariate beta estimates and produces an raseterized image (im for spatstat)
  #that is a linear combination of covariates*betas
  #time.fac provides a factor to modify the abundance by. This is 1 at MPA implementation and then increases from then
  
  #Extract site level covariates
  mpa.site.cov <- bbc[bbc$site == "BB1", ]
  ref.site.cov <- bbc[bbc$site == "BB4", ]
  
  #Extract covariate data from list based on mc sample
  cov.pars <- pars[[1]]
  
  
  #Extract intercpet term - this represent year 5
  int.start <- cov.pars["Intercept"]  
  #Calculate a vector of intercepts for the 25 year time-series
  int.all.mpa <- log(exp(int.start)*brf.mpa.N)
  int.all.ref <- rep((int.start + pars[[1]]['BB4']), 25)
  
  
  #Assign the intercept and environmental covariate coefs
  beta1 <- cov.pars["Depth2"] #depth^2 coefficient
  beta2 <- cov.pars["bpi_2_10"] #bpi_2_10 coefficient
  beta3 <- cov.pars["bpi_3_15"] #bp1 3_15 coefficient
  beta4 <- cov.pars["bpi_5_25"] #bpi 5_25 coefficient
  beta5 <- cov.pars["bpi_10_50"] #bpi 10_50 coefficient
  beta6 <- cov.pars["curvature"] #curvature coefficient
  beta7 <- cov.pars["slope"] #slope coefficient
  beta8 <- cov.pars["vrm_5"] #vrm_5 coefficient
  beta9 <- cov.pars["vrm_15"] #vrm_15 coefficient
  beta10 <- cov.pars["eastness"] #eastness coefficient
  beta11 <- cov.pars["northness"] #northness coefficient
  beta12 <- cov.pars["dist_hard"] #dist_hard coefficient
  beta13 <- cov.pars["rdmv"] #rdmv coefficient
  beta14 <- cov.pars["SF_char_Indicator_of_mixed"] #mixed coefficient
  beta15 <- cov.pars["SF_char_Indicator_of_soft"] #soft coefficient
  beta16 <- cov.pars["BB1"] #BB1 coefficient
  beta17 <- cov.pars["BB4"] #BB4 coefficient
  
  #Create image that has cells that are the linear combination
  mpa.site.cov[, "Depth2"] <- mpa.site.cov[, "Depth2"]*beta1
  mpa.site.cov[, "bpi_2_10"] <- mpa.site.cov[, "bpi_2_10"]*beta2
  mpa.site.cov[, "bpi_3_15"] <- mpa.site.cov[, "bpi_3_15"]*beta3
  mpa.site.cov[, "bpi_5_25"] <- mpa.site.cov[, "bpi_5_25"]*beta4
  mpa.site.cov[, "bpi_10_50"] <- mpa.site.cov[, "bpi_10_50"]*beta5
  mpa.site.cov[, "curvature"] <- mpa.site.cov[, "curvature"]*beta6
  mpa.site.cov[, "slope"] <- mpa.site.cov[, "slope"]*beta7
  mpa.site.cov[, "vrm_5"] <- mpa.site.cov[, "vrm_5"]*beta8
  mpa.site.cov[, "vrm_15"] <- mpa.site.cov[, "vrm_15"]*beta9
  mpa.site.cov[, "eastness"] <- mpa.site.cov[, "eastness"]*beta10
  mpa.site.cov[, "northness"] <- mpa.site.cov[, "northness"]*beta11
  mpa.site.cov[, "dist_hard"] <- mpa.site.cov[, "dist_hard"]*beta12
  mpa.site.cov[, "rdmv"] <- mpa.site.cov[, "rdmv"]*beta13
  mpa.site.cov[, "SF_char_Indicator_of_mixed"] <- mpa.site.cov[, "SF_char_Indicator_of_mixed"]*beta14
  mpa.site.cov[, "SF_char_Indicator_of_soft"] <- mpa.site.cov[, "SF_char_Indicator_of_soft"]*beta15
  
  ref.site.cov[, "Depth2"] <- ref.site.cov[, "Depth2"]*beta1
  ref.site.cov[, "bpi_2_10"] <- ref.site.cov[, "bpi_2_10"]*beta2
  ref.site.cov[, "bpi_3_15"] <- ref.site.cov[, "bpi_3_15"]*beta3
  ref.site.cov[, "bpi_5_25"] <- ref.site.cov[, "bpi_5_25"]*beta4
  ref.site.cov[, "bpi_10_50"] <- ref.site.cov[, "bpi_10_50"]*beta5
  ref.site.cov[, "curvature"] <- ref.site.cov[, "curvature"]*beta6
  ref.site.cov[, "slope"] <- ref.site.cov[, "slope"]*beta7
  ref.site.cov[, "vrm_5"] <- ref.site.cov[, "vrm_5"]*beta8
  ref.site.cov[, "vrm_15"] <- ref.site.cov[, "vrm_15"]*beta9
  ref.site.cov[, "eastness"] <- ref.site.cov[, "eastness"]*beta10
  ref.site.cov[, "northness"] <- ref.site.cov[, "northness"]*beta11
  ref.site.cov[, "dist_hard"] <- ref.site.cov[, "dist_hard"]*beta12
  ref.site.cov[, "rdmv"] <- ref.site.cov[, "rdmv"]*beta13
  ref.site.cov[, "SF_char_Indicator_of_mixed"] <- ref.site.cov[, "SF_char_Indicator_of_mixed"]*beta14
  ref.site.cov[, "SF_char_Indicator_of_soft"] <- ref.site.cov[, "SF_char_Indicator_of_soft"]*beta15
  
  
  #Create a list of "images" that are based on a linear combination of covariates with
  #a changing intercept through time
  l.im.mpa <- lapply(int.all.mpa, function(x){ 
      beta0 <- x - log(4) #intercept adjusted to 4m^2 grid cell size offset
      
      lin.df.mpa <- data.frame(mpa.site.cov[,c("x", "y")], (mpa.site.cov[,"Depth2"] + mpa.site.cov[, "bpi_2_10"] + 
                                                    mpa.site.cov[, "bpi_3_15"] +mpa.site.cov[, "bpi_5_25"] + mpa.site.cov[, "bpi_10_50"] +
                                                    mpa.site.cov[, "curvature"] + mpa.site.cov[, "slope"] + mpa.site.cov[, "vrm_5"] +
                                                    mpa.site.cov[, "vrm_15"] + mpa.site.cov[, "eastness"] + mpa.site.cov[, "northness"] + 
                                                    mpa.site.cov[, "dist_hard"] + mpa.site.cov[, "rdmv"] + 
                                                    mpa.site.cov[, "SF_char_Indicator_of_soft"] + mpa.site.cov[, "SF_char_Indicator_of_mixed"] +
                                                    beta0 + beta16))
    im.mpa <- as.im(lin.df.mpa)
    return(im.mpa)
    })
  
  l.im.ref <- lapply(int.all.ref, function(x){ 
    beta0 <- x - log(4) #intercept adjusted to 4m^2 grid cell size offset
    
    lin.df.ref <- data.frame(ref.site.cov[,c("x", "y")], (ref.site.cov[,"Depth2"] + ref.site.cov[, "bpi_2_10"] + 
                                                            ref.site.cov[, "bpi_3_15"] +ref.site.cov[, "bpi_5_25"] + ref.site.cov[, "bpi_10_50"] +
                                                            ref.site.cov[, "curvature"] + ref.site.cov[, "slope"] + ref.site.cov[, "vrm_5"] +
                                                            ref.site.cov[, "vrm_15"] + ref.site.cov[, "eastness"] + ref.site.cov[, "northness"] + 
                                                            ref.site.cov[, "dist_hard"] + ref.site.cov[, "rdmv"] + 
                                                            ref.site.cov[, "SF_char_Indicator_of_soft"] + ref.site.cov[, "SF_char_Indicator_of_mixed"] +
                                                            beta0 + beta17))
    im.ref <- as.im(lin.df.ref)
    return(im.ref)
  }) 
  all.imgs <- list(l.im.ref, l.im.mpa)
  return(all.imgs)
  
}

######################################################################################################
# Function that takes a point pattern timeline and assigns marks 

mark.func <- function (fish.xy, time.pt, mpa.fac) {
  #Extract total number of fish for the time point
  fish.tot <- nrow(fish.xy)
  
  #Divide total observed fish into size classes appropriate for the time period
  if(mpa.fac){
  size.propns <- brf.mpa.size[, time.pt]#extract appropriate row
  fish.bins <- round(size.propns*fish.tot, 0)
  diff.tots <- sum(fish.bins) - fish.tot #difference between total number of points and size bin allocations
  #Subtract or add fish to the smallest size class (i.e. 1) so that totals match
  fish.bins[1] <- fish.bins[1] - diff.tots} else {
    size.propns <- brf.ref.size
    fish.bins <- round(size.propns*fish.tot, 0)
    diff.tots <- sum(fish.bins) - fish.tot #difference between total number of points and size bin allocations
    #Subtract or add fish to the smallest size class (i.e. 1) so that totals match
    fish.bins[1] <- fish.bins[1] - diff.tots
    
  }
  
  #Create a new matrix that has sizes and counts
  fish.length.dist <- cbind(as.numeric(names(fish.bins)), fish.bins)
  #Create a vector of lengths
  fish.lengths <- unlist(apply(fish.length.dist, 1, function (x)  rep(x[1], x[2])))
  
  #Random assignment of lengths to observations
  samp.ord <- sample(1:fish.tot, fish.tot, replace = FALSE)
  fish.out <- as.data.frame(cbind(fish.xy, fish.lengths[samp.ord]))
  colnames(fish.out) <- c("x", "y", "length")
  return (fish.out)
}



#############################################################################################################################
####LGCP FUNCTION###
#This function takes the linear combination "image" created in the img.func and the spatial paramaters (phi and sigma^2 estimates)
#and creates a random realisation of an LGCP point process for the given image

lgcp.func <- function(sim.imgs, all.pars){
  
  #Create time-series of images
  ref.imgs <- sim.imgs[[1]]
  mpa.imgs <- sim.imgs[[2]]
  
  #Extract corresponding spatial parameters
  sig.par.mpa <- all.pars$sig.pars["BB1"]
  sig.par.ref <- all.pars$sig.pars["BB4"]
  phi.par.mpa <- all.pars$phi.pars["BB1"]
  phi.par.ref <- all.pars$phi.pars["BB4"]
  
  #Run rLGCP to create time-series data for the MPA
  mpa.xy.fish <- sapply(mpa.imgs, function (x){
    lg.s <- rLGCP('exp', x,
                var=sig.par.mpa, scale=phi.par.mpa, win = as.owin(x))
    xy <- cbind(lg.s$x, lg.s$y)
    return (xy)
  })
  
  
  
  ref.xy.fish <- sapply(ref.imgs, function (x){
    lg.s <- rLGCP('exp', x,
                  var=sig.par.ref, scale=phi.par.ref, win = as.owin(x))
    xy <- cbind(lg.s$x, lg.s$y)
    return (xy)
  }) 
  
  
  #Add marks
  mpa.output <- lapply(1:25, function (n){
    mark.func(mpa.xy.fish[[n]], mpa.fac = TRUE, time.pt = n)
  } )
  
  ref.output <- lapply(1:25, function (n){
    mark.func(ref.xy.fish[[n]], mpa.fac = FALSE, time.pt = n)
  } )

  
  output <- list(mpa.output, ref.output)
  return(output)
  
}

##########################################################################################################
###TRANSECT FUNCTION###
#Function to simulate a given number of transect lines over a given site
#Random starting point depends on number of transects across area
#Starting point is picked randomly in first section and then spaced systematically
#sample.site needs to be the owin object read in from the shape file (e.g. bb.1, bb.4)
trans.func <- function (num.trans, sample.site) 
  
{
  #Function to simulate a given number of transect lines over a given site
  #Random starting point depends on number of transects across area
  #Starting point is picked randomly in first section and then spaced systematically
  #sample.site needs to be the owin object read in from the shape file (e.g. bb.1, bb.4)
  
  #First find which boundary points define the bottom and top apexes
  bdrys <- cbind(sample.site$bdry[[1]]$x, sample.site$bdry[[1]]$y) #boundary points in a matrix
  cent <- unlist(centroid.owin(sample.site)) #centre point of site
  
  se.1 <- which(bdrys[,1] > cent[1])
  se.1 <- bdrys[se.1,] #rhs short edge
  se.2 <- which(bdrys[,1] < cent[1])
  se.2 <- bdrys[se.2,] #lhs short edge
  
  tops <- rbind(se.1[which.max(se.1[,2]),], se.2[which.max(se.2[,2]),])
  bots <- rbind(se.1[which.min(se.1[,2]),], se.2[which.min(se.2[,2]),])
  
  #Calculate site width for bottom of site
  sw <- sqrt((tops[1,1] - tops[2,1])^2 + (tops[1,2] - tops[2,2])^2) #width of site
  
  #Fix size.secs to be along bottom boundary
  size.secs <- sw/num.trans #size of each section 
  dist.rand <- runif(1, 0, size.secs)
  xy1 <- c(((1 - dist.rand/sw)*bots[1,1] + dist.rand/sw*bots[2,1]), ((1 - dist.rand/sw)*bots[1,2] + dist.rand/sw*bots[2,2]))
  #Points along a line equation taken from: https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point
  xystarts <- sapply(seq(2:num.trans), function (x) c(((1 - size.secs*x/(sw-dist.rand))*xy1[1] + size.secs*x/(sw-dist.rand)*bots[2,1]), ((1 - size.secs*x/(sw-dist.rand))*xy1[2] + size.secs*x/(sw-dist.rand)*bots[2,2]))) 
  xystarts <- cbind(xy1, xystarts)
  xystarts <- t(xystarts)
  
  #Fix size.secs to be along top boundary
  xy2 <- c(((1 - dist.rand/sw)*tops[1,1] + dist.rand/sw*tops[2,1]), ((1 - dist.rand/sw)*tops[1,2] + dist.rand/sw*tops[2,2]))
  #Points along a line equation taken from: https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point
  xyends <- sapply(seq(2:num.trans), function (x) c(((1 - size.secs*x/(sw-dist.rand))*xy2[1] + size.secs*x/(sw-dist.rand)*tops[2,1]), ((1 - size.secs*x/(sw-dist.rand))*xy2[2] + size.secs*x/(sw-dist.rand)*tops[2,2]))) 
  xyends <- cbind(xy2, xyends)
  xyends <- t(xyends)
  
  
  lines <- list()
  
  for(i in 1:nrow(xystarts)){
    #sequence of y values moving 10 cm forward at a time
    rw <- sapply(seq(1:5001), function (x) c(((1 - 0.1*x/500)*xystarts[i,1] + 0.1*x/500*xyends[i,1]), ((1 - 0.1*x/500)*xystarts[i,2] + 0.1*x/500*xyends[i,2]))) 
    rw <- t(rw)
    
    lines[[i]] <- rw
  }
  
  lines2 <- lapply(lines, spLines)
  
  return(lines2)
  
}

#####################################################################################################################
###EXTRACT FUNCTION###
# Function that takes simulated transect lines at a site and returns a list containing the "observed" fish locations
# and the associated covariates as well as the associated fish lengths

extract.func <- function (t.lines, lgcp.out, mpa.fac){
  spg.ras <- if (mpa.fac) {spg.bb1.ras} else {spg.bb4.ras}
  out <- lapply(t.lines, function (x) extract(spg.ras, x, cellnumbers = TRUE))
  cov.df <- unlist(out, recursive = FALSE)
  cov.df <- as.data.frame(do.call(rbind, cov.df))
  cov.df <- na.omit(cov.df)
  
  fish.cells <- extract(spg.ras, lgcp.out[, c("x", "y")], cellnumbers = TRUE)
  fish.cells <- fish.cells[, 'cells']
  
  m <- match(fish.cells, cov.df$cell)
  m <- m[which(!is.na(m))]
  
  cov.df$count <-  rep(0, nrow(cov.df))
  
  for(i in 1:length(m)){
    x <- m[i]
    cov.df$count[x] <- cov.df$count[x] + 1
  }
  
  #Match length data to observations
  lgcp.out$cell <- fish.cells
  lengths.out <- lgcp.out$length[!is.na(match(fish.cells, cov.df$cell))]
  
  final.out <- list(cov.df, lengths.out)
  return(final.out)
  
}




###SIMULATION###
#Below is code to conduct the simulations using the above functions
#Note that the simulations are computationally intensive and will require some time to run
#If possible it is recommended that the simulations be run using parallel processing
#The loop below will save the results of the simulations in separate R data files in the working directory
#for the MPA and reference site examples, with the simulation number noted

B <- 10 #number of simulations - 10 used as an example

#Specify the number of transects
no.trans <- c(8, 12, 16, 20)

#Generate the time-series of "images" that are linear combinations of covariates for each site across which simulations will be conducted
all.imgs <- img.func(pars = all.pars)

#Loop for simulations
#
for (bb in 1:B){
  #Simulate distributions for MPA and REF
  all.dists <- lgcp.func(sim.imgs = all.imgs, all.pars = all.pars)
  
  mpa.res <- lapply(1:25, function(tt){
    mpa.f <- all.dists[[1]][[tt]]
    mpa.trans <- sapply(no.trans, trans.func, sample.site = bb.1)
    mpa.out <- lapply(mpa.trans, extract.func, lgcp.out = mpa.f, mpa.fac = TRUE)
    names(mpa.out) <- c("8_trans", "12_trans", "16_trans", "20_trans")
    return(mpa.out)
  })
  
  for (tt in 1:25){
  
    mpa.f <- all.dists[[1]][[tt]]
    ref.f <- all.dists[[2]][[tt]]
    
    #Simulate transects
    mpa.trans <- sapply(no.trans, trans.func, sample.site = bb.1)
    ref.trans <- sapply(no.trans, trans.func, sample.site = bb.4)
    
    #Get results
    mpa.out <- lapply(mpa.trans, extract.func, lgcp.out = mpa.f, mpa.fac = TRUE)
    names(mpa.out) <- c("8_trans", "12_trans", "16_trans", "20_trans")
    ref.out <- lapply(ref.trans, extract.func, lgcp.out = ref.f, ref.fac = TRUE)
    names(ref.out) <- c("8_trans", "12_trans", "16_trans", "20_trans")
    
    
  }
  save(mpa.out, file = paste(paste("MPA_output_brf_BB_time", tt, "sim", bb, sep = "_"), "RData", sep = "."))
  save(ref.out, file = paste(paste("REF_output_brf_BB_time", tt, "sim", bb, sep = "_"), "RData", sep = "."))
}
