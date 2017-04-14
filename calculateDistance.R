library(raster)
library(sp)
library(foreach)
library(parallel)
library(doParallel)
library(SDMTools)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

wd.dir <- paste(args[1],"/",sep="")
rm(args)

#wd.dir <- "/Users/zhang/Documents/Teaching/ETHZ/data/"
setwd(wd.dir)

current.filename <- paste(wd.dir,"current_climate/","currentclimate_1971-2000.grd",sep="")
future.files <- list.files(paste(wd.dir,"future_climate",sep=""),pattern="*.grd",recursive = T, full=T)
dummy1 <- strsplit(future.files, "/")
len.dummy <- length(dummy1[[1]])
mat  <- matrix(unlist(dummy1), ncol=len.dummy, byrow=TRUE)
future.filenames <- mat[,len.dummy]

#t=0.1 for Temperature, t= 5 for Precipitation
calSpatialGradient <- function(present,future,t){
  #####
  #Crop a subregion
  #present <- crop(present,extent(0,1500000,-480000,600000))
  #future <- crop(future,extent(0,1500000,-480000,600000))
  
  #mask out NA pixels in both raster
  sp.mask <- is.na(present) | is.na(future)
  sp.mask[sp.mask==1] <- NA
  future <- raster::mask(future,sp.mask)
  present <- raster::mask(present,sp.mask)
  
  #set threshold
  t <- 1/(t*2)
  
  current.xyz <- rasterToPoints(present)
  future.xyz  <- rasterToPoints(future)
  
  x <- current.xyz[,1]
  y <- current.xyz[,2]
  
  p <- round(current.xyz[,3]*t)/t     # vector of rounded present climate values
  f <- round(future.xyz[,3]*t)/t      # vector of rounded future climate values
  d <- vector(length = length(p))     # empty vector to write distance to climate
  
  u     <- base::unique(p)[order(base::unique(p))]    # list of unique climate values in p
  match <- function(u){c(which(u==f))}    # function finding climate matches of u with f
  m     <- sapply(u, match)               # list of climate matches for unique values
  
  for(i in 1:length(p)){                  # loop for all grid cells of p
      mi   <- m[[which(u==p[i])]]          # recalls list of climate matches for p
      d[i] <- sqrt(min((x[i]-x[mi])^2 + (y[i]-y[mi])^2))    # distance to closest match
    
  }
  
  d <- d/1000 #to km
  d[d==Inf] <- 100000  # sets no analogue to 10,000km
  
  sp.ra <- rasterFromXYZ(cbind(x,y,d))
  crs(sp.ra) <- crs(present)
  return(sp.ra)
}

no_cores <- detectCores()
cl <- makeCluster(no_cores/2)
registerDoParallel(cl)
strt <- Sys.time()

foreach(i=1:length(future.filenames),
        .combine='c',
        .packages=c("raster","SDMTools")) %dopar% {
          current.grd <- stack(current.filename)
          current.grd.tas <- current.grd$tas/100 
          current.grd.prcp <- current.grd$prcp
          
          future.grd <- stack(future.files[i])
          future.grd.tas <- future.grd$tas
          future.grd.prcp <- future.grd$prcp*10  #convert cm to mm
          
          output.grd.tas  <- calSpatialGradient(current.grd.tas,future.grd.tas,0.25)
          output.grd.prcp <- calSpatialGradient(current.grd.prcp,future.grd.prcp,5)
          writeRaster(output.grd.tas,paste("tas_",future.filenames[i],sep=""))
          writeRaster(output.grd.prcp,paste("prcp_",future.filenames[i],sep=""))
          
}

print(Sys.time() - strt)
stopCluster(cl)


