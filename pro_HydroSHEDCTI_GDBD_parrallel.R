#Calculate TOPMODEL parameter based on Toby Marthews's CTI product
#Zhen Zhang, Sep. 7th
library(Hmisc)
library(raster)
library(maptools)
library(ncdf)
library(fBasics)
library(foreach)
library(doParallel)
library(parallel)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
 
wd.dir <- paste(args[1],"/",sep="")
 
rm(args)

#wd.dir <- "/Users/zhang/workspace/CTI_cal/"
setwd(wd.dir)
###########################
cal_fmax <- function(x){
  x.mean <- mean(x)
  #calculate Fmax
  ecd.x <- ecdf(x)
  Fmax <- 1- ecd.x(x.mean)
  return(Fmax)
}  
###########################
cal_f <- function(x){
  x.mean <- mean(x)
  
  #calculate Fmax
  ecd.x <- ecdf(x)
  Fmax <- 1- ecd.x(x.mean)     
  
  #calculate f 
  x.95 <- quantile(x,0.95)  # x.max will cause potential error ,maybe 95% (95 or 98?) value is a good choice 
  x.98 <- quantile(x,0.98)
  z.height <- 2
  f.95 <- (x.95-x.mean)/ z.height
  #f.98 <- (x.98-x.mean)/ z.height
  return(f.95)
}
###########################
calCs <- function(x,...){
  # handling the potential errors
  tryCatch({
#    x <- x[x>0]
    x.mean <- mean(x)
    #calculate Fmax
    ecd.x <- ecdf(x)
    Fmax <- 1- ecd.x(x.mean)     
    
    #Cs can be derived by fitting the exponential function to the CDF of CTI which is larger than mean CTI
    x.ecd <- Ecdf(x, what='1-F', pl=F)
    x.fit <- x.ecd$x[x.ecd$x > x.mean]
    y.fit <- x.ecd$y[x.ecd$x > x.mean]
    data.fit <- data.frame(x.fit,y.fit)
    
    fit <- nls(y.fit ~ Fmax *exp(-Cs*(x.fit - x.mean)), data = data.fit, start = list(Cs = 0.5))      
    fit.result <- summary(fit)
    Cs <- fit.result$parameters[1]
    return(Cs)
  }, error = function(err){
    #print(err)
  }, finally = {
    return(0.5)
  }) #End try catch
}

###################################################################
###################################################################
#calculate mean CTI for each basin using our basin dataset
#read in the CTI dataset at 15 arcsecond
#con_list <- c("af","as","au","ca","eu","na","sa")
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
strt <- Sys.time()
cti <- raster("ga2_Global_15s_geotiff/ga2.tif")
#extract values from cti for shapefile
con_list <- c("af","as","eu","na","oc","sa") #GDMD

foreach(i=1:length(con_list),
                   #.export = c("con_list","cti"),
                   .combine='c',
                   .packages=c("Hmisc","raster","maptools","ncdf","fBasics")) %dopar% {
    basin_shape <- shapefile(paste("Global_Drainage/",con_list[i],"/",con_list[i],"_basins_for_zhang.shp",sep=""))
    basin_shape <- spTransform(basin_shape, proj4string(cti))
    cti_sub <- crop(cti,extent(basin_shape))
    #Test
    # cti_sub <- crop(cti,extent(-5,5,35,40))
    # basin_shape_sub <- crop(basin_shape,extent(-5,5,35,40))
    # mcti <- extract(cti_sub,basin_shape_sub,fun=mean,na.rm=TRUE,df=TRUE,small=TRUE)
    # cs <- extract(cti_sub,basin_shape_sub,fun=calCs,na.rm=TRUE,df=TRUE,small=TRUE)
    #end Test
    mcti <- extract(cti_sub,basin_shape,fun=mean,na.rm=TRUE,df=TRUE,small=TRUE)
    cs <- extract(cti_sub,basin_shape,fun=calCs,na.rm=TRUE,df=TRUE,small=TRUE)
    write.csv(mcti,file=paste("mcti_",con_list[i],".csv",sep=""))
    write.csv(cs,file=paste("cs_",con_list[i],".csv",sep=""))
    rm(cti_sub)
    gc()
}
print(Sys.time() - strt)
stopCluster(cl)

