#Calculate TOPMODEL parameter based on Toby Marthews's CTI product
#Zhen Zhang, Sep. 7th
setwd("/mnt/lustrefs/store/zhen.zhang/data/CTI/")
require(Hmisc)
require(raster)
require(maptools)
require(ncdf)
require(fBasics)
library(doParallel)
require(parallel)
library(foreach)
library(rgdal)
library(sp)
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
cal_cs <- function(x,...){
  x.mean <- mean(x)
  
  #calculate Fmax
  ecd.x <- ecdf(x)
  Fmax <- 1- ecd.x(x.mean)     
  
  #Cs can be derived by fitting the exponential function to the CDF of CTI which is larger than mean CTI
  x.ecd <- Ecdf(x, what='1-F', pl=F)
  x.fit <- x.ecd$x[x.ecd$x > x.mean]
  y.fit <- x.ecd$y[x.ecd$x > x.mean]
  data.fit <- data.frame(x.fit,y.fit)
  # handling the potential errors
  tryCatch({
    fit <- nls(y.fit ~ Fmax *exp(-Cs*(x.fit - x.mean)), data = data.fit, start = list(Cs = 0.5))      
    fit.result <- summary(fit)
    Cs <- fit.result$parameters[1]
  }, error = function(err){
    print(err)
    #Cs <- 0.5      
    return(0.5)
  }, finally = {
    
  }) #End try catch
  return(Cs)
}
###################################################################
#Calculate the number of cores
#no_cores <- detectCores() - 1
#base <- 3
##Initiate cluster
#cl <- makeCluster(no_cores)
#clusterExport(cl,"base")
#clusterEvalQ(cl,library(pryr))
##Run
#test <- parLapply(cl,
#          2:4,
#          function(exponent) base^exponent)
#test <- parSapply(cl,as.character(2:4),
#                  function(exponent){
#                    x <- as.numeric(exponent)
#                    c(base = base^x, self=x^x)
#                  })
#parSapply(cl,x = 1:10, function(x) address(x)) == address(x)
##Finish
#stopCluster(cl)
####################################################################
##Using foreach package for parallel
no_cores <- detectCores()
print(paste("no_cores=",no_cores,sep=""))
cl <- makeCluster(2)
registerDoParallel(cl)
strt <- Sys.time()
cti <- raster("ga2_Global_15s_geotiff/ga2.tif")
#extract values from cti for shapefile
con_list <- c("af","as","eu","na","oc","sa") #GDMD

results <- foreach(i=1:length(con_list),
#.export = c("con_list","cti"),
                   .packages=c("Hmisc","raster","maptools","ncdf",     "fBasics")
                   ) %dopar% {
    basin_shape <- shapefile(paste("Global_Drainage/",con_list[i],"/",con_list[i],"_basins_for_zhang.shp",sep=""))
    basin_shape <- spTransform(basin_shape, proj4string(cti))
    cti_sub <- crop(cti,extent(basin_shape))
    mcti <- extract(cti_sub,basin_shape,fun=mean,na.rm=TRUE,df=TRUE)
    cs <- extract(cti_sub,basin_shape,fun=cal_cs,na.rm=TRUE,df=TRUE)
    write.csv(mcti,file=paste("mcti_",con_list[i],".csv",sep=", "))
    write.csv(cs,file=paste("cs_",con_list[i],".csv",sep=","))
    print(con_list[i])
}
print(Sys.time() - strt)
stopCluster(cl)
###################################################################
#calculate mean CTI for each basin using our basin dataset
#read in the CTI dataset at 15 arcsecond
#cti <- raster("ga2_Global_15s_geotiff/ga2.tif")
##Merge with HYDRO1K CTI dataset to cover the whole globle
#
##extract values from cti for shapefile
#con_list <- c("af","as","eu","na","oc","sa") #GDMD
##con_list <- c("af","as","au","ca","eu","na","sa")
#for(i in 1:length(con_list)){
#  basin_shape <- shapefile(paste("Global_Drainage/",con_list[i],"/",con_list[i],"_basins_for_zhang.shp",sep=""))
#  basin_shape <- spTransform(basin_shape, proj4string(cti))
#  cti_sub <- crop(cti,extent(basin_shape))
#  mcti <- extract(cti_sub,basin_shape,fun=mean,na.rm=TRUE,df=TRUE)
#  cs <- extract(cti_sub,basin_shape,fun=cal_cs,na.rm=TRUE,df=TRUE)
#  write.table(mcti,file=paste("mcti_",con_list[i],".csv",sep=","))
#  write.table(cs,file=paste("cs_",con_list[i],".csv",sep=""))
#  print(con_list[i])
#}
#

#Calculate basin versioned mean CTI based on HYDRO1k basin using HYDRO1k cti for region higher than 60N
# extent_cti <- c(-180,180,60,84)
# basin <- raster("/Users/zhang/workspace/CTI_cal/HYDRO/hydro1k_bas_filled_onedeg.nc",varname="BAS")
# mask <- raster("/Users/zhang/workspace/wetland_frac/input/hydro/cti.mean.grd")
# mask <- crop(mask, extent_cti)
# basin <- crop(basin,extent_cti)
# basin <- resample(basin,mask,method='ngb')
# output_ra <- basin
# 
# cti_new <- crop(cti,extent_cti)
# basin_new <- resample(basin,cti_new,method='ngb')
# cti_zonal <- zonal(cti_new, basin_new,fun='mean')
# 
# list_zonal <- cti_zonal[,2]
# name_zonal <- names(list_zonal)  
# min_basinid <- min(getValues(basin_new),na.rm = T)
# max_basinid <- max(getValues(basin_new),na.rm = T)
# for(j in min_basinid:max_basinid){
#   print(j)
#   if(length(name_zonal[name_zonal==j]) ==0){
#     next;
#   }else{
#     cti_mean <- as.numeric(list_zonal[names(list_zonal)== j])
#     if(is.nan(cti_mean)){
#       output_ra[output_ra==j] <- NA
#     }else{
#       output_ra[output_ra==j] <- cti_mean
#     }
#     
#   }
# }
