#diagnostics for CMIP5 GYE downscaling
#ZZ Sep. 4th 2015
#Check ensemble time series
require(ncdf)
require(field)
require(raster)

rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

ResizeMat <- function(mat, ndim=dim(mat)){
  if(!require(fields)) stop("`fields` required.")
  
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  
  # interpolation
  ans[ncord] <- interp.surface(obj, loc)
  
  ans
}

fdir <- "/mnt/lustrefs/store/poulterlab/Climate/Processed/Regional/GYE/CMIP5/tas/"
fdir.daymet <- "/mnt/lustrefs/store/poulterlab/Climate/Processed/Regional/GYE/Daymet2/"

# rcpList <- c("rcp26","rcp45","rcp60","rcp85")
rcpList <- c("rcp26")
nmonths <- 2868
years <- sort(rep(1861:2099,12))

#Calculate DAYMET
daymet.monthly <- 0
daymet.file <- open.ncdf(paste(fdir.daymet, "tmean_gye_1980-2012_daily_latlon_monthly-gye_shift.nc", sep=""))
for(month in 1:396){
  #Get the data by year
  daymet.data <- get.var.ncdf(daymet.file, "tmean", count=c(-1,-1,1), start=c(1,1,month))
  
  #Calculate monthly global temperature
  daymet.monthly[month] <- mean(daymet.data, na.rm=T)
  
  #Keep track on loop
  print(paste("daymet",month,sep=" "))
}
daymet.annual <- aggregate(daymet.monthly, list(sort(rep(1980:2012,12))), mean)[,2]

#get model name list
filenames <- dir(paste(fdir,rcpList[1],sep = ""))
modelList <- character(0)
for(i in 1:length(filenames)){
  tmp <- strsplit(filenames[i],"_")
  modelList[i] <- tmp[[1]][3]
}
modelList <- unique(modelList)

for(i in 1:length(modelList)){
  modelname <- modelList[i]  
  
  #Read in the original dataset 1st run using ratio approach
  rcp_ori.monthly <- matrix(0,nmonths, length(rcpList))
  for(rcp in rcpList){    
    filepath <- paste(fdir,rcp,"/","tas_Amon_",modelname,"_",rcp,"_r1i1p1_1861-2099-gye.nc",sep="")
    if(!file.exists(filepath)){
      next
    }
    rcp_ori.file <- open.ncdf(filepath)
    
    print(modelname)
    for(month in 1:nmonths){
      #Get the data by year
      rcp_ori.data <- get.var.ncdf(rcp_ori.file, "tas", count=c(-1,-1,1), start=c(1,1,month))
      
      #Remove NaN data
      #rcp_ori.data[rcp_ori.data < 0] <- NA
      
      #Calculate land weights      
      #resample CRU mask to RCP's resolution
      #land.mask <- ResizeMat(cru.mask,c(nlon,nlat))
      
      #rcp_ori.weights <- lat.area/sum(lat.area,na.rm=T)*land.mask
      
      #Calculate monthly mean temperature
      rcp_ori.monthly[month, which(rcp==rcpList)] <- mean(rcp_ori.data, na.rm=T)
      
      #Keep track on loop
      print(paste(rcp,"original", month, sep=" "))
    }
    
  }
  write.table(rcp_ori.monthly,paste(modelname,"r1i1p1","orginal",".txt",sep="_"),row.names=FALSE, col.names=FALSE)
}

