#December 2014
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

fdir <- "/mnt/lustrefs/store/poulterlab/Climate/Processed/Regional/GYE/CMIP5/pr/"
fdir.daymet <- "/mnt/lustrefs/store/poulterlab/Climate/Processed/Regional/GYE/Daymet2/"

#rcpList <- c("rcp26","rcp45","rcp60","rcp85")
rcpList <- c("rcp26")
nmonths <- 2868
years <- sort(rep(1861:2099,12))

#Calculate DAYMET
#daymet.monthly <- 0
#daymet.file <- open.ncdf(paste(fdir.daymet, "tmean_gye_1980-2012_daily_latlon_monthly-gye_shift.nc", sep=""))
#for(month in 1:396){
  #Get the data by year
#  daymet.data <- get.var.ncdf(daymet.file, "tmean", count=c(-1,-1,1), start=c(1,1,month))
  
  #Calculate monthly global temperature
#  daymet.monthly[month] <- mean(daymet.data, na.rm=T)
  
  #Keep track on loop
#  print(paste("daymet",month,sep=" "))
#}
#daymet.annual <- aggregate(daymet.monthly, list(sort(rep(1980:2012,12))), mean)[,2]

#get model name list
filenames <- dir(paste(fdir,rcpList[1],"-wsl",sep = ""))
modelList <- character(0)
for(i in 1:length(filenames)){
  tmp <- strsplit(filenames[i],"_")
  modelList[i] <- tmp[[1]][3]
}
modelList <- unique(modelList)

for(i in 1:length(modelList)){
  modelname <- modelList[i]  
  #Read and calculate the summary for origional, wsl, ratio respectively 
  
  #################################################
  
  #Read in the RCP r1i1p1 using WSL approach
  rcp_wsl.monthly <- matrix(0,nmonths, length(rcpList))
  for(rcp in rcpList){
    filepath <- paste(fdir,rcp,"-wsl/","pr_Amon_",modelname,"_",rcp,"_1861-2099_r1i1p1_global-cor.nc",sep="")
    if(!file.exists(filepath)){
      next
    }
    rcp_wsl.file <- open.ncdf(filepath)
    
    for(month in 1:nmonths){
      #Get the data by year
      rcp_wsl.data <- get.var.ncdf(rcp_wsl.file, "pr", count=c(-1,-1,1), start=c(1,1,month))
      
      #Remove NaN data
      rcp_wsl.data[rcp_wsl.data < 0] <- NA
      
      #Calculate land weights
      #                rcp.mask <- rcp.data
      #               rcp.mask[!is.na(rcp.mask)] <- 1
      
      #Calculate monthly global temperature (weight by area)
      rcp_wsl.monthly[month, which(rcp==rcpList)] <- mean(rcp_wsl.data, na.rm=T)
      
      #Keep track on loop
      print(paste(rcp,"wsl", month, sep=" "))
    }
    
  }
  write.table(rcp_wsl.monthly,paste("pr",modelname,"r1i1p1","gye",".txt",sep="_"),row.names=FALSE, col.names=FALSE)
  
}

