#Process ECV Soil Moisture
#Zhen Zhang Mar 4th, 2016
require(raster)
require(ncdf)
# Data span from 1978 to 2014

startYear <- 1978
endYear   <- 2014
years <- seq(startYear,endYear)
naval <- -99999
#inpath <- "/Users/zhang/Data/ECV Soil Moisture/data"
outpath <- "/mnt/lustrefs/work/zhen.zhang/ecv_sm/ecv_sm.nc"
inpath <- "/mnt/lustrefs/work/zhen.zhang/ecv_sm/daily_files/combined"

getMonthlyMax <- function(iyear,nmonth){
  filelist <- dir(paste(inpath,'/',iyear,'/',sep=''),pattern=paste(iyear,sprintf("%02d",nmonth),sep=''),full.names=T)
  tmp_sm <- brick(filelist[1],varname='sm')
  for(i in 2:length(filelist)){
    tmp <- raster(filelist[i],varname='sm')
    tmp_sm <- raster::addLayer(tmp_sm, tmp)
    tmp <- NULL
    tmp_sm[is.na(tmp_sm)] <- -1     
  }
  tmp_max <- calc(tmp_sm,max)   
  return(as.matrix(tmp_max))
}

#Set up output NC file
#Match lon / lat sequences
latseq <- rev(seq(-89.875, 89.875, 0.25))
lonseq <- seq(-179.875, 179.875, 0.25)
months <- seq(as.Date(paste(startYear, "-11-01", sep="")), as.Date(paste(endYear, "-12-01", sep="")), by="months")
nmonths <- 12

#Set up ncdf dims
x <- dim.def.ncdf( "lon", "degrees_east", lonseq)
y <- dim.def.ncdf( "lat", "degrees_north", latseq)
tmonths <- dim.def.ncdf("time", paste("months since ", startYear, "-11-01", sep=""), 0:(length(months)-1))
sm_var <- var.def.ncdf("sm", "m3 m-3", list(x,y,tmonths), naval, longname="Volumetric Soil Moisture", prec="double")
ncnew <- create.ncdf( outpath, sm_var)

imonth <- 1
for(iyear in 1:length(years)){
  if(iyear == 1){
    for(nmonth in 11:nmonths){
      
      varArray <- getMonthlyMax(years[iyear],nmonth)
      put.var.ncdf( ncnew, sm_var, start=c(1,1,imonth), count=c(-1,-1,1) , t(varArray))
      imonth <- imonth +1
      print(imonth)
    }
  }else{
    for(nmonth in 1:nmonths){
      
      varArray <- getMonthlyMax(years[iyear],nmonth)
      put.var.ncdf( ncnew, sm_var, start=c(1,1,imonth), count=c(-1,-1,1) , t(varArray))
      imonth <- imonth + 1 
      print(imonth)
    }
  }
}

att.put.ncdf(ncnew, 0, "project", "MAIOLICA-II")
att.put.ncdf(ncnew, 0, "satellite", "ECV")
att.put.ncdf(ncnew, 0, "contact", "zhen.zhang@wsl.ch")
att.put.ncdf(ncnew, 0, "institution", "WSL&MSU")
close(ncnew)

