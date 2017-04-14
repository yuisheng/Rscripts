#Convert Nonfrozen days dataset from NTSG
#Read in Bin file and convert into Netcdf file

require(rgdal)
require(tis)
require(raster)

dpath <- "/mnt/lustrefs/poulterlab/Benchmark/Original/AMSR_E_2"  #/2002/AMSRU_Mland_2002152A.fw
prjmap <- raster("/mnt/lustrefs/poulterlab/Benchmark/Original/AMSR_E_2/AMSRE_36V_AM_FT_2002_day170.tif")
mkgrid <- function(gridin){
  ease_r <- readBin(paste(dpath,"/sample_reading_code/anci/globland_r",sep=""),integer(),size=2,209091)
  ease_s <- readBin(paste(dpath,"/sample_reading_code/anci/globland_c",sep=""),integer(),size=2,209091)
  gridout <- matrix(NA,586,1383)
  
  for(i in 1:209091){
    gridout[ease_r[i]+1,ease_s[i]+1] <- gridin[i];    
  }
  return(gridout);
}


fname <- "Fw-DailyA-AMSRE.nc"
lon <- 1383
lat <- 586
startYear <- 2002
endYear <- 2014
outNyears <- length(startYear:endYear)
days <- seq(as.Date(paste(startYear, "-06-01", sep="")), as.Date(paste(endYear, "-06-30", sep="")), by="days")
months <- seq(as.Date(paste(startYear, "-06-01", sep="")), as.Date(paste(endYear, "-06-01", sep="")), by="months")
years <- seq(as.Date(paste(startYear, "-01-01", sep="")), as.Date(paste(endYear, "-01-01", sep="")), by="years")
ntstep <- 12
nmonths <- 12
nyear <- endYear - startYear + 1
naval <- -99999


#Match lon / lat sequences
latseq <- rev(seq(-89.75, 89.75, 0.5))
lonseq <- seq(-179.75, 179.75, 0.5)
lonmatch <- match(outLon, lonseq)
latmatch <- match(outLat, latseq)

#Set up ncdf dims
x <- dim.def.ncdf( "lon", "degrees_east", lonseq)
y <- dim.def.ncdf( "lat", "degrees_north", latseq)
tdays <- dim.def.ncdf("time", paste("days since", startYear, "-06-01", sep=""), (0:(length(days)-1)))

ft_var <- var.def.ncdf("fw", "fraction", list(x,y,tdays), naval, longname="Fraction of Water", prec="double")
ncnew <- create.ncdf( fname, ft_var)    
start_i <- 1

for(year in startYear:endYear){
  if(isLeapYear(year)){
    len_day=366
  }else{
    len_day=365
  }
  
  if(year == 2002){
    startday=152
  }else{
    startday=1
  }
  
  if(year == 2014){
    len_day= 181
  }
  for(index in startday:len_day){
    fw_bin <- readBin(paste(dpath,"/2002/AMSRU_Mland_",year,index,"A.fw",sep=""),double(),n=209091)
    fw_bin[fw_bin <= 0.0] <- -1
    fw_array <- mkgrid(gridin = fw_bin)
    ft_tmp <- raster(fw_array,crs=crs(prjmap))
    extent(ft_tmp) <- extent(prjmap)
    
    ra_tmp <- projectRaster(from=ft_tmp, res=0.5,crs=crs("+proj=longlat +ellps=WGS84 +no_def +towgs=0,0,0"),method='ngb')
    ra_tmp <- extend(ra_tmp_a,c(-180,180,-90,90))
    
    varArray <- getValues(ra_tmp, format="matrix")
    varArray[is.na(varArray)] <- naval
    varArray[is.infinite(varArray)] <- naval
    put.var.ncdf( ncnew, ft_var, start=c(1,1,start_i), count=c(-1,-1,1) , t(varArray))
    start_i <- start_i +1
    print(paste(year,index))
  }
}
att.put.ncdf(ncnew, 0, "project", "MAIOLICA-II")
att.put.ncdf(ncnew, 0, "daypass", "ascending")
att.put.ncdf(ncnew, 0, "description", "Daily open water fraction dataset from NTSG")
att.put.ncdf(ncnew, 0, "run_number", "01")
att.put.ncdf(ncnew, 0, "contact", "zhen.zhang@wsl.ch")
att.put.ncdf(ncnew, 0, "institution", "WSL&MSU")
close(ncnew)

