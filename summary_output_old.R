#April 2015
#For summary the output by Transcom regions
# 4 input parameters
#usage: varname rcp caltype isareal
#1. regional CH4 emission Transcom
#2. format: col1 region1, col2 region2 ....
library(sp)
library(raster)

options(echo=TRUE) #if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

varname <- args[1]
rcp <- args[2]
caltype <- args[3]
isareal <- (args[4] > 0)  # 0 or 1
rm(args)

#Calculate area weights
zonedir <- "/mnt/lustrefs/work/zhen.zhang/data/transcom_semiarid_mask.nc"
zones <- raster(zonedir)
cellarea <- raster::area(zones)
transcomreg <- c("Boreal North America","Temperate North America", "Tropical South America", "Temperate South America", "North Africa", "South Africa", "Boreal Eurasia", "Temperate Eurasia", "Tropical Asia", "Australia", "Europe", "North Africa semi-arid", "South Africa semi-arid")

#nmonths <- 139 *12
#years <- sort(rep(1961:2099,12))

#rcpList <- c("rcp26","rcp45","rcp60","rcp85")
for(i in 1:length(rcp)){
  cmip5_path <- paste("/mnt/lustrefs/work/zhen.zhang/output/",rcp[i],"/", sep="")
  
  driving_path <- paste("/mnt/lustrefs/store/poulterlab/Climate/Processed/Global/CMIP5/tas/",rcp[i],"/",sep="")
  
  model_list <- dir(cmip5_path)
  for(j in 1:length(model_list)){
    #fpath.lpj <- paste(cmip5_path,model_list[j],"/merge/LPJ_",varname,".nc", sep="")    
    fpath.lpj <- paste(driving_path,varname,"_Amon_",model_list[j],"_",rcp[i],"_1861-2099_r1i1p1_global-cor.nc", sep="")    
    
    varnc <- brick(fpath.lpj)
    if(isareal){
      varnc <- varnc * cellarea  
    }
    varnc_zonal <- zonal(varnc, zones,fun=caltype)
    
    #write output file
    #check if stats directory exist, if not, create it
    outpath <- paste(cmip5_path,model_list[j],"/stats",sep='')
    if(!file.exists(outpath)){
      dir.create(outpath)
    }
    foutpath <- paste(outpath,"/",varname,"_",caltype,".txt",sep='')
    write(varnc_zonal, foutpath, ncolumns=length(transcomreg), sep=",")
    #read method: read.table("sample.txt",header=T,sep=",")
    
    print(model_list[j])
  }
}
