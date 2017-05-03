#Read in SMAP Fw binary file and convert it to NCDF
library(raster)
library(sp)

inpath <- "/mnt/lustrefs/work/zhen.zhang/data/SMAP_FW/"
ease_lats <- readBin(paste(inpath,"lats.bin",sep=""),double(),size=4,3856*1624)
ease_lons <- readBin(paste(inpath,"lons.bin",sep=""),double(),size=4,3856*1624)

filespath <- list.files(paste(inpath,sep=""),pattern='SMAP')
ra.mask <- raster("/mnt/lustrefs/work/zhen.zhang/data/mask/global_mask.nc")
ra.mask[ra.mask > 0 ] <- 0
ra.mask <- disaggregate(ra.mask,5)

ra.output <- stack(ra.mask)
output.ra <- raster(nrow=1800,ncols=3600,res=0.1)
for(i in 1:length(filespath)){
  input <- data.frame(readBin(paste(inpath,filespath[i],sep=""),double(),size=4,3856*1624))
  coordinates(input) <- cbind(ease_lons,ease_lats)

  sample.ra <- raster(nrow=1800,ncols=3600,res=0.1)
  tmp1 <- rasterize(input,sample.ra)
  tmp2 <- tmp1[[2]]
  tmp2[tmp2 <0|is.na(tmp2)] <- 0
  
  test <- mask(tmp2,ra.mask)
  ra.output <- addLayer(ra.output,test)
}

dropLayer(ra.output,1)
writeRaster(ra.output,paste(inpath,"SMAP_Fw.nc",sep=""),format="CDF",varname='Fw')

