#Oct 2014
#Assemble LPJ netcdfs from monthly LPJ variables
#Submit as: qsub -q long runR.sh
#Load libraries
print(date())
library(sp)
library(ncdf)
library(hash)

#Set parameters
fpath.lpj <- "/mnt/lustrefs/store/zhen.zhang/output/CRU/CRU_2015_USDA_MLIT_PERMAFROST/merge/"
startYear <- 1901
endYear <- 2015

#Fixed parameters
conversionFactor <- 1
outNyears <- length(startYear:endYear)
ntstep <- 12
nmonths <- 12
nyear <- endYear - startYear + 1
naval <- 9.99999e14
ncell.tot <- file.info(paste(fpath.lpj,'grid.out',sep=''))$size/4

dayofmonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)

#Parameters
#For CMIP5
# startYear <- 1961
# endYear <- 2099
# outNyears <- length(startYear:endYear)
# ntstep <- 12
# nmonths <- 12
# nyear <- endYear - startYear + 1
# naval <- -99999
# ncell.tot <- 62482

#For CRU

#Set up the LPJ output read function
readFile <- function(infile, cell.tot, nmonths, year){
  fileName <- file(infile, 'rb')
  seek(fileName, 4*cell.tot*nmonths*(year-1))
  binOut <- readBin(fileName, double(), size=4, n = cell.tot*nmonths)
  close(fileName)
  return(matrix(binOut, ncell.tot, nmonths))
}


#Read outgrid
outGrid <- readBin(paste(fpath.lpj,"grid.out", sep=''), integer(), size=2, n=ncell.tot*2)/100
outLon <- outGrid[seq(1,ncell.tot*2,2)]
outLat <- outGrid[seq(2,ncell.tot*2,2)]

#Match lon / lat sequences
latseq <- seq(-89.75, 89.75, 0.5)
lonseq <- seq(-179.75, 179.75, 0.5)
lonmatch <- match(outLon, lonseq)
latmatch <- match(outLat, latseq)

#Set up ncdf dims
x <- dim.def.ncdf( "lon", "degrees_east", lonseq)
y <- dim.def.ncdf( "lat", "degrees_north", latseq)
nsoil <- dim.def.ncdf( "soil", paste("soil id ", 1:2, sep=""), 0:1)
npft <- dim.def.ncdf( "pft", paste("pft id ", 1:9, sep=""), 0:8)
tmonth <- dim.def.ncdf( "time", paste("months since ", startYear, "-01-01", sep=""), 0:(outNyears*nmonths-1))
tyear <- dim.def.ncdf( "time", paste("years since ", startYear, "-01-01", sep=""), 0:(outNyears-1))


#Make ncdf
varNC <- var.def.ncdf("emco2nep", "kg C m^-2 s^-1", list(x,y,tmonth), naval, "Net Ecosystem Exchange", prec="double")
ncnew <- create.ncdf( paste(fpath.lpj, "LPJ_NEE.nc", sep=''), varNC)    

#Loop through years
for(year in 1:(endYear - startYear +1)){
  
  varData_rh <- readFile(paste(fpath.lpj, "mrh.bin", sep=''), ncell.tot, nmonths, year) / conversionFactor
  varData_npp <- readFile(paste(fpath.lpj, "mnpp.bin", sep=''), ncell.tot, nmonths, year) / conversionFactor
  varData_firec <- readFile(paste(fpath.lpj, "firec.bin", sep=''), ncell.tot, 1, year) / conversionFactor
  varData_luc <- readFile(paste(fpath.lpj, "flux_luc.bin", sep=''), ncell.tot, 1, year) / conversionFactor
  varData_harvest <- readFile(paste(fpath.lpj, "flux_harvest.bin", sep=''), ncell.tot, 1, year) / conversionFactor
  varData_estab <- readFile(paste(fpath.lpj, "flux_estab.bin", sep=''), ncell.tot, 1, year) / conversionFactor

  varData_rh[is.na(varData_rh)] <- naval
  varData_rh[is.infinite(varData_rh)] <- naval
  varData_npp[is.na(varData_npp)] <- naval
  varData_npp[is.infinite(varData_npp)] <- naval
  varData_firec[is.na(varData_firec)] <- naval
  varData_firec[is.infinite(varData_firec)] <- naval
  varData_luc[is.na(varData_luc)] <- naval
  varData_luc[is.infinite(varData_luc)] <- naval
  varData_harvest[is.na(varData_harvest)] <- naval
  varData_harvest[is.infinite(varData_harvest)] <- naval
  varData_estab[is.na(varData_estab)] <- naval
  varData_estab[is.infinite(varData_estab)] <- naval

  varArray <- array(naval, c(length(lonseq),length(latseq),nmonths))
  for(month in 1:nmonths){
    for(cell in 1:ncell.tot){
      varArray[lonmatch[cell], latmatch[cell], month] <- (varData_rh[cell,month] - varData_npp[cell,month] + varData_firec[cell,1]/12 + varData_luc[cell,1]/12 + varData_harvest[cell,1]/12 - varData_estab[cell,1]/12)/(dayofmonth[month]*24*60*60*1000) #convert to kg/s
    }
    
    put.var.ncdf( ncnew, varNC, start=c(1,1,month+(year-1)*nmonths), count=c(-1,-1,1) , varArray[,,month])
  }
  
  print(year)
}
att.put.ncdf(ncnew, 0, "project", "MAIOLICA-II")
att.put.ncdf(ncnew, 0, "model", "LPJwsl")
att.put.ncdf(ncnew, 0, "run_number", "01")
att.put.ncdf(ncnew, 0, "contact", "zhen.zhang@wsl.ch")
att.put.ncdf(ncnew, 0, "institution", "WSL&MSU")
close(ncnew)

print(fpath.lpj)
print(date())
