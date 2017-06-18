#June 16th, 2016
#Generatedaily LPJ NEE netcdfs from daily LPJ variables and yearly variables(carbon fluxes from disturbances)
# HAVE NOT DONE YET
#Load libraries
print(date())
library(sp)
library(ncdf)
library(hash)

#Set parameters
fpath.lpj <- "/mnt/lustrefs/store/zhen.zhang/output/MERRA2/MERRA2_2016_USDA_DLIT_PERMAFROST/merge/"
startYear <- 2014
endYear <- 2016

#Fixed parameters
conversionFactor <- 1
outNyears <- length(startYear:endYear)
ntstep <- 12
ndaily <- 365
nmonths <- 12
nyear <- endYear - startYear + 1
naval <- -99999
ncell.tot <- file.info(paste(fpath.lpj,'grid.out',sep=''))$size/4

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
varNC <- var.def.ncdf("mnee", "g C m-2 month-1", list(x,y,tmonth), naval, "Net Ecosystem Exchange", prec="double")
ncnew <- create.ncdf( paste(fpath.lpj, "LPJ_NEE.nc", sep=''), varNC)    

#Loop through years
for(year in 1:(endYear - startYear +1)){
  
  varData_rh <- readFile(paste(fpath.lpj, "mrh.bin", sep=''), ncell.tot, nmonths, year) / conversionFactor
  varData_npp <- readFile(paste(fpath.lpj, "mnpp.bin", sep=''), ncell.tot, nmonths, year) / conversionFactor
  varData_firec <- readFile(paste(fpath.lpj, "firec.bin", sep=''), ncell.tot, 1, year) / conversionFactor
  varData_luc <- readFile(paste(fpath.lpj, "flux_luc.bin", sep=''), ncell.tot, 1, year) / conversionFactor
  varData_harvest <- readFile(paste(fpath.lpj, "flux_harvest.bin", sep=''), ncell.tot, 1, year) / conversionFactor
  varData_estab <- readFile(paste(fpath.lpj, "flux_estab.bin", sep=''), ncell.tot, 1, year) / conversionFactor

  varData_rh[is.na(varData_rh)] <- -99999
  varData_rh[is.infinite(varData_rh)] <- -99999
  varData_npp[is.na(varData_npp)] <- -99999
  varData_npp[is.infinite(varData_npp)] <- -99999
  varData_firec[is.na(varData_firec)] <- -99999
  varData_firec[is.infinite(varData_firec)] <- -99999
  varData_luc[is.na(varData_luc)] <- -99999
  varData_luc[is.infinite(varData_luc)] <- -99999
  varData_harvest[is.na(varData_harvest)] <- -99999
  varData_harvest[is.infinite(varData_harvest)] <- -99999
  varData_estab[is.na(varData_estab)] <- -99999
  varData_estab[is.infinite(varData_estab)] <- -99999

  varArray <- array(naval, c(length(lonseq),length(latseq),nmonths))
  for(month in 1:nmonths){
    for(cell in 1:ncell.tot){
      varArray[lonmatch[cell], latmatch[cell], month] <- varData_rh[cell,month] - varData_npp[cell,month] + varData_firec[cell,1]/12 + varData_luc[cell,1]/12 + varData_harvest[cell,1]/12 - varData_estab[cell,1]/12
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
