#June 17th 2016
# Generate GEOS input of NEP and CH4E from LPJ outputs
# GEOS inputs could be either daily or 3-hourly 
# List of LPJ output:
# dch4e
# dnep = drh - dnpp + (firec + flux_luc + flux_harvest - flux_estab)/365
print(date())
library(sp)
library(ncdf)
library(hash)

is.leapyear <- function(year){
      #http://en.wikipedia.org/wiki/Leap_year
      return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}


#Set parameters
fpath.lpj <- "/Users/zhang/Research/data/output/JRA55/JRA55_2016_USDA_DLIT_PERMAFROST/merge/"
#This is start year of binary output of LPJ
startYear <- 1901
endYear <- 2016

# Start year of daily output
dstartYear <- 1980
dendYear <- 2016
dlength <- length(seq(ISOdate(dstartYear,1,1),ISOdate(dendYear,12,31),'days'))

#Fixed parameters
conversionFactor <- 1*1000*86400 # convert from g/m2/day to kg/m2/s
outNyears <- length(startYear:endYear)
ntstep <- 12
ndays <- 365
nmonths <- 12
nyear <- endYear - startYear + 1
naval <- 1.e15
ncell.tot <- file.info(paste(fpath.lpj,'grid.out',sep=''))$size/4

dayofmonth <- c(31,28,31,30,31,30,31,31,30,31,30,31) # without leap day

#Set up the LPJ output read function
readFile <- function(infile, cell.tot, nlength, year){
  fileName <- file(infile, 'rb')
  seek(fileName, 4*cell.tot*nlength*(year-1))
  binOut <- readBin(fileName, double(), size=4, n = cell.tot*nlength)
  close(fileName)
  return(matrix(binOut, ncell.tot, nlength))
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
tday <- dim.def.ncdf("time",paste("days since ", dstartYear, "-01-01 00:00:00",sep=""), 0:(dlength-1))  # including leap days
tmonth <- dim.def.ncdf( "time", paste("months since ", startYear, "-01-01 00:00:00", sep=""), 0:(outNyears*nmonths-1))
tyear <- dim.def.ncdf( "time", paste("years since ", startYear, "-01-01 00:00:00", sep=""), 0:(outNyears-1))


#Make ncdf
varCH4ENC <- var.def.ncdf("CH4tot", "kg C m^-2 s^-1", list(x,y,tday), naval, "Total Wetland CH4 emission", prec="double")
ncnew.ch4e <- create.ncdf( paste(fpath.lpj, "LPJ_ch4e_daily_GEOS_",dstartYear,"-",dendYear,".nc", sep=''), varCH4ENC)

#Loop through years
for(year in dstartYear:dendYear){
  yearID  <- year - startYear + 1 
  dyearID <- year - dstartYear + 1
  #count.leaps <- length(which(is.leapyear(seq(dstartYear,dstartYear+dyearID-2)))) #count leap years between dstartyear and year
  
  #read in daily outputs
  varData.dch4e <- readFile(paste(fpath.lpj, "dch4e.bin", sep=''), ncell.tot, ndays, dyearID)/conversionFactor
  varData.dch4e[is.na(varData.dch4e)] <- naval
  varData.dch4e[is.infinite(varData.dch4e)] <- naval
  varArray.dch4e <- array(naval, c(length(lonseq),length(latseq),ndays))

  if(is.leapyear(year)){
    for(day in 1:(ndays+1)){
      print(day)
      if(day < 60){
        for(cell in 1:ncell.tot){
          varArray.dch4e[lonmatch[cell], latmatch[cell], day] <- varData.dch4e[cell,day]
        }
        put.var.ncdf( ncnew.ch4e, varCH4ENC, start=c(1,1,day+(dyearID-1)*ndays), count=c(-1,-1,1) , varArray.dch4e[,,day])
        
      }else if(day == 60){
        put.var.ncdf( ncnew.ch4e, varCH4ENC, start=c(1,1,day+(dyearID-1)*ndays), count=c(-1,-1,1) , varArray.dch4e[,,59])
        
      }else{
        for(cell in 1:ncell.tot){
          varArray.dch4e[lonmatch[cell], latmatch[cell], day-1] <- varData.dch4e[cell,day-1]
        }
        put.var.ncdf( ncnew.ch4e, varCH4ENC, start=c(1,1,day+(dyearID-1)*ndays), count=c(-1,-1,1) , varArray.dch4e[,,day-1])
      }
    }
  }else{
    for(day in 1:(ndays)){
      print(day)
      for(cell in 1:ncell.tot){
        varArray.dch4e[lonmatch[cell], latmatch[cell], day] <- varData.dch4e[cell,day]
      }
      put.var.ncdf( ncnew.ch4e, varCH4ENC, start=c(1,1,day+(dyearID-1)*ndays), count=c(-1,-1,1) , varArray.dch4e[,,day])
      
    }
  }
  
  print(year)
}
att.put.ncdf(ncnew.ch4e, 0, "project", "MAIOLICA-II")
att.put.ncdf(ncnew.ch4e, 0, "model", "LPJwsl")
att.put.ncdf(ncnew.ch4e, 0, "run_number", "01")
att.put.ncdf(ncnew.ch4e, 0, "contact", "zhen.zhang@wsl.ch")
att.put.ncdf(ncnew.ch4e, 0, "institution", "WSL&MSU")

close(ncnew.ch4e)

print(fpath.lpj)
print(date())
