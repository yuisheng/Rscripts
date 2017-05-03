#Oct 2014
#Assemble LPJ netcdfs from monthly LPJ variables
#Submit as: qsub -q long runR.sh
#Load libraries
print(date())
library(sp)
library(ncdf)
library(hash)

#Set parameters
fpath.lpj <- "/Users/zhang/Research/data/output/MERRA2/MERRA2_2016_USDA_DLIT_PERMAFROST/"
varName <- "dch4e"
tstep <- 1   # for determine the output time step (1: daily;2: monthly;3: yearly)
startYear <- 2015
endYear <- 2016

#Fixed parameters
conversionFactor <- 1
#caltype <- "average" # "average"   #determine summary type when desired output is annually
filename <- paste(varName, ".bin",sep="")

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

#For historical
outNyears <- length(startYear:endYear)
ntstep <- 12
ndays <- 365
nmonths <- 12
nyear <- endYear - startYear + 1
naval <- -99999
ncell.tot <- file.info(paste(fpath.lpj,'grid.out',sep=''))$size/4

#generate corresponding hashmap for varName and its attributes
unit_hash <- hash()
outname_hash <- hash()
longname_hash <- hash()
#Now we have 62 outputs
#List of variable in names and out names to loop through
#LPJ output names to be read
varNameDaily <- c("dch4e","dgpp","dnpp","drh","dsm1","dsm2")
varNameMonthly <- c("mswc1","mswc2","mrunoff","mdischarge","mevap","mtransp",
                    "mgpp","mra","mnpp","mrh",
                    "msoiltemp","mtemp_soil","mtsoil_0","mtsoil_10","mtsoil_25","mtsoil_50","mtsoil_100","mtsoil_150","mtsoil_200",
                    "msnowpack","msnowdepth","frozen_days",
                    "mthaw_depth","mFwater","mFice","mice_frac1","mice_frac2",
                    "wtd","wet_frac","mch4e","ch4o"
)
varNameAnnual <- c("vegc","litc","soilc","firec","flux_luc") #, "flux_estab","flux_harvest","firef")

varNameAnnualPFT <- c("pft_lai", "pft_gpp","pft_gc","pft_transp","fpc")
varNameMonthlyPFT <- c("mpft_ci","mpft_gc","mpft_gpp","mpft_lai","mpft_transp")
varNameDeprecated <- c("minterc","mirrig_wd","mpet",    # monthly output
                       "sdate","waterstress",         # ?
                       "pft_bimonfpar","pft_fO3uptake","pft_harvest","pft_maxphenday","pft_mort","pft_nind","pft_npp","pft_rharvest","pft_vegc"  #PFT OUTPUT
)


#Official output variable names MUST KEEP ORDER SAME AS BELOW FOR VARxxx
outNameDaily   <- c("dch4e","dgpp","dnpp","drh","dsm1","dsm2")
outNameMonthly <-    c("mswc1","mswc2","mrunoff","mdischarge","mevap","mtransp",
                       "mgpp","mra","mnpp","mrh",
                       "msoiltemp","mtemp_soil","mtsoil_0","mtsoil_10","mtsoil_25","mtsoil_50","mtsoil_100","mtsoil_150","mtsoil_200",
                       "msnowpack","msnowdepth","frozen_days",
                       "mthaw_depth","mFwater","mFice","mice_frac1","mice_frac2",
                       "wtd","wet_frac","mch4e","ch4o")
outNameAnnual <-     c("cVeg","cLitter","cSoil","fFire","fLuc")  #,"fGrazing","burntArea")
outNameAnnualPFT <-  c("pft_alai", "pft_agpp","pft_agc","pft_atransp","landCoverFrac")
outNameMonthlyPFT <- c("mpft_ci","mpft_gc","mpft_gpp","mpft_lai","mpft_transp")
outNameDeprecated <- c("minterc","mirrig_wd","mpet",    # monthly output
                       "sdate","waterstress",         # ?
                       "pft_bimonfpar","pft_fO3uptake","pft_harvest","pft_maxphenday","pft_mort","pft_nind","pft_npp","pft_rharvest","pft_vegc"  #PFT OUTPUT
)

unitNameDaily  <- c("g CH4 /m2 /day","kg C m-2 month-1","kg C m-2 month-1","kg C m-2 month-1","fraction","fraction")
unitNameMonthly <- c("fraction","fraction","mm month-1","mm month-1","mm month-1","mm h2o m-2 month-1",
                    "kg C m-2 month-1","kg C m-2 month-1","kg C m-2 month-1","kg C m-2 month-1",
                    "degreesC","degreesC","degreesC","degreesC","degreesC","degreesC","degreesC","degreesC","degreesC",
                    "mm month-1","mm month-1","days",
                    "mm month-1","fraction","fraction","fraction","fraction",
                    "meter month-1","fraction"," g CH4 /m2 /month","g CH4 /m2 /month"
)

unitNameAnnual <- c("kg C m-2 yr-1","kg C m-2 yr-1","kg C m-2 yr-1","kg C m-2 yr-1","kg C m-2 yr-1") #, "flux_estab","flux_harvest","firef")

unitNameAnnualPFT <- c("-", "kg C m-2 month-1","kg h2o m-2 month-1","kg h2o m-2 month-1","-")
unitNameMonthlyPFT <- c("kg C m-2 month-1","kg C m-2 month-1","kg C m-2 month-1","-","kg h2o m-2 month-1")
unitNameDeprecated <- c("NA","NA","NA",    # monthly output
                       "NA","NA",         # ?
                       "NA","NA","NA","NA","NA","NA","NA","NA","NA"  #PFT OUTPUT
)

longNameDaily <- c("CH4 emission","Gross Primary Production","Net Primary Production","Heterotrophic Respiration","Relative Soil Moisture in Upper Layer","Relative Soil Moisture in Lower Layer")
longNameMonthly <- c("Relative Soil Moisture in Upper Layer","Relative Soil Moisture in Lower Layer","Total Runoff","Total Discharge","Total Evapo-Transpiration","Transpiration",
                    "Gross Primary Production","Autotrophic respiration","Net Primary Production","Heterotrophic Respiration",
                    "Mean Temperature of Soil","Mean Temperature of Soil (Permafrost)","Soil Temperature at 0 cm","Soil Temperature at 10 cm","Soil Temperature at 25 cm","Soil Temperature at 50 cm","Soil Temperature at 100 cm","Soil Temperature at 150 cm","Soil Temperature at 200 cm",
                    "Snowpack Depth","Snowpack Depth (Permafrost)","Frozen Days",
                    "Thaw Depth","Water Fraction in 10 cm soil","Ice Fraction in 10 cm soil","Ice Fraction in Upper Layer","Ice Fraction in Lower Layer",
                    "Water Table Depth","Wetland Fraction","CH4 emission","Oxidated CH4 Production"
)
longNameAnnual <- c("Carbon in Vegetation","Carbon in Litter","Carbon in Soil","CO2 emissions from fire","CO2 Flux from Land Use Change") #, "flux_estab","flux_harvest","firef")

longNameAnnualPFT <- c("Leaf Area Index", "Gross Primary Production","Stomatal Conductance","Transpiration","Fractional Land Cover of PFT")
longNameMonthlyPFT <- c("mpft_ci","mpft_gc","mpft_gpp","mpft_lai","mpft_transp")
longNameDeprecated <- c("minterc","mirrig_wd","mpet",    # monthly output
                       "sdate","waterstress",         # ?
                       "pft_bimonfpar","pft_fO3uptake","pft_harvest","pft_maxphenday","pft_mort","pft_nind","pft_npp","pft_rharvest","pft_vegc"  #PFT OUTPUT
)

for(i in 1:length(varNameDaily)){
  unit_hash[varNameDaily[i]] <- unitNameDaily[i]
  outname_hash[varNameDaily[i]] <- outNameDaily[i]
  longname_hash[varNameDaily[i]] <- longNameDaily[i]
}

for(i in 1:length(varNameMonthly)){
  unit_hash[varNameMonthly[i]] <- unitNameMonthly[i]
  outname_hash[varNameMonthly[i]] <- outNameMonthly[i]
  longname_hash[varNameMonthly[i]] <- longNameMonthly[i]
}
for(i in 1:length(varNameAnnual)){
  unit_hash[varNameAnnual[i]] <- unitNameAnnual[i]
  outname_hash[varNameAnnual[i]] <- outNameAnnual[i]
  longname_hash[varNameAnnual[i]] <- longNameAnnual[i]
}
for(i in 1:length(varNameAnnualPFT)){
  unit_hash[varNameAnnualPFT[i]] <- unitNameAnnualPFT[i]
  outname_hash[varNameAnnualPFT[i]] <- outNameAnnualPFT[i]
  longname_hash[varNameAnnualPFT[i]] <- longNameAnnualPFT[i]
}
for(i in 1:length(varNameMonthlyPFT)){
  unit_hash[varNameMonthlyPFT[i]] <- unitNameMonthlyPFT[i]
  outname_hash[varNameMonthlyPFT[i]] <- outNameMonthlyPFT[i]
  longname_hash[varNameMonthlyPFT[i]] <- longNameMonthlyPFT[i]
}
for(i in 1:length(varNameDeprecated)){
  unit_hash[varNameDeprecated[i]] <- unitNameDeprecated[i]
  outname_hash[varNameDeprecated[i]] <- outNameDeprecated[i]
  longname_hash[varNameDeprecated[i]] <- longNameDeprecated[i]
}





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
tday <- dim.def.ncdf("time",paste("days since", startYear, "-01-01",sep=""), 0:(outNyears*ndays-1))
tmonth <- dim.def.ncdf( "time", paste("months since ", startYear, "-01-01", sep=""), 0:(outNyears*nmonths-1))
tyear <- dim.def.ncdf( "time", paste("years since ", startYear, "-01-01", sep=""), 0:(outNyears-1))


#Make ncdf
if(tstep==1){
  varNC <- var.def.ncdf(varName, hash::values(unit_hash[varName]), list(x,y,tday), naval, longname=hash::values(longname_hash[varName]), prec="double")
  ncnew <- create.ncdf( paste(fpath.lpj, "LPJ_d", hash::values(outname_hash[varName]),".nc", sep=''), varNC)    
}else if(tstep==2){
  varNC <- var.def.ncdf(varName, hash::values(unit_hash[varName]), list(x,y,tmonth), naval, longname=hash::values(longname_hash[varName]), prec="double")
  ncnew <- create.ncdf( paste(fpath.lpj, "LPJ_m", hash::values(outname_hash[varName]),".nc", sep=''), varNC)    
}else{
  varNC <- var.def.ncdf(varName, hash::values(unit_hash[varName]), list(x,y,tyear), naval, longname=hash::values(longname_hash[varName]), prec="double")
  ncnew <- create.ncdf( paste(fpath.lpj, "LPJ_a", hash::values(outname_hash[varName]),".nc", sep=''), varNC)    
}

#Loop through years
for(year in (startYear:endYear)-(startYear-1)){
  yearID <- which(year == (startYear:endYear)-(startYear-1))

  #Set up empty array
  if(tstep==1){
    #For daily outputs
    varData <- readFile(paste(fpath.lpj, filename, sep=''), ncell.tot, ndays, year) / conversionFactor
    varData[is.na(varData)] <- -99999
    varData[is.infinite(varData)] <- -99999
    
    varArray <- array(naval, c(length(lonseq),length(latseq),ndays))
    for(day in 1:ndays){
      for(cell in 1:ncell.tot){
        varArray[lonmatch[cell], latmatch[cell], day] <- varData[cell,day]
      }
      put.var.ncdf( ncnew, varNC, start=c(1,1,day+(yearID-1)*ndays), count=c(-1,-1,1) , varArray[,,day])
    }
  }else if(tstep==2){
    #For monthly outputs
    varData <- readFile(paste(fpath.lpj, filename, sep=''), ncell.tot, nmonths, year) / conversionFactor
    varData[is.na(varData)] <- -99999
    varData[is.infinite(varData)] <- -99999
    
    varArray <- array(naval, c(length(lonseq),length(latseq),nmonths))
    for(month in 1:nmonths){
      for(cell in 1:ncell.tot){
        varArray[lonmatch[cell], latmatch[cell], month] <- varData[cell,month]
      }
      
      put.var.ncdf( ncnew, varNC, start=c(1,1,month+(yearID-1)*nmonths), count=c(-1,-1,1) , varArray[,,month])
    }
  }else{
    #For yearly outputs
    varData <- readFile(paste(fpath.lpj, filename, sep=''), ncell.tot, 1, year) / conversionFactor
    varData[is.na(varData)] <- -99999
    varData[is.infinite(varData)] <- -99999
    
    varArray <- array(naval, c(length(lonseq),length(latseq),1))
    for(cell in 1:ncell.tot){
      varArray[lonmatch[cell], latmatch[cell], 1] <- varData[cell,1]
    }
    #varArray[varArray == -99999] = NA
    put.var.ncdf( ncnew, varNC, start=c(1,1,yearID), count=c(-1,-1,1) , varArray[,,1])
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
