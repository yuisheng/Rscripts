#Oct 2014
#Assemble LPJ netcdfs from monthly LPJ variables
#Usage: varname rcp ismonthly
#Load libraries
library(sp)
library(ncdf)
library(hash)

options(echo=TRUE) #if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
varname <- args[1]
rcp <- args[2]
ismonthly <- (args[3] > 0)
rm(args)

#rcp <- c("rcp26","rcp45","rcp60","rcp85")
for (i in 1:length(rcp)){
  cmip5_path <- paste("/mnt/lustrefs/work/zhen.zhang/output/", rcp[i],"/", sep='') 
  model_list <- dir(cmip5_path)
  for (j in 1:length(model_list)){
    fpath.lpj <- paste(cmip5_path,model_list[j],"/merge/", sep="")
    
    #set input parameters
    varName <- varname
    conversionFactor <- 1
    #caltype <- "average" # "average"   #determine summary type when desired output is annually
    filename <- paste(varName, ".bin",sep="")
    
    #Parameters
    #For CMIP5
    startYear <- 1961
    endYear <- 2099
    outNyears <- length(startYear:endYear)
    ntstep <- 12
    nmonths <- 12
    nyear <- endYear - startYear + 1
    naval <- -99999
    ncell.tot <- 62482
    
    #For CRU
    # startYear <- 1901
    # endYear <- 2013
    # outNyears <- length(startYear:endYear)
    # ntstep <- 12
    # nmonths <- 12
    # nyear <- 113
    # naval <- -99999
    # ncell.tot <- 62482
    
    #generate corresponding hashmap for varName and its attributes
    unit_hash <- hash()
    outname_hash <- hash()
    longname_hash <- hash()
    #Now we have 62 outputs
    #List of variable in names and out names to loop through
    #LPJ output names to be read
    varNameMonthly <- c("mswc1","mswc2","mrunoff","mdischarge","mevap","mtransp",
                        "mgpp","mra","mnpp","mrh",
                        "msoiltemp","mtemp_soil","mtsoil_0","mtsoil_10","mtsoil_25","mtsoil_50","mtsoil_100","mtsoil_150","mtsoil_200",
                        "msnowpack","msnowdepth","frozen_days",
                        "mthaw_depth","mFwater","mFice","mice_frac1","mice_frac2",
                        "wtd","wet_frac","ch4e","ch4o"
    )
    varNameAnnual <- c("vegc","litc","soilc","firec","flux_luc") #, "flux_estab","flux_harvest","firef")
    
    varNameAnnualPFT <- c("pft_lai", "pft_gpp","pft_gc","pft_transp","fpc")
    varNameMonthlyPFT <- c("mpft_ci","mpft_gc","mpft_gpp","mpft_lai","mpft_transp")
    varNameDeprecated <- c("minterc","mirrig_wd","mpet",    # monthly output
                           "sdate","waterstress",         # ?
                           "pft_bimonfpar","pft_fO3uptake","pft_harvest","pft_maxphenday","pft_mort","pft_nind","pft_npp","pft_rharvest","pft_vegc"  #PFT OUTPUT
    )
    
    
    #Official output variable names MUST KEEP ORDER SAME AS BELOW FOR VARxxx
    outNameMonthly <-    c("mswc1","mswc2","mrunoff","mdischarge","mevap","mtransp",
                           "mgpp","mra","mnpp","mrh",
                           "msoiltemp","mtemp_soil","mtsoil_0","mtsoil_10","mtsoil_25","mtsoil_50","mtsoil_100","mtsoil_150","mtsoil_200",
                           "msnowpack","msnowdepth","frozen_days",
                           "mthaw_depth","mFwater","mFice","mice_frac1","mice_frac2",
                           "wtd","wet_frac","ch4e","ch4o")
    outNameAnnual <-     c("cVeg","cLitter","cSoil","fFire","fLuc")  #,"fGrazing","burntArea")
    outNameAnnualPFT <-  c("pft_alai", "pft_agpp","pft_agc","pft_atransp","landCoverFrac")
    outNameMonthlyPFT <- c("mpft_ci","mpft_gc","mpft_gpp","mpft_lai","mpft_transp")
    outNameDeprecated <- c("minterc","mirrig_wd","mpet",    # monthly output
                           "sdate","waterstress",         # ?
                           "pft_bimonfpar","pft_fO3uptake","pft_harvest","pft_maxphenday","pft_mort","pft_nind","pft_npp","pft_rharvest","pft_vegc"  #PFT OUTPUT
    )
    
    unitNameMonthly <- c("fraction","fraction","mm month-1","mm month-1","mm month-1","mm h2o m-2 month-1",
                         "g C m-2 month-1","g C m-2 month-1","g C m-2 month-1","g C m-2 month-1",
                         "degreesC","degreesC","degreesC","degreesC","degreesC","degreesC","degreesC","degreesC","degreesC",
                         "mm month-1","mm month-1","days",
                         "mm month-1","fraction","fraction","fraction","fraction",
                         "meter month-1","fraction","10-1 g CH4 /m2 /month"," 10-1 g CH4 /m2 /month"
    )
    
    unitNameAnnual <- c("g C m-2 yr-1","g C m-2 yr-1","g C m-2 yr-1","g C m-2 yr-1","g C m-2 yr-1") #, "flux_estab","flux_harvest","firef")
    
    unitNameAnnualPFT <- c("-", "kg C m-2 month-1","kg h2o m-2 month-1","kg h2o m-2 month-1","-")
    unitNameMonthlyPFT <- c("kg C m-2 month-1","kg C m-2 month-1","kg C m-2 month-1","-","kg h2o m-2 month-1")
    unitNameDeprecated <- c("NA","NA","NA",    # monthly output
                            "NA","NA",         # ?
                            "NA","NA","NA","NA","NA","NA","NA","NA","NA"  #PFT OUTPUT
    )
    
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
    latseq <- rev(seq(-89.75, 89.75, 0.5))
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
    if(ismonthly){
      varNC <- var.def.ncdf(varName, values(unit_hash[varName]), list(x,y,tmonth), naval, longname=values(longname_hash[varName]), prec="double")
      ncnew <- create.ncdf( paste(fpath.lpj, "LPJ_", values(outname_hash[varName]),".nc", sep=''), varNC)    
    }else{
      varNC <- var.def.ncdf(varName, values(unit_hash[varName]), list(x,y,tyear), naval, longname=values(longname_hash[varName]), prec="double")
      ncnew <- create.ncdf( paste(fpath.lpj, "LPJ_a", values(outname_hash[varName]),".nc", sep=''), varNC)    
    }
    
    #Loop through years
    for(year in (startYear:endYear)-(startYear-1)){
      yearID <- which(year == (startYear:endYear)-(startYear-1))
      #Set up empty array
      if(ismonthly){
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
        varData <- readFile(paste(fpath.lpj, filename, sep=''), ncell.tot, 1, year) / conversionFactor
        varData[is.na(varData)] <- -99999
        varData[is.infinite(varData)] <- -99999

        varArray <- array(naval, c(length(lonseq),length(latseq),1))
        for(cell in 1:ncell.tot){
          varArray[lonmatch[cell], latmatch[cell], 1] <- varData[cell,1]
        }
        put.var.ncdf( ncnew, varNC, start=c(1,1,yearID), count=c(-1,-1,1) , varArray[,,1])
        }
      
      #print(year)
    }
    att.put.ncdf(ncnew, 0, "project", "MAIOLICA-II")
    att.put.ncdf(ncnew, 0, "model", "LPJwsl")
    att.put.ncdf(ncnew, 0, "run_number", "01")
    att.put.ncdf(ncnew, 0, "contact", "zhen.zhang@wsl.ch")
    att.put.ncdf(ncnew, 0, "institution", "WSL&MSU")
    close(ncnew)
    
    print(fpath.lpj)
    
  }
}


