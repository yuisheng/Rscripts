#Summerize Time Series of LPJ outputs and export it as R object
#Zhen 21st Sep. 2016
#minterc.bin, pft_* and firec were not summerized

#Set Parameters
dirs.lpj <- c("/mnt/lustrefs/store/zhen.zhang/output/CRU/MERRA2_2015_USDA_DLIT/merge/",
              "/mnt/lustrefs/store/zhen.zhang/output/CRU/CRU_2015_USDA_DLIT_NOFIRE/merge/")

varlist <- c("mnpp","mrh","firec","flux_luc","flux_harvest","flux_estab","mgpp","litc","soilc","vegc",
             "msnowpack","mra","mpet","msoiltemp","mswc1","mswc2","mevap","mdischarge","mrunoff","mtransp"
             )
ismonthly <- c(T,T,F,F,F,F,T,F,F,F,
               T,T,T,T,T,T,T,T,T,T)
isareal <- c(T,T,T,T,T,T,T,T,T,T,
             F,F,F,F,F,F,F,F,F,F)
# divided by varfactor
varfacter <- c(10^15,10^15,10^15,10^15,10^15,10^15,10^15,10^15,10^15,10^15,
               1,1,1,1,1,1,1,1,1,1)
startyear <- 1901
endyear <- 2015
ncells <- c()

#Fixed parameters
for(i in 1:length(dirs.lpj)){
  ncells[i] <- file.info(paste(dirs.lpj[i],'grid.out',sep=''))$size/4  
}
nmonths <- 12
nyears <- endyear - startyear + 1
nvars <- length(varlist)  #currently evalute 11 carbon-related variables


#Set up the LPJ output read function
readYFile <- function(infile, ncell, nmonths, year){
  fileName <- file(infile, 'rb')
  seek(fileName, 4*ncell*nmonths*(year-1))
  binOut <- readBin(fileName, double(), size=4, n = ncell*nmonths)
  close(fileName)
  return(matrix(binOut, ncell, nmonths))
}
readMFile <- function(infile, ncell, nmonth, year){
   fileName <- file(infile, 'rb')
   seek(fileName, 4*ncell*(12*(year-1)+nmonth-1))
   binOut <- readBin(fileName, double(), size=4, n = ncell)
   close(fileName)
   return(matrix(binOut, ncell, 1))
}


for(dir in dirs.lpj){
  dirID <- which(dir==dirs.lpj)
  
  #Create R list for storing output
  lpjout.gl <- list()
  lpjout.n  <- list()
  lpjout.tr <- list()
  lpjout.s  <- list()
  lpjout.gl[['time']] <- c(startyear,endyear)
  lpjout.n[['time']] <- c(startyear,endyear)
  lpjout.tr[['time']] <- c(startyear,endyear)
  lpjout.s[['time']] <- c(startyear,endyear)
  
  
  outGrid <- readBin(paste(dirs.lpj[dirID], "grid.out", sep=''), integer(), size=2, n=ncells[dirID]*2)/100
  outLon <- outGrid[seq(1,ncells[dirID]*2,2)]
  outLat <- outGrid[seq(2,ncells[dirID]*2,2)]
  deg2m <- (111*111*10^6*0.5*0.5*cos(outLat*(pi/180)))
  
  for(ivar in 1:nvars){
    if(ismonthly[ivar]){
      for(year in 1:nyears){
        for(imonth in 1:12){
          varRead <- readMFile(paste(dirs.lpj[dirID],varlist[ivar],".bin", sep=""),ncells[dirID],imonth,year)
          #check if there is NAN or Infinite
          varRead[which(is.na(varRead)==T)] <- 0
          varRead[which(is.infinite(varRead)==T)] <- 0
          if(isareal[ivar]){
            lpjout.gl[[varlist[ivar]]][(year-1)*12+imonth] <- sum(rowSums(varRead)*deg2m)/varfacter[ivar]  
            lpjout.n[[varlist[ivar]]][(year-1)*12+imonth] <- sum(varRead[which(outLat >= 30),]*deg2m[which(outLat >= 30)])/varfacter[ivar]
            lpjout.tr[[varlist[ivar]]][(year-1)*12+imonth] <- sum(varRead[which(outLat > -30 & outLat < 30),]*deg2m[which(outLat > -30 & outLat < 30)])/varfacter[ivar]
            lpjout.s[[varlist[ivar]]][(year-1)*12+imonth] <- sum(varRead[which(outLat <= -30),]*deg2m[which(outLat <= -30)])/varfacter[ivar]
            
          }else{
            lpjout.gl[[varlist[ivar]]][(year-1)*12+imonth] <- mean(rowMeans(varRead))/varfacter[ivar]
            lpjout.n[[varlist[ivar]]][(year-1)*12+imonth] <- mean(varRead[which(outLat >= 30),])/varfacter[ivar]
            lpjout.tr[[varlist[ivar]]][(year-1)*12+imonth] <- mean(varRead[which(outLat > -30 & outLat < 30),])/varfacter[ivar]
            lpjout.s[[varlist[ivar]]][(year-1)*12+imonth] <- mean(varRead[which(outLat <= -30),])/varfacter[ivar]
            
          }
          
        }
      }
    }else{
      for(year in 1:nyears){
        varRead <- readYFile(paste(dirs.lpj[dirID],varlist[ivar],".bin", sep=""),ncells[dirID],1,year)
        #check if there is NAN or Infinite
        varRead[which(is.na(varRead)==T)] <- 0
        varRead[which(is.infinite(varRead)==T)] <- 0
        if(isareal[ivar]){
          lpjout.gl[[varlist[ivar]]][year] <- sum(rowSums(varRead)*deg2m)/varfacter[ivar]  
          lpjout.n[[varlist[ivar]]][year]  <- sum(varRead[which(outLat >= 30),]*deg2m[which(outLat >= 30)])/varfacter[ivar]
          lpjout.tr[[varlist[ivar]]][year] <- sum(varRead[which(outLat > -30 & outLat < 30),]*deg2m[which(outLat > -30 & outLat < 30)])/varfacter[ivar]
          lpjout.s[[varlist[ivar]]][year] <- sum(varRead[which(outLat <= -30),]*deg2m[which(outLat <= -30)])/varfacter[ivar]
        }else{
          lpjout.gl[[varlist[ivar]]][year] <- mean(rowMeans(varRead))/varfacter[ivar]
          lpjout.n[[varlist[ivar]]][year] <- mean(varRead[which(outLat >= 30),])/varfacter[ivar]
          lpjout.tr[[varlist[ivar]]][year] <- mean(varRead[which(outLat > -30 & outLat < 30),])/varfacter[ivar]
          lpjout.s[[varlist[ivar]]][year] <- mean(varRead[which(outLat <= -30),])/varfacter[ivar]
          
        }
        
        
      }
    }
  print(paste(dir,varlist[ivar],sep=" "))  
  }
  #write R object
  save(lpjout.gl,lpjout.n,lpjout.tr,lpjout.s,file=paste(dirs.lpj[dirID],"Output.RData",sep=''))
  lpjout.gl <- NULL
  lpjout.n  <- NULL
  lpjout.tr <- NULL
  lpjout.s  <- NULL
}

