#Check fluxes July 2015

#Set parameters for each model run
fpath.out <- "/home/zhen.zhang/lpj-transient2015_wflux_MERRA2_vs_CRU_version2.pdf"
dirs.lpj <- c("/mnt/lustrefs/store/zhen.zhang/output/MERRA2/MERRA2_2015_USDA_MLIT_PERMAFROST/merge/","/mnt/lustrefs/store/zhen.zhang/output/CRU/CRU_2015_USDA_MLIT_PERMAFROST/merge/")
legend.label <- c("MERRA2_PERMAFROST","CRU")
main.label <- c("Soil Moisture top","Soil Moisture bot","Snow Pack","PET","Soil Temp","Water Stress","EVAP","Runoff","TRANSP", "Discharge","INTERC")
varlist <- c("mswc1","mswc2","msnowpack","mpet","msoiltemp","waterstress","mevap","mrunoff","mtransp","mdischarge","minterc")
ismonthly <- c(T,T,T,T,T,T,T,T,T,T,T)
startyear <- 1901
endyear <- 2015

#Constant parameters
ncells <- c(file.info(paste(dirs.lpj[1],'grid.out',sep=''))$size/4,file.info(paste(dirs.lpj[2],'grid.out',sep=''))$size/4)
nmonths <- 12
nyears <- endyear - startyear + 1
nvars <- length(varlist)  #currently evalute 11 carbon-related variables

#Set up the LPJ output read function
readFile <- function(infile, ncell, nmonths, year){
	fileName <- file(infile, 'rb')
	seek(fileName, 4*ncell*nmonths*(year-1))
	binOut <- readBin(fileName, double(), size=4, n = ncell*nmonths)
	close(fileName)
	return(matrix(binOut, ncell, nmonths))
}


#Read years
varAnn <- array(0,    c(nyears, length(varlist), length(dirs.lpj)))
varAnn.n <- array(0,  c(nyears, length(varlist), length(dirs.lpj)))
varAnn.tr <- array(0, c(nyears, length(varlist), length(dirs.lpj)))
varAnn.s <- array(0,  c(nyears, length(varlist), length(dirs.lpj)))
for(dir in dirs.lpj){
  dirID <- which(dir==dirs.lpj)

  outGrid <- readBin(paste(dirs.lpj[dirID], "grid.out", sep=''), integer(), size=2, n=ncells[dirID]*2)/100
  outLon <- outGrid[seq(1,ncells[dirID]*2,2)]
  outLat <- outGrid[seq(2,ncells[dirID]*2,2)]
  deg2m <- (111*111*10^6*0.5*0.5*cos(outLat*(pi/180)))
  
  for(year in 1:nyears){
    for (ivar in 1:nvars){
      
      if(ismonthly[ivar]){
        varRead <- readFile(paste(dir,varlist[ivar],".bin", sep=""),ncells[dirID],nmonths,year)
        #check if there is NAN or Infinite
        varRead[which(is.na(varRead)==T)] <- 0
        varRead[which(is.infinite(varRead)==T)] <- 0
        varAnn[year,ivar,dirID] <- mean(varRead)
        varAnn.n[year,ivar,dirID] <- mean(varRead[which(outLat >= 30)])
        varAnn.tr[year,ivar,dirID] <- mean(varRead[which(outLat > -30 & outLat < 30)])
        varAnn.s[year,ivar,dirID] <- mean(varRead[which(outLat <- -30 & outLat > 30)])
      }else{
        varRead <- readFile(paste(dir,varlist[ivar],".bin", sep=""),ncells[dirID],1,year)
        varAnn[year,ivar,dirID] <- mean(varRead)     
        varAnn.n[year,ivar,dirID] <-  mean(varRead[which(outLat >= 30)])
        varAnn.tr[year,ivar,dirID] <-  mean(varRead[which(outLat > -30 & outLat < 30)])
        varAnn.s[year,ivar,dirID] <- mean(varRead[which(outLat <- -30 & outLat > 30)])
      }
      
    }
    
    print(year)
  }
  
}

pdf(paste(fpath.out, sep=""))
par(mfrow=c(4,3),mar=c(2,2,2,2),xpd=F)

for(ii in 1:(nvars)){
    plot(startyear:endyear,varAnn[,ii,1], type='l',main=main.label[ii],
         ylab="PgC",xlab="Year",
         ylim=c(floor(min(varAnn[,ii,1],varAnn[,ii,2])),ceiling(max(varAnn[,ii,1],varAnn[,ii,2]))))
    lines(startyear:endyear,varAnn[,ii,2], col=2)  
    
    
}

legend(1902,4.5,legend.label,col=c(1,2),lwd=2)
par(xpd=NA)

dev.off()


# write.table(cbind(startyear:endyear, -1*((varAnn[,3]+varAnn[,4]+varAnn[,5]+varAnn[,6])-(varAnn[,2]+varAnn[,7])),
# 	-1*((varAnn[,3]+varAnn[,4]+varAnn[,5]+varAnn[,6])-(varAnn[,2]+varAnn[,7]))),
# 	paste(fpath.out,"lpj_nbp.txt",sep=""), row.names=F,col.names=c("Year","S2","S3"))
# 
# write.table(cbind(startyear:endyear, -1*((varAnn[,3]+varAnn[,4]+varAnn[,5]+varAnn[,6])-(varAnn[,2]+varAnn[,7])),
# 								-1*((varAnn.n[,3]+varAnn.n[,4]+varAnn.n[,5]+varAnn.n[,6])-(varAnn.n[,2]+varAnn.n[,7])),
# 								-1*((varAnn.tr[,3]+varAnn.tr[,4]+varAnn.tr[,5]+varAnn.tr[,6])-(varAnn.tr[,2]+varAnn.tr[,7])),
# 								-1*((varAnn.s[,3]+varAnn.s[,4]+varAnn.s[,5]+varAnn.s[,6])-(varAnn.s[,2]+varAnn.s[,7]))),
# 	paste(fpath.out,"lpj_nbp_regions.txt",sep=""), row.names=F,col.names=c("Year","Globe","N30","TR","S30"))
