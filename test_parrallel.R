#test foreach
library(foreach)
library(numbers)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
wd.dir <- args[1]
rm(args)

setwd(wd.dir)

max.eig <- function(N, sigma) {
d <- matrix(rnorm(N**2, sd = sigma), nrow = N)
E <- eigen(d)$values
abs(E)[[1]]
}

system.time(
            E <- sapply(1:10000,function(n) { sapply(1:5,function(m) {max.eig(m,n)})})
            )

system.time(
            F <- foreach(n=1:100) %:% when(isPrime(n)) %:% foreach(m=1:5) %do% max.eig(n,m)
            )
dir.create(paste(wd.dir,"/test",sep=""))

write(E,paste(wd.dir,"/test/a.rda",sep=""))

