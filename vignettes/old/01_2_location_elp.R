rm(list=ls())
setwd("~/hdeeg")
library(eeguana)
library(eegusta)


files=dir("./data_bdf",full.names = TRUE,pattern = "elp$")

# for(i in 1:length(files)){
#   print(i)
#   res=read_elp(files[i])
# }


locs=lapply(files,read_elp)

LOC=locs[[1]]
for(i in 2:length(locs)){
  LOC[,5:7]=LOC[,5:7]+locs[[i]][,5:7]
}

LOC[,5:7]=LOC[,5:7]/length(locs)

D=c()
for(i in 1:length(locs)){
  D[i]=sum((LOC[,5:7]-locs[[i]][,5:7])^2)
}

plot(D)
plot(sort(D))
summary(D)

dists=rowSums(locs[[1]][,5:7]-LOC[,5:7])^2
summary(dists)
plot(dists)

locs[[1]][which.max(dists),]
LOC[which.max(dists),]


save(LOC, file="locations_eeguana.Rdata")
