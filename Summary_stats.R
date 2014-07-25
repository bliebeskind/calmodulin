files <- list(DeltaG="EF_DeltaG.csv",Flexibility="EF_Flex.csv",Isoelectric="EF_IsoEl.csv",Volume="EF_volumes.csv")

## Iterate over Biophysical parameters

for(name in names(files)){
	src <- read.csv(files[[name]],header=TRUE)

	EF <- c("EF1","EF2","EF3","EF4")

## Write csv files of means and standard deviations for each for each
## biophysical property at each site, for each EF hand. 
## Creates 8 files total (mean and stdev for each of 4 properties).

	means <- list()
	stdevs <- list()

	for(i in EF){
		ef <- src[src$EF==i,]
		means[[i]] <- sapply(ef[2:13],'mean')
		stdevs[[i]] <- sapply(ef[2:13],'sd')
	}

	Site <- paste(1:12)
	means_frame <- cbind(Site,data.frame(means))
	stdev_frame <- cbind(Site,data.frame(stdevs))
	write.csv(means_frame,paste(name,"_means.csv",sep=''),row.names=FALSE)
	write.csv(stdev_frame,paste(name,"_stdevs.csv",sep=''),row.names=FALSE)
}