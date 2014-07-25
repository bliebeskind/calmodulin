files <- list(DeltaG="EF_DeltaG_cleaned.csv",Flexibility="EF_Flex_cleaned.csv",Isoelectric="EF_IsoEl_cleaned.csv",Volume="EF_volumes_cleaned.csv")

for(name in names(files)){
	src <- read.csv(files[[name]],header=TRUE)
	std_vals <- as.data.frame(scale(src[2:13])) # Standardize variables
	pca <- prcomp(std_vals)
	outfile <- paste(name,"_scree.pdf")
	pdf(outfile)
	screeplot(pca,type="barplot",main=name)
	dev.off()
}