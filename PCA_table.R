library(ggbiplot)
library(grid)
library(gridExtra)

files <- list(DeltaG="EF_DeltaG_cleaned.csv",Flexibility="EF_Flex_cleaned.csv",Isoelectric="EF_IsoEl_cleaned.csv",Volume="EF_volumes_cleaned.csv")

## Iterate over Biophysical parameters
## Perform PCA analysis

plot_list <- list()

for(name in names(files)){
	src <- read.csv(files[[name]],header=TRUE)

	EF <- c("EF1","EF2","EF3","EF4")

	std_vals <- as.data.frame(scale(src[2:13])) # Standardize variables
	pca <- prcomp(std_vals)

	g <- ggbiplot(pca,choices=1:2,obs.scale=1,var.scale=1,groups=src$EF,ellipse=TRUE,var.axes=FALSE)
	g <- g + scale_color_brewer(palette='PuOr',name='')
	g <- g + theme(legend.direction="horizontal",legend.position='top')
	g <- g + ggtitle(name)
	plot_list[[name]] <- g
}

## Function for extracting legend to share between plots
## https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

## Function for making a grid
## http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


ef_legend <- g_legend(plot_list[["DeltaG"]])
ef_legend <- editGrob(ef_legend,vp=viewport(x=.75,y=.5))
layout <- matrix(1:4,ncol=2,byrow=TRUE)


multiplot(plot_list[["DeltaG"]] + theme(legend.position="none"),
		plot_list[["Flexibility"]] + theme(legend.position="none"),
		plot_list[["Isoelectric"]] + theme(legend.position="none"),
		plot_list[["Volume"]] + theme(legend.position="none"),
		layout=layout)

grid.draw(ef_legend)
