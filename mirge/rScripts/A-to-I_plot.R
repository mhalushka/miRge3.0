library(ggplot2)
options(warn=-1)
args <- commandArgs(TRUE)
cancerCellsData <- read.csv(file=args[1], header=TRUE, sep=",")
midValue <- median(cancerCellsData$log2RPM, na.rm=TRUE)

sample1 <- unique(cancerCellsData$sample)
miRNA2 <-  unique(cancerCellsData$miRNA.position)
miRNA1 <- rev(miRNA2)

pdf(args[2], width=8,height=8)

ggplot(cancerCellsData, aes(y = factor(cancerCellsData$miRNA.position, levels=miRNA1), x = factor(cancerCellsData$sample, levels=sample1)))+       ## global aes
  ## to get the rect filled
geom_point(aes(colour = log2RPM, size =A.to.I.percentage))+        ## geom_point for circle illusion
scale_colour_gradient2(low='lightyellow',mid='orange',high='red',midpoint=midValue)+     ## color of the corresponding aes
scale_size(limits=c(0,100), range = c(1, 8), breaks=c(0,20,40,60,80,100))+             ## to tune the size of circles
theme_bw()+
labs(x = "", y = "")+
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="bottom",
      legend.direction='horizontal',
      legend.box='horizontal',
      legend.box.just='bottom',
      axis.text.x=element_text(angle=45, hjust=1))
garbage <- dev.off()



