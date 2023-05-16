
#==========================================
# This code will produce
#   Figure 4:  If region='North'
#   Figure S6: If region='South'
#==========================================
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

source("tools/ecov_utils.R")

# Load functions
library(ggplot2)

traits <- c('Grain yield'="yield",
            'Days-to-anthesis'="anthesis",
            'Days-to-silking'="silking",
            'ASI'="ASI")

region <- c("North","South")[2] # Select a region

outfolder <- "4_ecov_WAS/output"

# Get results for all traits
dat <- c()
for(i in 1:length(traits)){
  filename <- paste0(outfolder,"/",region,"/ecov_WAS_results_",traits[i],".RData")
  if(file.exists(filename)){
    load(filename, verbose=T)
    dat <- rbind(dat,out)
  }
}
head(dat)

colnames(dat)[grep("^Pr",colnames(dat))] <- "Pvalue"
range(dat$Pvalue)

index <- match(unique(dat$ecov),dat$ecov)
ecov_annotation <- dat[index,"Type",drop=FALSE]
rownames(ecov_annotation) <- dat[index,"ecov"]

plot1 <- plot_EWAS(data=dat, ecov_annotation=ecov_annotation,
           select=traits[traits%in%unique(dat$trait)],
           color_pal=c("azure2","lightsteelblue"),
           alpha=0.5, text.size=3, ylab=expression(-log[10]~'(P)'),
           expand.y=0.01,annotation_height=0.07)

# Figure 4 for North, and Figure S6 for South
tmp <- ifelse(region=="North",4,"S6")
png(paste0(outfolder,"/Figure_",tmp,".png"),
    res=300, units="in", width=8.5, height=5.5)
plot(plot1)
dev.off()
