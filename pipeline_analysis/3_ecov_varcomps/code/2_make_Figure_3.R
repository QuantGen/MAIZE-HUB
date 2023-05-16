
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

library(ggplot2)
library(RColorBrewer)
source("tools/ecov_utils.R")

outfolder <- "3_ecov_varcomps/output"

load(paste0(outfolder,"/ecov_varcomps_all_regions.RData"))

comp_names <- c('Region'="Region",
                'Location'="Location",
                'Year(Location)'="Year_Location")

index <- match(unique(OUT$ecov),OUT$ecov)
ecov_annotation <- OUT[index,"Type",drop=FALSE]
rownames(ecov_annotation) <- OUT[index,"ecov"]

regions <- unique(OUT$region)
names(regions) <- gsub("_"," ",regions)

comp_mean <- lapply(regions,function(x)"Location")
comp_mean[["Across regions"]] <- c("Region","Location")

comp_color <- brewer.pal(8, "Paired")[c(2,1,3)]
comp_color <- c("royalblue","lightskyblue","darkseagreen1")

plot1 <- plot_varcomp(data=OUT, comp_names, ylab="Proportion of variance",
             ecov_annotation=ecov_annotation,
             comp_mean=comp_mean, select=regions[],
             color_pal=NULL, bar.color=NA, comp_color=comp_color)

png(paste0(outfolder,"/Figure_3.png"),
   res=300, units="in", width=8.5, height=6.5)
print(plot1)
dev.off()
