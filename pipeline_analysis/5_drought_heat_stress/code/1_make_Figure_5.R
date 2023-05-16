
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

library(ggplot2)
library(reshape2)
source("tools/ecov_utils.R")

pheno <- read.csv("../data/PHENO.csv")
ECOV0 <- read.csv("../data/ECOV_layered.csv",row.names=1)
ECOV <- read.csv("../data/ECOV.csv",row.names=1)
ECOV_KL <- read.csv("../data/ECOV_KL.csv")

ECOV_info <- get_ec_info(ECOV)

# Phenological stages
index <- match(levels(ECOV_info$period), ECOV_info$period)
ECOV_info0 <- ECOV_info[index,c("period","start","end")]
ECOV_info0 <- ECOV_info0[-c(1,9),]

# Get SDR
tmp <- data.frame(ECOV[,grep("Eo_|CoverGreen_",colnames(ECOV))],ECOV0)
SDR <- get_SDR(EC=tmp, EC_KL=ECOV_KL)

out <- c()
names_col <- c("year_loc","year","location","region")
for(i in 1:nrow(ECOV_info0)){
  stg <- ECOV_info0[i,]

  sdr <- SDR[,paste0("SDR_",stg$period)]

  info <- pheno[match(rownames(SDR), pheno$year_loc),names_col]

  tmp <- data.frame(x=i, info, stg[rep(1,nrow(SDR)),],SDR=sdr)
  out <- rbind(out, tmp)
}
out$year <- factor(as.character(out$year))
# head(out)

dat1 <- data.frame(x=c(seq(nrow(ECOV_info0))-0.5,nrow(ECOV_info0)+0.5),
                   label=c(ECOV_info0$start,ECOV_info0$end[nrow(ECOV_info0)]))
rgx <- range(out$x)

maxSDR <- 4

out$SDR2 <- out$SDR
out[out$SDR>maxSDR,'SDR2'] <- maxSDR

# Data for annotation
dat2 <- do.call(rbind,lapply(split(out,out$region),function(x){
  tt <- x[which.min(x$x),]
  tt$x <- tt$x-0.5  #-0.125*rgx[2]
  tt$y <- maxSDR# + 0.085*maxSDR
  tt
}))
dat2$label <- ifelse(dat2$region=="North","A","B")

theme0 <- theme(
      legend.position=c(0.01,0.01), legend.justification=c(0,0),
      axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
      legend.margin=margin(t=-0.15,r=0.2,b=0.15,l=0.2,unit='line'),
      legend.box.background = element_rect(colour = "black",linewidth=0.75),
      legend.key.height=unit(0.69,"line"),
      legend.text = element_text(size=8),
      strip.text.x = element_text(size=10, margin=margin(t=2,b=2)),
      panel.spacing.x = unit(1, "lines")
    )

pp <- ggplot(out, aes(x, SDR2)) +
    geom_line(aes(group=year_loc, color=year),linewidth=0.27) +
    geom_point(aes(group=year_loc, color=year),size=0) +
    labs(x="Phenological stage", y="SDR", color=NULL) +
    theme_bw() + theme0 +
    geom_text(dat=dat2,aes(x=x, y=y, label=label),
              vjust=-0.45,hjust=1.3,size=5,fontface = "bold") +
    coord_cartesian(xlim=c(rgx[1]-0.5,rgx[2]+0.5), ylim=c(0,maxSDR), clip="off") +
    facet_wrap(~region, ncol=2) +
    scale_x_continuous(breaks=dat1$x,labels=dat1$label,expand=c(0,0)) +
    scale_y_continuous(expand=c(0.01,0.02)) +
    guides(color = guide_legend(override.aes=list(size = 2.5)))

# Save outputs
outfolder <- "5_drought_heat_stress/output"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

png(paste0(outfolder,"/Figure_5.png"), res=300, units="in", width=8.0, height=4.3)
print(pp)
dev.off()
