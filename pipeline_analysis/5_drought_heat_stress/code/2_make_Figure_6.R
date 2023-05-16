
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

library(ggplot2)
source("tools/ecov_utils.R")

pheno <- read.csv("../data/PHENO.csv")
ECOV0 <- read.csv("../data/ECOV_layered.csv",row.names=1)
ECOV <- read.csv("../data/ECOV.csv",row.names=1)
ECOV_KL <- read.csv("../data/ECOV_KL.csv")

# Get a period around flowering
period <- c("FlaFlw","FlwStG","StGEnG")

# HI30 index
HI30 <- rowSums(ECOV[,paste0("HI30_",period)])

# SDR index
tmp <- data.frame(ECOV[,grep("Eo_|CoverGreen_",colnames(ECOV))], ECOV0)
SI <- get_SDR(EC=tmp, EC_KL=ECOV_KL)
SDR <- rowMeans(SI[,paste0("SDR_",period)])

dat <- data.frame(year_loc=names(SDR),SDR,HI30)
dat$region <- factor(pheno[match(dat$year_loc, pheno$year_loc),"region"])

thrHI <- 36
thrSDR <- 0.54

levelsHI <- paste0("HI30",c("<=",">"),thrHI)
levelsSDR <- paste0("SDR",c("<=",">"),thrSDR)

dat$HI_group <- factor(ifelse(dat$HI30<=thrHI,levelsHI[1],levelsHI[2]),levels=levelsHI)
dat$SDR_group <- factor(ifelse(dat$SDR<=thrSDR,levelsSDR[1],levelsSDR[2]),levels=levelsSDR)

dat2 <- do.call(rbind,lapply(split(dat,paste(dat$SDR_group,dat$HI_group)),function(x){
  tt <- table(x$region)
  data.frame(x[rep(1,length(tt)),c("HI_group","SDR_group")],
                    region=names(tt),nYL=as.vector(tt))
}))

rx <- range(dat$HI30)
ry <- range(dat$SDR)
dy <- diff(ry)
dat2$hjust <- ifelse(dat2$HI_group==levelsHI[2],1.2,-0.5)
dat2$x <- ifelse(dat2$HI_group==levelsHI[2],rx[2]-0.015*diff(rx),rx[1]+0.015*diff(rx))
dat2$y <- ifelse(dat2$SDR_group==levelsSDR[2],thrSDR+0.03*dy,thrSDR-0.03*dy)

dat3 <- rbind(   # For text
  data.frame(x=rx[1]+0.17*diff(rx),y=thrSDR,label=paste0("SDR=",sprintf('%.2f',thrSDR))),
  data.frame(x=thrHI,y=ry[1]+0.82*diff(ry),label=paste0("HI30=",sprintf('%.0f',thrHI)))
)

pp <- ggplot(dat, aes(HI30, SDR)) + theme_bw()+
   geom_point(aes(fill=region),shape=21, size=3) + labs(fill=NULL) +
   geom_hline(yintercept=thrSDR, linetype="dashed", color="blue") +
   geom_vline(xintercept=thrHI, linetype="dashed", color="blue") +
   geom_label(data=dat2[dat2$region=="North",], aes(x=x,y=y,label=nYL,fill=region),
              hjust=1, size=3.3,
              label.padding=unit(2.0,"pt"),label.r=unit(0,"pt")) +
   geom_label(data=dat2[dat2$region=="South",], aes(x=x,y=y,label=nYL,fill=region),
              hjust=0, size=3.3,
              label.padding=unit(2.0,"pt"),label.r=unit(0,"pt")) +
   geom_text(data=dat3[1,],aes(x=x,y=y,label=label),
             color="blue",vjust=-0.5) +
   geom_text(data=dat3[2,],aes(x=x,y=y,label=label),angle=90,
             color="blue",vjust=-0.5) +
   labs(x="HI30 (days)") +
   theme(legend.position=c(0.99,0.99), legend.justification=c(1,1)) +
   guides(fill = guide_legend(override.aes = aes(label = "")))

# Save outputs
outfolder <- "5_drought_heat_stress/output"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

png(paste0(outfolder,"/Figure_6.png"),res=200,units="in",width=5.6,height=5.3)
print(pp)
dev.off()
