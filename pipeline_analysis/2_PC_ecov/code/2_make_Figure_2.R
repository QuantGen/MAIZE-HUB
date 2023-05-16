
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

library(ggplot2)
source("tools/ecov_utils.R")

pheno <- read.csv("../data/PHENO.csv")

outfolder <- "2_PC_ecov/output"

load(paste0(outfolder,"/W_EVD.RData"), verbose=TRUE)

names0 <- c('location','latitude','longitude','city','region')

theme0 <- theme(
              legend.position=c(0.03,0.03), legend.justification=c(0,0),
              axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
              panel.background = element_blank())

#==================================================
# Panel A: PC1 vs PC2 plot
#==================================================
PC <- sweep(W_EVD$vectors[,1:2], 2, sqrt(W_EVD$values[1:2]),FUN='*')
PC[1:10,]

info <- pheno[match(rownames(PC),pheno$year_loc),names0]

varPC <- W_EVD$values/sum(W_EVD$values)

res <- rotate_PC(PC, latitude=info$latitude,
                     longitude=info$longitude, signPC=c(1,1))

tmp <- paste0("PC",c(1,2)," (",sprintf('%.1f',100*varPC)[c(1,2)],"%)")
dat0 <- data.frame(PC=tmp,
                  rbind(res$x.axis[1,], res$y.axis[1,]),
                  rbind(res$x.axis[2,], res$y.axis[2,]))

dat <- data.frame(res$PC, res$xy, info)
dat$state <- substr(dat$location,1,2)
head(dat)

datm <- do.call(rbind,lapply(split(dat, dat$location),function(x){
  data.frame(x[1,c("region","state","location")],
             x=median(x$x),y=median(x$y))
}))

pt1 <- ggplot(dat, aes(x,y)) +
  geom_segment(data=dat0, aes(x=x,y=y,xend=x.1,yend=y.1),color="gray30",
               arrow = arrow(length = unit(6, "pt"), ends="last")) +
  geom_point(aes(fill=region), color="gray20", size=1.1, shape=21) +
  geom_text(data=dat0[1,],aes(x=x.1,y=y.1,label=PC), angle=180+res$shift,
            size=2.6,vjust=-0.5, hjust=-0.3) +
  geom_text(data=dat0[2,],aes(x=x.1,y=y.1,label=PC), angle=90+res$shift,
            size=2.6,vjust=-0.5, hjust=1.3) +
  geom_label(data=datm, aes(fill=region, label=state), size=2.1,
             label.padding = unit(1.66, "pt")) +
  labs(x=NULL, y=NULL, fill=NULL) + theme0 +
  theme(legend.key.height = unit(0.8, 'pt'),
        legend.key.width = unit(1.5, 'pt')) +
  guides(fill = guide_legend(override.aes = aes(label = "")))

#==================================================
# Panel B: Eigen values plot
#==================================================
dat2 <- data.frame(d=c(0,W_EVD$values)/sum(W_EVD$values),x=0:length(W_EVD$values))
dat2$PCcumvar <- cumsum(dat2$d)

thr <- 0.95
dat3 <- data.frame(x=0,y=thr,label=paste0(100*thr,"%"))
tmp <- dat2[which.min(abs(dat2$PCcumvar-thr)),'x']

tmp <- dat2[which.min(abs(dat2$PCcumvar-thr)),'x']
dat4 <- data.frame(x=tmp,y=0,label=paste0(tmp," PCs"))

pt2 <- ggplot(dat2, aes(x,PCcumvar)) +
  geom_line() +
  geom_point(color="red", shape=1) +
  geom_hline(yintercept=c(0.5,thr), color="gray70", linetype="dashed") +
  geom_vline(xintercept=tmp, color="gray70", linetype="dashed") +
  scale_x_continuous(breaks=c(0,25,50, 75, 100, 125)) +
  geom_text(data=dat3, aes(x=x,y=y,label=label), size=3.1, vjust=-0.5,color="blue") +
  geom_text(data=dat4, aes(x=x,y=y,label=label),angle=90,
            size=3.1, vjust=-0.5,hjust=0.1,color="blue") +
  theme_bw() +
  labs(x="Number of eigenvectors", y="Proportion of variance")

pp <- ggpubr::ggarrange(pt1,pt2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, widths=c(0.495, 0.505))

png(paste0(outfolder,"/Figure_2.png"),res=300,units="in",width=8.8,height=4)
print(pp)
dev.off()
