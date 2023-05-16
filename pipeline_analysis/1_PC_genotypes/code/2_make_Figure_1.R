
rm(list=ls())
library(ggplot2)

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

outfolder <- "1_PC_genotypes/output"

load(paste0(outfolder,"/G_EVD.RData"), verbose=TRUE)

ID <- rownames(G_EVD$vectors)

tester <- unlist(lapply(strsplit(ID,"/"),function(x)x[length(x)]))

tester_list <- list(
    'CG'=c("C103","CG102","CG108","CG110","CG111","CG120","CG123","CG44","CG60",
           "CGR01","CGR03"),
    'DK3IIH6'=c("DK3IIH6"),
    'LH'=c("LH145","LH162","LH210","LH38","LH51"),
    'LH82'=c("LH82"),
    'LH198'=c("LH198"),
    'LH185'=c("LH185"),
    'LH195'=c("LH195"),
    'PH'=c("PB80","PH207","PHAJ0","PHG29","PHG35","PHG47",
                "PHG83","PHHV4","PHJ65",
                "PHM49","PHM57","PHN11","PHN47","PHN82",
                "PHR03","PHR25","PHR55","PHR63","PHRE1",
                "PHW03","PHW30","PHW53","PHW65"),
    'PHB'=c("PHB47"),
    'PHK'=c("PHK05","PHK56","PHK76"),
    'PHP'=c("PHP02","PHP60"),
    'PHT'=c("PHT69","PHTD5"),
    'PHZ'=c("PHZ51"),
    'TX'=c("TX205","TX6252","TX714","TX775","TX777","TX779"),
    'Other'=c("MO17","NK787","OH43","OH7B","H95","Q381","S8324")
)

INFO <- data.frame(ID,tester)
INFO$tester_group <- rep(NA,nrow(INFO))
for(i in 1:length(tester_list)){
  INFO$tester_group[tester %in% tester_list[[i]]] <- names(tester_list)[i]
}
INFO$tester_group <- factor(as.character(INFO$tester_group), levels=names(tester_list))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color0 <- gg_color_hue(5)
shape0 <- c(21,22,24)

n0 <- nlevels(INFO$tester_group)
color2 <- rep(color0,ceiling(n0/length(color0)))[seq(n0)]
shape2 <- rep(shape0,each=ceiling(n0/length(shape0)))[seq(n0)]

names(color2) <- names(shape2) <- levels(INFO$tester_group)

#==================================================
# Panel A: PC1 vs PC2 plot
#==================================================
PC <- sweep(G_EVD$vectors[,1:2], 2, sqrt(G_EVD$values[1:2]),FUN='*')

dat <- data.frame(PC,INFO[match(rownames(PC),INFO$ID),])
pVar <- paste0("PC",1:2," (",sprintf('%.1f',100*G_EVD$values[1:2]/sum(G_EVD$values)),"%)")

pt1 <- ggplot(dat,aes(X1,X2)) + theme_bw() +
       geom_point(aes(shape=tester_group,fill=tester_group),
                  size=0.78, color="gray20") +
       labs(x=pVar[1],y=pVar[2],shape=NULL,fill=NULL) +
       scale_fill_manual(values=color2) +
       scale_shape_manual(values=shape2) +
       theme(
             legend.position=c(0.01,0.99),legend.justification=c(0,1),
             legend.margin=margin(t=-0.2,b=0.15,l=0.15,r=0.15,unit='line'),
             legend.key.height=unit(0.5,"line"),
             legend.text=element_text(size=7.5),
             plot.margin = margin(t=5,r=5,b=5,l=5, "pt")) +
       guides(fill = guide_legend(override.aes = aes(size = 2.2)))


#==================================================
# Panel B: Eigen values plot
#==================================================
dat2 <- data.frame(d=c(0,G_EVD$values)/sum(G_EVD$values),x=0:length(G_EVD$values))
dat2$PCcumvar <- cumsum(dat2$d)

thr <- 0.95
dat3 <- data.frame(x=0,y=thr,label=paste0(100*thr,"%"))

tmp <- dat2[which.min(abs(dat2$PCcumvar-thr)),'x']
dat4 <- data.frame(x=tmp,y=0,label=paste0(tmp," PCs"))

pt2 <- ggplot(dat2, aes(x,PCcumvar)) +
  geom_line() +
  geom_point(color="red", shape=1) +
  geom_hline(yintercept=c(0.5,thr), color="gray70", linetype="dashed") +
  geom_vline(xintercept=tmp, color="gray70", linetype="dashed") +
  geom_text(data=dat3, aes(x=x,y=y,label=label),
            size=3.1, vjust=-0.5,color="blue") +
  geom_text(data=dat4, aes(x=x,y=y,label=label),angle=90,
            size=3.1, vjust=-0.5,hjust=0.1,color="blue") +
  theme_bw() +
  theme(plot.margin = margin(t=5,r=5,b=5,l=12, "pt")) +
  labs(x="Number of eigenvectors", y="Proportion of variance")

pp <- ggpubr::ggarrange(pt1,pt2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, widths=c(0.495, 0.505))

png(paste0(outfolder,"/Figure_1.png"),res=300,units="in",width=8.8,height=4)
print(pp)
dev.off()
