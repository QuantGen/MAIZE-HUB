
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

# Load functions
library(ggplot2)
library(ggpubr)

source("../tools/ecov_utils.R")

traits <- c('Grain yield (ton/ha)'="yield",
            'Anthesis (days)'="anthesis",
            'Silking (days)'="silking",
            'ASI (days)'="ASI")[c(2,3,4,1)]
regions <- c("Across_regions","North","South")[]

outfolder <- "5_drought_heat_stress/output"

dat <- c()
for(rn in 1:length(regions)){
    region <- regions[rn]
    message("### Region =",region)

    filename <- paste0(outfolder,"/mean_effects_results_",region,"_all_traits.RData")
    if(file.exists(filename)){
      load(filename, verbose=T)

      out <- data.frame(trait2=names(traits)[match(out$trait,traits)], out)
      dat <- rbind(dat, out)
   }
}
head(dat)

dat$trait <- factor(as.character(dat$trait), levels=traits)
dat$trait2 <- factor(as.character(dat$trait2), levels=names(traits))
dat$region <- factor(as.character(dat$region),
                     levels=regions[regions%in%unique(dat$region)])

ft <- 0.3
color0 <- list(
  'HI30' = RColorBrewer::brewer.pal(8,"Spectral")[c(4,2)],
  'SDR' = RColorBrewer::brewer.pal(8,"Paired")[c(2,1)],
  'SI' = RColorBrewer::brewer.pal(8,"Set2")[c(1,2)]
)
#color0 <- RColorBrewer::brewer.pal(8,"Set2")[c(1,2)]

PP <- vector("list",length(levelsSI))
names(PP) <- names(levelsSI)

for(k in seq_along(levelsSI)){
  SI_type0 <- names(levelsSI)[k]
  dat0 <- dat[dat$SI_type==SI_type0,]

  levels0 <- levelsSI[[k]]
  dat0$group <- factor(as.character(dat0$group), levels=levels0)

  dat1 <- do.call(rbind,lapply(split(dat0, dat0$trait),function(x){
    tt <- min(x$total-x$SE.gp) - 0.4*(max(x$total+x$SE.gp) - min(x$total-x$SE.gp))
    x$ymin=ifelse(tt < 0, 0, tt)
    x
  }))

  dat1$x1 <- as.numeric(dat1$region)
  dat1$x2 <- dat1$x1 + ft*((-1)^as.numeric(dat1$group))

  tmp <- split(dat1, paste(dat1$region,dat1$trait))
  out0 <- do.call(rbind,lapply(tmp,function(x){
    data.frame(x[x$group==levels0[2],],ypos1=min(x$total-x$SE.gp),ypos2=max(x$total+x$SE.gp))
  }))
  out0$label <- sprintf('%.2f',out0$x)
  out0$label <- ifelse(out0$x>0, paste0("+",out0$label),out0$label)
  rownames(out0) <- rownames(dat0) <- NULL

  pp <- ggplot(dat1, aes(group=group,fill=group)) +
      geom_rect(aes(ymin=ymin,ymax=total,xmin=x1,xmax=x2)) +
      geom_errorbar(aes(x=(x1+x2)/2, y=total,ymin=total-SE.gp, ymax=total+SE.gp), width=0.1) +
      labs(x=NULL,y=NULL,fill=NULL,title=NULL) +
      geom_text(data=out0,aes(x=x1, y=ypos1, label=paste0("(",sprintf('%.3f',p.value.LRT),")")),
                size=2.8, vjust=1.5, color="gray5") +
      geom_text(data=out0,aes(x=x1,y=ypos2, label=label), size=2.9, vjust=-0.5, color="blue") +
      scale_fill_manual(values=color0[[SI_type0]]) +
      scale_x_continuous(breaks=seq_along(levels(dat0$region)),
                         labels=gsub("_"," ",levels(dat0$region))) +
      scale_y_continuous(expand=expansion(mult = c(0, 0.1))) +
      facet_wrap(~trait2, scales="free_y") + theme_bw() +
      theme(
       strip.text.x = element_text(size=9, margin=margin(t=2,b=2)),
       legend.margin=margin(t=-0.2,b=0.15,l=0.15,r=0,unit='lines'),
       legend.key.size=unit(0.6,"lines"),
       legend.text=element_text(size=8, margin=margin(l=-4)),
       legend.position=c(0.99,0.995), legend.justification=c(1,1),
       plot.margin = margin(t=5,r=5,b=5,l=12, "pt")
      )
  PP[[k]] <- pp

}

# Save Figure 7
pp <- ggarrange(PP[['SDR']],PP[['HI30']],
                labels = c("A", "B"),
                ncol = 2, nrow = 1, widths=c(0.5, 0.5))

pp <- annotate_figure(pp, bottom=text_grob("Region", vjust=-0.3, size=10),
                      left=text_grob(expression(Mean %+-% SE), vjust=1.5, rot=90, size=10))

png(paste0(outfolder,"/Figure_7.png"), res=300,units="in",width=9.2,height=4)
print(pp)
dev.off()

# # Save Figure S7: Combining SDR and HI
pp <- annotate_figure(PP[['SI']], bottom=text_grob("Region", vjust=-0.3, size=10),
                      left=text_grob(expression(Mean %+-% SE), vjust=1.5, rot=90, size=10))

png(paste0(outfolder,"/Figure_S7.png"),res=300,units="in",width=5.2,height=4)
print(pp)
dev.off()
