
#==========================================
# This code will produce
#   Figure 8:  If region='North'
#   Figure S8: If region='South'
#==========================================

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

library(ggplot2)
library(reshape2)
library(ggpubr)
library(RColorBrewer)

source("tools/Functions.R")

traits <- c('Grain yield'="yield",
            'Days-to-anthesis'="anthesis",
            'Days-to-silking'="silking",
            'ASI'="ASI")

models <- c('Random Effects'="YEAR+LOC+YL+H+HL",
            'Random Effects+PC'="YEAR+LOC+YL+H+HL+GPC+WPC",
            'SNP+EC'="SNP+EC+SNPxEC")

region <- c("North","South")[2]  # Select a region

outfolder <- "6_genomic_prediction/output"

# Function to get within year_loc correlation
get_corr <- function(dat, y_name="y", yHat_name="yHat", by="year_loc"){
  res <- do.call(rbind,lapply(split(dat, as.character(dat[,by])),function(x){
         x$y <- x[,y_name]
         x$yHat <- x[,yHat_name]
         cond1=(!is.na(x$y))
         cond2=(!is.na(x$yHat))
         nRecords=sum(cond1&cond2)
         correlation=cor(x$y,x$yHat,use='pairwise.complete')
         tmp <- (x$y-x$yHat)^2
         mse <- sum(tmp[cond1&cond2])/nRecords
         data.frame(trait,region,model,
                    x[1,c(by),drop=F], correlation, mse, nRecords,
                    SE_cor=sqrt((1-correlation^2)/nRecords)
         )}))
 rownames(res) <- NULL
 res
}

#======================================
# Collect 10F-CV results
#======================================
folds <- c(1:10)
RES1 <- c()
for(tr in 1:length(traits)){
  trait <- as.vector(traits[tr])
  cat("Trait ",trait, ". Region ",region,"\n")

  for(mm in 1:length(models)){
    model <- as.vector(models[mm])
    out2 <- c()
    for(ff in 1:length(folds)){
      fold <- folds[ff]
      filename <- paste0(outfolder,"/10F_CV/",model,"/",trait,"/",region,
                    "/results_fold_",fold,".RData")

      if(file.exists(filename)){
        load(filename)
        out2 <- rbind(out2, data.frame(fold=ff,out))
      }else{
        cat("File '",filename,"'does not exists\n")
      }
    }

    # Within year correlation
    res1 <- get_corr(out2, by="year_loc")
    RES1 <- rbind(RES1, data.frame(CV="10F-CV",res1))
  }
}
head(RES1)

#======================================
# Collect LYO-CV results
#======================================
years <- c(2014:2021)
RES2 <- c()
for(tr in 1:length(traits)){
  trait <- traits[tr]
  cat("Trait ",trait, ". Region ",region,"\n")

  for(mm in 1:length(models)){
    model <- as.vector(models[mm])
    for(yr in 1:length(years)){
      year <- years[yr]
      filename <- paste0(outfolder,"/LYO_CV/",model,"/",trait,"/",region,
                    "/results_",year,".RData")

      if(file.exists(filename)){
        load(filename)
        res1 <- get_corr(out, by="year_loc")
        RES2 <- rbind(RES2, data.frame(CV="LYO-CV",res1))
      }else{
        cat("File '",filename,"'does not exists\n")
      }
    }
  }
}
head(RES2)

out <- rbind(RES1, RES2)   # Within year-locations correlation

table(out$CV)
out$CV <- factor(as.character(out$CV), levels=c("10F-CV","LYO-CV"))
out$year <- factor(as.character(out$year))

out$trait <- names(traits)[match(out$trait, traits)]
out$trait <- factor(capitalize(out$trait), levels=names(traits))

out$model2 <- names(models)[match(out$model, models)]
out$model2 <- factor(out$model2, levels=names(models))

wd <- 0.23   # Separation between boxplots within group
out$x1 <- as.numeric(out$CV)
out$x2 <- out$x1 -wd*(nlevels(out$model2)-1)/2 + wd*(as.numeric(out$model2)-1)

namesX <- levels(out$CV)

color0 <- c("#56B4E9","palegreen1","#E69F00")
names(color0) <- names(models)

datm <- do.call(rbind,lapply(split(out, paste(out$trait,out$CV,out$model2)),
 function(x){
    z <- 0.5*log((1+x$correlation)/(1-x$correlation))  # fisher transformation
    ni <- x$nRecords
    SE <- 1/sqrt(ni-3)
    VAR <- sum((ni-1)*(SE^2))/(sum(ni)-nrow(x))
    MEAN <- sum(ni*z)/sum(ni)
    prec <- 1/(x$SE^2)
    data.frame(x[1,c("region","trait","CV","model2","x2")],
               correlation=sum((x$correlation*prec))/(sum(prec)),
               SE=sqrt(1/sum(prec)),
               ypos1=min(out$correlation)
               )
}))
datm$labelSD <- paste0("(",sprintf("%.3f",datm$SE),")")
rownames(datm) <- NULL

pp <- ggplot(out, aes(x=x2, y=correlation, group=CV:model2)) +
      stat_boxplot(geom="errorbar", width=wd/2, position=position_dodge(width =wd)) +
      geom_boxplot(aes(fill=model2), width=wd*0.9, outlier.shape=21,
                   outlier.size=1, position=position_dodge(width=wd)) +
      facet_wrap(~trait) + theme_bw() +
      geom_label(data=datm, aes(x2,y=ypos1,label=sprintf("%.2f",correlation)), size=3,
                vjust=-0.08, color="gray40", label.padding=unit(2.0,"pt"), label.size=NA) +
      geom_label(data=datm, aes(x2,y=ypos1,label=labelSD), size=2.5,
                vjust=0.9, color="gray40", label.padding=unit(2.0,"pt"), label.size=NA) +
      scale_x_continuous(breaks=seq_along(levels(out$CV)), labels=levels(out$CV)) +
      scale_y_continuous(expand=expansion(0.04,0.04)) +
      scale_fill_manual(values=color0) +
      labs(title=NULL, x="Cross-validation", fill=NULL,
           y='Within year-location correlation') +
      theme(plot.title=element_text(hjust=0.5),
            legend.key.height=unit(0.7,"line"),
            legend.text=element_text(size=7),
            strip.text.x = element_text(size=9.0, margin=margin(t=2,b=2)),
            legend.position=c(0.99,0.47), legend.justification=c(1,1),
            legend.margin=margin(t=-0.20,b=0.17,l=0.17,r=0.17,unit='line')
           )

# Figure 8 for North, and Figure S8 for South
tmp <- ifelse(region=="North",8,"S8")
png(paste0(outfolder,"/Figure_",tmp,".png"), res=300, units="in", width=6.3, height=5.3)
print(pp)
dev.off()
