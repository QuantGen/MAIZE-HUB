
#==========================================
# This code will produce
#   Table 3:  If region='North'
#   Table S4: If region='South'
#==========================================
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

library(reshape2)

traits <- c("yield","anthesis","silking","ASI")
models <- c('Random Effects'="YEAR+LOC+YL+H+HL",
            'SNP+EC'="SNP+EC+SNPxEC")

region <- c("North","South")[2] # Select a region

sources <- c("YEAR","LOC","YL","(Total YL)","EC",
             "H","SNP","HL","SNPxEC","Error")

format_VC <- function(VC, models=models, sources=sources){
  out <- do.call(cbind, lapply(split(VC, VC[,"model"]),function(x){
    tt <- paste0(sprintf('%.3f',x$var),"\n(",sprintf('%.3f',x$sd),")")
    ttt <- rep("--",length(sources)); names(ttt) <- sources
    ttt[match(x$source, sources)] <- tt
    ttt
  }))
  out <- out[,models]
  if(!is.null(names(models))){
    index <- which(names(models)!="")
    colnames(out)[index] <- names(models)[index]
  }

  if(!is.null(names(sources))){
    index <- which(names(sources) != "")
    rownames(out)[index] <- names(sources)[index]
  }
  out <- data.frame('source/model'=rownames(out),out, check.names=F)
  rownames(out) <- NULL
  out
}

merge_traits <- function(VC){
  OUT <- VC[[1]][,1, drop=FALSE]
  for(k in 1:length(VC)){
    tmp <- VC[[k]][,-1]
    colnames(tmp) <- paste0(names(VC)[k],". ", colnames(tmp))
    OUT <- data.frame(OUT,tmp, check.names=FALSE)
  }
  return(OUT)
}

capitalize <- function(string){
  substr(string,1,1) <- toupper(substr(string,1,1))
  string
}

outfolder <- "6_genomic_prediction/output"

### All traits togethe (by region)
VC2 <- vector('list', length(traits))

for(tr in 1:length(traits)){
  trait <- traits[tr]

  VC0 <- c()
  for(mod in 1:length(models)){
    model <- as.vector(models[mod])
    filename <- paste0(outfolder,"/ANOVA/",model,"/",trait,"/",region,"/VC.RData")
    if(file.exists(filename)){
      load(filename)
      VC0 <- rbind(VC0,data.frame(model, VC))
    }
  }
  VC0$model <- factor(VC0$model, levels=models)

  vc0 <- format_VC(VC0, models=models, sources=sources)
  VC2[[tr]] <- vc0
  names(VC2)[tr] <- capitalize(trait)

  cat("\n")
}

OUT <- merge_traits(VC2)

# Table 3 for North, and Table S4 for South
tmp <- ifelse(region=="North",3,"S4")
write.table(OUT, file=paste0(outfolder,"/Table_",tmp,".csv"),
            row.names=FALSE, sep=",")
