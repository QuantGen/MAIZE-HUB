
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

source("tools/ecov_utils.R")

# Load functions
library(lme4)

regions <- c("Across_regions","North","South")

pheno <- read.csv("../data/PHENO.csv")
ECOV <- read.csv("../data/ECOV.csv", row.names=1)

ECOV <- ECOV[,-grep("HI30_",colnames(ECOV))]

comp_names <- c("Region","Location","Year_Location")
OUT <- c()

for(reg in seq_along(regions))
{
  region <- regions[reg]

  if(region == "Across_regions"){
    pheno0 <- pheno[]
  }else{
    pheno0 <- pheno[pheno$region %in% region,]
  }

  ECOV0 <- ECOV[rownames(ECOV) %in% pheno0$year_loc, ]
  W <- scale(ECOV0)

  info <- pheno0[match(rownames(W),pheno0$year_loc),c("year_loc","year","location","region")]
  Region <- factor(as.character(info$region))
  Location <- factor(as.character(info$location))

  out <- c()
  for(j in 1:ncol(W)){
    if(all(!is.na(W[,j]))){
      if(region == "Across_regions"){
        fm <- aov(W[,j] ~ Region + Location)
      }else{
        fm <- aov(W[,j] ~ Location)
      }
      res <- as.data.frame(summary(fm)[[1]])
      rownames(res) <- gsub(" ","",rownames(res))
      rownames(res)[grep("Residuals",rownames(res))] <- "Year_Location"

      tmp <- res[comp_names,'Sum Sq']
    }else{
      tmp <- rep(NA,length(comp_names))
    }
    names(tmp) <- comp_names
    tmp <- data.frame(region,ecov=colnames(W)[j], t(tmp))

    out <- rbind(out, tmp)

  }
  OUT <- rbind(OUT, out)
  cat("Region",region,"\n")
}

str(OUT)

# Add ecov type
ECOV_info <- get_ec_info(ECOV)
OUT$Type <- ECOV_info[match(OUT$ecov, ECOV_info$name),"variable"]

head(OUT)

outfolder <- "3_ecov_varcomps/output"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

save(OUT, file=paste0(outfolder,"/ecov_varcomps_all_regions.RData"))
