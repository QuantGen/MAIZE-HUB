#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=LDprune
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=08:59:00
#SBATCH --mem=180G
#SBATCH --constraint=intel16|intel18

#=========================================
# LD-pruning and imputation
#=========================================
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_ecov")

source("tools/LD_prune.R")
library(data.table)

#==========================================
# MAF and freqNA thresholds.
# SNPs are already filterd by MAF and freqNA;
# however we can apply a second round of filter
#==========================================
MAF.min <- 0.03
pNA.max <- 0.1
#==========================================

pheno_folder <- "1_phenotypes"
geno_folder <- "2_genotypes"

cat("Reading data ...\n")
pheno <- read.csv(paste0(pheno_folder,"/output/PHENO.csv"))
prefix <- "geno_numeric_MAF_3perc_filtered"
geno <- fread(paste0(geno_folder,"/output/",prefix,".vcf.012"),
              check.names=FALSE, data.table=FALSE)
geno <- as.matrix(geno[,-1])  # First column are row numbers
ID <- fread(paste0(geno_folder,"/output/",prefix,".vcf.012.indv"),
                   header=FALSE, data.table=FALSE)
MAP <- fread(paste0(geno_folder,"/output/",prefix,".vcf.012.pos"),
             data.table=FALSE)
colnames(MAP) <- c("chr","pos")
MAP <- data.frame(name=paste0("S",MAP$chr,"_",MAP$pos), MAP)
dimnames(geno) <- list(toupper(ID[,1]), MAP$name)

pheno[1:10,1:16]
geno[1:10,1:16]
dim(geno)

#=====================================
cat("Matching genotypic and phenotypic data ...\n")
#=====================================
ID <- intersect(rownames(geno), unique(pheno$genotype))

pheno <- pheno[pheno$genotype %in% ID,]
dim(pheno)

# Replace PHENO file
write.csv(pheno, paste0(pheno_folder,"/output/PHENO.csv"),
          row.names=FALSE)

geno <- geno[rownames(geno) %in% ID,]
dim(geno)

#=====================================
cat("MAF filter ...\n")
#=====================================
geno[geno == -1] <- NA

MAP$MAF <- apply(geno, 2, function(x) mean(x,na.rm=T)/2)
range(MAP$MAF)
drop <- which(MAP$MAF < MAF.min | MAP$MAF > 1-MAF.min)
if(length(drop) > 0){
  geno <- geno[,-drop]
  MAP <- MAP[-drop,]
  cat(length(drop),"markers with MAF<",MAF.min," were dropped\n")
}
dim(geno)
# [1]   4178 317135

#=====================================
cat("Missing values filter ...\n")
#=====================================
MAP$freqNA <- apply(geno, 2, function(x) mean(is.na(x)))
range(MAP$freqNA)
drop <- which(MAP$freqNA > pNA.max)
if(length(drop) > 0){
   geno <- geno[,-drop]
   MAP <- MAP[-drop,]
   cat(length(drop),"markers with pNA>",pNA.max," were dropped\n")
}
dim(geno)
# [1]   4178 316297

#=====================================
cat("LD pruning ...\n")
#=====================================
mc.cores <- 1
threshold <- c(0.75,0.80,0.85)
out <- LD_prune(geno, MAP, threshold=threshold,
                 mc.cores=mc.cores, d.max=1E6)

out2 <- data.frame(matrix(FALSE, nrow=nrow(MAP), ncol=length(threshold)))
colnames(out2) <- paste0("pruneIn",100*threshold)
out2[MAP$name %in% out$pruneIn[[1]]$NAME, 1] <- TRUE
out2[MAP$name %in% out$pruneIn[[2]]$NAME, 2] <- TRUE
out2[MAP$name %in% out$pruneIn[[3]]$NAME, 3] <- TRUE
MAP <- data.frame(MAP, out2)

#=====================================
cat("Naive imputation ...\n")
#=====================================
indexNA <- which(MAP$freqNA>0)
tmp <- seq(0.1, 1, by=0.1)
aa <- round(tmp*length(indexNA)); names(aa) <- paste0(tmp*100,"%")
if(length(indexNA) > 0){
  for(j in 1:length(indexNA)){
    tmp <- as.vector(geno[,indexNA[j]])
    geno[,indexNA[j]] <- ifelse(is.na(tmp),mean(tmp,na.rm=T),tmp)

    if(any(aa==j)) cat("  Imputation:",names(which(aa==j)),"done\n")
  }
}

geno <- geno[,out2$pruneIn85]
MAP <- MAP[out2$pruneIn85,]
dim(geno)
dim(MAP)
stopifnot(all(colnames(geno) == MAP$name))

#=====================================
cat("Saving output ...\n")
#=====================================
write.csv(geno, paste0(geno_folder,"/output/GENO.csv"), row.names=T)
write.csv(MAP, paste0(geno_folder,"/output/MAP.csv"), row.names=F)
