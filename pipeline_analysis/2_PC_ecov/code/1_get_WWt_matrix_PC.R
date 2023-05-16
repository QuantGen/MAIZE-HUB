#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=W_PC
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

pheno <- read.csv("../data/PHENO.csv")
ECOV <- read.csv("../data/ECOV.csv", row.names=1)

# Remove HI30 ecov
ECOV <- ECOV[,-grep("HI30_",colnames(ECOV))]

outfolder <- "2_PC_ecov/output"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

WWt <- tcrossprod(scale(ECOV))
WWt <- WWt/mean(diag(WWt))

W_EVD <- eigen(WWt)
W_EVD$vectors[1:10,1:5]
rownames(W_EVD$vectors) <- rownames(WWt)
index <- which(W_EVD$values>1E-8)
W_PC <- sweep(W_EVD$vectors[,index], 2, sqrt(W_EVD$values[index]),FUN='*')
W_PC[1:10,1:6]
dim(W_PC)
# 136 132

save(WWt, file=paste0(outfolder,"/WWt.RData"))
save(W_PC, file=paste0(outfolder,"/W_PC.RData"))
save(W_EVD, file=paste0(outfolder,"/W_EVD.RData"))

# Stratified within region
regions <- unique(pheno$region)

for(k in 1:length(regions)){
  region <- regions[k]
  message("Region=",region)
  pheno0 <- pheno[pheno$region %in% region,]
  pheno0$year_loc <- factor(as.character(pheno0$year_loc))
  WWt0 <- WWt[levels(pheno0$year_loc),levels(pheno0$year_loc)]
  WWt0 <- WWt0/mean(diag(WWt0))

  W_EVD <- eigen(WWt0)
  rownames(W_EVD$vectors) <- rownames(WWt0)
  index <- which(W_EVD$values>1E-8)
  W_PC <- sweep(W_EVD$vectors[,index], 2, sqrt(W_EVD$values[index]),FUN='*')

  message("  nPC=",ncol(W_PC))
  save(W_PC, file=paste0(outfolder,"/W_PC_",region,".RData"))
  save(W_EVD, file=paste0(outfolder,"/W_EVD_",region,".RData"))
}
