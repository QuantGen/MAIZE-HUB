#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=G_PC
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G

rm(list=ls())
library(data.table)

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

pheno <- read.csv("../data/PHENO.csv")
geno <- fread("../data/GENO.csv", data.table=FALSE)
ID <- geno[,1]
geno <- as.matrix(geno[,-1])
rownames(geno) <- ID
geno[1:10,1:6]
dim(geno)

outfolder <- "1_PC_genotypes/output"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

X <- scale(geno,center=TRUE,scale=FALSE)
G <- tcrossprod(X)
G <- G/mean(diag(G))

G_EVD <- eigen(G)
G_EVD$vectors[1:10,1:5]
rownames(G_EVD$vectors) <- rownames(G)
index <- which(G_EVD$values>1E-8)
G_PC <- sweep(G_EVD$vectors[,index], 2, sqrt(G_EVD$values[index]),FUN='*')
G_PC[1:10,1:6]
dim(G_PC)
# 4372 4371

save(G, file=paste0(outfolder,"/G.RData"))
save(G_PC, file=paste0(outfolder,"/G_PC.RData"))
save(G_EVD, file=paste0(outfolder,"/G_EVD.RData"))

# Stratified within region
regions <- unique(pheno$region)

for(k in 1:length(regions)){
  region <- regions[k]
  message("Region=",region)
  pheno0 <- pheno[pheno$region %in% region,]
  pheno0$genotype <- factor(as.character(pheno0$genotype))
  G0 <- G[levels(pheno0$genotype),levels(pheno0$genotype)]
  G0 <- G0/mean(diag(G0))

  G_EVD <- eigen(G0)
  rownames(G_EVD$vectors) <- rownames(G0)
  index <- which(G_EVD$values>1E-8)
  G_PC <- sweep(G_EVD$vectors[,index], 2, sqrt(G_EVD$values[index]),FUN='*')

  message("  nPC=",ncol(G_PC))
  save(G_PC, file=paste0(outfolder,"/G_PC_",region,".RData"))
  save(G_EVD, file=paste0(outfolder,"/G_EVD_",region,".RData"))
}
