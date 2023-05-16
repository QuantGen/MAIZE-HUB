#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=model_comps
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=160G
#SBATCH --array=1-2

# -----------------------------------------------------
# This scripts will produce inputs used for later analyses
# -----------------------------------------------------

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

# Load functions
source("tools/TENSOR_EVD.R")
source("tools/ecov_utils.R")

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "2")) # By default: South region

regions <- c("North","South")

region <- regions[job]

pheno <- read.csv("../data/PHENO.csv")
ECOV <- read.csv("../data/ECOV.csv",row.names=1)

ECOV <- ECOV[,-grep("HI30_",colnames(ECOV))]

# Folders to read genotypes and ecov PCs
GPC_folder <- "1_PC_genotypes/output"
WPC_folder <- "2_PC_ecov/output"

outfolder <- paste0("6_genomic_prediction/output/model_comps/",region)
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

pheno <- pheno[pheno$region %in% region,]

pheno$year <- factor(as.character(pheno$year))
pheno$location <- factor(as.character(pheno$location))
pheno$year_loc <- factor(as.character(pheno$year_loc))
pheno$genotype <- factor(as.character(pheno$genotype))
pheno$gen_loc <- factor(paste0(pheno$genotype,":",pheno$location))

#---------------------------------------------
message("SNP data ...")
# Read previously calculated G_PC (module 1)
load(paste0(GPC_folder,"/G_PC_",region,".RData"), verbose=T)
dim(G_PC)
SNP <- G_PC[as.character(pheno$genotype),]
SNP[1:10,1:4]
save(SNP, file=paste0(outfolder,"/SNP.RData"))

#---------------------------------------------
message("ECOV data ...")
# Read previously calculated W_PC (module 2)
load(paste0(WPC_folder,"/W_PC_",region,".RData"), verbose=T)
EC <- W_PC[as.character(pheno$year_loc),]
EC[1:10,1:4]
save(EC, file=paste0(outfolder,"/EC.RData"))

#---------------------------------------------
Z.YEAR=scale(as.matrix(model.matrix(~ year-1,data=pheno)),center=TRUE,scale=FALSE)
Z.LOC=scale(as.matrix(model.matrix(~ location-1,data=pheno)),center=TRUE,scale=FALSE)
Z.YL=scale(as.matrix(model.matrix(~ year_loc-1,data=pheno)),center=TRUE,scale=FALSE)
Z.H=scale(as.matrix(model.matrix(~ genotype-1,data=pheno)),center=TRUE,scale=FALSE)

save(Z.YEAR, Z.LOC, Z.YL, Z.H, file=paste0(outfolder,"/Z_YL.RData"))

#---------------------------------------------

# Interaction terms
message("HxLOC terms ...")
Z.HL=scale(as.matrix(model.matrix(~ gen_loc-1,data=pheno)),center=TRUE,scale=FALSE)

save(Z.HL, file=paste0(outfolder,"/Z_HL.RData"))

#---------------------------------------------
message("SNPxW terms ...")
# Load G and WWt first
load(paste0(GPC_folder,"/G.RData"), verbose=T)
load(paste0(WPC_folder,"/WWt.RData"), verbose=T)

G <- G[levels(pheno$genotype),levels(pheno$genotype)]
G <- G/mean(diag(G))

WWt <- WWt[levels(pheno$year_loc),levels(pheno$year_loc)]
WWt <- WWt/mean(diag(WWt))

tmp <- TENSOR.EVD(G, WWt, pheno$genotype, pheno$year_loc, 0.975,
                  path=paste0(getwd(),"/tools"),verbose=T)
SNPxEC <- sweep(tmp$vectors,2,sqrt(tmp$values),FUN="*")
dim(SNPxEC)
save(SNPxEC, file=paste0(outfolder,"/SNPxEC.RData"))
