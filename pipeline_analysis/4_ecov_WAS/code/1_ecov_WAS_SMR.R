#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=EWAS
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --array=1-4

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

source("tools/ecov_utils.R")

# Load functions
library(poolr)
library(lme4)

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

traits <- c("yield",
            "anthesis",
            "silking",
            "ASI")
region <- c("Across_regions","North","South")[3]  # Select a region
model <- c("genotype+year_loc") # Model (1) in manuscript
nPC <- 5

trait <- traits[job]

comps <- unlist(strsplit(model,"\\+"))

GPC_folder <- "1_PC_genotypes/output"
WPC_folder <- "2_PC_ecov/output"

pheno <- read.csv("../data/PHENO.csv")
ECOV <- read.csv("../data/ECOV.csv", row.names=1)

# Remove HI30 index
ECOV <- ECOV[,-grep("HI30_",colnames(ECOV))]

if(region != "Across_regions"){
  pheno <- pheno[pheno$region %in% region,]
}

length(unique(pheno$year_loc))

pheno$year <- factor(as.character(pheno$year))
pheno$location <- factor(as.character(pheno$location))
pheno$year_loc <- factor(as.character(pheno$year_loc))
pheno$genotype <- factor(as.character(pheno$genotype))

ECOV <- ECOV[levels(pheno$year_loc),]

# PC for number of independent test
load(paste0(WPC_folder,"/W_EVD.RData"), verbose=T)
index <- which(W_EVD$values>1E-8)
nt <- meff(eigen=W_EVD$values[index], method="galwey")

# Add PCs from genotypes
load(paste0(GPC_folder,"/G_PC_",region,".RData"), verbose=T)
GPC <- G_PC[as.character(pheno$genotype),1:nPC]  #Expand and select top PCs
colnames(GPC) <- paste0("GPC",1:nPC)

# Add PCs from environmental covariates
load(paste0(WPC_folder,"/W_PC_",region,".RData"), verbose=T)
WPC <- W_PC[as.character(pheno$year_loc),1:nPC]  # Expand and select top PCs
colnames(WPC) <- paste0("WPC",1:nPC)

formula0 <- paste("y~",paste(paste0("(1|",comps,")"),collapse="+"))
formula0 <- paste(c(formula0,colnames(GPC),colnames(WPC)),collapse="+")

W <- ECOV[as.character(pheno$year_loc),]
W <- scale(W)

dat <- data.frame(y=pheno[,trait],pheno[,comps],GPC,WPC,W)

out <- c()
fm0 <- lmer(formula(formula0), data=dat, REML=FALSE)  # Reduced model: without EC
for(j in 1:ncol(W)){
   formula1 <- paste0(formula0,"+",colnames(W)[j])
   fm1 <- lmer(formula(formula1), data=dat, REML=FALSE) # Full model: with EC[,j]
   res <- anova(fm0,fm1)
   tmp <- data.frame(trait, ecov=colnames(W)[j],data.frame(res, check.names=F)[2,6:8],
              bonferroni=0.05/nt,check.names=FALSE)
   rownames(tmp) <- NULL
   out <- rbind(out, tmp)

   message("ECOV ",j,". ",colnames(W)[j],": -log10(Pr)=",round(-log10(tmp[,5]),6))
}

# Add ecov type
ECOV_info <- get_ec_info(ECOV)
out$Type <- ECOV_info[match(out$ecov, ECOV_info$name),"variable"]

outfolder <- paste0("4_ecov_WAS/output/",region)
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

save(out,file=paste0(outfolder,"/ecov_WAS_results_",trait,".RData"))
