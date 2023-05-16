#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=ANOVA
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --array=1-8

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

# Load functions
library(BGLR)
source("tools/get_variance.R")

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

traits <- c("yield","anthesis","silking","ASI")
models <- c("YEAR+LOC+YL+H+HL",
            "SNP+EC+SNPxEC")

region <- c("North","South")[2] # Select a region

JOBS <- expand.grid(trait=traits,model=models)
model <- as.vector(JOBS[job,"model"])
trait <- as.vector(JOBS[job,"trait"])

pheno <- read.csv("../data/PHENO.csv")

pheno <- pheno[pheno$region %in% region,]

length(unique(pheno$year_loc))

outfolder <- "6_genomic_prediction/output"

if(model == "YEAR+LOC+YL+H+HL"){
  load(paste0(outfolder,"/model_comps/",region,"/Z_YL.RData"), verbose=T)
  load(paste0(outfolder,"/model_comps/",region,"/Z_HL.RData"), verbose=T)

  ETA <- list(
           YEAR=list(X=Z.YEAR,model='FIXED',saveEffects=TRUE),
           LOC=list(X=Z.LOC,model='BRR',saveEffects=TRUE),
           YL=list(X=Z.YL,model='BRR',saveEffects=TRUE),
           H=list(X=Z.H,model='BRR',saveEffects=TRUE),
           HL=list(X=Z.HL,model='BRR',saveEffects=TRUE)
         )
}
if(model=="SNP+EC+SNPxEC"){
  load(paste0(outfolder,"/model_comps/",region,"/SNP.RData"), verbose=T)
  load(paste0(outfolder,"/model_comps/",region,"/EC.RData"), verbose=T)
  load(paste0(outfolder,"/model_comps/",region,"/SNPxEC.RData"), verbose=T)

  ETA <- list(
           SNP=list(X=SNP,model='BRR',saveEffects=TRUE),
           EC=list(X=EC,model='BRR',saveEffects=TRUE),
           SNPxEC=list(X=SNPxEC,model='BRR',saveEffects=TRUE)
         )
}

lapply(ETA, function(x)dim(x$X))

# In the manuscript we used nIter=35000 and burnIn=5000
nIter <- 2000
burnIn <- 500

y <- scale(pheno[,trait])

outfolder2 <- paste0(outfolder,"/ANOVA/",model,"/",trait,"/",region)
if(!file.exists(outfolder2)) dir.create(outfolder2,recursive=TRUE)

fm <- BLRXy(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn, saveAt=paste0(outfolder2,"/"))

names(fm$ETA) <- names(ETA)

save(fm, file=paste0(outfolder2,"/fm.RData"))

comps <- unlist(strsplit(model,"\\+"))
compNames <- ifelse(comps %in% c('YEAR','LOC','YL','H','HL'),paste0("Z.",comps),comps)

names(comps) <- compNames

out <- get_variance_components(fm, comps, nIter=nIter, burnIn=burnIn)
VC <- out$VC
save(VC, file=paste0(outfolder2,"/VC.RData"))
