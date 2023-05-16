#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=LYO_CV
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --array=1-96

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

# Load functions
library(BGLR)

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

years <- as.character(2014:2021)
traits <- c("yield","anthesis","silking","ASI")

models <- c("YEAR+LOC+YL+H+HL",
            "SNP+EC+SNPxEC",
            "YEAR+LOC+YL+H+HL+GPC+WPC")

region <- c("North","South")[2]  # Select a region
nPC <- 5

JOBS <- expand.grid(year=years,trait=traits,model=models)
model <- as.vector(JOBS[job,"model"])
trait <- as.vector(JOBS[job,"trait"])
year <- as.vector(JOBS[job,"year"])

pheno <- read.csv("../data/PHENO.csv")

table(pheno$year)
pheno <- pheno[pheno$region %in% region,]

if(year %in% unique(pheno$year)){
  tst <- which(pheno$year == year)
  trn <- which(pheno$year != year)
}else{
  stop(" Year ",year," was not found in region ",region)
}
length(trn); length(tst)

comps <- unlist(strsplit(model,"\\+"))

outfolder <- "6_genomic_prediction/output"

if(model %in% c("YEAR+LOC+YL+H+HL","YEAR+LOC+YL+H+HL+GPC+WPC")){
  load(paste0(outfolder,"/model_comps/",region,"/Z_YL.RData"), verbose=T)
  load(paste0(outfolder,"/model_comps/",region,"/Z_HL.RData"), verbose=T)

  ETA <- list(
           YEAR=list(X=Z.YEAR[trn,],model='FIXED',saveEffects=FALSE),
           LOC=list(X=Z.LOC[trn,],model='BRR',saveEffects=FALSE),
           YL=list(X=Z.YL[trn,],model='BRR',saveEffects=FALSE),
           H=list(X=Z.H[trn,],model='BRR',saveEffects=FALSE),
           HL=list(X=Z.HL[trn,],model='BRR',saveEffects=FALSE)
         )

  if("GPC" %in% comps){ # Add PCs
    load(paste0(outfolder,"/model_comps/",region,"/SNP.RData"), verbose=T)
    GPC <- SNP[,1:nPC]
    ETA$GPC <- list(X=GPC[trn,],model='FIXED',saveEffects=FALSE)
  }
  if("WPC" %in% comps){ # Add PCs
    load(paste0(outfolder,"/model_comps/",region,"/EC.RData"), verbose=T)
    WPC <- EC[,1:nPC]
    ETA$WPC <- list(X=WPC[trn,],model='FIXED',saveEffects=FALSE)
  }
}
if(model=="SNP+EC+SNPxEC"){
  load(paste0(outfolder,"/model_comps/",region,"/SNP.RData"), verbose=T)
  load(paste0(outfolder,"/model_comps/",region,"/EC.RData"), verbose=T)
  load(paste0(outfolder,"/model_comps/",region,"/SNPxEC.RData"), verbose=T)

  ETA <- list(
           SNP=list(X=SNP[trn,],model='BRR',saveEffects=FALSE),
           EC=list(X=EC[trn,],model='BRR',saveEffects=FALSE),
           SNPxEC=list(X=SNPxEC[trn,],model='BRR',saveEffects=FALSE)
         )
}
lapply(ETA, function(x)dim(x$X))

# In the manuscript we used nIter=12000 and burnIn=4000
nIter <- 1000
burnIn <- 400

outfolder2 <- paste0(outfolder,"/LYO_CV/",model,"/",trait,"/",region)
if(!file.exists(outfolder2)) dir.create(outfolder2,recursive=TRUE)

y <- pheno[,trait]

tmp <- matrix(NA,ncol=length(comps)+1,nrow=length(tst))
colnames(tmp) <- c("INTERCEPT",comps)
out <- data.frame(pheno[tst,c("year_loc","year","location","genotype")],
                  y=y[tst],tmp, check.names=FALSE)

cat("Starting the model ...\n")
fm <- BLRXy(y=y[trn], ETA=ETA, nIter=nIter, burnIn=burnIn,
            saveAt=paste0(outfolder2,"/",year,"_"))

names(fm$ETA) <- names(ETA)

compNames <- ifelse(comps %in% c('YEAR','LOC','YL','H','HL'),paste0("Z.",comps),comps)

# Compute Effects
out[,"INTERCEPT"] <- rep(fm$mu,length(tst))
for(k in seq_along(comps)){
  comp <- comps[k]
  b <- fm$ETA[[comp]]$b
  run <- paste0("eff <- ",compNames[k], "[tst,] %*% b")
  eval(parse(text=run))
  out[,comp] <- as.vector(eff)
}

comps2 <- comps[!comps %in% c("YEAR")]
out$yHat <- apply(out[,c("INTERCEPT",comps2)],1,sum)

cat("YEAR=",year,"\n")

# cor(out$y, out$yHat)
# plot(out$y, out$yHat)

save(out, file=paste0(outfolder2,"/results_",year,".RData"))
unlink(paste0(outfolder2,"/",year,"_*.dat"))
