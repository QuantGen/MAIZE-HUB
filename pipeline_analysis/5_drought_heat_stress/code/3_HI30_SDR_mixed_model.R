#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=HI_effect
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-3

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_analysis")

# Load functions
library(lme4)
library(lmtest)
source("tools/ecov_utils.R")

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "3"))  # By default run region "South"

traits <- c("yield","anthesis","silking","ASI")
regions <- c("Across_regions","North","South")

region <- regions[job]

pheno <- read.csv("../data/PHENO.csv")
ECOV0 <- read.csv("../data/ECOV_layered.csv",row.names=1)
ECOV <- read.csv("../data/ECOV.csv",row.names=1)
ECOV_KL <- read.csv("../data/ECOV_KL.csv")

# Get a period around flowering
period <- c("FlaFlw","FlwStG","StGEnG")

# HI30 index
HI30 <- rowSums(ECOV[,paste0("HI30_",period)])

# SDR index
tmp <- data.frame(ECOV[,grep("Eo_|CoverGreen_",colnames(ECOV))], ECOV0)
SI <- get_SDR(EC=tmp, EC_KL=ECOV_KL)
SDR <- rowMeans(SI[,paste0("SDR_",period)])

dat <- data.frame(year_loc=names(SDR),SDR,HI30)
dat$region <- factor(pheno[match(dat$year_loc, pheno$year_loc),"region"])

thrHI <- 36
thrSDR <- 0.54

# Levels for the indices indicating 'No stress'/'Stress'
levHI30 <- paste0("HI30",c("<=",">"),thrHI)
levSDR <- paste0("SDR",c(">","<="),thrSDR)
levSI <- c("No-Stress","Stress")     # For the combined index HI30&SDR

dat$HI30_group <- factor(paste0("HI30",ifelse(dat$HI30<=thrHI,"<=",">"),thrHI),levels=levHI30)
dat$SDR_group <- factor(paste0("SDR",ifelse(dat$SDR<=thrSDR,"<=",">"),thrSDR),levels=levSDR)

index <- dat$HI30_group==levHI30[2] & dat$SDR_group==levSDR[2]  # Both HI30&SDR stresses
dat$SI_group <- factor(ifelse(index,levSI[2],levSI[1]), levels=levSI)

if(region != "Across_regions"){
  pheno <- pheno[pheno$region %in% region,]
}else{
  pheno <- pheno[]
}

pheno$year <- factor(as.character(pheno$year))
pheno$location <- factor(as.character(pheno$location))
pheno$year_loc <- factor(as.character(pheno$year_loc))
pheno$genotype <- factor(as.character(pheno$genotype))
pheno$gen_loc <- factor(paste0(pheno$genotype,":",pheno$location))

ECOV <- ECOV[levels(pheno$year_loc),]

#formula1 <- "(1|year)+(1|location)+(1|year_loc)+(1|genotype)+(1|gen_loc)"
formula1 <- "(1|year_loc)+(1|genotype)+(1|gen_loc)"   # Simpler model

# Stress type
SI_type <- c("HI30","SDR","SI")

# k <- 4
out <- c()
levelsSI <- sapply(SI_type,function(x)NULL)
for(k in seq_along(SI_type)){
  message("Stress type: ",SI_type[k])
  pheno$SI_group <- dat[match(pheno$year_loc, dat$year_loc),paste0(SI_type[k],"_group")]
  SI_levels <- levels(dat[,paste0(SI_type[k],"_group")])
  levelsSI[[SI_type[k]]] <- SI_levels
  si0 <- ifelse(SI_type[k]=="HI30",thrHI,
                ifelse(SI_type[k]=="SDR",thrSDR,paste0(thrHI,"&",thrSDR)))

  nYL <- length(unique(pheno$year_loc))
  tmp <- table(as.character(pheno$year_loc), pheno$SI_group)
  nYL2 <- apply(tmp,2,function(x)sum(x>0))
  names(nYL2) <- paste0("nYL_",names(nYL2))

  for(tr in seq_along(traits))
  {
      trait <- traits[tr]
      fm <- lmer(formula=paste(trait,"~ SI_group+",formula1), data=pheno)
      rm <- lmer(formula=paste(trait,"~",formula1), data=pheno)

      LRT <- lrtest(rm, fm)
      Chisq <- LRT[2,"Chisq"]
      p.value.LRT <- LRT[2,"Pr(>Chisq)"]
      B <- summary(fm)$coefficients
      rownames(B) <- gsub("SI_group","",rownames(B))

      ff <- B[2,]; names(ff)[] <- c(paste0("x"),"SE","t.value")
      Intercept <- B["(Intercept)","Estimate"]
      p.value <- pnorm(abs(as.vector(ff['t.value'])), lower.tail=FALSE)

      total <- rep(Intercept,2); names(total) <- SI_levels
      total[SI_levels[2]] <- total[SI_levels[2]] + ff['x']

      tmp <- drop(sqrt(c(1,0)%*%vcov(fm)%*%c(1,0)))
      SE.gp <- rep(tmp,2); names(SE.gp) <- SI_levels
      SE.gp[SI_levels[2]] <- drop(sqrt(c(1,1)%*%vcov(fm)%*%c(1,1)))

      out0 <- data.frame(t(ff),p.value,Chisq,p.value.LRT,check.names=F)
      tmp <- out0; tmp[] <- NA
      out0 <- rbind(tmp, out0)

      out0 <- data.frame(SI_type=SI_type[k],si=si0, region,trait,group=SI_levels,
                    nYL,nYL_gp=nYL2,
                    Intercept,out0, total, SE.gp,check.names=F)
      rownames(out0) <- NULL
      message("   Trait ",trait)

    out <- rbind(out, out0)
  }
}
head(out)
levelsSI

outfolder <- "5_drought_heat_stress/output"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

save(out, levelsSI, file=paste0(outfolder,"/mean_effects_results_",region,"_all_traits.RData"))
