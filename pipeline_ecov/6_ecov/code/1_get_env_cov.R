#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=ECOV
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --mem-per-cpu=32G

rm(list=ls())

library(lubridate)

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_ecov")

source("tools/ecov_utils.R")

pheno_folder <- "1_phenotypes"
met_folder <- "4_weather"
sim_folder <- "5_APSIM"
ecov_folder <- "6_ecov"

files <- list.files(paste0(sim_folder,"/output/simulations"), pattern = ".csv")
year_loc <- gsub(".csv","",files)

DATA <- list()
for(i in seq_along(year_loc)){
  env <- read.csv(paste0(sim_folder,"/output/simulations/",files[i]), check.names=F)
  DATA[[length(DATA)+1]] <- data.frame(year_loc=year_loc[i],env, check.names=F)
}
DATA <- do.call(rbind, DATA)
DATA[1:10,1:5]

variable <- colnames(DATA)[c(4:ncol(DATA))]

DATA[,variable] <- apply(DATA[,variable],2,as.numeric)

stages <- c(Ger="Germination", Eme="Emergence", EnJ="EndJuvenile",
            Flo="FloralInitiation", Fla="FlagLeaf", Flw="Flowering",
            StG="StartGrainFill", EnG="EndGrainFill", Mat="Maturity",
            Har="HarvestRipe")

periods <- list("Ger", "Eme", c('EnJ','Flo'), 'Flo', 'Fla',
                'Flw', 'StG', c('EnG','Mat','Har'),
                c('Mat','Har'), 'Har')

period_list <- lapply(2:length(periods),function(k){
  a1 <- periods[[k-1]]
  a2 <- periods[[k]]
  list(period=paste0(a1[1],a2[1]),start=a1,end=a2)
})

# ---------------------------------------------------
datn <- c()
datKL <- c()
result <- list()
rr <- 0  # 0/1 to finish the period 0 or 1 days prior to the start to the other period
for(i in seq_along(year_loc)){
  cat("---------- Year-location(",i,"): ",year_loc[i],"----------\n")
  apsim_file <- paste0(year_loc[i],".apsimx")

  datai <- DATA[DATA$year_loc%in%year_loc[i],]
  datai$Date = as.Date(datai$Date, tryFormats = c("%Y-%m-%d", "%m/%d/%Y"))

  # Make sure the date is actually a date
  if(class(datai$Date) != 'Date') stop('Date is not actually a date')

  # Get dates for the period between emergence and end of juvenile (or floral initiation, what occurs first)
  index <- lapply(period_list,function(x){
    a1 <- min(which(datai$Stage %in% stages[x$start]))
    a2 <- min(which(datai$Stage %in% stages[x$end]))
    out <- a1
    if(a2>a1){ out <- a1:a2 }
    out
  })
  names(index) <- unlist(lapply(period_list,function(x)x$period))
  datn <- rbind(datn, data.frame(year_loc=year_loc[i],
                                 t(unlist(lapply(index, length)))))

  tmp <- unlist(lapply(index,length))
  if(any(tmp < 1)){
    message(" Some periods were not found in location=",year_loc[i],":")
    message("  ",paste(names(tmp)[tmp<1], collapse=","))
  }

  tmp <- apsimx::inspect_apsimx(apsim_file,paste0(sim_folder,"/output/apsim_model"),
                                node="Soil", soil.child="Physical")
  KL <- unlist(apsimx::extract_values_apsimx(apsim_file,paste0(sim_folder,"/output/apsim_model"),
                                             paste0(tmp,".MaizeSoil.KL")))
  LL <- unlist(apsimx::extract_values_apsimx(apsim_file,paste0(sim_folder,"/output/apsim_model"),
                                             paste0(tmp,".MaizeSoil.LL")))
  datKL <- rbind(datKL, data.frame(year_loc=year_loc[i],layer=1:length(KL),
                                   KL=KL,LL=LL))

  W0 <- datai[,variable]
  W0[,'yield'] <- W0[,'yield']*10/1000   # convert to ton/ha

  # Get heat index
  met <- apsimx::read_apsim_met(paste0(year_loc[i],".met"),
                                paste0(met_folder,"/output"),
                                verbose=FALSE)
  hi <- rep(NA,nrow(datai))
  for(k in 1:nrow(datai)){
    tmp <- which(met$year==year(datai[k,"Date"]) & met$day==yday(datai[k,"Date"]))
    met0 <- drop(as.matrix(met[tmp,]))
    hi[k] <- FtoC(get_HI(CtoF(met0['maxt']), met0['rh']))
  }
  HI <- cbind(HI30=1*(hi>30))

  # Get cumulated EC within each period
  dat <- c()
  for(j in 1:ncol(HI)){
    tmp <- unlist(lapply(index,function(ii)sum(HI[ii, j])))
    dat <- rbind(dat,  data.frame(t(tmp)))
  }
  dat <- rbind(dat, t(apply(dat,1,cumsum)))

  for(j in 1:ncol(W0)){
    tmp <- unlist(lapply(index,function(ii){
      sum(W0[ii, j])/length(ii)
    }))
    dat <- rbind(dat,  data.frame(t(tmp)))
  }
  rownames(dat) <- c(colnames(HI),paste0("Cum",colnames(HI)),colnames(W0))

  result[[length(result) + 1]] <- data.frame(year_loc=year_loc[i],
                                             ecov=rownames(dat), dat)
}
result <- do.call(rbind,result)
result[1:10,1:5]
datn[,1:5]
datKL[1:10,]

# Transpose results
ecov_names <- unique(result$ecov)
out <- do.call(rbind,lapply(split(result, result$year_loc),function(x){
  tt <- x[,-c(1:2)]
  rownames(tt) <- x$ecov
  data.frame(year_loc=x[1,"year_loc"],period=colnames(tt), t(tt[ecov_names,]), check.names=F)
}))
rownames(out) <- NULL
out[1:10,1:5]

# Reshaping output
result2 <- reshape2::melt(out, id=c("year_loc","period"))
tmp <- strsplit(as.character(result2$variable),"\\(|\\)")
result2$variable <- unlist(lapply(tmp, function(x)x[1]))
result2$layer <- unlist(lapply(tmp, function(x)x[2]))
result2$name <- paste0(result2$variable,"_",result2$period,
                ifelse(is.na(result2$layer),"",paste0("_",result2$layer)))
ecov_names <- unique(result2$name)
ECOV <- do.call(rbind, lapply(split(result2, result2$year_loc),function(x){
  rownames(x) <- x$name
  data.frame(t(x[ecov_names,"value", drop=F]))
}))
ECOV[1:10,1:5]

# Rename some ECOVs
colnames(ECOV) <- gsub("AccumulatedTT","CumTT", colnames(ECOV))

# Delete repeated variables
  # LeachNO3 is FlowNO3(10)
  # Drainage is Flux(10)
  # Water(i) is SWmm(i) and SW(i)
  # PAW(i) is PAWmm(i)
drop <- c( grep("LeachNO3",colnames(ECOV),value=T),
           grep("Drainage",colnames(ECOV),value=T),
           grep("^Water_",colnames(ECOV),value=T),
           grep("^SWmm",colnames(ECOV),value=T),
           grep("^PAW_",colnames(ECOV),value=T),
           grep("^PAWmm_",colnames(ECOV),value=T)
         )
ECOV <- ECOV[, !colnames(ECOV) %in% drop]

# Save ECOV that are layered
ecov_info <- get_ec_info(ECOV)
ECOV0 <- ECOV[,ecov_info[!is.na(ecov_info$layer),'name']]

# Aggregate soil ECOV across layers
ECOV1 <- aggregate_ec(ECOV)
dim(ECOV1)

# Quality control ECs based on coeficient of variation (CV=SD/mean)
# This will be done within region as well
pheno <- read.csv(paste0(pheno_folder,"/output/PHENO.csv"))

get_CV <- function(ecov){
  abs(apply(ecov,2,sd,na.rm=T) / apply(ecov,2,mean,na.rm=T))
}

tmp <- lapply(split(pheno,pheno$region),function(x)unique(x$year_loc))

CV <- do.call(cbind,lapply(tmp,function(x) get_CV(ECOV1[x,])))
drop <- apply(CV,1,function(x)sum(x<1E-4))
drop <- which(ifelse(is.na(drop),TRUE,drop>0))
if(length(drop)>0){
  message(length(drop)," of ",ncol(ECOV1)," ECs with a |cv|<1E-4 were removed!")
  ECOV1 <- ECOV1[,-as.vector(drop)]
}

cv <- get_CV(ECOV1)  # CV across all ECOVs
drop <- which(ifelse(is.na(cv),TRUE,ifelse(cv<1E-4,TRUE,FALSE)))
if(length(drop)>0){
  message(length(drop)," of ",ncol(ECOV1)," ECs with a |cv|<1E-4 were removed!")
  ECOV1 <- ECOV1[,-as.vector(drop)]
}

# Trim the Layered ECOVs
ecov_names <- unique(get_ec_info(ECOV1)$variable)
index <- c()
for(j in 1:length(ecov_names)){
  index <- c(index ,grep(paste0("^",ecov_names[j],"_"),colnames(ECOV0),value=T))
}
ECOV0 <- ECOV0[,index]

# Save outputs
outfolder <- paste0(ecov_folder,"/output/")
if(!file.exists(outfolder)) dir.create(outfolder)

stopifnot(all(rownames(ECOV0) == datn[,1]))
stopifnot(all(rownames(ECOV1) == datn[,1]))
write.csv(datKL, paste0(outfolder,"/ECOV_KL.csv"), row.names=FALSE)
write.csv(ECOV0, paste0(outfolder,"/ECOV_layered.csv"))
write.csv(ECOV1, paste0(outfolder,"/ECOV.csv"))

# Save period note
period_note <- do.call(rbind,lapply(period_list,function(x){
  data.frame(period=x[[1]], start=stages[x[[2]][1]], end=stages[x[[3]][1]])
}))

dat <- reshape2::melt(datn, id="year_loc")
colnames(dat)[2:3] <- c("period","ndays")

dat <- data.frame(dat,
          period_note[match(dat$period, period_note$period),c("start","end")])

write.csv(dat, paste0(outfolder,"/ECOV_period.csv"), row.names=FALSE)
