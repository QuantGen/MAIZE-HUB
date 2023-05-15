
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_ecov")

pheno_folder <- "1_phenotypes"
sim_folder <- "5_APSIM"

pheno <- data.table::fread(paste0(pheno_folder,"/output/PHENO.csv"),
                           data.table=FALSE)

# Get GDD for flowering traits
# Thermal time will be read from APSIM outputs
traits <- c("anthesis","silking")

YEAR_LOC <- unique(as.character(pheno$year_loc))

# Object to store GDD traits
out <- matrix(NA, ncol=length(traits), nrow=nrow(pheno))
colnames(out) <- traits

for(i in 1:length(YEAR_LOC)){
 sim <- data.table::fread(paste0(sim_folder,"/output/simulations/",YEAR_LOC[i],".csv"),
                          data.table=FALSE)

 index <- which(pheno$year_loc==YEAR_LOC[i])
 for(j in 1:length(traits)){
   sowing <- pheno[index,"date_plant"]
   trait <- pheno[index,paste0("date_",traits[j])]
   index_sowing <- match(sowing, sim$Date)
   index_trait <- match(trait, sim$Date)
   stopifnot(all(!is.na(index_sowing)))
   stopifnot(all(!is.na(index_trait)))

   gdd <- unlist(lapply(1:length(index),function(k){
     sum(sim[seq(index_sowing[k],index_trait[k]),"TT"])
   }))
   out[index,j] <- gdd
   #stopifnot(all(sim[index_sowing,"Date"]==sowing))
   #stopifnot(all(sim[index_trait,"Date"]==trait))
 }
 message("[",i,"] Year-loc: ", YEAR_LOC[i])
}
head(out)
stopifnot(all(!is.na(out)))
apply(out,2,function(x)any(is.na(x)))

out <- data.frame(out,ASI=out[,"silking"]-out[,"anthesis"])
colnames(out) <- paste0(colnames(out),"_GDD")
# plot(pheno$anthesis,out$anthesis_GDD)
# plot(pheno$silking,out$silking_GDD)
# plot(pheno$ASI,out$ASI_GDD)

# Add the GDD flowering traits after ASI (days)
tmp <- which(colnames(pheno)=="ASI")
pheno <- data.frame(pheno[,1:tmp],out,pheno[,(tmp+1):ncol(pheno)])

# Replace PHENO file
write.csv(pheno, paste0(pheno_folder,"/output/PHENO.csv"),
          row.names=FALSE)
