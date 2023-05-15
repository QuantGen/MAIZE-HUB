#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=pheno
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --mem-per-cpu=32G

rm(list=ls())

library(tidyverse)
#library(lubridate)

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_ecov")

source('tools/Functions.R')
source('tools/read_phenotype.R')
source('tools/read_metadata.R')

outfolder <- "1_phenotypes"

## Combine phenotype with metadata
pheno <- read_phenotype(
                pheno_folder="source/Phenotype",
                factors_list=NULL,
                dates_list=NULL,
                traits_list=NULL
         )

head(pheno)

meta <- read_metadata(
              metadata_folder="source/Metadata",
              agronomic_folder="source/Agronomic_information",
              factors_list=NULL,
              latlon_list=NULL,       # field latitude and longitude names
              latlon_field_list=NULL, # latitude and longitude names
              location_changes=list(MOH1="MOH1-rep[1-2]"),
              verbose=2
        )
head(meta)

# Change_names
pheno$genotype <- toupper(gsub(" $","",pheno$genotype))
pheno$genotype <- gsub("^\\?+","",pheno$genotype)

tmp <- intersect(pheno$year, meta$year)
pheno <- pheno[pheno$year %in% tmp,]
meta <- meta[meta$year %in% tmp,]

# Year-locations that are not in metadata file
unique(pheno$year_loc[!pheno$year_loc %in% meta$year_loc])

# Make some changes in names
pheno$location <- gsub("NYS1","NYH1",pheno$location)
pheno$location <- gsub("IAH1[a-z]","IAH1",pheno$location) # Make experiments a,b,c in IAH1 a single one

pheno$year_loc <- paste0(pheno$year,"-",pheno$location)

# Add flowering traits
pheno$anthesis <- as.numeric(pheno$date_anthesis - pheno$date_plant)
pheno$silking <- as.numeric(pheno$date_silking - pheno$date_plant)
pheno$ASI <- pheno$silking-pheno$anthesis

YEAR_LOC <- unique(pheno$year_loc)
length(YEAR_LOC)

# Metadata object to append to pheno
meta0 <- data.frame(matrix(NA,ncol=5,nrow=nrow(pheno)))
colnames(meta0) <- c("city","latitude","longitude","irrigated","region")

for(i in 1:length(YEAR_LOC)){
  index <- names(unlist(sapply(meta$year_loc,function(x)grep(x, YEAR_LOC[i]))))
  if(length(index)==0){
    cat("i=",i,". Year-loc=",YEAR_LOC[i],"could not be found in metadata\n")
  }else{
    stopifnot(length(index)==1)
    meta0[which(pheno$year_loc==YEAR_LOC[i]), ] <-
            meta[which(meta$year_loc==index), colnames(meta0)]
  }
}
any(is.na(meta0))
pheno <- data.frame(pheno, meta0)

## NA check in date_related columns
sum(is.na(pheno$anthesis))
sum(is.na(pheno$silking))

# number of year_loc in the raw phenotype
length(unique(pheno$year_loc))

## NA removal date_related columns
## remove NA observations in any of four date columns
dim(pheno)
pheno <- pheno%>%
  drop_na(date_plant, date_harvest, date_anthesis, date_silking)
dim(pheno)

## Date order check
pheno <- pheno%>%
  mutate(h_p = date_harvest - date_plant,
         a_p = date_anthesis - date_plant,
         s_p = date_silking - date_plant,
         h_a = date_harvest - date_anthesis,
         h_s = date_harvest - date_plant)%>%
  filter(h_p > 0)%>%
  filter(s_p > 0)%>%
  filter(h_s > 0)
dim(pheno)

## Planting/harvesting interval check
##check if planting/harvesting was done in a short interval
YEAR_LOC <- unique(pheno$year_loc)
summary_date <- data.frame(year_loc=YEAR_LOC, interval_plant=NA, interval_harvest=NA)

for (i in 1:length(YEAR_LOC)){
  data <- pheno%>%
    filter(year_loc == summary_date$year_loc[i])
  summary_date[summary_date$year_loc==YEAR_LOC[i],"interval_plant"] = max(data$date_plant)-min(data$date_plant)
  summary_date[summary_date$year_loc==YEAR_LOC[i],"interval_harvest"] = max(data$date_harvest)-min(data$date_harvest)
}

problematic_year_loc <- summary_date%>%
  filter(interval_plant>10 | interval_harvest >10)
print(problematic_year_loc)

# take a look at problematic loc-year raw data
prob_year_loc_pheno <- pheno%>%
  filter(year_loc%in%problematic_year_loc$year_loc)
dim(prob_year_loc_pheno)

## Problematic loc-year

### Effect of planting date
# check problem in date_plant
table(prob_year_loc_pheno$year_loc, prob_year_loc_pheno$date_plant)

### Effect of harvesting date
table(prob_year_loc_pheno$year_loc, prob_year_loc_pheno$date_harvest)

# **Conclusion: need to separate 2019-NEH1**
# **Conclusion: the third harvesting date of 2016-NYH2 can be separated or excluded**
# **Conclusion: need to separate 2021-SCH1**
# **Conclusion: exclude plots harvested in 2021-8-13**

## Separate problematic loc-year

# separated 2019-NEH2
pheno[pheno$year_loc=="2019-NEH2" & pheno$date_plant=="2019-06-08","location"] <- "NEH2L"

#separate 2021-SCH1
pheno[pheno$year_loc=="2021-SCH1" & pheno$date_plant=="2021-05-02","location"] <- "SCH1L"

#separate 2016-NYH2
pheno[pheno$year_loc=="2016-NYH2" & pheno$date_harvest=="2016-12-08","location"] <- "NYH2L"

# exclude plot harvest in 2021-08-13 for 2021-TXH3
pheno <- pheno[-which(pheno$year_loc=="2021-TXH3"&pheno$date_harvest=="2021-08-13"),]

pheno$year_loc <- paste0(pheno$year,"-",pheno$location)
length(unique(pheno$year_loc))

## Yield outlier removal
# create a function to calculate lower (Q1-2.5IQR) and upper boundary (Q3+2.5IQR)
lower <- function(i){
  iqr <- quantile(i,0.75,na.rm = TRUE)-quantile(i,0.25,na.rm = TRUE)
  quantile(i,0.25,na.rm = TRUE)-2.5*iqr
}
upper <- function(i,na.rm=TRUE){
  iqr <- quantile(i,0.75,na.rm = TRUE)-quantile(i,0.25,na.rm = TRUE)
  quantile(i,0.75,na.rm = TRUE)+2.5*iqr
}

# yield outlier removal
table_range <- data.frame(matrix(ncol =3, nrow = length(unique(pheno$year_loc))))
colnames(table_range) <- c('year_loc', 'low', 'high')
table_range$year_loc=unique(pheno$year_loc)
for(i in 1:nrow(table_range)){
  df <- pheno%>%
    filter(year_loc==table_range$year_loc[i])%>%
    select(yield)%>%
    mutate(yield = as.numeric(yield))%>%
    as.vector()
  table_range[i, "low"]=lower(df)
  table_range[i,"high"]=upper(df)
}
dim(table_range)

pheno$outlier <- "no"
for(i in 1:nrow(table_range)){
  tmp <- table_range[i, ]
  index <- which(pheno$year_loc == tmp$year_loc)
  outliers <- which(pheno[index,"yield"]<tmp$low | pheno[index, "yield"]>tmp$high)
  if(length(outliers)>0){
    pheno[index[outliers],"outlier"] <- "yes"
  }
}
table(pheno$outlier)
range(pheno[pheno$outlier=="no","yield"])
range(pheno[pheno$outlier=="yes","yield"])

pheno <- pheno[which(pheno$outlier == "no"),]
pheno <- pheno[which(pheno$irrigated=="no"),]
dim(pheno)

pheno_clean <- pheno%>%
  mutate(plot_area_m2 = plot_area*0.0929)%>%
  select(-h_p, -s_p, -h_s, -a_p, -h_a, -plot_area, -outlier)

#summary(pheno_clean)
dim(pheno_clean)

# QC using coeficicient of variation
# Remove year-locations for which a trait is a constant
traits <- c("yield","anthesis","silking","ASI")
CV <- do.call(rbind,lapply(split(pheno_clean, pheno_clean$year_loc),function(x){
  abs(apply(x[,traits],2, sd)/apply(x[,traits],2, mean))
}))
drop <- which(apply(CV,1,function(x)any(is.na(x) | x<1E-5)))
if(length(drop)>0){
  cat("Year-location removed:\n")
  cat(paste(names(drop),collapse=","),"\n")
  pheno_clean <- pheno_clean[!pheno_clean$year_loc %in% names(drop),]
}

# QC on ASI
pheno_clean <- pheno_clean[pheno_clean$ASI > -15 & pheno_clean$ASI<15,]

write.csv(pheno_clean, paste0(outfolder,"/output/PHENO.csv"),
          row.names = FALSE)
