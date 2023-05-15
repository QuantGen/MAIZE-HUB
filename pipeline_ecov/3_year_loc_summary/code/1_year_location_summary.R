
rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_ecov")

pheno_folder <- "1_phenotypes"
info_folder <- "3_year_loc_summary"

pheno <- data.table::fread(paste0(pheno_folder,"/output/PHENO.csv"),
                           data.table=FALSE)

factors <- c("year_loc","year","location","city","irrigated","region")
dates <- c("date_plant","date_anthesis","date_silking","date_harvest")
traits <- c("yield","anthesis","silking","ASI")
latlon <- c("latitude","longitude")

year_loc_sum <- do.call(rbind,lapply(split(pheno, pheno$year_loc),
  function(x){
    dat <- lapply(seq_along(dates),function(i)data.frame(mean(x[,dates[i]])))
    dat <- do.call(cbind,dat)
    names(dat) <- dates
    data.frame(x[1,c(factors,latlon)],
          dat, t(apply(x[,traits],2, mean)),
          plant_density=round(mean(x$plants_stand, na.rm=T)/mean(x$plot_area_m2),0)
      )
}))
rownames(year_loc_sum) <- NULL

# check which year-loc do not have the stand info and impute with mean
year_loc_sum$plant_density <- ifelse(is.na(year_loc_sum$plant_density),
                                     mean(year_loc_sum$plant_density,na.rm=T),
                                     year_loc_sum$plant_density)

head(year_loc_sum)

write.csv(year_loc_sum, paste0(info_folder,"/output/year_location_info.csv"),
          row.names = FALSE)
