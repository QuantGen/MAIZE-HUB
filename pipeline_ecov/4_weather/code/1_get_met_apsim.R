#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=met
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-136              # Array that goes from 1 to nrow(info)

rm(list=ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_ecov")

info <- data.table::fread("3_year_loc_summary/output/year_location_info.csv",
                          data.table=FALSE)

suppressPackageStartupMessages(library(apsimx))
apsimx_options(
  exe.path = "Models",
  examples.path = paste0(Sys.getenv("APSIM_HOME"), "/lib/Examples/")
)

library(apsimx)

extra <- 120  # extra days before sowing and after harvesting

met_folder <- "4_weather"

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

date_plant <- info$date_plant[job]
date_harvest <- info$date_harvest[job]
sim_period <- c(date_plant-extra, date_harvest+extra)
if(sim_period[2] > Sys.Date()) sim_period[2] <- Sys.Date()
site_coord <- c(info$longitude[job], info$latitude[job])
weather <- get_power_apsim_met(lonlat = site_coord, dates = sim_period)

# Impute if missing
if(any(is.na(weather))){
  weather <- impute_apsim_met(weather)
}

# Save weather data in .met file
write_apsim_met(weather,
                wrt.dir=paste0(met_folder,"/output"),
                filename=paste0(info$year_loc[job],".met"))
