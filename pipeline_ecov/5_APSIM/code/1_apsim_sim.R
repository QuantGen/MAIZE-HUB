
rm(list = ls())

setwd("/mnt/research/quantgen/projects/G2F/resource_paper/pipeline_ecov")

suppressPackageStartupMessages(library(apsimx))
apsimx_options(
  exe.path = "Models",
  examples.path = paste0(Sys.getenv("APSIM_HOME"), "/lib/Examples/"),
  warn.find.apsimx = FALSE
)

# Get the functions that are needed
source('tools/APSIM_functions.R')
library(apsimx)

info <- data.table::fread("3_year_loc_summary/output/year_location_info.csv",
                          data.table=FALSE)

sim_folder <- "5_APSIM"

args=(commandArgs(TRUE))
if(length(args)==0){
  stop('No args provided')
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

dat <- info[job,]

# Parameters
parameters <- list(
    sim_dir = paste0(sim_folder,'/output/apsim_model'), #  Will be created if it does not exist
    sim_name = paste0(dat$year_loc),#'mysite_example', # This is the main file name name
    lat = dat$latitude,    # latitude for the location
    lon = dat$longitude,   # longitude for the location
    plant_date = format(strptime(dat$date_plant, format="%Y-%m-%d"), usetz=FALSE),    # planting date
    flower_date = format(strptime(dat$date_silking, format="%Y-%m-%d"), usetz=FALSE), # either provide this (the function will iterate to reach this date) or provide ggd_juv and ggd_gf
    weather_file = paste0(dat$year_loc,'.met'), # if this file does not exist, it will be downloaded
    weather_dir = "4_weather/output/",
    fertilize_NO3N = 200, # This is kg/ha of nitrate nitrogen
    plant_population = dat$plant_density # This is plant population per square meter
    # it is possible to add:
    # weather_days_before = 120 # this is the number of days before sowing to start the clock and get weather data
    # weather_dir = '' # By default is the same as sim_dir but it can be changed
    # ggd_juv = 200 # Growing degree days from emergence to end of juvenile stage (flowering initiation)
    # ggd_gf = 550 # Growing degree days from start of grain filling to end of grain filling
)

# Run the simulation
tmp <- apsim_simulate(parameters, trim=FALSE)

outfolder <- paste0(sim_folder,"/output/simulations")
if(!file.exists(outfolder)) dir.create(outfolder)

write.csv(tmp, paste0(outfolder,"/",dat$year_loc,".csv"),
          row.names=FALSE)
