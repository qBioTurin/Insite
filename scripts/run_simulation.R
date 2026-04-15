#source("scripts/install_libraries.R")
source("scripts/libraries.R")
source("scripts/Utils.R")
source("scripts/Population.R")
source("scripts/Local_Params.R")
source("scripts/Population_with_size_nmut.R")
source("scripts/Simulation_Algorithm.R")

option_list<-list(
  make_option(
    c("--seed"),
    type="integer",
    default = as.integer(Sys.time()),
    help = "Seed for the first simulation"
  ),
  make_option(
    c("--Nexp"),
    type="integer",
    default = 1,
    help = "simulation number"
  ),
  make_option(
    c("--params"),
    type="character",
    default = "params.json",
    help = "json file with parameters"
  ),
  make_option(
    c("--dir"),
    type="character",
    default = "raw",
    help = "directory for the files"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

json_data <- fromJSON(file=opt$params)
path<-opt$dir

if(!dir.exists(path)){dir.create(path)}

import_json_par(json_data,path)

load(paste(path,"/Parameters.RData",sep=""))

write(opt$seed,file=paste(path,"/seed.txt",sep=""))

Nexp<-opt$Nexp

simulazione(Nexp=Nexp,
            seed=opt$seed,
            path = path,
            starting_gen = starting_gen,
            starting_fun_eff=starting_fun_eff,
            tmax=tmax,
            Ncellsmax = Ncellsmax,
            Ncells_start = Ncells_start,
            parameters = parameters)

