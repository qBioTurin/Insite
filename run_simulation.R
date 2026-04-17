library(Insite)
if (!require("optparse")) install.packages("optparse",repos = "https://cloud.r-project.org")
if (!require("rjson")) install.packages("rjson",repos = "https://cloud.r-project.org")

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

simulation(Nexp=Nexp,
            seed=opt$seed,
            path = path,
            starting_gen = starting_gen,
            starting_fun_eff=starting_fun_eff,
            tmax=tmax,
            Ncellsmax = Ncellsmax,
            Ncells_start = Ncells_start,
            parameters = parameters,
            epsilon_rel=10^(-3))

