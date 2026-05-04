library(Insite)
library(tibble)
library(dplyr)
library(stringr)
if (!require("optparse")) {
  install.packages("optparse",repos = "https://cloud.r-project.org")
  library(rjson)
}

option_list<-list(
  make_option(
    c("--sim_dir"),
    type="character",
    default = "raw/sim1",
    help = "path to the folder in which the simulation outputs is stored"
  ),
  make_option(
    c("--path_out"),
    type="character",
    default = "raw",
    help = "folder in which the output plot files are stored"
  ),
  make_option(
    c("--seq_day"),
    type="numeric",
    default = Inf,
    help = "Day in which the sequencing is to be performed. Default is the end of the simulation"
  ))

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_sim <- opt$sim_dir
seq_day<-opt$seq_day
path_out <- opt$path_out

sim_files <- list.files(path_sim)
sim_files <- sim_files[grepl("Zprovv", sim_files)]

index_files <- as.numeric(str_remove(str_remove(sim_files, "Zprovv"), ".RData"))

if (seq_day == Inf) {
  Zprovv_file <- sim_files[which.max(index_files)]
  seq_day_name<-"Final"
} else {
  Zprovv_file <- sim_files[
    index_files == which.min(abs(parameters@print_time - seq_day))
  ]
  seq_day_name<-paste0("d",seq_day)
}
load(file.path(path_sim, Zprovv_file))

nD_indices<-get_nD_index(Zprovv,10^(-2))%>%
  mutate(time=time_provv)


write.table(nD_indices,
            file = paste0(path_out,"/nD_indices_",seq_day_name,".txt"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
