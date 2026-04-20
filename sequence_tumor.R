library(Insite)
if (!require("optparse")) {
  install.packages("optparse",repos = "https://cloud.r-project.org")
  library(rjson)
}
if (!require("stringr")) {
  install.packages("stringr",repos = "https://cloud.r-project.org")
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
    c("--params"),
    type="character",
    default = "raw/Parameters.RData",
    help = "path of the .RData file with the elaborated parameters"
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
  ),
  make_option(
    c("--nregions"),
    type="numeric",
    default = 1,
    help = "Number of regions of the sequencing, if multiregional or bluk (in such case 1)"
  ),
  make_option(
    c("--ncells"),
    type="numeric",
    default = 1000,
    help = "Number of cells to be sequenced in each region, use Inf for sequencing the whole mass"
  ),
  make_option(
    c("--neighborhood"),
    type="character",
    default = NULL,
    help = "If the neighborhood for the simulation has already been computed, indicate here the path to the file (default name Clones_ordered_day.RData).
    The most expensive part of the sequencing is the computation of these neighborhood, which once done is not required for repeating the sequencing"
  ),
  make_option(
    c("--seed"),
    type="numeric",
    default = NULL,
    help = "To retrieve the same sequencing, fix the seed. Default is not fixed."
  ),
  make_option(
    c("--repeat"),
    type="numeric",
    default = 1,
    help = "How many times the sequencing should be repeated (produce a file per repetition)"
  ),
  make_option(
    c("--dens"),
    type="character",
    default = "Data/dens.RData",
    help = "path to the .RData file with the DP probability density to be used in the sequencing for coverage"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_sim <- opt$sim_dir
path_params <- opt$params
path_out <- opt$path_out
seq_day<-opt$seq_day
n_regions<-opt$nregions
n_seq_cells<-opt$ncells
Clones_ordered_path<-opt$neighborhood
seed<-opt$seed
Nrep<-opt$`repeat`
dens_path<-opt$dens

load(path_params)

sim_files <- list.files(path_sim)
sim_files <- sim_files[grepl("Zprovv", sim_files)]

index_files <- as.numeric(str_remove(str_remove(sim_files, "Zprovv"), ".RData"))

if (seq_day == Inf) {
  Zprovv_file <- sim_files[which.max(index_files)]
} else {
  Zprovv_file <- sim_files[
    index_files == which.min(abs(parameters@print_time - seq_day))
  ]
}

load(file.path(path_sim, Zprovv_file))

load(dens_path)

if(is.null(Clones_ordered_path)){
  Clones_df<-Insite:::get_ordered_clones_sequencing(Zprovv)
  save(Clones_df,file = paste0(path_out,"/Clones_ordered_",seq_day,".RData"))
}else{
  load(Clones_ordered_path)
}

vcf_list <- sequencing(
  Zprovv = Zprovv,
  Clones_df = Clones_df,
  parameters = parameters,
  seed = seed,
  n_regions = n_regions,
  n_seq_cells = n_seq_cells,
  Nrep = Nrep,
  dens = dens
)

lapply(seq_along(vcf_list), function(i) {
  write.table(
    vcf_list[[i]][[1]],
    row.names = FALSE,
    file = file.path(path_out, paste0("vcf", i, ".txt"))
  )
})
