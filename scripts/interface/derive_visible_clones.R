source("scripts/libraries.R")
source("scripts/Utils.R")
source("scripts/Population.R")
source("scripts/Local_Params.R")
source("scripts/Population_with_size_nmut.R")

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
    help = "folder in which the output obs_tumor.RData and label_color.json file is stored"
  ),
  make_option(
    c("--depth"),
    type="numeric",
    default = 3,
    help = "only clones reaching at least 1/10^depth within their lifetime prevalence are retrieved"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_sim <- opt$sim_dir
path_params <- opt$params
path_out <- opt$path_out

depth<-10^{-as.numeric(opt$depth)}

load(path_params)

obs_tum<-get_obs_tum(path_sim = path_sim,depth = depth,parameters = parameters)

save(obs_tumor,file = paste(path_out,"/obs_tumor.RData",sep=""))

obs_pop_ordered<-obs_tum$obs_tumor_tibble%>%
  group_by(Population_ID)%>%
  summarise(max_ncells=max(Ncells))%>%
  merge(obs_pop_id)%>%
  arrange(desc(max_ncells))%>%
  pull(Populations)

fun_eff_labels<-unique(sapply(obs_pop_ordered,
                              function(pop){
                                paste(
                                  names(
                                    parameters@functional_effects[
                                      sort(
                                        phenotype(pop)
                                      )
                                    ]
                                  ),
                                  collapse = ", ")
                              }
)
)

json_data <- lapply(fun_eff_labels,
                    function(fun_eff_label) {
                      list(label = fun_eff_label, color = "")
                    }
)

write(toJSON(
  json_data
),paste(path_out,"/label_color.json",sep=""))

