source("scripts/libraries.R")
source("scripts/Utils.R")
source("scripts/Population.R")
source("scripts/Local_Params.R")
source("scripts/Population_with_size_nmut.R")

option_list<-list(
  make_option(
    c("--obs_tum"),
    type="character",
    default = "raw/obs_tumor.RData",
    help = "path to the .RData file with the observable clones"
  ),
  make_option(
    c("--params"),
    type="character",
    default = "raw/Parameters.RData",
    help = "path to the .RData file with the elaborated parameters"
  ),
  make_option(
    c("--json_palette_file"),
    type="character",
    default = "raw/label_color.json",
    help = "path to the json file with colors for populations"
  ),
  make_option(
    c("--path_out"),
    type="character",
    default = "output",
    help = "path to the folder in which the plots and tibbles are going to be saved"
  ),
  make_option(
    c("--path_out"),
    type="character",
    default = "output",
    help = "path to the folder in which the plots and tibbles are going to be saved"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_obs_tum<-opt$obs_tum
path_params<-opt$params

path_out<-opt$path_out
if(!dir.exists(path_out)){dir.create(path_out)}

json_palette_file<-opt$json_palette_file

load(path_obs_tum)
load(path_params)
json_palette<-fromJSON(file=json_palette_file)
palette<-sapply(json_palette,function(el){el$color})

if(all(palette=="")){
    hues <- seq(0, 360, length.out = length(palette) + 1)[-1]  
    palette <- hsv(h = hues / 360, s = 0.7, v = 0.9)  
}
names(palette)<-sapply(json_palette,function(el){el$label})


Clones_df_absolute<-get_muller_plot_info(obs_Pop_ID = obs_tumor$obs_Pop_ID,
                     obs_tumor_tibble = obs_tumor$obs_tumor_tibble,
                     freq = FALSE,
                     functional_effects = parameters@functional_effects)

Clones_df_relative<-get_muller_plot_info(obs_Pop_ID = obs_tumor$obs_Pop_ID,
                                         obs_tumor_tibble = obs_tumor$obs_tumor_tibble,
                                         freq = TRUE,
                                         functional_effects = parameters@functional_effects)

plot_show_absolute<-get_muller_plot_show(Clones_df = Clones_df_absolute,
                                freq = FALSE,
                                palette = palette)

save(plot_show_absolute,Clones_df_absolute,file=paste(path_out,"Clones_df_absolute.RData",sep="/"))

ggsave(plot_show_absolute,device = "png",
       path = path_out,
       width = 9,height = 5,
       filename="plot_show_absolute.png")

max_size_reached<-max(obs_tumor$obs_tumor_tibble%>%group_by(time)%>%summarise(Ncells=sum(Ncells)))
little_label_k<-scales::scientific(max_size_reached, digits = 1)
little_label_k<-str_split(little_label_k,"e",simplify = TRUE)
little_label_k<-as.numeric(little_label_k)
names(little_label_k)<-c("base","exponent")

write(toJSON(little_label_k),file = paste(path_out,"little_label_k.json",sep="/"))


plot_show_relative<-get_muller_plot_show(Clones_df = Clones_df_relative,
                                         freq = TRUE,
                                         palette = palette)

ggsave(plot_show_relative,device = "png",
       path = path_out,
       width = 9,height = 5,
       filename="plot_show_relative.png")

plot_download_absolute<-get_muller_plot_download(Clones_df = Clones_df_absolute,
                                                 freq = FALSE,
                                                 palette = palette)

ggsave(plot_download_absolute,device = "pdf",
       path = path_out,
       width = 9,height = 5,
       filename="plot_download_absolute.pdf")
plot_download_absolute
plot_download_relative<-get_muller_plot_download(Clones_df = Clones_df_relative,
                                                 freq = TRUE,
                                                 palette = palette)

ggsave(plot_download_relative,device = "pdf",
       path = path_out,
       width = 9,height = 5,
       filename="plot_download_relative.pdf")

