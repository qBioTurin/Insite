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
    help = "folder in which the output plot files are stored"
  ),
  make_option(
    c("--depth"),
    type="numeric",
    default = 3,
    help = "only clones reaching at least 1/10^depth within their lifetime prevalence are retrieved"
  ),
  make_option(
    c("--relative"),
    type="logical",
    default = FALSE,
    help = "Should the plot be drawn with relative (TRUE) or absolute (FALSE) aboundance as y-axis"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_sim <- opt$sim_dir
path_params <- opt$params
path_out <- opt$path_out
relative<-opt$relative

if(!dir.exists(path_out)){dir.create(path_out)}

depth<-10^{-as.numeric(opt$depth)}

load(path_params)

obs_tumor<-get_obs_tum(path_sim = path_sim,depth = depth,parameters = parameters)


Clones_df<-get_muller_plot_info(obs_Pop_ID = obs_tumor$obs_Pop_ID,
                                obs_tumor_tibble = obs_tumor$obs_tumor_tibble,
                                freq = relative,
                                functional_effects = parameters@functional_effects)


hues <- seq(0, 360, length.out = length(unique(Clones_df$fun_eff)) + 1)[-1]  
col_palette <- hsv(h = hues / 360, s = 0.3, v = 0.7) 
names(col_palette)<-unique(Clones_df$fun_eff)

p<-get_muller_plot_download(Clones_df = Clones_df,
                         freq = relative,
                         palette = col_palette
                         )

p <- p + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

plot_name<-paste(ifelse(relative,"rel","abs"),"_aboundance_muller_plot.pdf",sep="")

save(Clones_df,p,file = paste(path_out,"/",ifelse(relative,"rel","abs"),"_aboundance_muller_plot.RData",sep=""))
ggsave(filename = plot_name,plot = p,path = path_out,device = "pdf",scale = 2)


fun_eff<-parameters@functional_effects
base_colors <- c(
  "growth"      = "#C2D4DD",
  "mutation"    = "#E4DDC3",
  "competition" = "#E6D1D1",
  "space"       = "#AEAEC2",
  "passenger"   = "#B6CBBC"
)

jitter_color <- function(hex_color) {
  rgb_obj <- hex2RGB(hex_color)
  hls_obj <- as(rgb_obj, "HLS")
  
  coords <- hls_obj@coords
  coords[1, "L"] <- max(0, min(1, coords[1, "L"] + runif(1, -0.06, 0.06)))
  coords[1, "S"] <- max(0, min(1, coords[1, "S"] + runif(1, -0.04, 0.04)))
  hls_mod <- HLS(coords[1, "H"], coords[1, "L"], coords[1, "S"])
  return(hex(hls_mod))
}

clamp <- function(x, min, max) pmin(pmax(x, min), max)

palette_tree <- sapply(fun_eff, function(x) jitter_color(base_colors[x]))

names(palette_tree) <- names(fun_eff)

tree_plot<-get_tree_plot_download(path_sim = path_sim,
                                  depth=depth,palette = palette_tree,parameters = parameters)
ggsave(filename = "tree_plot.pdf",plot = tree_plot,path = path_out,device = "pdf",scale = 1.5)
