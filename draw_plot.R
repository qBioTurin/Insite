library(Insite)
if (!require("optparse")) install.packages("optparse",repos = "https://cloud.r-project.org")
if (!require("ggplot2")) install.packages("ggplot2",repos = "https://cloud.r-project.org")
#if (!require("patchwork")) install.packages("patchwork",repos = "https://cloud.r-project.org")
#if (!require("colorspace")) install.packages("colorspace")

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

p <- get_muller_plot(path_sim=path_sim,depth = depth,parameters = parameters)

plot_name<-paste(ifelse(relative,"rel","abs"),"_aboundance_muller_plot.pdf",sep="")

ggsave(filename = plot_name,plot = p,path = path_out,device = "pdf",scale = 2)

tree_plot<-get_tree_plot(path_sim = path_sim,
                                  depth=depth,
                         parameters = parameters)
ggsave(filename = "tree_plot.pdf",plot = tree_plot,path = path_out,device = "pdf",scale = 1.5)
