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
    help = "path to the .RData file with the elaborated parameters"
  ),
  make_option(
    c("--clones_ordered"),
    type="character",
    default = "output/Clones_ordered_muller.RData",
    help = "path to the .RData file with ordered clones"
  ),
  make_option(
    c("--plot_data"),
    type="character",
    default = "output/Clones_df_absolute.RData",
    help = "path to the .RData file with the dataset to build the figure"
  ),
  make_option(
    c("--num_seq"),
    type="numeric",
    default = 1,
    help = "simulation step to be sequenced"
  ),
  make_option(
    c("--path_out"),
    type="character",
    default = "output",
    help = "path to the folder in which the plot and vcf is going to be saved"
  ),
  make_option(
    c("--json_palette_file"),
    type="character",
    default = "raw/label_color.json",
    help = "json file with colors for populations (path)"
  ),
  make_option(
    c("--seed"),
    type="character",
    default = "raw/seed_seq.txt",
    help = "path to the .txt file with the seed to be used"
  ),
  make_option(
    c("--dens"),
    type="character",
    default = "dens.RData",
    help = "path to the .RData file with the DP probability density to be used in the sequencing for coverage"
  ),
  make_option(
    c("--mut_names_tbl"),
    type="character",
    default = "output/mut_names_tbl.RData",
    help = "path to the .RData file with the mutation names to be used in the sequencing"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_sim <- opt$sim_dir
path_params <- opt$params
path_out<-opt$path_out
path_seed<-opt$seed
path_ord_clones<-opt$clones_ordered
path_plot_data<-opt$plot_data
path_dens<-opt$dens
path_mut_names_tbl<-opt$mut_names_tbl

num_seq<-opt$num_seq
json_palette_file<-opt$json_palette_file

if(!file.exists(path_seed)|is.null(path_seed)){
  runif(1)
  seed_selected<-as.integer(Sys.time())
  write(seed_selected,file=path_seed)
}else{
  seed_selected<- as.numeric(read.table(path_seed))
}
set.seed(seed_selected)

load(path_params)

json_palette<-fromJSON(file=json_palette_file)
palette<-sapply(json_palette,function(el){el$color})
if(all(palette=="")){
  hues <- seq(0, 360, length.out = length(palette) + 1)[-1]  
  palette <- hsv(h = hues / 360, s = 0.7, v = 0.9)  
}
names(palette)<-sapply(json_palette,function(el){el$label})

load(paste(path_sim,"/Zprovv",num_seq,".RData",sep=""))

time_provv<-parameters@print_time[which.min(abs(time_provv-parameters@print_time))]

pop<-lapply(Zprovv,Population)
Pop_ID<-1:length(pop)
ncells<-sapply(Zprovv,Ncells)

gen<-lapply(pop,genotype)
fun_eff<-lapply(pop,functional_effect)
fun_eff_label<-lapply(fun_eff,function(f){names(parameters@functional_effects)[f]})
unique_mut_id<-lapply(gen,function(g){
  unique_id_mut<-vector()
  for(i in 1:length(g)){
    unique_id_mut<-c(unique_id_mut,paste(g[1:i],collapse="_"))
  }
  return(unique_id_mut)
})

load(path_ord_clones)

n_seq_cells<-round(0.1*sum(ncells))
if(n_seq_cells<1){n_seq_cells<-sum(ncells)}
seq_min_y<-runif(1,min(Clones_df$y_lower),max(Clones_df$y_upper)-n_seq_cells)
seq_max_y<-seq_min_y+n_seq_cells

Clones_df_seq<-Clones_df%>%
  filter(y_lower<seq_max_y,y_upper>seq_min_y)%>%
  merge(tibble(clone=Pop_ID,
               mut=sapply(unique_mut_id,function(muts){muts[length(muts)]}),
               fun_eff=sapply(fun_eff_label,function(funcs){funcs[length(funcs)]})))%>%
  rowwise()%>%
  mutate(Ncells_seq=round(min(y_upper,seq_max_y)-max(y_lower,seq_min_y)),
         prob=Ncells_seq/(2*(seq_max_y-seq_min_y)))%>%
  ungroup()%>%
  dplyr::select(mut,fun_eff,Ncells_seq,prob)

load(path_dens)

sample_DP<-round(sample(x = dens_coverage$x, nrow(Clones_df_seq), prob = dens_coverage$y, replace=TRUE) + rnorm(1, 0, dens_coverage$bw))
sample_DP[sample_DP<0]<-0
sample_AD<-mapply(rbinom,prob=Clones_df_seq$prob,size=sample_DP,MoreArgs = list(n=1))


vcf_sample<-Clones_df_seq%>%
  bind_cols(sample_DP=sample_DP,sample_AD=sample_AD)%>%
  filter(sample_AD>0)%>%
  mutate(VAF=sample_AD/sample_DP)%>%
  dplyr::select(-c(prob,Ncells_seq))

if(file.exists(path_mut_names_tbl)){
  load(paste(path_out,"mut_names_tbl.RData",sep="/"))
  vcf_sample<-vcf_sample%>%
    merge(mut_names_tbl)%>%
    dplyr::select(-mut)%>%
    rename("mut"="names")
}


write(jsonlite::toJSON(vcf_sample,auto_unbox = FALSE),file=paste(path_out,"vcf_sampled.json",sep="/"))

load(path_plot_data)

range_plot_zoom_x<-unique(Clones_df_absolute$time)[which(sort(unique(Clones_df_absolute$time))==time_provv)+c(-1,1)]
range_plot_zoom_y<-c(min(Clones_df_absolute$y_lower[Clones_df_absolute$time==time_provv]),
                     max(Clones_df_absolute$y_upper[Clones_df_absolute$time==time_provv]))
xmin_rect<-time_provv-diff(range_plot_zoom_x)/50
xmax_rect<-time_provv+diff(range_plot_zoom_x)/50
y_trasl<-min(Clones_df_absolute$y_lower[Clones_df_absolute$time==time_provv])

p<-plot_show_absolute+
  coord_cartesian(xlim =range_plot_zoom_x,
                  ylim = range_plot_zoom_y)+
  scale_fill_manual(values=palette)+
  geom_vline(xintercept = time_provv,color="white",alpha=0.4)+
  geom_rect(aes(xmin = xmin_rect,
                xmax = xmax_rect,
                ymin = seq_min_y+y_trasl,
                ymax=seq_max_y+y_trasl),
            fill="transparent",
            color="black",
            linetype = 2,
            linewidth = 0.5)
p
vcf_sample

ggsave(plot=p,filename = "zoom_sequence_plot.png",device = "png",width = 5,height = 5,path = path_out)
