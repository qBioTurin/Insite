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
    c("--num_seq"),
    type="numeric",
    default = 1,
    help = "simulation step to be sequenced"
  ),
  make_option(
    c("--path_out"),
    type="character",
    default = "output",
    help = "path to the folder in which the tibble is going to be saved"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_sim <- opt$sim_dir
path_params <- opt$params
path_out<-opt$path_out
num_seq<-opt$num_seq

if(!dir.exists(path_out)){dir.create(path_out)}

load(path_params)


load(paste(path_sim,"/Zprovv",num_seq,".RData",sep=""))
time_provv<-parameters@print_time[which.min(abs(time_provv-parameters@print_time))]


pop<-lapply(Zprovv,Population)
gen<-lapply(pop,genotype)
fun_eff_num<-lapply(pop,functional_effect)
fun_eff<-lapply(fun_eff_num,function(n){names(parameters@functional_effects)[n]})
ncells<-sapply(Zprovv,Ncells)
tot_ncells<-sum(ncells)
unique_mut_id<-lapply(gen,function(g){
  unique_id_mut<-vector()
  for(i in 1:length(g)){
    unique_id_mut<-c(unique_id_mut,paste(g[1:i],collapse="_"))
  }
  return(unique_id_mut)
})
nmut_in_gen<-lengths(gen)

all_mut<-unique(unlist(unique_mut_id))
used_nums<-as.numeric(gsub(pattern = "Mut",x = names(mut_names),replacement = ""))
mut_nums<-1:length(all_mut)
not_used_mut_nums<-mut_nums[!mut_nums%in%used_nums]
mut_names_tbl<-tibble(mut=c(mut_names,all_mut[!all_mut%in%mut_names]),
                      names=c(names(mut_names),paste0("Mut",not_used_mut_nums,sep="",recycle0 = TRUE)))

save(mut_names_tbl,file=paste(path_out,"mut_names_tbl.RData",sep="/"))
pop_nmut<-tibble(ncells=sapply(split(ncells,nmut_in_gen), sum),
       npop=table(nmut_in_gen))

write(toJSON(pop_nmut),file = paste(path_out,"seq_barplot_df.json",sep="/"))

parents<-lapply(unique_mut_id,lag)
mut_generation<-sapply(gen, function(g){1:length(g)})

composition<-tibble(mut=unique_mut_id,parents,fun_eff,ncells,mut_generation)%>%
  unnest(c(mut,parents,fun_eff,mut_generation))%>%
  group_by(mut,parents,fun_eff,mut_generation)%>%
  summarise(ncells=sum(ncells),
            frequency=ncells/tot_ncells)%>%
  ungroup()%>%
  arrange(desc(frequency))%>%
  merge(mut_names_tbl)


hist_df<-composition%>%dplyr::select(fun_eff,frequency)
write(toJSON(hist_df),file = paste(path_out,"seq_hist_df.json",sep="/"))
hist_plot<-ggplot(hist_df)+
  geom_histogram(aes(x=frequency),fill="#AFABAB")+
  theme_void()+
  xlim(-0.05,1.05)+
  theme(axis.text.x = element_text(color="black"))

ggsave(hist_plot,device = "png",
       path = path_out,
       height = 3,
       width = 4,
       filename="hist_plot.png")

hist_plot_fun_eff<-ggplot(hist_df)+
  geom_histogram(aes(x=frequency,fill=fun_eff))+
  theme_void()+
  scale_fill_manual(values=palette)+
  xlim(-0.05,1.05)+
  theme(legend.position = "none",
        axis.text.x = element_text(color="black"))
ggsave(hist_plot_fun_eff,device = "png",
       path = path_out,
       height = 3,
       width = 4,
       filename="hist_plot_fun_eff.png")

roots<-composition$mut[is.na(composition$parents)]

if(length(roots)>1){
  composition<-composition%>%
    mutate(parents=ifelse(is.na(parents),"0",parents))%>%
    bind_rows(tibble(mut="0"))
}

p<-get_tree_plot_app(composition,palette)

p
ggsave(p,device = "png",
       path = path_out,
       filename="plot_tree_sequenced.png")


table_pops<-tibble(pop_id=1:length(pop),mut=unique_mut_id,fun_eff)%>%
  unnest(c(mut,fun_eff))%>%
  merge(mut_names_tbl)%>%
  group_by(pop_id)%>%
  summarise(
    mut_names=list(names),
    fun_eff=list(fun_eff),
  )%>%
  bind_cols(tibble(ncells,nmut=lengths(gen)))%>%
  arrange(desc(ncells))%>%
  group_by(nmut)%>%
  summarise(pop_names=list(mut_names),fun_eff=list(fun_eff),ncells=list(ncells))

write(jsonlite::toJSON(table_pops,auto_unbox = FALSE),file = paste(path_out,"table_pops.json",sep="/"))

