library(Insite)
library(stringr)
library(colorspace)

tum_types<-c("AML","breast","kidney","lung","mesothelioma","uveal")

col_ranges<-list(
  c("#F0A18D","#e76f51","#902C14"),
  c("#FBDEC6","#f4a261","#D2660F"),
  c("#F7EACA","#e9c46a","#A07918"),
  c("#BAD1B3","#8ab17d","#47663D"),
  c("#9EE5DD","#2a9d8f","#15514A"),
  c("#4D8EA8","#264653","#13242A"))
names(col_ranges)<-tum_types

for(tum_type in tum_types){
  path<-paste("Data/best_fit/",tum_type,sep="")
  
  best_sim<-as.character(read.table(paste(path,"n_sim_centroid.txt",sep="/")))
  
  if(!best_sim%in%list.dirs(path,recursive = FALSE,full.names = FALSE)|
     length(list.files(paste(path,best_sim,sep="/")))==0){
    sim_n<-as.numeric(str_remove(best_sim,pattern = "sim"))
    seed<-as.numeric(read.table(paste(path,"seed.txt",sep="/")))
    
    system(paste("Rscript scripts/run_simulation.R --seed",seed,
                 "--Nexp",sim_n,
                 "--params",paste0(path,"/params.json"),
                 "--dir",path))
  }
  
  load(paste(path,"Parameters.RData",sep="/"))
  obs_tum<-Insite:::get_obs_tum(path_sim = paste(path,best_sim,sep="/"),
                                depth = 10^(-3),
                                parameters = parameters)
  
  Clones_df<-Insite:::get_muller_plot_info(obs_Pop_ID = obs_tum$obs_Pop_ID,
                                obs_tumor_tibble = obs_tum$obs_tumor_tibble,
                                functional_effects = parameters@functional_effects,
                                freq = FALSE)
  
  n<-length(unique(Clones_df$fun_eff))
  palette<-colorRampPalette(col_ranges[[tum_type]])(n)
  palette_dark<-darken(as.character(palette),amount = 0.3)
  names(palette)<-unique(Clones_df$fun_eff)
  
  Insite:::get_muller_plot_download(Clones_df = Clones_df,
                              freq = FALSE,
                              palette = palette)
  
  get_tree_plot(path_sim = paste(path,best_sim,sep="/"),depth = 10^(-3),parameters = parameters)
}


