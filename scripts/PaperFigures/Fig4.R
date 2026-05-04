library(Insite)
library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)

dir<-"Data/nD_exp/Sim_expanded"
folders<-list.dirs(dir,recursive = FALSE)

nD_all<-tibble()
for(folder in folders){
  
  json_params<-paste(folder,"params.json",sep="/")
  seed<-read.table(paste(folder,"seed.txt",sep="/"))$V1
  sim_n<-read.table(paste(folder,"sim_n.txt",sep="/"))$V1
  
  if(!dir.exists(paste(folder,paste0("sim",sim_n),sep="/"))){
    system(paste("Rscript scripts/run_simulation.R --seed",seed,
                 "--Nexp",sim_n,
                 "--params",json_params,
                 "--dir",folder)
    )
  }
  
  load(paste(folder,"Parameters.RData",sep="/"))
  obs_tum<-Insite:::get_obs_tum(path_sim = paste(folder,paste0("sim",sim_n),sep="/"),
                                depth = 10^(-3),
                                parameters = parameters)
  Clones_df<-Insite:::get_muller_plot_info(obs_Pop_ID = obs_tum$obs_Pop_ID,
                                           obs_tumor_tibble = obs_tum$obs_tumor_tibble,
                                           functional_effects = parameters@functional_effects,
                                           freq = FALSE)
  fun_eff_nmut<-stringr::str_count(unique(Clones_df$fun_eff), ",")+1
  colors<-vector()
  colors[fun_eff_nmut==1]<-"#C7D66D"
  colors[fun_eff_nmut==2]<-"#D58936"
  colors[fun_eff_nmut==3]<-"#AF939F"
  colors[fun_eff_nmut==4]<-"#e5d352"
  colors[fun_eff_nmut==5]<-"#28666E"
  names(colors)<-unique(Clones_df$fun_eff)
  
  p<-Insite:::get_muller_plot_download(Clones_df = Clones_df,freq = FALSE,palette = colors)
  
  print(p+
        ggtitle(str_remove(folder,"Data/nD_exp/Sim_expanded/"))+
          theme(plot.title = element_text(face="bold",size=14)))
  
  system(paste("Rscript scripts/derive_nD_indices.R --sim_dir",paste(folder,paste0("sim",sim_n),sep="/"),
               "--path_out",folder,
               "--seq_day",Inf)
  )
  
  nD_tum<-read.table(paste(folder,"nD_indices_Final.txt",sep="/"),header = TRUE)%>%
    mutate(par_set=str_remove(folder,"Data/nD_exp/Sim_expanded/"))
  nD_all<-bind_rows(nD_all,nD_tum)
}

nD_all<-nD_all%>%mutate(FE=c("competition","competition","null","growth"))
colors_palette<-c("Proliferation"="#DDA77B",
                  "CompetitionMild"="#C6A5A5",
                  "CompetitionStrong"="#714747",
                  "Passenger"="#A8BCDF"
)

area_grey<-tibble(
  n=seq(1,2,by=0.01)
)%>%
  mutate(D=(2-n)^(-2))%>%
  filter(D<ceiling(max(nD_all$D)+1))

mean_branch<-tibble(
  n=seq(1,ceiling(max(nD_all$n)+1),by=0.01)
)%>%
  mutate(D=9*(2*n-1)/8)

lin_tree<-tibble(
  n=seq(1,ceiling(max(nD_all$n)+1),by=0.01)
)%>%
  mutate(
    frac_n=n-floor(n),
    D=((1-frac_n)^2+frac_n^2)^(-1))

ggplot(nD_all)+
  geom_ribbon(data=area_grey,aes(x=n,ymin=D,ymax = ceiling(max(nD_all$D)+1)),fill="#DFE2E0")+
  geom_line(data=lin_tree,aes(x=n,y=D),color="#DFE2E0")+
  geom_point(aes(x=n,y=D,fill=par_set,shape=FE),
              color="white",stroke = 0.2,size=6,alpha=1)+
  scale_shape_manual(values = 21:24)+
  scale_fill_manual(values = colors_palette)+
  scale_y_continuous(limits = c(1,ceiling(max(nD_all$D)+1)))+
  scale_x_continuous(limits = c(1,ceiling(max(nD_all$n)+1)),breaks = 1:ceiling(max(nD_all$n)+1))+
  theme_minimal()+
  theme(legend.position = "none",
        panel.spacing.x = unit(30,"pt"),
        panel.grid.major = element_line(color="#DFE2E0",linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face="bold"))
