library(Insite)
if (!require("ggplot2")) install.packages("ggplot2",repos = "https://cloud.r-project.org")
if (!require("stringr")) install.packages("stringr",repos = "https://cloud.r-project.org")
if (!require("colorspace")) install.packages("colorspace",repos = "https://cloud.r-project.org")

path<-"Data/Clonal_Architectures"
par_json<-list.files(path)
par_json<-par_json[grepl(x = par_json,pattern = ".json")]
seed<-read.table(paste(path,"seed.txt",sep="/"))
colnames(seed)<-c("exp","n")

for(json_file in par_json){
  dir<-str_remove(json_file,".json")
  if(!dir.exists(paste(path,dir,sep="/"))){
    dir.create(paste(path,dir,sep="/"))
  }
  if(!dir.exists(paste(path,dir,"sim1",sep = "/"))){
      system(paste(
        "Rscript scripts/run_simulation.R",
        "--seed",seed$n[seed$exp==dir],
        "--Nexp",1,
        "--params",paste(path,json_file,sep="/"),
        "--dir",paste(path,dir,sep="/")
      ))
    }
  
  system(paste("Rscript scripts/draw_plot.R --sim_dir",paste(path,dir,"sim1",sep="/"),
               "--params",paste(path,dir,"Parameters.RData",sep="/"),
               "--path_out",paste(path,dir,sep="/"),
               "--depth",3
  ))
  
  load(paste(path,dir,"abs_aboundance_muller_plot.RData",sep="/"))
  
  fun_eff_nmut<-stringr::str_count(unique(Clones_df$fun_eff), ",")+1
  colors<-vector()
  colors[fun_eff_nmut==1]<-"#C7D66D"
  colors[fun_eff_nmut==2]<-"#D58936"
  colors[fun_eff_nmut==3]<-"#AF939F"
  colors[fun_eff_nmut==4]<-"#e5d352"
  colors[fun_eff_nmut==5]<-"#28666E"
  names(colors)<-unique(Clones_df$fun_eff)
  
  colors_dark<-darken(colors,amount = 0.3)
  names(colors_dark)<-names(colors)
  muller_plot<-p[[1]]+
    scale_fill_manual(values = colors)+
    scale_color_manual(values = colors_dark)+
    theme(legend.position = "none")
  
  ggsave(filename = paste(path,dir,"abs_aboundance_muller_plot.pdf",sep="/"),
         plot = muller_plot,device = "pdf",
         width = 9,height = 5)

}

