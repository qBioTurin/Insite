if (!require("rjson")) {
  install.packages("rjson",repos = "https://cloud.r-project.org")
  library(rjson)
}
library(Insite)
library(tibble)
library(stringr)
library(dplyr)

load("Data/best_fit/Noble_tum_data.RData")

dir_results<-"OutputSim/TumorFit/results"
sim_folders<-list.dirs(dir_results,recursive = FALSE)
sim_folders<-sim_folders[grepl("params",sim_folders)]
nD_all<-lapply(sim_folders, function(sim_folder){
  sims<-list.dirs(sim_folder,recursive = FALSE,full.names = FALSE)
  sims<-sims[grepl("sim",sims)]
  nD_indices<-lapply(sims,function(sim){
    Zprovv_file <- list.files(paste(sim_folder,sim,sep="/"))
    load(paste(sim_folder,sim,Zprovv_file,sep="/"))
    nD_indices<-get_nD_index(Zprovv,10^(-2))%>%
      mutate(time_end=time_provv,sim=sim)
  })%>%
    bind_rows()%>%
    mutate(
      par_set=str_remove(string = sim_folder,pattern = paste0(dir_results,"/"))
    )
})%>%
  bind_rows()

print("computed nD_all")
dir_json<-"OutputSim/TumorFit/params"
par_sets<-unique(nD_all$par_set)

cell_life<-sapply(
  par_sets,
  function(par_set){
    json_file<-paste0(dir_json,"/",par_set,".json")
    json_data <- fromJSON(file=json_file)
    return(json_data$cellLife)
  }
)
ending_size<-sapply(
  par_sets,
  function(par_set){
    json_file<-paste0(dir_json,"/",par_set,".json")
    json_data <- fromJSON(file=json_file)
    return(json_data$endingSize)
  }
)

info_tum<-tibble(tum_type=c("AML","breast", "kidney", "lung", "mesothelioma", "uveal"),
                 min_time=c(0,0,7300,0,7300,0),
                 max_time=c(3650,3650,18250,3650,18250,5475),
                 cell_life=c(4,4,100,4,30,365),
                 ending_size=c(10^9,10^8,10^8,10^8,10^8,10^7))

sim_data<-tibble()
for(tum_type in c("AML","breast", "kidney", "mesothelioma", "lung", "uveal")){
  
  print(paste("start",tum_type))
  
  real_points_tum<-real_points[grepl(tum_type,real_points$dataset),]
  
  min_time_tum<-info_tum$min_time[info_tum$tum_type==tum_type]
  max_time_tum<-info_tum$max_time[info_tum$tum_type==tum_type]
  cell_life_tum<-info_tum$cell_life[info_tum$tum_type==tum_type]
  ending_size_tum<-info_tum$ending_size[info_tum$tum_type==tum_type]
  nD_tum<-nD_all%>%
    filter(par_set%in%par_sets[(ending_size==ending_size_tum)&(cell_life==cell_life_tum)],
           time_end>min_time_tum,
           time_end<max_time_tum)
  nD_tum_nested<-nD_tum%>%
    nest_by(par_set)
  nD_tum_list<-nD_tum_nested$data
  names(nD_tum_list)<-nD_tum_nested$par_set
  Likelihood<-sapply(nD_tum_list,
                     function(last_sim_points){
                       if(mean(last_sim_points$time_end)>max_time_tum|
                          mean(last_sim_points$time_end)<min_time_tum){return(NA)}
                       else{
                         dens<-MASS::kde2d(last_sim_points$n,last_sim_points$D,
                                           lims=c(min(real_points_tum$n),
                                                  max(real_points_tum$n),
                                                  min(real_points_tum$D),
                                                  max(real_points_tum$D)),
                                           n=100,h=1)
                         Likelihood<-sum(mapply(
                           function(n_real,D_real){
                             return(log(dens$z[which.min(abs(n_real-dens$x)),which.min(abs(D_real-dens$y))]))
                           },
                           n_real=real_points_tum$n,
                           D_real=real_points_tum$D
                         ))
                         return(Likelihood)
                       }
                     })
  
  best_par_set<-names(which.max(Likelihood))
  
  centroid <- colMeans(nD_tum[nD_tum$par_set==best_par_set, c("D","n")])
  distances <- sqrt((nD_tum$D[nD_tum$par_set==best_par_set] - centroid[1])^2 + (nD_tum$n[nD_tum$par_set==best_par_set] - centroid[2])^2)
  most_central_idx <- which.min(distances)
  representative_sim<-nD_tum$sim[nD_tum$par_set==best_par_set][most_central_idx]
  
  out_dir<-paste0("Data/best_fit/",tum_type)
  
  seed_tbl<-read.table("OutputSim/TumorFit/seeds.txt")
  write.table(seed_tbl$V2[seed_tbl$V1==best_par_set],file = paste0(out_dir,"/seed.txt"),
              append = FALSE,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  file.copy(
    from = paste0(dir_json,"/",best_par_set,".json"),
    to   = paste0(out_dir,"/params.json"),
    overwrite = TRUE
  )
  write.table(x = representative_sim,file = paste0(out_dir,"/n_sim_centroid.txt"),row.names = FALSE,col.names = FALSE)
  
  sim_data<-bind_rows(sim_data,nD_tum%>%
    filter(par_set==best_par_set)%>%
    select(sim,D,n,time_end)%>%
    mutate(dataset=tum_type,centroid=(sim==representative_sim)))
  
  
  Insite:::import_json_par(json_data = rjson::fromJSON(file=paste0(dir_json,"/",best_par_set,".json")),
                           path = paste0(dir_results,"/best_sim_",tum_type))
  
  print(paste("Params.RData created for best_sim",tum_type))
  load(paste0(dir_results,"/best_sim_",tum_type,"/Parameters.RData"))
  sims<-list.dirs(paste(dir_results,best_par_set,sep="/"))
  sims<-sims[grepl("sim",sims)]
  vcf_multiregional_all<-lapply(
    sims,
    function(sim){
      load(paste(sim,list.files(sim),sep="/"))
      vcf<-sequencing(
        Zprovv = Zprovv,
        parameters = parameters,
        seed = 1,n_regions = 10,n_seq_cells = 10^4
      )[[1]][[1]]
      
      return(vcf%>%mutate(sim=str_remove(sim,paste(dir_results,best_par_set,"",sep="/"))))
    }
  )%>%
    bind_rows()
  
  print(paste("sequence",tum_type))
  unlink(paste0(dir_results,"/best_sim_",tum_type),recursive = TRUE)
  
  save(vcf_multiregional_all,paste0(out_dir,"/vcf_multiregional_allsims.RData"))
  }

save(sim_data,file = "Data/best_fit/nD_points_bestfit.RData")

print("FINISHED")

