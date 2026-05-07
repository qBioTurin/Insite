library(Insite)
library(dplyr)
library(tidyr)
library(stringr)

dir<-"OutputSim/nD_exp"
folder_sim<-paste(dir,"results",sep="/")
exp_folders<-list.dirs(folder_sim,full.names = TRUE,recursive = FALSE)
folder_out<-"Data/nD_exp"

if(!dir.exists(folder_out)){
  dir.create(folder_out,recursive = TRUE)
}

lapply(exp_folders,
       function(exp_folder){
         read.table(paste(exp_folder,"seed.txt",sep="/"),header = TRUE)%>%
           pivot_longer(all_of(1),names_to = "par_set",values_to = "seed")%>%
           return()
       })%>%
  bind_rows()%>%
  write.table(file = paste(dir,"seeds.txt",sep="/"),
              row.names = FALSE,col.names = TRUE)

lapply(exp_folders,
       function(exp_folder){
          time_tbl<-read.table(paste(exp_folder,"time_report.txt",sep="/"),header = FALSE)
          colnames(time_tbl)<-c("sim","time","unit")
          time_tbl<-time_tbl%>%bind_cols(par_set=str_remove(exp_folder,".*/"))
          return(time_tbl)
       })%>%
  bind_rows()%>%
  write.table(file = paste(folder_out,"time_report.txt",sep="/"),
              row.names = FALSE,col.names = TRUE)


nD_pts<-lapply(exp_folders,
               function(exp_folder){
                 json_file<-paste(dir,"params",paste0(str_remove(exp_folder,folder_sim),".json"),sep="/")
                 endingSize<-jsonlite::fromJSON(json_file)$endingSize
                 sims<-list.files(exp_folder)
                 sims<-sims[grepl("sim",sims)]
                 lapply(sims,function(sim){
                   Zprovv_file<-list.files(paste(exp_folder,sim,sep="/"))
                   load(paste(exp_folder,sim,Zprovv_file,sep="/"))
                   Noble_index<-get_nD_index(Zprovv,10^(-2))%>%
                     mutate(sim=sim,time_end=time_provv,det_size_reached=(sum(sapply(Zprovv,Ncells))>=endingSize))
                   return(Noble_index)
                 })%>%
                   bind_rows()%>%
                   mutate(par_set=str_remove(exp_folder,".*/"))%>%
                   return()
               })%>%
  bind_rows()

save(nD_pts,file = paste0(folder_out,"/nD_pts.RData"))
file.copy(from = paste0(dir,"/seeds.txt"), to = paste0(folder_out,"/seeds.txt"))
