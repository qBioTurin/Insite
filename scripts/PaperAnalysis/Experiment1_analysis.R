library(Insite)
library(tibble)
library(jsonlite)
library(dplyr)
library(stringr)

dir<-"OutputSim/Experiment1"
folders<-list.dirs(dir,recursive = FALSE)

for(folder in folders){
  json_files<-list.files(paste(folder,"params",sep="/"))
  json_files<-json_files[grepl(x=json_files,".json")]
  
  sim_folder<-"params_passenger"
  sims<-list.files(paste(folder,"results",sim_folder,sep="/"))
  passenger_VAF<-tibble(VAF=sapply(sims, function(sim){
    Zprovv_file<-list.files(paste(folder,"results",sim_folder,sim,sep="/"))
    load(paste(folder,"results",sim_folder,sim,Zprovv_file,sep="/"))
    ncells<-sapply(Zprovv, Ncells)
    pop<-sapply(Zprovv, Population)
    gen<-sapply(pop, genotype)
    VAF<-ncells[lengths(gen)==2]/(2*sum(ncells))
    if(length(VAF)==0){VAF<-0}
    return(VAF)
  }))
  passenger_prob_surv<-sum(passenger_VAF$VAF>0)/nrow(passenger_VAF)
  
  growth_json_files<-json_files[grepl(x=json_files,"growth")]
  growth_VAF<-tibble()
  growth_sampled_VAF<-tibble()
  for(file in growth_json_files){
    sel_adv<-fromJSON(paste(folder,"params",file,sep="/"))$functionalEvents$params$proliferativeAdvantage[2]
    sim_folder<-gsub(x = file,pattern = ".json",replacement = "")
    sims<-list.files(paste(folder,"results",sim_folder,sep="/"))
    if(length(sims)>0){
      VAFs<-sapply(sims, function(sim){
        Zprovv_file<-list.files(paste(folder,"results",sim_folder,sim,sep="/"))
        load(paste(folder,"results",sim_folder,sim,Zprovv_file,sep="/"))
        ncells<-sapply(Zprovv, Ncells)
        pop<-sapply(Zprovv, Population)
        gen<-sapply(pop, genotype)
        VAF<-ncells[lengths(gen)==2]/(2*sum(ncells))
        if(length(VAF)==0){VAF<-0}
        return(VAF)
      })%>%unname()
      
      growth_VAF<-bind_rows(growth_VAF, tibble(sel_adv=sel_adv,VAF=VAFs))
    }
  }
  
  growth_prob_surv<-growth_VAF%>%
    group_by(sel_adv)%>%
    summarise(prob_surv_postK=sum(VAF>0)/length(VAF))
  
  
  space_json_files<-json_files[grepl(x=json_files,"space")]
  space_VAF<-tibble()
  for(file in space_json_files){
    add_space<-fromJSON(paste(folder,"params",file,sep="/"))$functionalEvents$params$additionalSpace[2]
    
    sim_folder<-gsub(x = file,pattern = ".json",replacement = "")
    sims<-list.files(paste(folder,"results",sim_folder,sep="/"))
    VAFs<-sapply(sims, function(sim){
      Zprovv_file<-list.files(paste(folder,"results",sim_folder,sim,sep="/"))
      load(paste(folder,"results",sim_folder,sim,Zprovv_file,sep="/"))
      ncells<-sapply(Zprovv, Ncells)
      pop<-sapply(Zprovv, Population)
      gen<-sapply(pop, genotype)
      VAF<-ncells[lengths(gen)==2]/(2*sum(ncells))
      if(length(VAF)==0){VAF<-0}
      return(VAF)
    })
    
    space_VAF<-bind_rows(space_VAF, tibble(add_space=add_space,VAF=VAFs))
  }
  
  space_prob_surv<-space_VAF%>%
    group_by(add_space)%>%
    summarise(prob_surv_postK=sum(VAF>0)/length(VAF))
  
  
  comp_json_files<-json_files[grepl(x=json_files,"comp")]
  
  comp_VAF<-tibble()
  for(file in comp_json_files){
    sus<-fromJSON(paste(folder,"params",file,sep="/"))$functionalEvents$params$susceptibility[2]
    off<-fromJSON(paste(folder,"params",file,sep="/"))$functionalEvents$params$offensiveScore[2]
    
    sim_folder<-gsub(x = file,pattern = ".json",replacement = "")
    sims<-list.files(paste(folder,"results",sim_folder,sep="/"))
    sims<-sims[grepl(x=sims,pattern = "sim")]
    if(length(sims)>0){
      VAFs<-sapply(sims, function(sim){
        Zprovv_file<-list.files(paste(folder,"results",sim_folder,sim,sep="/"))
        if(length(Zprovv_file)>0){
          load(paste(folder,"results",sim_folder,sim,Zprovv_file,sep="/"))
          ncells<-sapply(Zprovv, Ncells)
          pop<-sapply(Zprovv, Population)
          gen<-sapply(pop, genotype)
          VAF<-ncells[lengths(gen)==2]/(2*sum(ncells))
          if(length(VAF)==0){VAF<-0}
          return(VAF)
        }else{return(NA)}
      })
      
      comp_VAF<-bind_rows(comp_VAF, tibble(sus, off,VAF=VAFs))
    }
    
  }
  
  comp_prob_surv<-comp_VAF%>%
    group_by(sus,off)%>%
    summarise(prob_surv_postK=sum(VAF>0,na.rm = TRUE)/length(VAF))
  
  VAFs<-list(passenger_VAF,growth_VAF,comp_VAF,space_VAF)
  
  dir.create(paste0("Data/Experiment1",str_remove(folder,dir)))
  save(VAFs,file=paste0("Data/Experiment1",str_remove(folder,dir),"/VAF_FunEff.RData"))
  
}
