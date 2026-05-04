library(Insite)
library(tibble)
library(dplyr)
library(stringr)

dir<-"OutputSim/Experiment2"
json_files<-list.files(paste(dir,"params",sep="/"))
json_files<-json_files[grepl(x=json_files,".json")]
params_names<-list.files(paste(dir,"params",sep="/"))
params_names<-str_remove(params_names[grepl(".json",params_names)],".json")
prop_cells_seq<-0.01

# PASSENGER
par_set<-"params_passenger"
sims<-list.files(paste(dir,"results",par_set,sep="/"))
sims<-sims[grepl("sim",sims)]
set.seed(1)
vcf_list<-lapply(sims,
       function(sim){
         Zprovv_file<-list.files(paste(dir,"results",par_set,sim,sep="/"))
         load(paste(dir,"results",par_set,sim,Zprovv_file,sep="/"))
         
         load(paste(dir,"results",par_set,"Parameters.RData",sep="/"))
         starting_mut<-sapply(starting_gen, paste0,collapse = "_")
         
         vcf_pure<-Insite:::theoretical_sequencing(Zprovv=Zprovv,
                                         parameters = parameters,
                                         starting_mut = starting_mut)
         
         n_seq_cells<-round(prop_cells_seq*sum(sapply(Zprovv,Ncells)))
         vcf_sampled<-sequencing(Zprovv = Zprovv,
                                 parameters = parameters,
                                 n_regions = 1,
                                 Nrep = 1,
                                 n_seq_cells = n_seq_cells)[[1]][[1]]%>%
           filter(mut!=starting_mut)
         return(list(vcf_pure,vcf_sampled))
       })
vcf_list_1<-lapply(vcf_list,function(two_elem){return(two_elem[[1]])})
vcf_list_2<-lapply(vcf_list,function(two_elem){return(two_elem[[2]])})

passenger_synth_vcf<-bind_rows(vcf_list_1)
passenger_synth_vcf_sampled<-bind_rows(vcf_list_2)

save(passenger_synth_vcf_sampled,passenger_synth_vcf,
     file = paste("Data/Experiment2/passenger_vcf_synth.RData",sep=""))

# GROWTH
par_sets<-params_names[grepl("growth",params_names)]
growth_synth_vcf_sampled<-tibble()
growth_synth_vcf<-tibble()
set.seed(1)
for(par_set in par_sets){
    file<-paste0(dir,"/params/",par_set,".json")
    sel_adv<-jsonlite::fromJSON(txt =  file)$functionalEvents$params$proliferativeAdvantage[2]
    sims<-list.files(paste(dir,"results",par_set,sep="/"))
    sims<-sims[grepl("sim",sims)]
    vcf_list<-lapply(sims,
                     function(sim){
                       Zprovv_file<-list.files(paste(dir,"results",par_set,sim,sep="/"))
                       load(paste(dir,"results",par_set,sim,Zprovv_file,sep="/"))
                       
                       load(paste(dir,"results",par_set,"Parameters.RData",sep="/"))
                       starting_mut<-sapply(starting_gen, paste0,collapse = "_")
                       
                       vcf_pure<-Insite:::theoretical_sequencing(Zprovv=Zprovv,
                                                                 parameters = parameters,
                                                                 starting_mut = starting_mut)
                       
                       n_seq_cells<-round(prop_cells_seq*sum(sapply(Zprovv,Ncells)))
                       vcf_sampled<-sequencing(Zprovv = Zprovv,
                                               parameters = parameters,
                                               n_regions = 1,
                                               Nrep = 1,
                                               n_seq_cells = n_seq_cells)[[1]][[1]]%>%
                         filter(mut!=starting_mut)
                       return(list(vcf_pure,vcf_sampled))
                     })
    vcf_list_1<-lapply(vcf_list,function(two_elem){return(two_elem[[1]])})
    vcf_list_2<-lapply(vcf_list,function(two_elem){return(two_elem[[2]])})
    
    synth_vcf<-bind_rows(vcf_list_1)%>%mutate(sel_adv=as.numeric(sel_adv))
    synth_vcf_sampled<-bind_rows(vcf_list_2)%>%mutate(sel_adv=as.numeric(sel_adv))
    
    growth_synth_vcf_sampled<-bind_rows(growth_synth_vcf_sampled,synth_vcf_sampled)
    growth_synth_vcf<-bind_rows(growth_synth_vcf,synth_vcf)
  }

save(growth_synth_vcf,growth_synth_vcf_sampled,
     file = paste("Data/Experiment2/growth_vcf_synth.RData",sep=""))


# COMPETITION

par_sets<-params_names[grepl("comp",params_names)]
comp_synth_vcf_sampled<-tibble()
comp_synth_vcf<-tibble()
set.seed(1)
for(par_set in par_sets){
  file<-paste0(dir,"/params/",par_set,".json")
  sus<-jsonlite::fromJSON(file)$functionalEvents$params$susceptibility[2]
  off<-jsonlite::fromJSON(file)$functionalEvents$params$offensiveScore[2]
  
  sims<-list.files(paste(dir,"results",par_set,sep="/"))
  sims<-sims[grepl("sim",sims)]
  vcf_list<-lapply(sims,
                   function(sim){
                     Zprovv_file<-list.files(paste(dir,"results",par_set,sim,sep="/"))
                     load(paste(dir,"results",par_set,sim,Zprovv_file,sep="/"))
                     
                     load(paste(dir,"results",par_set,"Parameters.RData",sep="/"))
                     starting_mut<-sapply(starting_gen, paste0,collapse = "_")
                     
                     vcf_pure<-Insite:::theoretical_sequencing(Zprovv=Zprovv,
                                                               parameters = parameters,
                                                               starting_mut = starting_mut)
                     
                     n_seq_cells<-round(prop_cells_seq*sum(sapply(Zprovv,Ncells)))
                     vcf_sampled<-sequencing(Zprovv = Zprovv,
                                             parameters = parameters,
                                             n_regions = 1,
                                             Nrep = 1,
                                             n_seq_cells = n_seq_cells)[[1]][[1]]%>%
                       filter(mut!=starting_mut)
                     return(list(vcf_pure,vcf_sampled))
                   })
  vcf_list_1<-lapply(vcf_list,function(two_elem){return(two_elem[[1]])})
  vcf_list_2<-lapply(vcf_list,function(two_elem){return(two_elem[[2]])})
  
  synth_vcf<-bind_rows(vcf_list_1)%>%mutate(sel_adv=as.numeric(sel_adv))
  synth_vcf_sampled<-bind_rows(vcf_list_2)%>%mutate(sel_adv=as.numeric(sel_adv))
  
  comp_synth_vcf_sampled<-bind_rows(comp_synth_vcf_sampled,synth_vcf_sampled)
  comp_synth_vcf<-bind_rows(comp_synth_vcf,synth_vcf)
}


save(comp_synth_vcf,comp_synth_vcf_sampled,
     file = paste("Data/Experiment2/comp_vcf_synth.RData",sep=""))

# SPACE

par_sets<-params_names[grepl("space",params_names)]
space_synth_vcf_sampled<-tibble()
space_synth_vcf<-tibble()
set.seed(1)
for(par_set in par_sets){
  file<-paste0(dir,"/params/",par_set,".json")
  add_space<-jsonlite::fromJSON(txt =  file)$functionalEvents$params$additionalSpace[2]
  sims<-list.files(paste(dir,"results",par_set,sep="/"))
  sims<-sims[grepl("sim",sims)]
  vcf_list<-lapply(sims,
                   function(sim){
                     Zprovv_file<-list.files(paste(dir,"results",par_set,sim,sep="/"))
                     load(paste(dir,"results",par_set,sim,Zprovv_file,sep="/"))
                     
                     load(paste(dir,"results",par_set,"Parameters.RData",sep="/"))
                     starting_mut<-sapply(starting_gen, paste0,collapse = "_")
                     
                     vcf_pure<-Insite:::theoretical_sequencing(Zprovv=Zprovv,
                                                               parameters = parameters,
                                                               starting_mut = starting_mut)
                     
                     n_seq_cells<-round(prop_cells_seq*sum(sapply(Zprovv,Ncells)))
                     vcf_sampled<-sequencing(Zprovv = Zprovv,
                                             parameters = parameters,
                                             n_regions = 1,
                                             Nrep = 1,
                                             n_seq_cells = n_seq_cells)[[1]][[1]]%>%
                       filter(mut!=starting_mut)
                     return(list(vcf_pure,vcf_sampled))
                   })
  vcf_list_1<-lapply(vcf_list,function(two_elem){return(two_elem[[1]])})
  vcf_list_2<-lapply(vcf_list,function(two_elem){return(two_elem[[2]])})
  
  synth_vcf<-bind_rows(vcf_list_1)%>%mutate(add_space=as.numeric(add_space))
  synth_vcf_sampled<-bind_rows(vcf_list_2)%>%mutate(add_space=as.numeric(add_space))
  
  space_synth_vcf_sampled<-bind_rows(space_synth_vcf_sampled,synth_vcf_sampled)
  space_synth_vcf<-bind_rows(space_synth_vcf,synth_vcf)
}


save(space_synth_vcf,space_synth_vcf_sampled,
     file = paste("Data/Experiment2/space_vcf_synth.RData",sep=""))


# PASSENGER num bases varied

par_sets<-params_names[grepl("passenger",params_names)]
passenger_synth_vcf_sampled<-tibble()
passenger_synth_vcf<-tibble()
set.seed(1)
for(par_set in par_sets){
  file<-paste0(dir,"/params/",par_set,".json")
  nmut<-jsonlite::fromJSON(txt =  file)$mutableBases
  sims<-list.files(paste(dir,"results",par_set,sep="/"))
  sims<-sims[grepl("sim",sims)]
  vcf_list<-lapply(sims,
                   function(sim){
                     Zprovv_file<-list.files(paste(dir,"results",par_set,sim,sep="/"))
                     load(paste(dir,"results",par_set,sim,Zprovv_file,sep="/"))
                     
                     load(paste(dir,"results",par_set,"Parameters.RData",sep="/"))
                     starting_mut<-sapply(starting_gen, paste0,collapse = "_")
                     
                     vcf_pure<-Insite:::theoretical_sequencing(Zprovv=Zprovv,
                                                               parameters = parameters,
                                                               starting_mut = starting_mut)
                     
                     n_seq_cells<-round(prop_cells_seq*sum(sapply(Zprovv,Ncells)))
                     vcf_sampled<-sequencing(Zprovv = Zprovv,
                                             parameters = parameters,
                                             n_regions = 1,
                                             Nrep = 1,
                                             n_seq_cells = n_seq_cells)[[1]][[1]]%>%
                       filter(mut!=starting_mut)
                     return(list(vcf_pure,vcf_sampled))
                   })
  vcf_list_1<-lapply(vcf_list,function(two_elem){return(two_elem[[1]])})
  vcf_list_2<-lapply(vcf_list,function(two_elem){return(two_elem[[2]])})
  
  synth_vcf<-bind_rows(vcf_list_1)%>%mutate(nmut=as.numeric(nmut))
  synth_vcf_sampled<-bind_rows(vcf_list_2)%>%mutate(nmut=as.numeric(nmut))
  
  passenger_synth_vcf_sampled<-bind_rows(passenger_synth_vcf_sampled,synth_vcf_sampled)
  passenger_synth_vcf<-bind_rows(passenger_synth_vcf,synth_vcf)
}


save(passenger_synth_vcf,passenger_synth_vcf_sampled,
     file = paste("Data/Experiment2/passenger_vcf_synth_nmut.RData",sep=""))

