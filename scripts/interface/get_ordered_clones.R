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

load(path_params)

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

daughters<-lapply(pop,
                  function(p){
                    return(which(
                      sapply(pop,
                             function(p1){
                               is_parent(p1,p)
                             }))
                    )
                  })

clones<-lapply(Pop_ID,function(i){c(i,daughters[[i]])})

ncells_clone=sapply(clones,function(c){
  sum(ncells[c])
})

root<-setdiff(Pop_ID,unique(unlist(daughters[lengths(daughters)>0])))
pop_with_sons<-which(lengths(daughters)>0)

if(length(root)>1){
  Pop_ID<-c(0,Pop_ID)
  ncells_clone<-c(sum(ncells_clone[root]),ncells_clone)
  ncells<-c(0,ncells)
  pop_with_sons<-c(0,pop_with_sons)
  daughters<-append(list(root),daughters)
  root<-0
}

Clones_df<-tibble(clone=Pop_ID,Ncells_clone=ncells_clone,Ncells=ncells)%>%
  filter(clone%in%root)%>%
  mutate(y_lower=-Ncells_clone/2,
         y_upper=Ncells_clone/2)

for(p in pop_with_sons){
  daughters_p<-daughters[Pop_ID==p][[1]]
  n_siblings_d<-length(daughters_p)
  
  for(d in daughters_p){
    sibling_number<-which(daughters_p==d)
    older_siblings<-daughters_p[0:(sibling_number-1)]
    
    prop_y_start<-sibling_number/(n_siblings_d+1)
    
    center_parent<-Clones_df%>%
      filter(clone==p)%>%
      mutate(center=(y_upper+y_lower)/2)%>%
      pull(center)
    
    Ncells_d<-ncells[Pop_ID==d]
    
    Ncells_parent<-ncells[Pop_ID==p]
    
    Ncells_parent_clone<-ncells_clone[Pop_ID==p]
    
    Ncells_older_siblings<-ncells_clone[Pop_ID%in%older_siblings]%>%
      sum()
    
    d_df<-tibble(
      clone=d,
      Ncells_clone=ncells_clone[Pop_ID==d],
      Ncells=Ncells_d,
      Ncells_p=Ncells_parent,
      Ncells_p_clone=Ncells_parent_clone,
      Ncells_s=Ncells_older_siblings,
      center=center_parent
    )%>%
      mutate(y_lower=sum(center,prop_y_start*Ncells_p,-Ncells_p_clone/2,Ncells_s,na.rm=TRUE),
             y_upper=y_lower+Ncells_clone)%>%
      dplyr::select(clone,Ncells_clone,y_lower,y_upper,Ncells)
    
    Clones_df<-rbind(Clones_df,d_df)
  }
}
Clones_df<-Clones_df%>%
  filter(clone!=0)

if(root==0){
  ncells_clone<-ncells_clone[-1]
  ncells<-ncells[-1]
  pop_with_sons<-pop_with_sons[-1]
  root<-daughters[[1]]
  daughters<-daughters[-1]
  Pop_ID<-Pop_ID[-1]
  }

trasl<-min(Clones_df$y_lower)
Clones_df$y_lower<-Clones_df$y_lower-trasl
Clones_df$y_upper<-Clones_df$y_upper-trasl
Clones_df$clone<-factor(Clones_df$clone,levels=Pop_ID)

Clones_df<-distinct(Clones_df)

save(Clones_df,file = paste(path_out,"Clones_ordered_muller.RData",sep="/"))
