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
    c("--path_out"),
    type="character",
    default = "raw",
    help = "folder in which the output plot files are stored"
  ),
  make_option(
    c("--seq_day"),
    type="numeric",
    default = Inf,
    help = "Day in which the sequencing is to be performed. Default is the end of the simulation"
  ),
  make_option(
    c("--nregions"),
    type="numeric",
    default = 1,
    help = "Number of regions of the sequencing, if multiregional or bluk (in such case 1)"
  ),
  make_option(
    c("--ncells"),
    type="numeric",
    default = 1000,
    help = "Number of cells to be sequenced in each region, use Inf for sequencing the whole mass"
  ),
  make_option(
    c("--neighborhood"),
    type="character",
    default = NULL,
    help = "If the neighborhood for the simulation has already been computed, indicate here the path to the file (default name Clones_ordered_day.RData).
    The most expensive part of the sequencing is the computation of these neighborhood, which once done is not required for repeating the sequencing"
  ),
  make_option(
    c("--seed"),
    type="numeric",
    default = NULL,
    help = "To retrieve the same sequencing, fix the seed. Default is not fixed."
  ),
  make_option(
    c("--repeat"),
    type="numeric",
    default = 1,
    help = "How many times the sequencing should be repeated (produce a file per repetition)"
  ),
  make_option(
    c("--dens"),
    type="character",
    default = "Data/dens.RData",
    help = "path to the .RData file with the DP probability density to be used in the sequencing for coverage"
  )
)

opt_parser<-OptionParser(option_list = option_list)
opt<-parse_args(opt_parser)

path_sim <- opt$sim_dir
path_params <- opt$params
path_out <- opt$path_out
seq_day<-opt$seq_day
n_regions<-opt$nregions
n_seq_cells<-opt$ncells
Clones_ordered_path<-opt$neighborhood
seed<-opt$seed
Nrep<-opt$`repeat`
dens_path<-opt$dens


sim_files<-list.files(path_sim)
sim_files<-sim_files[grepl("Zprovv",sim_files)]

index_files<-as.numeric(str_remove(str_remove(string = sim_files,pattern = "Zprovv"),".RData"))

load(path_params)

if(seq_day==Inf){
  n_sim<-which.max(index_files)
  Zprovv_file<-sim_files[n_sim]
}else{
  Zprovv_file<-sim_files[index_files==(which.min(abs(parameters@print_time-seq_day)))]
}
load(paste(path_sim,Zprovv_file,sep = "/"))

if(is.null(Clones_ordered_path)){
  
  Clones_df<-get_ordered_clones_sequencing(Zprovv)
  
  save(Clones_df,file = paste(path_out,"/Clones_ordered_",seq_day,".RData",sep=""))
  
}else{
  load(Clones_ordered_path)
}

if(is.null(seed)){
    runif(1)
    seed_selected<-as.integer(Sys.time())
  }else{
    seed_selected<-seed
  }

set.seed(seed_selected)

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

lapply(1:Nrep,function(nrep){
if(n_regions*n_seq_cells>=sum(ncells)){
  seq_min_y_list<-0
  n_seq_cells<-sum(ncells)
  seq_min_y_list<-min(Clones_df$y_lower)
  print("The number of sequenced cells exeded the total number of cells of the mass, hence the whole mass has been sequenced")
  
}else{
  L <- (max(Clones_df$y_upper)-min(Clones_df$y_lower))-(n_regions - 1)*n_seq_cells
  u <- sort(runif(n_regions,0,L))
  seq_min_y_list <- min(Clones_df$y_lower)+u+(0:(n_regions-1))*n_seq_cells
}

load(dens_path)
vcf_multiregional<-bind_rows(lapply(seq_min_y_list,function(seq_min_y){
  seq_max_y<-seq_min_y+n_seq_cells
  
  Clones_df_seq<-Clones_df%>%
    filter(y_lower<seq_max_y,y_upper>seq_min_y)
  
  mut_fe_to_add<-lapply(Clones_df_seq$clone,function(c){
    g<-gen[[c]]
    lg<-length(g)
    if(lg==1){
      return(tibble())
    }
    f<-fun_eff[[c]]
    
    mta<-vector()
    feta<-vector()
    for(i in 1:(lg-1)){
      newg<-head(g,-i)
      newf<-names(parameters@functional_effects)[f[lg-i]]
      if(any(sapply(gen[Clones_df_seq$clone],function(g1){identical(g1,newg)}))){
        next
      }else{
        mta<-c(mta,paste(newg,collapse="_"))
        feta<-c(feta,paste(newf,collapse="_"))
      }
    }
    return(tibble(mut=as.character(mta),
                  fun_eff=as.character(feta),
                  y_lower=seq_min_y,y_upper=seq_max_y))
  })%>%bind_rows()%>%distinct()
  
  Clones_df_seq<-merge(Clones_df_seq,tibble(clone=Pop_ID,
                                            mut=sapply(unique_mut_id,function(muts){muts[length(muts)]}),
                                            fun_eff=sapply(fun_eff_label,function(funcs){funcs[length(funcs)]})))%>%
    bind_rows(mut_fe_to_add)%>%
    rowwise()%>%
    mutate(fun_eff_category=parameters@functional_effects[fun_eff],
           y2=ifelse(y_upper<seq_max_y,y_upper,seq_max_y),
           y1=ifelse(y_lower>seq_min_y,y_lower,seq_min_y),
           Ncells_seq=round(y2-y1),
           prob=Ncells_seq/(2*(seq_max_y-seq_min_y)))%>%
    ungroup()%>%
    dplyr::select(mut,fun_eff,fun_eff_category,Ncells_seq,prob)
  
  
  
  sample_DP<-round(sample(x = dens_coverage$x, nrow(Clones_df_seq), prob = dens_coverage$y, replace=TRUE) + rnorm(1, 0, dens_coverage$bw))
  sample_DP[sample_DP<0]<-0
  sample_AD<-mapply(rbinom,prob=Clones_df_seq$prob,size=sample_DP,MoreArgs = list(n=1))
  
  
  vcf_sample<-Clones_df_seq%>%
    bind_cols(sample_DP=sample_DP,sample_AD=sample_AD)%>%
    filter(sample_AD>0)%>%
    mutate(VAF=sample_AD/sample_DP)%>%
    dplyr::select(-c(prob,Ncells_seq))
  
  return(vcf_sample)
  }),.id	="region")

  write.table(vcf_multiregional,row.names = FALSE,
              file = paste(path_out,"/vcf",nrep,".txt",sep=""))
  
})
