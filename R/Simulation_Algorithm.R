#' Run a simulation experiment.
#'
#' Runs one stochastic simulation experiment of the population dynamics model,
#' starting from the initial populations and saving intermediate states to disk.
#'
#' @param Nexp Numeric index of the simulation experiment.
#'  If more than one simulation is wanted, it is possible to cycle or parallelize 
#'  over this number letting it vary from 1 to a wanted Nexp_max.
#' @param seed Numeric seed used to initialize the random number generator
#'  the actual seed set by the simulator is seed+Nexp to let multiple simulation
#'  runs being replicable but avoid producing the same output repeatedly
#' @param path Character string giving the output directory.
#' @param starting_gen List of starting genotypes.
#' @param starting_fun_eff List of starting functional effects associated with
#'   the starting genotypes.
#' @param Ncells_start Numeric vector giving the initial number of cells for
#'   each starting population.
#' @param parameters A [`Parameters-class`] object containing global simulation
#'   parameters.
#' @param tmax Maximum simulation time.
#' @param Ncellsmax Optional maximum population size used as stopping criterion.
#' @param epsilon_rel Maximum error ammitted at each simulative step (relatively to the simulation parameters)
#' @param only_last logical. If TRUE only the last time stamp will we saved, 
#'    representing the ending status of the simulation.
#' @return
#' This creates a simulation output directory, saves intermediate `.RData` files
#' containing  `Zprovv`, which are snapshot of tumor evolution at different times and
#' their corresponding time (`time_provv`). It writes an `error.log` file if an error occurs.
#'
#' @details
#' The simulation is restarted until a non-empty trajectory is produced. During
#' each run, the function repeatedly updates local parameters, samples
#' birth-death outcomes, applies mutation events, advances time, and stores
#' snapshots when checkpoint times are reached.
#'
#' @name simulazione
#' @rdname simulazione
#' @export
simulation<-function(
    Nexp,
    seed,
    path,
    starting_gen,
    starting_fun_eff,
    Ncells_start,
    parameters,
    tmax,
    Ncellsmax,
    epsilon_rel,
    only_last=FALSE){
  tryCatch(expr = {
    set.seed(seed+Nexp)
    if(length(starting_fun_eff)!=length(starting_gen)){stop("Check genotype and functional effect association")}
    path_sim<-paste(path,"/sim",Nexp,sep="")
    unlink(path_sim,recursive =TRUE)
    dir.create(path_sim,recursive =TRUE)
    
    gc()
    
    starting_pop<-mapply(genotype=starting_gen,
                         functional_effect=starting_fun_eff,
                         phenotype=lapply(starting_fun_eff,unique),
                         new,MoreArgs = list(Class="Population")
    )
    last_mut<-sapply(starting_gen,function(gen){gen[length(gen)]})
    Nmut<-sapply(starting_gen,
                 function(gen){
                   if(length(gen)<max(lengths(starting_gen))){
                     max(last_mut[lengths(starting_gen)==length(gen)+1])
                   }else{0}
                 })
    
    Zprovv<-mapply(Population=starting_pop,
                   Ncells=Ncells_start,
                   Nmut=Nmut,
                   new,
                   MoreArgs = list(Class="Population_with_size_nmut")
    )
    
    time_provv<-0
    
    
    check_cond_end<-TRUE
    while(check_cond_end | length(Zprovv)==0 ||sum(sapply(Zprovv,Ncells))==0){
      
      count<-0
      Zprovv<-mapply(Population=starting_pop,
                     Ncells=Ncells_start,
                     Nmut=Nmut,
                     new,
                     MoreArgs = list(Class="Population_with_size_nmut")
      )
      
      time_provv<-0
      unlink(path_sim,recursive =TRUE)
      dir.create(path_sim,recursive =TRUE)
      if(!only_last){save(list = c("Zprovv","time_provv"),
           file=paste(path_sim,"/Zprovv",count,".RData",sep=""))}
      check_cond_end<-TRUE
      
      
      while (check_cond_end) {
        local_params<-get_local_params(Parameters = parameters,
                                       list_of_pops_with_size_nmut =Zprovv,
                                       count=count,
                                       time_provv=time_provv,
                                       epsilon = epsilon_rel)

        p<-get_p(local_params)
        
        W<-lapply(Zprovv,get_bd,p)
        
        
        Zprovv<-mapply(get_mut,
                       Zprovv,
                       W,
                       MoreArgs =list(local_params),SIMPLIFY = TRUE)%>%unlist()
        
        time_provv<-time_provv+local_params@Delta*parameters@av_lifespan
        print(time_provv)
        if(is.null(Zprovv)|length(Zprovv)==0||sum(sapply(Zprovv,Ncells))==0){
          check_cond_end<-FALSE
        }else if(!is.null(Ncellsmax)){
          check_cond_end<-(sum(sapply(Zprovv,Ncells))<Ncellsmax)&(time_provv<max(parameters@print_time))
        }else{
          check_cond_end<-(time_provv<tmax)
        }
        
        if(time_provv>=parameters@print_time[count+1]){
          if(!only_last){
            save(list = c("Zprovv","time_provv"),
                 file=paste(path_sim,"/Zprovv",count,".RData",sep=""))
          }
          count<-count+1
        }
      }
      
      if(only_last){
        save(list = c("Zprovv","time_provv"),
             file=paste(path_sim,"/ZprovvFinal.RData",sep=""))
      }
    }
    },
  error=function(cond) {
    write.table(paste(conditionMessage(cond),"\n",seed),
                paste(path,"error.log",sep="/"))
      }
    )
}