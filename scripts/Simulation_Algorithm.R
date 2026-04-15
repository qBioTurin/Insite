epsilon_rel<-10^(-3) # errore massimo concesso a ogni step

simulazione<-function(Nexp,seed,path,starting_gen,starting_fun_eff,Ncells_start,parameters,tmax,Ncellsmax){
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
      save(list = c("Zprovv","time_provv"),
           file=paste(path_sim,"/Zprovv",count,".RData",sep=""))
      check_cond_end<-TRUE
      
      
      while (check_cond_end) {
        
        tot_Ncells<-sum(sapply(Zprovv,Ncells))
        epsilon<-tot_Ncells*epsilon_rel
        
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
          save(list = c("Zprovv","time_provv"),
               file=paste(path_sim,"/Zprovv",count,".RData",sep=""))
          count<-count+1
          #print(time_provv)
        }
      }
    }
    },
  error=function(cond) {
    write.table(paste(conditionMessage(cond),"\n",seed),
                paste(path,"error.log",sep="/"))
      }
    )
}
