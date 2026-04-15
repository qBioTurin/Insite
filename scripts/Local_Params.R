setClass("Local_Params",
         slots=c(phenotypes_local="list",
                 a_local="list",
                 b_local="list",
                 mu_local="list",
                 rel_freq_local="list",
                 Delta="numeric")
)

setGeneric("get_local_params",function(Parameters,list_of_pops_with_size_nmut,count,time_provv,epsilon) standardGeneric("get_local_params"))
setMethod("get_local_params",
          signature(Parameters="Parameters",
                    list_of_pops_with_size_nmut="list",
                    count="numeric",
                    time_provv="numeric",
                    epsilon="numeric"),
          function(Parameters,list_of_pops_with_size_nmut,count,time_provv,epsilon){
            list_of_pops<-lapply(list_of_pops_with_size_nmut,Population)
            phenotypes<-lapply(list_of_pops,phenotype)
            Ncells_by_phen<-vec_split(sapply(list_of_pops_with_size_nmut,Ncells),phenotypes)$val
            Ncells_sum<-sapply(Ncells_by_phen,sum)
            phenotypes_local<-phenotypes%>%unique()
            nphen<-length(phenotypes_local)

            populations_conv<-sapply(phenotypes_local,conversion,Parameters@I)
            Competition_local<-Parameters@Competition[populations_conv,populations_conv]
            if(nphen>1){
              weighted_sums_local<-rowSums(t(t(Competition_local)*Ncells_sum))
            }else{weighted_sums_local<-Competition_local*Ncells_sum}
            rel_freq_norm<-lapply(phenotypes_local,function(phenotype){
              rel_freq<-diag(Parameters@influence)
              #rel_freq[phenotype]<-1
              #if(length(phenotype)>1){
              #  rel_freq<-matrixStats::colProds(Parameters@influence[phenotype,])*rel_freq
              #}
              #else{Parameters@influence[phenotype,]*rel_freq}
              return(rel_freq/sum(rel_freq))
            })
            lambda_local<-Parameters@lambda[populations_conv]*(1-weighted_sums_local/Parameters@K[populations_conv])
            b_local<-rep(1,length(populations_conv))
            a_local<-(lambda_local+b_local)/(Parameters@mu[populations_conv]+1)
            b_local[a_local<0]<--lambda_local[a_local<0]
            a_local[a_local<0]<-0
            mu_local<-Parameters@mu[populations_conv]*a_local
            mu_local<-(mu_local+abs(mu_local))/2

            
            
            Delta<-epsilon*sum(a_local+b_local+mu_local)*min(
              abs(
                t(1/rowSums(t(Competition_local)*(a_local-b_local)*Ncells_sum)*
                    (Parameters@mu[populations_conv]+1)/max(Parameters@mu[populations_conv],1)*
                    Parameters@K[populations_conv]/Parameters@lambda[populations_conv]
                  )
                )
              )
            
            
            
            if(time_provv+Delta*Parameters@av_lifespan>Parameters@print_time[count+1]){ # SE SUPERO PRINT TIME mi fermo lì
              Delta<-(Parameters@print_time[count+1]-time_provv)/Parameters@av_lifespan
            }
            Delta<-min(5,Delta)
            Delta<-min(abs(log(2)*(.Machine$double.max.exp-1)/lambda_local),Delta)

            local_Params<-new("Local_Params",
                              phenotypes_local=phenotypes_local,
                              a_local=as.list(a_local),
                              b_local=as.list(b_local),
                              mu_local=as.list(mu_local),
                              rel_freq_local=rel_freq_norm,
                              Delta=Delta)
          }
)

setGeneric("get_p",function(local_params) standardGeneric("get_p"))
setMethod("get_p",
          "Local_Params",
          definition = function(local_params){
            a<-local_params@a_local
            b<-local_params@b_local
            Delta<-local_params@Delta
            p<-mapply(function(a,b,Delta){
              if(a==b){
                p_new<-(a*Delta)/(1+a*Delta)
                p<-p_new
                m<-2
                start_time <- Sys.time()
                
                while(p_new>0){
                  p_new<-(a*Delta)^(m-2)/(1+a*Delta)^m
                  p<-c(p,p_new)
                  m<-m+1
                  
                  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
                  if (elapsed > 60) {
                    message("Warning: approximation made on the probability estimation")
                    break
                  }
                }
                p<-p[-length(p)]
              }else if(a==0){
                p<-c(1-exp(-b*Delta),exp(-b*Delta))
                }else if(b==0){
                p_new<-exp(-a*Delta)/(1-exp(-a*Delta))
                p<-p_new
                start_time <- Sys.time()
                
                while(p_new>0){
                  p_new<-p_new*(1-exp(-a*Delta))
                  p<-c(p,p_new)
                  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
                  if (elapsed > 60) {
                    message("Warning: approximation made on the probability estimation")
                    break
                  }
                }
              }else{
                lambda<-a-b
                
                alpha=(exp(lambda*Delta)-1)/(a*exp(lambda*Delta)/b-1)
                beta=(exp(lambda*Delta)-1)/(exp(lambda*Delta)-b/a)
                if(beta==1){
                  if(a>b){beta<-1-.Machine$double.eps}
                  else{beta<-1+.Machine$double.eps}
                }
                if(alpha==0){
                  alpha<-.Machine$double.xmin
                }
                p_new<-alpha
                p<-p_new
                m<-2
                start_time <- Sys.time()
                
                while(p_new>0){
                  p_new<-(1-alpha)*(1-beta)*beta^(m-2)
                  p<-c(p,p_new)
                  m<-m+1
                  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
                  if (elapsed > 60) {
                    message("Warning: approximation made in get_p")
                    break
                  }
                }
                p<-p[-length(p)]
              }},
              a,
              b,
              Delta,
              SIMPLIFY = FALSE)
            return(list(p=p,phenotypes_local=local_params@phenotypes_local))
          }
)
