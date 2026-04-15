setClass("Population_with_size_nmut",
         slots=c(Population="Population",
                 Ncells="numeric",
                 Nmut="numeric"),
         prototype=list(Population=new("Population"),
                        Ncells=numeric(),
                        Nmut=numeric())
)

setGeneric("Population",function(Population_with_size_nmut) standardGeneric("Population"))
setMethod("Population",
          "Population_with_size_nmut",
          function(Population_with_size_nmut){
            return(Population_with_size_nmut@Population)
          })

setGeneric("Ncells",function(Population_with_size_nmut) standardGeneric("Ncells"))
setMethod("Ncells",
          "Population_with_size_nmut",
          function(Population_with_size_nmut){
            return(Population_with_size_nmut@Ncells)
          })

setGeneric("Nmut",function(Population_with_size_nmut) standardGeneric("Nmut"))
setMethod("Nmut",
          "Population_with_size_nmut",
          function(Population_with_size_nmut){
            return(Population_with_size_nmut@Nmut)
          })

setGeneric("get_bd",function(Population_with_size_nmut,p) standardGeneric("get_bd"))
setMethod("get_bd",
          signature(Population_with_size_nmut="Population_with_size_nmut",
                    p="list"),
          function(Population_with_size_nmut,p){
            Population<-Population_with_size_nmut@Population
            p_local<-p$p[sapply(p$phenotypes_local,identical,phenotype(Population))][[1]]
            Pmax<-length(p_local)
            size<-Population_with_size_nmut@Ncells
            if(size<.Machine$integer.max){
              x<-rmultinom(1,size=size,prob=p_local)
              y<-0:(Pmax-1)
              Population_with_size_nmut@Ncells<-sum(x*y)
            }
            else{
              mu<-size*p_local[1:(Pmax-1)]
              Sigma<--t(matrix(size,nrow = Pmax-1,ncol = Pmax-1)*p_local[1:(Pmax-1)])*p_local[1:(Pmax-1)]
              diag(Sigma)<-size*p_local[1:(Pmax-1)]*(1-p_local[1:(Pmax-1)])
              x<-mvrnorm(1,mu,Sigma)
              y<-0:(Pmax-1)
              Population_with_size_nmut@Ncells<-sum(x*y[1:(Pmax-1)])+y[Pmax]*(size-sum(x))
            }
            return(Population_with_size_nmut)
          })

setGeneric("get_mut",function(Population_with_size_nmut_old,Population_with_size_nmut,local_params) standardGeneric("get_mut"))
setMethod("get_mut",
          signature(Population_with_size_nmut_old="Population_with_size_nmut",
                    Population_with_size_nmut="Population_with_size_nmut",
                    local_params="Local_Params"),
          function(Population_with_size_nmut_old,Population_with_size_nmut,local_params){
            Population<-Population_with_size_nmut@Population
            Population_old<-Population_with_size_nmut_old@Population
            nmut_old<-Population_with_size_nmut@Nmut
            #if(nmut_old>0){browser()}
            mu_local<-local_params@mu_local[sapply(local_params@phenotypes_local,identical,phenotype(Population_old))][[1]]
            influenced_prob<-local_params@rel_freq_local[sapply(local_params@phenotypes_local,identical,phenotype(Population_old))][[1]]
            Delta<-local_params@Delta
            sum_size<-sum(Population_with_size_nmut_old@Ncells,Population_with_size_nmut@Ncells)
            if(mu_local*(Delta/2)*sum_size>.Machine$integer.max){
              nmut<-rnorm(1,mu_local*(Delta/2)*sum_size,sqrt(mu_local*(Delta/2)*sum_size))
            }
            else{nmut<-rpois(1,mu_local*(Delta/2)*sum_size)}
            if(length(nmut)==0||nmut==0){
              new_pop_list<-NULL}
            else{
              Population_with_size_nmut@Nmut<-nmut_old+nmut
              new_genotypes<-lapply(nmut_old+(1:nmut),function(n){c(Population@genotype,n)})
              new_functional_events<-sample(nmut,x=length(influenced_prob),prob=influenced_prob,replace = TRUE)
              new_functional_effect<-lapply(new_functional_events,function(f){c(Population@functional_effect,f)})
              new_phenotypes<-lapply(new_functional_events,function(f){unique(c(Population@phenotype,f))})
              new_pop_list<-mapply(function(new_genotypes,
                                            new_functional_effect,
                                            new_phenotypes){new("Population_with_size_nmut",
                                                                Population=new("Population",
                                                                               genotype=new_genotypes,
                                                                               functional_effect=new_functional_effect,
                                                                               phenotype=new_phenotypes),
                                                                Ncells=1,
                                                                Nmut=0)},
                                   new_genotypes,
                                   new_functional_effect,
                                   new_phenotypes)
            }
            if(Population_with_size_nmut@Ncells>0){
              return(append(Population_with_size_nmut,new_pop_list))}
            else{return(new_pop_list)}
          }
)

setGeneric("get_integral_error",function(Population_with_size_nmut_old,Population_with_size_nmut,local_params) standardGeneric("get_integral_error"))
setMethod("get_integral_error",
          signature(Population_with_size_nmut_old="Population_with_size_nmut",
                    Population_with_size_nmut="Population_with_size_nmut",
                    local_params="Local_Params"),
          function(Population_with_size_nmut_old,Population_with_size_nmut,local_params){
            Population<-Population_with_size_nmut@Population
            xT<-as.double(Population_with_size_nmut@Ncells)
            Population_old<-Population_with_size_nmut_old@Population
            x0<-as.double(Population_with_size_nmut_old@Ncells)
            a_local<-local_params@a_local[sapply(local_params@phenotypes_local,identical,phenotype(Population_old))][[1]]
            b_local<-local_params@b_local[sapply(local_params@phenotypes_local,identical,phenotype(Population_old))][[1]]
            mu_local<-local_params@mu_local[sapply(local_params@phenotypes_local,identical,phenotype(Population_old))][[1]]
            lambda_local<-a_local-b_local
            #Delta<-local_params@Delta
            
            if(lambda_local==0|mu_local==0){
              error<-function(Delta){return(0)}
            }
            else{error<-function(Delta){
              abs(mu_local*(Delta^3/12)*(lambda_local/sinh(Delta*lambda_local/2))*
                (lambda_local/(2*sinh(Delta*lambda_local/2))*(min(x0,xT)+max(x0,xT)*cosh(lambda_local*Delta))-
                   ((a_local+b_local)/2+lambda_local/sinh(Delta*lambda_local/2)*sqrt(x0*xT))*cosh(lambda_local/2*Delta)))}}
            return(error)
          }
)

setGeneric("get_Delta",function(Population_with_size_nmut_old,Population_with_size_nmut,local_params,epsilon) standardGeneric("get_Delta"))
setMethod("get_Delta",
          signature(Population_with_size_nmut_old="Population_with_size_nmut",
                    Population_with_size_nmut="Population_with_size_nmut",
                    local_params="Local_Params",
                    epsilon="numeric"),
          function(Population_with_size_nmut_old,Population_with_size_nmut,local_params,epsilon){
            error<-get_integral_error(Population_with_size_nmut_old,
                                      Population_with_size_nmut,
                                      local_params)
            inv_error_function<-function (y) {
              uniroot((function (Delta) y-error(Delta)),
                      lower = 10^(-200),
                      upper = local_params@Delta)[1]}
            
            root<-inv_error_function(epsilon)
            return(root$root)
          }
)

