#' Population_with_size_nmut class
#'
#' S4 class representing a population together with its current size and the
#' number of mutations generated from it.
#'
#' @slot Population A [`Population-class`] object.
#' @slot Ncells Numeric value giving the number of cells in the population.
#' @slot Nmut Numeric value giving the number of mutations generated from the
#'   population.
#'
#' @name Population_with_size_nmut-class
#' @rdname Population_with_size_nmut-class
#' @exportClass Population_with_size_nmut
setClass("Population_with_size_nmut",
         slots=c(Population="Population",
                 Ncells="numeric",
                 Nmut="numeric"),
         prototype=list(Population=new("Population"),
                        Ncells=numeric(),
                        Nmut=numeric())
)


#' Get the population component
#'
#' Returns the [`Population-class`] object stored inside a
#' [`Population_with_size_nmut-class`] object.
#'
#' @param Population_with_size_nmut A [`Population_with_size_nmut-class`]
#'   object.
#' @return A [`Population-class`] object.
#'
#' @rdname Population
#' @export
setGeneric("Population", function(Population_with_size_nmut) standardGeneric("Population"))

setMethod("Population",
          "Population_with_size_nmut",
          function(Population_with_size_nmut){
            return(Population_with_size_nmut@Population)
          })


#' Get the number of cells
#'
#' Returns the number of cells stored in a
#' [`Population_with_size_nmut-class`] object.
#'
#' @param Population_with_size_nmut A [`Population_with_size_nmut-class`]
#'   object.
#' @return A numeric value.
#'
#' @rdname Ncells
#' @export
setGeneric("Ncells", function(Population_with_size_nmut) standardGeneric("Ncells"))

setMethod("Ncells",
          "Population_with_size_nmut",
          function(Population_with_size_nmut){
            return(Population_with_size_nmut@Ncells)
          })


#' Get the number of mutations
#'
#' Returns the mutation count stored in a
#' [`Population_with_size_nmut-class`] object.
#'
#' @param Population_with_size_nmut A [`Population_with_size_nmut-class`]
#'   object.
#' @return A numeric value.
#'
#' @rdname Nmut

setGeneric("Nmut", function(Population_with_size_nmut) standardGeneric("Nmut"))

setMethod("Nmut",
          "Population_with_size_nmut",
          function(Population_with_size_nmut){
            return(Population_with_size_nmut@Nmut)
          })


#' Update population size after birth-death sampling
#'
#' Updates the number of cells in a population by sampling from the local
#' offspring distribution.
#'
#' @param Population_with_size_nmut A [`Population_with_size_nmut-class`]
#'   object.
#' @param p A list containing local offspring probabilities and phenotype
#'   information.
#'
#' @return A [`Population_with_size_nmut-class`] object with updated `Ncells`.
#'
#' @noRd
#' @importFrom MASS mvrnorm
setGeneric("get_bd", function(Population_with_size_nmut, p) standardGeneric("get_bd"))

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


#' Generate mutations from a population
#'
#' Simulates mutation events from a population over a local time step and
#' returns the updated population together with any newly created populations.
#'
#' @param Population_with_size_nmut_old A
#'   [`Population_with_size_nmut-class`] object representing the previous state.
#' @param Population_with_size_nmut A [`Population_with_size_nmut-class`]
#'   object representing the updated state.
#' @param local_params A [`Local_Params-class`] object.
#'
#' @return Either a single [`Population_with_size_nmut-class`] object or a list
#'   containing the updated population and newly created mutant populations.
#'
#' @noRd

setGeneric("get_mut", function(Population_with_size_nmut_old, Population_with_size_nmut, local_params) standardGeneric("get_mut"))

setMethod("get_mut",
          signature(Population_with_size_nmut_old="Population_with_size_nmut",
                    Population_with_size_nmut="Population_with_size_nmut",
                    local_params="Local_Params"),
          function(Population_with_size_nmut_old,Population_with_size_nmut,local_params){
            Population<-Population_with_size_nmut@Population
            Population_old<-Population_with_size_nmut_old@Population
            nmut_old<-Population_with_size_nmut@Nmut
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


#' Build the integral error function
#'
#' Constructs the error function used to control the local time discretization.
#'
#' @param Population_with_size_nmut_old A
#'   [`Population_with_size_nmut-class`] object representing the previous state.
#' @param Population_with_size_nmut A [`Population_with_size_nmut-class`]
#'   object representing the updated state.
#' @param local_params A [`Local_Params-class`] object.
#'
#' @return A function of `Delta` returning the corresponding integral error.
#'
#' @noRd

setGeneric("get_integral_error", function(Population_with_size_nmut_old, Population_with_size_nmut, local_params) standardGeneric("get_integral_error"))

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


#' Compute an adapted local time step
#'
#' Computes the local value of `Delta` by inverting the integral error function
#' at a target tolerance.
#'
#' @param Population_with_size_nmut_old A
#'   [`Population_with_size_nmut-class`] object representing the previous state.
#' @param Population_with_size_nmut A [`Population_with_size_nmut-class`]
#'   object representing the updated state.
#' @param local_params A [`Local_Params-class`] object.
#' @param epsilon Numeric tolerance used to determine the local time step.
#'
#' @return A numeric value for the adapted local time step.
#'
#' @noRd

setGeneric("get_Delta", function(Population_with_size_nmut_old, Population_with_size_nmut, local_params, epsilon) standardGeneric("get_Delta"))

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