#' Population class
#'
#' S4 class representing a population, including its genotype, functional
#' effects, and phenotype.
#'
#' @slot genotype Integer vector describing the genotype of the population.
#' @slot functional_effect Integer vector describing the functional effects
#'   associated with the genotype.
#' @slot phenotype Integer vector describing the phenotype of the population.
#'
#' @name Population-class
#' @rdname Population-class
#' @exportClass Population
setClass("Population",
         slots=c(genotype="vector",
                 functional_effect="vector",
                 phenotype="vector"),
         prototype=list(genotype=vector(mode="integer"),
                        functional_effect=vector(mode="integer"),
                        phenotype=vector(mode="integer"))
)


#' Get the genotype of a population
#'
#' Returns the genotype stored in a [`Population-class`] object.
#'
#' @param Population A [`Population-class`] object.
#' @return A vector representing the genotype.
#'
#' @rdname genotype
#' @export
setGeneric("genotype", function(Population) standardGeneric("genotype"))

setMethod("genotype",
          "Population",
          function(Population){
            return(Population@genotype)
          })


#' Get the functional effects of a population
#'
#' Returns the functional effects stored in a [`Population-class`] object.
#'
#' @param Population A [`Population-class`] object.
#' @return A vector representing the functional effects.
#'
#' @rdname functional_effect
#' @export
setGeneric("functional_effect", function(Population) standardGeneric("functional_effect"))

setMethod("functional_effect",
          "Population",
          function(Population){
            return(Population@functional_effect)
          })


#' Get the phenotype of a population
#'
#' Returns the sorted phenotype stored in a [`Population-class`] object.
#'
#' @param Population A [`Population-class`] object.
#' @return A sorted vector representing the phenotype.
#'
#' @rdname phenotype
#' @export
setGeneric("phenotype", function(Population) standardGeneric("phenotype"))

setMethod("phenotype",
          "Population",
          function(Population){
            return(sort(Population@phenotype))
          })


#' Compare two populations for equality
#'
#' Checks whether two populations have the same genotype.
#'
#' @param Population1 A [`Population-class`] object.
#' @param Population2 A [`Population-class`] object.
#' @return `TRUE` if the two populations have identical genotypes, `FALSE`
#'   otherwise.
#'
#' @noRd
setGeneric("is_equal_pop", function(Population1, Population2) standardGeneric("is_equal_pop"))
setMethod("is_equal_pop",
          signature("Population",
                    "Population"),
          function(Population1,Population2){
            gen1<-genotype(Population1)
            gen2<-genotype(Population2)
            if(length(gen1)!=length(gen2)){return(FALSE)}
            else{return(all(gen1==gen2))}
          })


#' Get the parent population
#'
#' Returns the parent of a population by removing the last genotype and
#' functional effect entry.
#'
#' @param Population A [`Population-class`] object.
#' @return A [`Population-class`] object corresponding to the parent population.
#'
#' @rdname parent
#' @export
setGeneric("parent", function(Population) standardGeneric("parent"))

setMethod("parent",
          "Population",
          function(Population){
            parent_genotype<-Population@genotype[-length(Population@genotype)]
            parent_functional_effect<-Population@functional_effect[-length(Population@functional_effect)]
            parent_phenotype<-unique(parent_functional_effect)
            parent<-new("Population",
                        genotype=parent_genotype,
                        functional_effect=parent_functional_effect,
                        phenotype=parent_phenotype)
            return(parent)
          })


#' Check whether one population is the direct parent of another
#'
#' Tests whether `Population_older` is the direct parent of
#' `Population_younger`.
#'
#' @param Population_younger A [`Population-class`] object.
#' @param Population_older A [`Population-class`] object.
#' @return `TRUE` if `Population_older` is the direct parent of
#'   `Population_younger`, `FALSE` otherwise.
#'
#' @noRd
setGeneric("is_parent", function(Population_younger, Population_older) standardGeneric("is_parent"))

setMethod("is_parent",
          signature("Population",
                    "Population"),
          function(Population_younger,Population_older){
            younger_genotype<-Population_younger@genotype
            younger_generation<-length(younger_genotype)
            older_genotype<-Population_older@genotype
            older_generation<-length(older_genotype)
            if(younger_generation!=(older_generation+1)){return(FALSE)}
            else{return(all(genotype(parent(Population_younger))==genotype(Population_older)))}
          })


#' Check whether one population is a descendant of another
#'
#' Tests whether `Population_younger` descends from `Population_older`.
#'
#' @param Population_younger A [`Population-class`] object.
#' @param Population_older A [`Population-class`] object.
#' @return `TRUE` if `Population_younger` is a descendant of `Population_older`,
#'   `FALSE` otherwise.
#'
#' @noRd
setGeneric("is_descendant", function(Population_younger, Population_older) standardGeneric("is_descendant"))

setMethod("is_descendant",
          signature("Population",
                    "Population"),
          function(Population_younger,Population_older){
            younger_genotype<-Population_younger@genotype
            younger_generation<-length(younger_genotype)
            older_genotype<-Population_older@genotype
            older_generation<-length(older_genotype)
            
            if(older_generation>younger_generation){return(FALSE)}
            else{return(all(younger_genotype[1:older_generation]==older_genotype))}
          })


#' Measure descendant distance
#'
#' Returns how many generations separate `Population_younger` from
#' `Population_older` when the former is a descendant of the latter.
#'
#' @param Population_younger A [`Population-class`] object.
#' @param Population_older A [`Population-class`] object.
#' @return A non-negative integer if `Population_younger` is a descendant of
#'   `Population_older`; otherwise `Inf`.
#'
#' @noRd

setGeneric("how_old_descendant", function(Population_younger, Population_older) standardGeneric("how_old_descendant"))

setMethod("how_old_descendant",
          signature("Population",
                    "Population"),
          function(Population_younger,Population_older){
            
            younger_genotype<-Population_younger@genotype
            younger_generation<-length(younger_genotype)
            older_genotype<-Population_older@genotype
            older_generation<-length(older_genotype)
            
            if(older_generation>younger_generation){return(Inf)}
            
            is_desc<-all(younger_genotype[1:older_generation]==older_genotype)
            if(is_desc){
              return(younger_generation-older_generation)
            }else{return(Inf)}
          })


#' Check whether one population is an ancestor of another
#'
#' Tests whether `Population_older` is an ancestor of `Population_younger`.
#'
#' @param Population_older A [`Population-class`] object.
#' @param Population_younger A [`Population-class`] object.
#' @return `TRUE` if `Population_older` is an ancestor of `Population_younger`,
#'   `FALSE` otherwise.
#'
#' @noRd
setGeneric("is_ancestor", function(Population_older, Population_younger) standardGeneric("is_ancestor"))

setMethod("is_ancestor",
          signature("Population",
                    "Population"),
          function(Population_older,Population_younger){
            younger_genotype<-Population_younger@genotype
            younger_generation<-length(younger_genotype)
            older_genotype<-Population_older@genotype
            older_generation<-length(older_genotype)
            if(older_generation>younger_generation){return(FALSE)}
            else{return(all(younger_genotype[1:older_generation]==older_genotype))}
          })


#' Get phenotype labels
#'
#' Returns human-readable labels for the phenotype of a population based on a
#' vector of functional effects.
#'
#' @param pop A [`Population-class`] object.
#' @param functional_effects A vector describing functional effects.
#' @return A character vector of phenotype labels.
#'
#' @noRd

setGeneric("get_phenotype_label", function(pop, functional_effects) standardGeneric("get_phenotype_label"))

setMethod("get_phenotype_label",
          signature("Population",
                    "vector"),
          function(pop,functional_effects){
            fun_eff<-sort(phenotype(pop))
            if(!is.null(names(functional_effects))){
              return(names(functional_effects[fun_eff]))
            }else{
              fun_eff_name<-functional_effects[fun_eff]
              fun_eff_num_int<-vector()
              for(i in 1:length(fun_eff)){
                if(sum(functional_effects==fun_eff_name[i])==1){
                  fun_eff_num_int<-c(fun_eff_num_int,"")
                }
                else{
                  fun_eff_num_int<-c(fun_eff_num_int,sum((functional_effects==fun_eff_name[i])[1:fun_eff[i]]))
                }
              }
              return(paste(fun_eff_name,
                           fun_eff_num_int,
                           sep=""))
            }
            
          })


#' Get functional effect labels
#'
#' Returns human-readable labels for the functional effects of a population.
#'
#' @param pop A [`Population-class`] object.
#' @param functional_effects A vector describing functional effects.
#' @return A character vector of functional effect labels.
#'
#' @noRd
setGeneric("get_fun_eff_label", function(pop, functional_effects) standardGeneric("get_fun_eff_label"))

setMethod("get_fun_eff_label",
          signature("Population",
                    "vector"),
          function(pop,functional_effects){
            fun_eff<-sort(functional_effect(pop))
            
            if(!is.null(names(functional_effects))){
              return(names(functional_effects)[fun_eff])
            }else{
              fun_eff_name<-functional_effects[fun_eff]
              fun_eff_num_int<-vector()
              for(i in 1:length(fun_eff)){
                if(sum(functional_effects==fun_eff_name[i])==1){
                  fun_eff_num_int<-c(fun_eff_num_int,"")
                }
                else{
                  fun_eff_num_int<-c(fun_eff_num_int,sum((functional_effects==fun_eff_name[i])[1:fun_eff[i]]))
                }
              }
              return(paste(fun_eff_name,
                           fun_eff_num_int,
                           sep=""))}
            
          })