setClass("Population",
         slots=c(genotype="vector",
                 functional_effect="vector",
                 phenotype="vector"),
         prototype=list(genotype=vector(mode="integer"),
                        functional_effect=vector(mode="integer"),
                        phenotype=vector(mode="integer"))
)

setGeneric("genotype",function(Population) standardGeneric("genotype"))
setMethod("genotype",
          "Population",
          function(Population){
            return(Population@genotype)
          })

setGeneric("functional_effect",function(Population) standardGeneric("functional_effect"))
setMethod("functional_effect",
          "Population",
          function(Population){
            return(Population@functional_effect)
          })

setGeneric("phenotype",function(Population) standardGeneric("phenotype"))
setMethod("phenotype",
          "Population",
          function(Population){
            return(sort(Population@phenotype))
          })

setGeneric("is_equal_pop",function(Population1,Population2) standardGeneric("is_equal_pop"))
setMethod("is_equal_pop",
          signature("Population",
                    "Population"),
          function(Population1,Population2){
            gen1<-genotype(Population1)
            gen2<-genotype(Population2)
            if(length(gen1)!=length(gen2)){return(FALSE)}
            else{return(all(gen1==gen2))}
          })

setGeneric("parent",function(Population) stantardGeneric("parent"))
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

setGeneric("is_parent",function(Population_younger,Population_older) stantardGeneric("is_parent"))
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


setGeneric("is_descendant",function(Population_younger,Population_older) stantardGeneric("is_descendant"))
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

setGeneric("how_old_descendant",function(Population_younger,Population_older) stantardGeneric("is_descendant"))
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


setGeneric("is_ancestor",function(Population_older,Population_younger) stantardGeneric("is_ancestor"))
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

setGeneric("get_phenotype_label",function(pop,functional_effects) stantardGeneric("get_phenotype_label"))
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
setGeneric("get_fun_eff_label",function(pop,functional_effects) stantardGeneric("get_fun_eff_label"))
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
