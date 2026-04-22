#' Parameters class
#'
#' S4 class used to store simulation parameters and derived quantities for the
#' tumor evolution model.
#'
#' @slot functional_effects A vector describing the functional effect associated
#'   with each event.
#' @slot I Number of functional events.
#' @slot lambda Numeric vector of clone-specific proliferation advantages.
#' @slot mu Numeric vector of clone-specific mutation rates.
#' @slot Competition Competition matrix between clones.
#' @slot influence Matrix describing the relative frequencies or influence of
#'   functional events.
#' @slot K Numeric vector of clone-specific carrying capacities.
#' @slot print_time Numeric vector of simulation times to save or print.
#' @slot av_lifespan Average cell lifespan.
#'
#' @name Parameters-class
#' @rdname Parameters-class
setClass("Parameters",
         slots=c(functional_effects="vector",
                 I="numeric",
                 lambda="numeric",
                 mu="numeric",
                 Competition="matrix",
                 influence="matrix",
                 K="numeric",
                 print_time="numeric",
                 av_lifespan="numeric"
         )
)


#' Import simulation parameters from JSON data
#'
#' Builds a Parameters object from input JSON data, computes derived
#' quantities used by the simulation, and saves the resulting objects to a
#' `Parameters.RData` file.
#'
#' @param json_data Input JSON data containing simulation settings, functional
#'   events, initial populations, and mutation information.
#' @param path Character string giving the directory where `Parameters.RData`
#'   will be saved.
#'
#' @return Saves simulation objects to disk.
#'
#' @rdname import_json_par
#' @export
setGeneric("import_json_par", function(json_data, path) standardGeneric("import_json_par"))

setMethod("import_json_par",
          signature(json_data="list",
                    path="character"),
          function(json_data,path){
            length_panel<-json_data$mutableBases
            K_base<-json_data$carryingCapacity
            tmax<-json_data$endingTime
            Np<-json_data$savingCheckpoints
            av_lifespan<-json_data$cellLife
            mut_rate_base<-json_data$mutationRate
            save_day<-json_data$savingFrequency
            Ncellsmax<-json_data$endingSize
            
            functional_effects<-sapply(json_data$functionalEvents,
                                       function(event){event$type})
            names(functional_effects)<-sapply(json_data$functionalEvents,
                                              function(event){event$name})
            
            starting_gen<-json_data$populations$genotype
            starting_fun_eff<-json_data$populations$phenotype
            Ncells_start<-json_data$populations$numCells
            
            s<-vector()
            for(event in json_data$functionalEvents){
              if(event$type=="growth"){
                s<-c(s,event$params$proliferativeAdvantage)
              }
            }
            
            m<-vector()
            for(event in json_data$functionalEvents){
              if(event$type=="mutation"){
                m<-c(m,event$params$mutationalAmplificationFactor)
              }
            }
            
            k<-vector()
            for(event in json_data$functionalEvents){
              if(event$type=="space"){
                k<-c(k,event$params$additionalSpace)
              }
            }
            
            alpha<-list()
            for(event in json_data$functionalEvents){
              if(event$type=="competition"){
                alpha<-append(alpha,list(c(event$params$susceptibility,
                                           event$params$offensiveScore)))
              }
            }
            
            rel_freq<-sapply(json_data$functionalEvents,function(fun_ev){fun_ev$frequency})
            
            I<-length(functional_effects)
            L<-2^I-1
            binary_mat <- number2binary(1:L,I)
            
            mu<-rep(mut_rate_base*length_panel,L)
            if(length(m)>0){
              mu<-matrix(binary_mat[,functional_effects=="mutation"]*rep(m,each=L),ncol=length(m))
              mu[mu==0]<-1
              mu<-apply(mu, 1, prod)*mut_rate_base*length_panel
            }
            
            influence<-matrix(1,I,I)
            diag(influence)<-rel_freq
            
            sel_adv=rep(0,I)
            sel_adv[functional_effects=="growth"]<-s
            lambda=rowSums(t(t(binary_mat)*sel_adv))
            
            alpha_complete<-list()
            alpha_complete[functional_effects=="competition"]<-alpha
            alpha_complete[functional_effects!="competition"]<-rep(list(c(1,1)),sum(functional_effects!="competition"))
            
            vuln_weights <- sapply(alpha_complete, function(a) a[1])
            att_weights  <- sapply(alpha_complete, function(a) a[2])
            
            vuln <- binary_mat %*% diag(vuln_weights)
            att  <- binary_mat %*% diag(att_weights)
            
            if(I==0){
              Competition<-matrix(1,L,L)
            }else{
              Competition<-matrix(0,L,L)
              for(i in 1:L){
                for(j in 1:L){
                  i_vec <- binary_mat[i, ]
                  j_vec <- binary_mat[j, ]
                  i_notin_j_idx <- which(i_vec == 1 & j_vec == 0)
                  j_notin_i_idx <- which(j_vec == 1 & i_vec == 0)
                  Competition[i, j] <- sum(vuln_weights[i_notin_j_idx] - 1) +
                    sum(att_weights[j_notin_i_idx] - 1) + 1
                }
              }
            }
            
            K<-rep(K_base,L)
            if(length(k)>0){
              K<-matrix(binary_mat[,functional_effects=="space"]*rep(k,each=L),ncol=length(k))
              K<-rowSums(K)+rep(K_base,L)
            }
            
            if(is.null(tmax)){tmax<-100*365}
            if(is.null(Np)){
              print_time<-cumsum(rep(save_day,floor(tmax/save_day)))
            }else{
              print_time<-round(seq(from=0,to=tmax, length.out = Np+1), 4)[-1]
            }
            
            parameters<-new("Parameters",
                            functional_effects=functional_effects,
                            I=I,
                            lambda=lambda,
                            mu=mu,
                            Competition=Competition,
                            influence=influence,
                            K=K,
                            print_time=print_time,
                            av_lifespan=av_lifespan
            )
            
            palette<-sapply(json_data$functionalEvents,function(event){event$color})
            names(palette)<-names(functional_effects)
            
            mut_names<-unlist(json_data$mutations)
            
            if(!(dir.exists(path))){dir.create(path)}
            save(list=c("parameters",
                        "starting_gen",
                        "starting_fun_eff",
                        "Ncells_start",
                        "Ncellsmax",
                        "tmax",
                        "palette",
                        "mut_names"),
                 file = paste(path,"/Parameters.RData",sep="")
            )
          }
)


#' Convert an integer to binary representation
#'
#' Converts a non-negative integer to its binary representation as an integer
#' vector.
#'
#' @param decimal_num A non-negative numeric value to convert.
#' @return An integer vector containing the binary representation of
#'   `decimal_num`.
#'
#' @noRd
setGeneric("int_to_binary", function(decimal_num) standardGeneric("int_to_binary"))

setMethod("int_to_binary",
          signature(decimal_num="numeric"),
          function(decimal_num){
            if (decimal_num == 0) {
              return(c(0))
            }
            
            binary_vector <- integer(0)
            
            while (decimal_num > 0) {
              remainder <- decimal_num %% 2
              binary_vector <- c(remainder, binary_vector)
              decimal_num <- decimal_num %/% 2
            }
            
            return(as.integer(binary_vector))
          }
)


#' Convert integers to a binary matrix
#'
#' Converts a vector of decimal numbers to a binary matrix with a fixed number
#' of bits.
#'
#' @param decimal_nums Numeric vector of non-negative integers.
#' @param NofBits Number of bits to use for the binary representation.
#' @return A matrix where each row is the binary representation of the
#'   corresponding element of `decimal_nums`.
#'
#' @noRd
#' @importFrom gmp sizeinbase

setGeneric("number2binary", function(decimal_nums, NofBits) standardGeneric("number2binary"))

setMethod("number2binary",
          signature(decimal_nums="numeric",
                    NofBits="numeric"),
          function(decimal_nums, NofBits){
            binary_matrix<-matrix(0,length(decimal_nums),NofBits)
            for(i in 1:length(decimal_nums)){
              int<-as.integer(decimal_nums[i])
              binary_matrix[i,NofBits:(NofBits-sizeinbase(int,2)+1)]<-rev(int_to_binary(int))
            }
            return(binary_matrix)
          }
)


#' Convert binary representation to decimal numbers
#'
#' Converts binary vectors or matrices into their decimal representation.
#'
#' @param binary_nums A numeric binary vector or a matrix whose columns or rows
#'   encode binary values.
#' @return A numeric value or vector of decimal values corresponding to the
#'   binary input.
#'
#' @noRd

setGeneric("binary2number", function(binary_nums) standardGeneric("binary2number"))

setMethod("binary2number",
          signature(binary_nums="matrix"),
          function(binary_nums){
            decimal_nums<-binary2number(binary_nums[,1])
            if(ncol(binary_nums)>1){
              for(col in 2:ncol(binary_nums)){
                decimal_nums<-c(decimal_nums,binary2number(binary_nums[,col]))
              }
            }
            return(decimal_nums)
          }
)

setMethod("binary2number",
          signature(binary_nums="numeric"),
          function(binary_nums){
            return(sum(2^(0:(length(binary_nums)-1))*rev(binary_nums)))
          }
)


#' Convert phenotype indices to binary-encoded number
#'
#' Converts a phenotype represented by active positions into a binary vector and
#' then into its decimal encoding.
#'
#' @param phenotype Numeric vector of active phenotype indices.
#' @param I Total number of possible phenotype positions.
#' @return A numeric value corresponding to the decimal encoding of the binary
#'   phenotype vector.
#' @noRd

setGeneric("conversion", function(phenotype, I) standardGeneric("conversion"))

setMethod("conversion",
          signature(phenotype="numeric",
                    I="numeric"),
          function(phenotype,I){
            new_vector_bin<-rep(0,I)
            new_vector_bin[phenotype]<-1
            return(binary2number(new_vector_bin))
          }
)

#' Build observed tumor data from simulation output
#'
#' Reads simulation snapshots from disk, filters populations by detection depth,
#' and constructs observed tumor data structures for downstream plotting and
#' analysis.
#'
#' @param path_sim Character string giving the path to the simulation output
#'   directory.
#' @param depth Detection threshold used to retain populations above a minimum
#'   relative abundance.
#' @param parameters A Parameters object containing simulation times
#'   and related settings.
#' @return A list with observed tumor information.
#'
#' @noRd
#' @importFrom tibble tibble
setGeneric("get_obs_tum", function(path_sim, depth, parameters) standardGeneric("get_obs_tum"))

setMethod("get_obs_tum",
          signature(path_sim="character",
                    depth="numeric",
                    parameters="Parameters"),
          function(path_sim,depth,parameters){
            filenames <- list.files(path_sim, full.names = FALSE)
            filenames<-filenames[grepl("Zprovv",filenames)]
            
            tumor<-lapply(filenames,function(filename){
              load(paste(path_sim,filename,sep="/"))
              time_provv<-parameters@print_time[which.min(abs(time_provv-parameters@print_time))]
              setNames(object = list(Zprovv), time_provv)
            })
            tumor<-unlist(tumor,recursive = FALSE)
            tumor<-tumor[order(as.numeric(names(tumor)))]
            
            obs_pop<-unique(unlist(sapply(tumor,
                                          function(tum_time_fix){
                                            Ncells_time_fix<-sapply(tum_time_fix, Ncells)
                                            Fcells_time_fix<-Ncells_time_fix/sum(Ncells_time_fix)
                                            tum_time_fix_filt<-tum_time_fix[Fcells_time_fix>depth]
                                            obs_pop<-lapply(tum_time_fix_filt,Population)
                                            return(obs_pop)
                                          })))
            
            obs_pop_id<-tibble(Population_ID=1:length(obs_pop),Populations=obs_pop)
            
            obs_tumor_tibble_timefix<-lapply(1:length(tumor),function(count){
              tum_time_fix<-tumor[[count]]
              pop_time_fix<-lapply(tum_time_fix,Population)
              obs_tumor_tibble<-tibble(
                time=as.numeric(names(tumor[count])),
                Population_ID=obs_pop_id$Population_ID,
                Ncells=0)
              for(pop_w_s_nmut in tum_time_fix){
                pop<-Population(pop_w_s_nmut)
                ncells<-Ncells(pop_w_s_nmut)
                desc<-sapply(obs_pop,
                             is_descendant,
                             Population_younger=pop)
                
                which_assign<-which.min(sapply(obs_pop,
                                               how_old_descendant,
                                               Population_younger=pop))
                
                obs_tumor_tibble$Ncells[obs_tumor_tibble$Population_ID==which_assign]<-obs_tumor_tibble$Ncells[obs_tumor_tibble$Population_ID==which_assign]+ncells
              }
              return(obs_tumor_tibble%>%filter(Ncells>0))
            }
            )
            obs_tumor_tibble<-bind_rows(obs_tumor_tibble_timefix)
            obs_tumor<-list()
            obs_tumor$obs_tumor_tibble<-obs_tumor_tibble
            obs_tumor$obs_Pop_ID<-obs_pop_id
            
            return(obs_tumor)
          })


#' Compute data for Muller plots
#'
#' Builds a data frame containing clone geometry and metadata needed to draw
#' Muller plots from observed tumor data.
#'
#' @param obs_Pop_ID Tibble mapping observed population identifiers to
#'   populations.
#' @param obs_tumor_tibble Tibble with observed tumor abundances over time.
#' @param functional_effects Vector of functional effect labels.
#' @param freq Logical; if `TRUE`, compute relative abundances, otherwise
#'   absolute abundances.
#' @return A tibble with plotting information.
#'
#' @noRd
#' @import dplyr
#' @import tidyr
#' 
setGeneric("get_muller_plot_info", function(obs_Pop_ID, obs_tumor_tibble, functional_effects, freq) standardGeneric("get_muller_plot_info"))

setMethod("get_muller_plot_info",
          signature(obs_Pop_ID="tbl",
                    obs_tumor_tibble="tbl",
                    functional_effects="vector",
                    freq="logical"),
          function(obs_Pop_ID,obs_tumor_tibble,functional_effects,freq){
            
            
            tryCatch(expr = {
              
              time_of_appearance<-obs_tumor_tibble%>%
                group_by(Population_ID)%>%
                summarise(time_of_appearance=time[1])%>%
                ungroup()
              
              pop<-obs_Pop_ID$Populations
              
              Pop_ID<-obs_Pop_ID$Population_ID%>%unique()
              
              ancestors<-lapply(pop,
                                function(p){
                                  Pop_ID[
                                    sapply(pop,
                                           function(p1){
                                             is_descendant(p,p1)
                                           })
                                  ]
                                })
              
              
              fun_eff<-sapply(pop,function(p){
                paste(names(functional_effects[sort(p@phenotype)]),collapse = ", ")
              })
              
              
              parent<-sapply(pop,
                             function(p){
                               anc<-Pop_ID[
                                 sapply(pop,
                                        function(p1){
                                          is_descendant(p,p1)
                                        })
                               ]
                               
                               return(anc[length(anc)-1])
                             })
              
              daughters_ordered_tbl<-tibble::tibble(pop=Pop_ID,Population_ID=as.vector(parent))%>%
                filter(lengths(Population_ID)>0)%>%
                group_by(Population_ID)%>%
                summarise(daughters=list(pop))
              
              daughters_ordered_tbl$daughters<-lapply(daughters_ordered_tbl$daughters, function(d){
                time_of_appearance%>%
                  filter(Population_ID%in%d)%>%
                  arrange(time_of_appearance)%>%
                  pull(Population_ID)
              })
              
              pop_with_sons<-unlist(daughters_ordered_tbl$Population_ID)
              
              if(freq==FALSE){
                obs_tumor_tibble_clones<-tibble::tibble(Population_ID=Pop_ID,ancestors)%>%
                  tidyr::unnest(ancestors)%>%
                  dplyr::full_join(obs_tumor_tibble)%>%
                  dplyr::group_by(ancestors,time)%>%
                  dplyr::mutate(Ncells_clone=sum(Ncells))%>%
                  dplyr::ungroup()%>%
                  dplyr::select(ancestors,time,Ncells_clone)%>%
                  dplyr::distinct()%>%
                  dplyr::rename("clone"="ancestors")%>%
                  dplyr::left_join(tibble(clone=Pop_ID,fun_eff))
                
                root<-unlist(ancestors[lengths(ancestors)==1])
                
                if(length(root)>1){
                  obs_tumor_tibble_clones<-dplyr::bind_rows(
                    obs_tumor_tibble_clones%>%
                      dplyr::filter(clone%in%root)%>%
                      dplyr::group_by(time)%>%
                      dplyr::summarize(clone=0,
                                       Ncells_clone=sum(Ncells_clone),
                                       fun_eff=NA),
                    obs_tumor_tibble_clones)
                  pop_with_sons<-c(0,pop_with_sons)
                  daughters_ordered_tbl<-dplyr::bind_rows(tibble(Population_ID=0,
                                                                 daughters=list(root)),daughters_ordered_tbl)
                  root<-0
                }
                
                Clones_df<-obs_tumor_tibble_clones%>%
                  dplyr::filter(clone%in%root)%>%
                  dplyr::mutate(y_lower=-Ncells_clone/2,
                                y_upper=Ncells_clone/2,
                                time_appearance=0)
                
                for(p in pop_with_sons){
                  daughters_p<-daughters_ordered_tbl$daughters[daughters_ordered_tbl$Population_ID==p][[1]]
                  n_siblings_d<-length(daughters_p)
                  time_appearance_p<-time_of_appearance$time_of_appearance[time_of_appearance$Population_ID==p]
                  
                  
                  for(d in daughters_p){
                    sibling_number<-which(daughters_p==d)
                    older_siblings<-daughters_p[0:(sibling_number-1)]
                    
                    prop_y_start<-sibling_number/(n_siblings_d+1)
                    
                    time_appearance_d<-time_of_appearance$time_of_appearance[time_of_appearance$Population_ID==d]
                    
                    Ncells_parent_appearance<-obs_tumor_tibble$Ncells[obs_tumor_tibble$Population_ID==p&
                                                                        obs_tumor_tibble$time==time_appearance_d]
                    center_parent<-Clones_df%>%
                      dplyr::filter(clone==p)%>%
                      dplyr::mutate(center=(y_upper+y_lower)/2)%>%
                      dplyr::select(time,center)
                    
                    Ncells_parent<-obs_tumor_tibble[obs_tumor_tibble$Population_ID==p,c("time","Ncells")]%>%
                      rename("Ncells_p"="Ncells")
                    
                    Ncells_parent_clone<-obs_tumor_tibble_clones[obs_tumor_tibble_clones$clone==p,c("time","Ncells_clone")]%>%
                      rename("Ncells_p_clone"="Ncells_clone")
                    
                    Ncells_older_siblings<-obs_tumor_tibble_clones[obs_tumor_tibble_clones$clone%in%older_siblings,]%>%
                      dplyr::group_by(time)%>%
                      dplyr::summarize(Ncells_s=sum(Ncells_clone))
                    
                    d_df<-full_join(Ncells_parent,Ncells_parent_clone,by="time")%>%
                      dplyr::full_join(Ncells_older_siblings,by="time")%>%
                      dplyr::full_join(obs_tumor_tibble_clones%>%
                                         filter(clone==d),by="time")%>%
                      dplyr::full_join(center_parent,by="time")%>%
                      dplyr::filter(!is.na(clone))%>%
                      dplyr::rowwise()%>%
                      dplyr::mutate(y_lower=sum(center,prop_y_start*Ncells_p,-Ncells_p_clone/2,Ncells_s,na.rm=TRUE),
                                    y_upper=y_lower+Ncells_clone,
                                    time_appearance=time_appearance_d)%>%
                      dplyr::select(clone,time,Ncells_clone,"fun_eff",y_lower,y_upper,time_appearance)
                    
                    Clones_df<-rbind(Clones_df,d_df)
                    
                  }
                }
                Clones_df<-Clones_df%>%
                  dplyr::filter(clone!=0)
                
                trasl<-min(Clones_df$y_lower)
                Clones_df$y_lower<-Clones_df$y_lower-trasl
                Clones_df$y_upper<-Clones_df$y_upper-trasl
                Clones_df$time<-as.numeric(Clones_df$time)
                Clones_df$clone<-factor(Clones_df$clone,levels=Pop_ID)
                
                return(Clones_df)
              }
              else{
                obs_tumor_tibble<-obs_tumor_tibble%>%
                  group_by(time)%>%
                  mutate(
                    tot_Ncells=sum(Ncells),
                    Fcells=Ncells/tot_Ncells)
                obs_tumor_tibble_clones<-tibble(Population_ID=Pop_ID,ancestors)%>%
                  unnest(ancestors)%>%
                  full_join(obs_tumor_tibble)%>%
                  group_by(ancestors,time)%>%
                  mutate(Fcells_clone=sum(Fcells))%>%
                  ungroup()%>%
                  dplyr::select(ancestors,time,Fcells_clone)%>%
                  distinct()%>%
                  rename("clone"="ancestors")%>%
                  left_join(tibble(clone=Pop_ID,fun_eff))
                
                root<-unlist(ancestors[lengths(ancestors)==1])
                
                if(length(root)>1){
                  obs_tumor_tibble_clones<-bind_rows(
                    tibble(clone=0,
                           time=unique(obs_tumor_tibble_clones$time),
                           Fcells_clone=1,
                           fun_eff="competition"),
                    obs_tumor_tibble_clones)
                  pop_with_sons<-c(0,pop_with_sons)
                  daughters_ordered_tbl<-bind_rows(tibble(Population_ID=0,
                                                          daughters=list(root)),daughters_ordered_tbl)
                  root<-0
                }
                
                Clones_df<-obs_tumor_tibble_clones%>%
                  filter(clone==root)%>%
                  mutate(y_lower_frac=-Fcells_clone/2,
                         y_upper_frac=Fcells_clone/2)
                
                for(p in pop_with_sons){
                  daughters_p<-daughters_ordered_tbl$daughters[daughters_ordered_tbl$Population_ID==p][[1]]
                  n_siblings_d<-length(daughters_p)
                  
                  for(d in daughters_p){
                    sibling_number<-which(daughters_p==d)
                    older_siblings<-daughters_p[0:(sibling_number-1)]
                    
                    prop_y_start<-sibling_number/(n_siblings_d+1)
                    
                    time_appearance_d<-time_of_appearance$time_of_appearance[time_of_appearance$Population_ID==d]
                    
                    Fcells_parent_appearance<-obs_tumor_tibble$Fcells[obs_tumor_tibble$Population_ID==p&
                                                                        obs_tumor_tibble$time==time_appearance_d]
                    center_parent<-Clones_df%>%
                      filter(clone==p)%>%
                      mutate(center=(y_upper_frac+y_lower_frac)/2)%>%
                      dplyr::select(time,center)
                    
                    Fcells_parent<-obs_tumor_tibble[obs_tumor_tibble$Population_ID==p,c("time","Fcells")]%>%
                      rename("Fcells_p"="Fcells")
                    
                    Fcells_parent_clone<-obs_tumor_tibble_clones[obs_tumor_tibble_clones$clone==p,c("time","Fcells_clone")]%>%
                      rename("Fcells_p_clone"="Fcells_clone")
                    
                    Fcells_older_siblings<-obs_tumor_tibble_clones[obs_tumor_tibble_clones$clone%in%older_siblings,]%>%
                      group_by(time)%>%
                      summarize(Fcells_s=sum(Fcells_clone))
                    
                    d_df<-full_join(Fcells_parent,Fcells_parent_clone,by="time")%>%
                      full_join(Fcells_older_siblings,by="time")%>%
                      full_join(obs_tumor_tibble_clones%>%
                                  filter(clone==d),by="time")%>%
                      full_join(center_parent,by="time")%>%
                      filter(!is.na(clone))%>%
                      rowwise()%>%
                      mutate(y_lower_frac=sum(center,prop_y_start*Fcells_p,-Fcells_p_clone/2,Fcells_s,na.rm=TRUE),
                             y_upper_frac=y_lower_frac+Fcells_clone)%>%
                      dplyr::select(time,clone,Fcells_clone,fun_eff,y_lower_frac,y_upper_frac)
                    
                    Clones_df<-rbind(Clones_df,d_df)
                    
                  }
                }
                
                Clones_df<-Clones_df%>%
                  dplyr::filter(clone!=0)
                Clones_df$y_lower_frac<-Clones_df$y_lower_frac+0.5
                Clones_df$y_upper_frac<-Clones_df$y_upper_frac+0.5
                Clones_df$clone<-factor(Clones_df$clone,levels=Pop_ID)
                Clones_df$time<-as.numeric(Clones_df$time)
                
                return(Clones_df)
              }
              
            },
            error=function(cond) {
              message(conditionMessage(cond))
            })
            
          })


#' Create an on-screen Muller plot
#'
#' Produces a compact Muller plot for visualization in the application.
#'
#' @param Clones_df Tibble containing clone geometry and abundances.
#' @param freq Logical; if `TRUE`, plot relative abundances, otherwise absolute
#'   abundances.
#' @param palette Named vector of colors associated with functional effects.
#' @return A `ggplot2` object.
#'
#' @noRd
#' @import ggplot2
#' @import colorspace
setGeneric("get_muller_plot_show", function(Clones_df, freq, palette) standardGeneric("get_muller_plot_show"))

setMethod("get_muller_plot_show",
          signature(Clones_df="tbl",
                    freq="logical",
                    palette="vector"),
          function(Clones_df,freq,palette){
            
            tryCatch(expr = {
              
              fun_eff<-Clones_df$fun_eff
              
              palette<-palette[sort(unique(fun_eff))]
              
              palette_light<-colorspace::lighten(palette,amount = 0.1)
              names(palette_light)<-names(palette)
              palette_dark<-colorspace::darken(palette,amount = 0.3)
              names(palette_dark)<-names(palette)
              
              
              if(freq==FALSE){
                
                p<-ggplot()+
                  geom_ribbon(data=Clones_df,
                              aes(x=time,
                                  ymin = y_lower,
                                  ymax=y_upper,
                                  group=clone,
                                  color=fun_eff,
                                  fill=fun_eff))+
                  scale_fill_manual(values=palette)+
                  scale_color_manual(values=palette_dark,guide = "none")+
                  scale_x_continuous(breaks=round(seq(min(Clones_df$time),max(Clones_df$time),length.out=5)))+
                  theme_void()+
                  theme(
                    legend.position ="none",
                    plot.margin =  unit(c(-17.5,-33,-17.5,-33), "pt"),
                    axis.ticks.length = unit(0,"cm"),
                    axis.ticks.margin = unit(0,"cm")
                  )
              }
              else{
                p<-ggplot(Clones_df)+
                  geom_ribbon(aes(x=time,
                                  ymin = y_lower_frac,
                                  ymax=y_upper_frac,
                                  group=clone,
                                  fill=fun_eff,
                                  col=fun_eff))+
                  geom_rect(xmin = min(Clones_df$time),
                            xmax=max(Clones_df$time),
                            ymin=0,
                            ymax=1,
                            fill="transparent",
                            color="black")+
                  theme_void()+
                  scale_fill_manual(values=palette)+
                  scale_color_manual(values=palette_dark,guide = "none")+
                  guides(fill=guide_legend(title="Functional effect:",override.aes = list(color = palette_dark)))+
                  theme(
                    legend.position ="none",
                    plot.margin =  unit(c(-17,-32,-17,-32), "pt"),
                    axis.ticks.length = unit(0,"cm"),
                    axis.ticks.margin = unit(0,"cm")
                  )
                
              }
              return(p)
            },
            error=function(cond) {
              message(conditionMessage(cond))
              ggplot() +
                annotate("text",
                         x = 1,
                         y = 1,
                         label = "Error in the muller plot",
                         size = 6,
                         fontface = "bold") +
                theme_void()
            })
          })


#' Create a downloadable Muller plot
#'
#' Produces a Muller plot formatted for export or download.
#'
#' @param Clones_df Tibble containing clone geometry and abundances.
#' @param freq Logical; if `TRUE`, plot relative abundances, otherwise absolute
#'   abundances.
#' @param palette Named vector of colors associated with functional effects.
#' @return A `ggplot2` object ready for export.
#'
#' @noRd
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom scales scientific
#' @import stringr
setGeneric("get_muller_plot_download", function(Clones_df, freq, palette) standardGeneric("get_muller_plot_download"))

setMethod("get_muller_plot_download",
          signature(Clones_df="tbl",
                    freq="logical",
                    palette="vector"),
          function(Clones_df,freq,palette){
            
            tryCatch(expr = {
              
              fun_eff<-Clones_df$fun_eff
              
              palette<-palette[sort(unique(fun_eff))]
              
              palette_light<-colorspace::lighten(palette,amount = 0.1)
              names(palette_light)<-names(palette)
              palette_dark<-colorspace::darken(palette,amount = 0.3)
              names(palette_dark)<-names(palette)
              
              if(freq==FALSE){
                
                latex_label_k<-scales::scientific(max(Clones_df$Ncells_clone), digits = 1)
                latex_label_k<-stringr::str_split(latex_label_k,"e",simplify = TRUE)
                if(latex_label_k[1,1]=="1"){latex_label_k[1,1]<-""}
                if(grepl(x=latex_label_k[1,2],pattern="+")){
                  latex_label_k[1,2]<-gsub(x=latex_label_k[1,2],pattern="+","",fixed = TRUE)}
                latex_label_k[1,2]<-gsub(x=latex_label_k[1,2],pattern="^0+","")
                if(latex_label_k[1,2]==""&latex_label_k[1,1]!=""){
                  latex_label_k<-paste("$",latex_label_k[1,1],"$",sep="")
                }else if(latex_label_k[1,1]==""&latex_label_k[1,2]!=""){
                  latex_label_k<-paste("$10^{",latex_label_k[1,2],"}$",sep="")
                }else if(latex_label_k[1,1]!=""&latex_label_k[1,2]!=""){
                  latex_label_k<-paste(paste(paste(c("$","\\cdot 10^{"),latex_label_k,sep=""),collapse = ""),"}$",sep="")
                }else{
                  latex_label_k<-"$1$"
                }
                
                p<-ggplot()+
                  geom_ribbon(data=Clones_df,
                              aes(x=time,
                                  ymin = y_lower,
                                  ymax=y_upper,
                                  group=clone,
                                  color=fun_eff,
                                  fill=fun_eff))+
                  scale_fill_manual(values=palette)+
                  geom_segment(aes(x=max(Clones_df$time)*(1+1/20),
                                   xend = max(Clones_df$time)*(1+1/20),
                                   y=min(Clones_df$y_lower),
                                   yend=max(Clones_df$y_upper)
                  ),arrow = arrow(ends="both",length = unit(5, "points")))+
                  geom_label(aes(x=max(Clones_df$time)*(1+1/20),
                                 y=(min(Clones_df$y_lower)+max(Clones_df$y_upper))/2),
                             label = latex2exp::TeX(latex_label_k),
                             size=3,fill = "white",
                             label.size = NA)+
                  scale_color_manual(values=palette_dark,guide = "none")+
                  xlab("Days")+
                  ylab("Absolute Aboundance")+
                  scale_x_continuous(breaks=round(seq(min(Clones_df$time),max(Clones_df$time),length.out=5)))+
                  guides(fill=guide_legend(title="Phenotype:",override.aes = list(color = palette_dark)))+
                  theme_void()+
                  theme(
                    axis.title.x = element_text(size=14),
                    axis.title.y = element_text(size=14,angle=90,vjust=2),
                    axis.text.y = element_blank(),
                    plot.margin = unit(c(0,0,0,0.2), "cm"),
                    axis.text.x = element_text(size=12,vjust = 3),
                    legend.position ="bottom",
                    legend.box = "vertical"
                  )
                
                
                return(p)
                
              }
              else{
                
                p<-ggplot(Clones_df)+
                  geom_ribbon(aes(x=time,
                                  ymin = y_lower_frac,
                                  ymax=y_upper_frac,
                                  group=clone,
                                  fill=fun_eff,
                                  col=fun_eff))+
                  geom_rect(xmin = min(Clones_df$time),
                            xmax=max(Clones_df$time),
                            ymin=0,
                            ymax=1,
                            fill="transparent",
                            color="black")+
                  theme_void()+
                  scale_fill_manual(values=palette)+
                  scale_color_manual(values=palette_dark,guide = "none")+
                  labs(x ="Days",
                       y="Relative Aboundance")+
                  guides(fill=guide_legend(title="Phenotype:",override.aes = list(color = palette_dark)))+
                  theme(
                    axis.title.x = element_text(size=14),
                    axis.title.y = element_text(size=14,angle=90,vjust = 2.5),
                    axis.text.y = element_text(size=12,hjust = -5),
                    plot.margin = unit(c(0,0,0,0.2), "cm"),
                    axis.text.x = element_text(size=12,vjust = 3),
                    legend.position = "none"
                  )
                return(p)
                
              }
              
            },
            error=function(cond) {
              message(conditionMessage(cond))
              ggplot() +
                annotate("text",
                         x = 1,
                         y = 1,
                         label = "Error in the muller plot",
                         size = 6,
                         fontface = "bold") +
                theme_void()
            })
          })


#' Compute ordered clone geometry from sequencing data
#'
#' Derives clone hierarchy and plotting coordinates from a list of simulated
#' populations, for use in clone-based visualizations.
#'
#' @param Zprovv A list of population objects.
#' @return A tibble containing clone identifiers, clone sizes, and vertical
#'   positions for plotting.
#'
#' @noRd
#' @import dplyr
setGeneric("get_ordered_clones_sequencing", function(Zprovv) standardGeneric("get_ordered_clones_sequencing"))

setMethod("get_ordered_clones_sequencing",
          signature(Zprovv="list"),
          function(Zprovv){
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
            
            generations <- lengths(gen)
            by_gen <- split(seq_along(pop), generations)
            
            daughters <- vector("list", length(pop))
            names(daughters) <- as.character(seq_along(pop))
            
            for (g in names(by_gen)) {
              g_next <- as.character(as.integer(g) + 1)
              if (!g_next %in% names(by_gen)) next
              
              parents_idx  <- by_gen[[g]]
              children_idx <- by_gen[[g_next]]
              
              for (i in parents_idx) {
                gi <- gen[[i]]
                
                hits <- children_idx[
                  vapply(children_idx,
                         function(j) identical(gen[[j]][-length(gen[[j]])], gi),
                         logical(1))
                ]
                
                if (length(hits)) {
                  daughters[[as.character(i)]] <- hits
                }
              }
            }
            
            get_descendants <- function(p, daughters) {
              res   <- integer(0)
              stack <- p
              
              while (length(stack)) {
                cur <- stack[1]
                stack <- stack[-1]
                
                kids <- daughters[[as.character(cur)]]
                if (is.null(kids) || !length(kids)) next
                
                res   <- c(res, kids)
                stack <- c(stack, kids)
              }
              
              res
            }
            
            
            clones <- lapply(Pop_ID, function(p) {
              c(p, get_descendants(p, daughters))
            })
            
            
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
              daughters[["0"]]<-root
              root<-0
            }
            
            Clones_df<-tibble::tibble(clone=Pop_ID,Ncells_clone=ncells_clone,Ncells=ncells)%>%
              filter(clone%in%root)%>%
              mutate(y_lower=-Ncells_clone/2,
                     y_upper=Ncells_clone/2)
            
            for(p in pop_with_sons){
              daughters_p<-daughters[[as.character(p)]]
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
                
                d_df<-tibble::tibble(
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
              Pop_ID<-Pop_ID[-1]
            }
            
            trasl<-min(Clones_df$y_lower)
            Clones_df$y_lower<-Clones_df$y_lower-trasl
            Clones_df$y_upper<-Clones_df$y_upper-trasl
            Clones_df$clone<-factor(Clones_df$clone,levels=Pop_ID)
            
            Clones_df<-distinct(Clones_df)
            return(Clones_df)
          })

#' Elbowed link: right then up
#'
#' Creates a path connecting two points with a horizontal segment,
#' a rounded corner, and a vertical segment going upwards.
#'
#' @param x1,y1 Numeric. Coordinates of the starting point.
#' @param x2,y2 Numeric. Coordinates of the ending point.
#' @param r Numeric. Radius of the rounded corner.
#' @noRd
#' @return A data.frame with columns `x` and `y` representing the path.
#' @importFrom dplyr bind_rows
elbowed_link_right_up<-function(x1,y1,x2,y2,r){
  
  theta <- seq(-pi/2, 0, length.out = 20)  
  arc <- data.frame(
    x = x2 - r + r * cos(theta),
    y = y1 +r + r * sin(theta)
  )
  
  path <- bind_rows(
    data.frame(x = c(x1, x2-r), y = c(y1, y1)),  
    arc,
    data.frame(x = c(x2, x2), y = c(y1 + r, y2))   
  )
  return(path)
}

#' Elbowed link: right then down
#'
#' Creates a path connecting two points with a horizontal segment,
#' a rounded corner, and a vertical segment going downwards.
#'
#' @param x1,y1 Numeric. Coordinates of the starting point.
#' @param x2,y2 Numeric. Coordinates of the ending point.
#' @param r Numeric. Radius of the rounded corner.
#' @noRd
#' @return A data.frame with columns `x` and `y` representing the path.
#' @importFrom dplyr bind_rows
elbowed_link_right_down<-function(x1,y1,x2,y2,r){
  theta <- seq(pi/2, 0, length.out = 20)
  arc <- data.frame(
    x = x2 -r + r * cos(theta),
    y = y1 -r  + r * sin(theta)
  )
  
  path <- bind_rows(
    data.frame(x = c(x1, x2-r), y = c(y1, y1)),
    arc,
    data.frame(x = c(x2, x2), y = c(y1 -r, y2))
  )
  return(path)
}

#' Elbowed link: up then right
#'
#' Creates a path connecting two points with a vertical segment,
#' a rounded corner, and a horizontal segment to the right.
#'
#' @param x1,y1 Numeric. Coordinates of the starting point.
#' @param x2,y2 Numeric. Coordinates of the ending point.
#' @param r Numeric. Radius of the rounded corner.
#' @noRd
#' @return A data.frame with columns `x` and `y` representing the path.
#' @importFrom dplyr bind_rows
elbowed_link_up_right<-function(x1,y1,x2,y2,r){
  theta <- seq(pi, pi/2, length.out = 20)
  arc <- data.frame(
    x = x1 +r + r * cos(theta),
    y = y2 -r  + r * sin(theta)
  )
  
  path <- bind_rows(
    data.frame(x = c(x1, x1), y = c(y1, y2-r)),
    arc,
    data.frame(x = c(x1+r, x2), y = c(y2, y2))
  )
  return(path)
}

#' Elbowed link: down then right
#'
#' Creates a path connecting two points with a vertical segment,
#' a rounded corner, and a horizontal segment to the right.
#' The vertical segment goes downward.
#'
#' @param x1,y1 Numeric. Coordinates of the starting point.
#' @param x2,y2 Numeric. Coordinates of the ending point.
#' @param r Numeric. Radius of the rounded corner.
#' @noRd
#' @return A data.frame with columns `x` and `y` representing the path.
#' @importFrom dplyr bind_rows
elbowed_link_down_right<-function(x1,y1,x2,y2,r){
  theta <- seq(pi, 3*pi/2, length.out = 20)
  arc <- data.frame(
    x = x1 + r + r * cos(theta),
    y = y2 + r  + r * sin(theta)
  )
  
  path <- bind_rows(
    data.frame(x = c(x1, x1), y = c(y1, y2+r)),
    arc,
    data.frame(x = c(x1+r, x2), y = c(y2, y2))
  )
  return(path)
}

#' Elbowed link between two points
#'
#' Generates a piecewise path between two points using horizontal/vertical
#' segments connected by rounded corners. The path is constructed by
#' splitting at an intermediate x position.
#'
#' @param x1,y1 Numeric. Coordinates of the first point.
#' @param x2,y2 Numeric. Coordinates of the second point.
#' @param r Numeric. Radius of the rounded corners.
#' @param x_mean Numeric. X-coordinate of the intermediate splitting point.
#' @noRd
#' @return A tibble/data.frame with columns `x` and `y` representing the path.
#'
#' @details
#' The function determines the relative position of the two points and
#' selects the appropriate combination of elbowed segments.
#'
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
elbowed_link<-function(x1,y1,x2,y2,r,x_mean){
  x_left<-min(x1,x2)
  x_right<-max(x1,x2)
  y_left<-c(y1,y2)[which.min(c(x1,x2))]
  y_right<-c(y1,y2)[which.max(c(x1,x2))]
  
  x_mid<-x_mean
  y_mid<-mean(c(y1,y2))
  
  if(y_left<y_right){
    link<-bind_rows(elbowed_link_right_up(x_left,y_left,x_mid,y_mid,r),
                    elbowed_link_up_right(x_mid,y_mid,x_right,y_right,r))
  }else if(y_left>y_right){
    link<-bind_rows(elbowed_link_right_down(x_left,y_left,x_mid,y_mid,r),
                    elbowed_link_down_right(x_mid,y_mid,x_right,y_right,r))
  }else{
    link<-tibble(x=c(x_left,x_right),y=c(y_left,y_right))
  }
  return(link)
}

#' Generate a tree plot for mutations
#'
#' Creates a ggplot2 visualization of a mutation tree from a data frame.
#' The tree is rendered using a dendrogram layout, with nodes represented
#' as labeled boxes and edges drawn as elbowed connections.
#'
#' @param df Data frame containing mutation information. Must include at least:
#' \describe{
#'   \item{mut}{Unique mutation identifier}
#'   \item{parents}{Parent mutation identifier (NA for roots)}
#'   \item{names}{Label to display}
#'   \item{fun_eff}{Functional effect (used for coloring)}
#'   \item{mut_generation}{Integer indicating depth in the tree}
#' }
#' @param palette Named vector of colors for `fun_eff`.
#'
#' @return A \code{ggplot} object representing the mutation tree.
#'
#' @details
#' For large trees, only a subset of mutations is displayed. Hidden nodes
#' are aggregated and shown as "n more" labels.
#'
#' If no parent relationships are present, a simple vertical layout is used.
#'
#' Edges are drawn using elbowed connections with rounded corners.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom tibble tibble as_tibble
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph create_layout
#' @importFrom scales rescale
#' @importFrom grid unit
#' @noRd
get_tree_plot_app<-function(df,palette){
  if(all(is.na(df$parents))){
    wanted_mut<-df$mut
    df$y<-1:length(wanted_mut)
    
    plot<-ggplot(df)+
      geom_label(aes(x=0,y=y,label=names,
                     fill=fun_eff),
                 label=df$names,
                 color="white",
                 size=6,
                 label.r =unit(0.5,"lines"),
                 label.size = 0)+
      xlim(-0.1,0.1)+
      ylim(0.9,length(df$mut)+0.1)+
      coord_fixed()+
      scale_fill_manual(values=palette,na.value = "white")+
      theme_void()+
      theme(legend.position = "none")
  }
  else{
    if(length(df$mut)<10){
      wanted_mut<-df$mut
    }
    else{
      wanted_mut<-unname(
        unlist(
          lapply(split(df$mut,
                       df$parents),
                 function(v){
                   v_filt<-v[(v%in%unique(df$parents))]
                   v_filt<-c(v_filt,v[!v%in%v_filt][1])
                   return(v_filt)
                 }
          )
        )
      )
    }
    wanted_mut<-unique(c(wanted_mut,unique(df$mut[is.na(df$parents)])))
    
    mut_info<-df%>%
      mutate(wanted=ifelse(mut%in%wanted_mut,TRUE,FALSE))%>%
      group_by(parents)%>%
      summarise(ndaught_drawn=sum(wanted),
                ndaught_hidden=n()-ndaught_drawn)%>%
      filter(!is.na(parents)&ndaught_hidden>0)%>%
      rowwise()%>%
      mutate(mut=paste(parents,"n",sep="_"),
             label=paste(ndaught_hidden,"more"),
             fun_eff=as.character(NA),
             mut_generation=unique(df$mut_generation[df$mut==parents])+1)%>%
      dplyr::select(mut,parents,label,mut_generation,fun_eff)%>%
      ungroup()
    
    mut_info<-df%>%
      filter(mut%in%wanted_mut)%>%
      mutate(label=names)%>%
      dplyr::select(mut,parents,label,fun_eff,mut_generation)%>%
      bind_rows(mut_info)
    
    g <- graph_from_data_frame(mut_info%>%dplyr::select(parents,mut)%>%filter(!is.na(parents)), directed = TRUE)
    layout_df <- create_layout(g, layout = "dendrogram", circular = FALSE)
    
    x_grid<-length(unique(layout_df$y))
    y_grid<-sum(layout_df$leaf)
    x_range<-c(0,1.2)
    y_range<-c(0,1)
    
    layout_df<-layout_df%>%
      filter(name!=0)%>%
      rowwise()%>%
      mutate(y=x_grid-unique(mut_info$mut_generation[mut_info$mut==name]))%>%
      ungroup()
    
    nodes_coord<-tibble(x=rescale(-layout_df$y,to=x_range),
                        y=rescale(-layout_df$x,to=y_range),
                        mut=layout_df$name)%>%
      merge(mut_info%>%dplyr::select(mut,label,fun_eff,mut_generation))%>%
      group_by(mut_generation)%>%
      mutate(n_mut_layer=n())%>%
      ungroup()
    
    size_label<-min(40/max(nodes_coord$n_mut_layer),6)
    x_lim<-range(nodes_coord$x)+c(-0.1,0.1)
    y_lim<-range(nodes_coord$y)+c(-0.1,0.1)
    
    if(sum(!is.na(unique(layout_df$y)))==1){
      plot<-ggplot()+
        geom_label(data=nodes_coord,
                   aes(x=x,
                       y=y,
                       label=label,
                       fill=fun_eff,
                       color=fun_eff),
                   size=size_label,
                   label.r =unit(0.5,"lines"),
                   label.size = 0)+
        xlim(x_lim)+
        ylim(y_lim)+
        coord_fixed()+
        scale_color_manual(values=rep("white",length(parameters@functional_effects)),na.value = "black")+
        scale_fill_manual(values=palette,na.value = "white")+
        theme_void()+
        theme(legend.position = "none")
    }
    else{
      
      mid_pts<-unique(nodes_coord$x)+x_range[2]/x_grid
      
      edges_coord<-tibble(
        x1=sapply(mut_info$parents[!is.na(mut_info$parents)],
                  function(parent){nodes_coord$x[nodes_coord$mut==parent]}),
        y1=sapply(mut_info$parents[!is.na(mut_info$parents)],
                  function(parent){nodes_coord$y[nodes_coord$mut==parent]}),
        x2=sapply(mut_info$mut[!is.na(mut_info$parents)],
                  function(mut){nodes_coord$x[nodes_coord$mut==mut]}),
        y2=sapply(mut_info$mut[!is.na(mut_info$parents)],
                  function(mut){nodes_coord$y[nodes_coord$mut==mut]}),
        x_mid=sapply(mut_info$mut[!is.na(mut_info$parents)],
                     function(mut){mid_pts[mut_info$mut_generation[mut_info$mut==mut]-1]})
      )
      
      radius<-edges_coord%>%
        rowwise()%>%
        mutate(r=min(abs(x1-x2),abs(y1-y2))/2)%>%
        pull(r)
      
      r<-min(radius[radius>0])
      
      link<-edges_coord%>%
        rowwise()%>%
        mutate(link=list(as_tibble(elbowed_link(x1,y1,x2,y2,r,x_mid))))%>%
        ungroup()%>%
        pull(link)%>%bind_rows(.id="mut")
      
      plot<-ggplot()+
        geom_path(data=link,
                  aes(x=x,y=y,group=mut))+
        geom_label(data=nodes_coord,
                   aes(x=x,
                       y=y,
                       label=label,
                       color=fun_eff,
                       fill=fun_eff),
                   size=size_label,
                   label.r =unit(0.5,"lines"),
                   label.size = 0)+
        xlim(x_lim)+
        ylim(y_lim)+
        coord_fixed()+
        scale_color_manual(values=rep("white",length(parameters@functional_effects)),na.value = "black")+
        scale_fill_manual(values=palette,na.value = "white")+
        theme_void()+
        theme(legend.position = "none")
    }
  }
  return(plot)
}

#' Create split tiles for ggplot-style visualizations
#'
#' Generates polygon coordinates to represent split tiles (two triangular halves)
#' within each grid cell, typically for use with \code{ggplot2::geom_polygon()}.
#'
#' @param df A data frame containing the variables to plot.
#' @param x_col Column for x-axis grouping (unquoted).
#' @param y_col Optional column for y-axis grouping (unquoted). If missing,
#'   all tiles are placed on a single row.
#' @param value_col Column containing values associated with each tile (unquoted).
#'   Not directly used for geometry but typically mapped to aesthetics.
#' @param half_col Column indicating which half of the tile each observation
#'   belongs to (unquoted). Must contain at least two distinct values.
#' @param direction Character. Direction of the split:
#'   \itemize{
#'     \item \code{"bl_tr"}: diagonal from bottom-left to top-right
#'     \item \code{"tl_br"}: diagonal from top-left to bottom-right
#'   }
#' @param half_levels Optional character vector of length 2 specifying which values
#'   in \code{half_col} correspond to the two halves. If \code{NULL}, the first two
#'   unique values are used.
#'
#' @details
#' The function maps each (x, y) combination to a unit square and splits it into
#' two triangular polygons according to \code{half_col}. It returns coordinates
#' suitable for polygon plotting.
#'
#' If more than two values are present in \code{half_col}, only the first two are used.
#' Observations with other values are ignored.
#'
#' @return A tibble containing the original data plus:
#' \describe{
#'   \item{x_label, y_label}{Character labels for axes.}
#'   \item{x_pos, y_pos}{Numeric positions on the grid.}
#'   \item{tile_idx}{Unique identifier for each tile.}
#'   \item{px, py}{Polygon coordinates for plotting.}
#' }
#'
#' @importFrom dplyr mutate pull left_join row_number
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr pmap
#' @importFrom rlang enquo quo_name
#'
#' @examples
#' # Example usage (with ggplot2):
#' # df_poly <- make_split_tiles(df, x, y, value, half)
#' # ggplot(df_poly, aes(px, py, group = interaction(tile_idx, half))) +
#' #   geom_polygon(aes(fill = value))
#'
make_split_tiles <- function(df,
                             x_col, y_col,
                             value_col, half_col,
                             direction = c("bl_tr", "tl_br"),
                             half_levels = NULL) {
  
  direction <- match.arg(direction)
  x_q <- enquo(x_col)
  value_q <- enquo(value_col)
  half_q <- enquo(half_col)
  
  if (!missing(y_col)) {
    y_q <- enquo(y_col)
    
    df2 <- df %>%
      as_tibble() %>%
      mutate(
        x_label = as.character(!!x_q),
        y_label = as.character(!!y_q)
      )
  }else{
    df2 <- df %>%
      as_tibble() %>%
      mutate(
        x_label = as.character(!!x_q),
      )
  }
  
  x_levels <- sort(unique(df%>%pull(!!x_q)))
  if (!missing(y_col)){
    y_levels <- sort(unique(df%>%pull(!!y_q)))
  }
  
  if (!missing(y_col)){
    df2 <- df2 %>%
      mutate(
        x_pos = match(x_label, x_levels),
        y_pos = match(y_label, y_levels)
      )
  }else{
    df2 <- df2 %>%
      mutate(
        x_pos = match(x_label, x_levels),
        y_pos = 1
      )
  }
  
  half_vec <- pull(df2, !!half_q) %>% as.character()
  if (is.null(half_levels)) {
    half_levels <- unique(half_vec)
    if (length(half_levels) < 2)
      stop("`half_col` must contain at least two distinct values.")
    if (length(half_levels) > 2)
      warning("More than two half values found; using first two.")
    half_levels <- half_levels[1:2]
  }
  
  hl1 <- half_levels[1]
  hl2 <- half_levels[2]
  hw <- 0.5
  hh <- 0.5
  
  df2 <- df2 %>% mutate(tile_idx = row_number())
  
  coords_list <- pmap(
    list(x = df2$x_pos, y = df2$y_pos, half = df2[[quo_name(half_q)]], idx = df2$tile_idx),
    function(x, y, half, idx) {
      if (half == hl1) {
        if (direction == "bl_tr") {
          tibble(idx = idx,
                 px = c(x - hw, x + hw, x + hw),
                 py = c(y - hh, y - hh, y + hh))
        } else {
          tibble(idx = idx,
                 px = c(x - hw, x - hw, x + hw),
                 py = c(y - hh, y + hh, y + hh))
        }
      } else if (half == hl2) {
        if (direction == "bl_tr") {
          tibble(idx = idx,
                 px = c(x - hw, x - hw, x + hw),
                 py = c(y - hh, y + hh, y + hh))
        } else {
          tibble(idx = idx,
                 px = c(x - hw, x + hw, x + hw),
                 py = c(y - hh, y - hh, y + hh))
        }
      } else {
        NULL
      }
    }
  )
  
  coords_df <- bind_rows(coords_list)
  
  out <- df2 %>%
    left_join(coords_df, by = c("tile_idx" = "idx"))
  
  return(out)
}