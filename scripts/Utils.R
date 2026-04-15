setClassUnion("numericOrNULL", c("numeric", "NULL"))


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
                # tmax="numericOrNULL",
                # Ncellsmax="numericOrNULL"
                )
)

setGeneric("import_json_par",function(json_data,path) standardGeneric("import_json_par"))
setMethod("import_json_par",
          signature(json_data="character",
                    path="character"),
          import_json_par<-function(json_data,path){
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
              mu<-rowProds(mu)*mut_rate_base*length_panel
            }
            
            influence<-matrix(1,I,I)
            diag(influence)<-rel_freq
            
            sel_adv=rep(0,I)
            sel_adv[functional_effects=="growth"]<-s
            lambda=rowSums(t(t(binary_mat)*sel_adv)) # vantaggio selettivo di ogni clone (somma vantaggi sue driver)
            
            
            alpha_complete<-list()
            alpha_complete[functional_effects=="competition"]<-alpha
            alpha_complete[functional_effects!="competition"]<-rep(list(c(1,1)),sum(functional_effects!="competition"))
            
            vuln_weights <- sapply(alpha_complete, function(a) a[1])
            att_weights  <- sapply(alpha_complete, function(a) a[2])
            
            vuln <- binary_mat %*% diag(vuln_weights)
            att  <- binary_mat %*% diag(att_weights)  
            
            if(I==0){
              Competition<-matrix(1,L,L)
            }else{Competition<-matrix(0,L,L)
            for(i in 1:L){
              for(j in 1:L){
                i_vec <- binary_mat[i, ]
                j_vec <- binary_mat[j, ]
                i_notin_j_idx <- which(i_vec == 1 & j_vec == 0)
                j_notin_i_idx <- which(j_vec == 1 & i_vec == 0)
                Competition[i, j] <- sum(vuln_weights[i_notin_j_idx] - 1) +
                  sum(att_weights[j_notin_i_idx] - 1) + 1
              }
            }}
            
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


setGeneric("int_to_binary",function(decimal_num) standardGeneric("int_to_binary"))
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

setGeneric("number2binary",function(decimal_nums, NofBits) standardGeneric("number2binary"))
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

setGeneric("binary2number",function(binary_nums) standardGeneric("binary2number"))
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

setGeneric("conversion",function(phenotype,I) standardGeneric("conversion"))
setMethod("conversion",
          signature(phenotype="numeric",
                    I="numeric"),
          function(phenotype,I){
            new_vector_bin<-rep(0,I)
            new_vector_bin[phenotype]<-1
            return(binary2number(new_vector_bin))
          }
)

setGeneric("get_obs_tum",function(path_sim,depth,parameters) standardGeneric("get_obs_tum"))
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

setGeneric("get_muller_plot_info",function(obs_Pop_ID,obs_tumor_tibble,functional_effects,freq) standardGeneric("get_muller_plot_info"))
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
              
              daughters_ordered_tbl<-tibble(pop=Pop_ID,Population_ID=as.vector(parent))%>%
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
                obs_tumor_tibble_clones<-tibble(Population_ID=Pop_ID,ancestors)%>%
                  unnest(ancestors)%>%
                  full_join(obs_tumor_tibble)%>%
                  group_by(ancestors,time)%>%
                  mutate(Ncells_clone=sum(Ncells))%>%
                  ungroup()%>%
                  dplyr::select(ancestors,time,Ncells_clone)%>%
                  distinct()%>%
                  rename("clone"="ancestors")%>%
                  left_join(tibble(clone=Pop_ID,fun_eff))
                
                
                root<-unlist(ancestors[lengths(ancestors)==1])

                if(length(root)>1){
                  obs_tumor_tibble_clones<-bind_rows(
                    obs_tumor_tibble_clones%>%
                      filter(clone%in%root)%>%
                      group_by(time)%>%
                      summarize(clone=0,
                                Ncells_clone=sum(Ncells_clone),
                                fun_eff=NA),
                    obs_tumor_tibble_clones)
                  pop_with_sons<-c(0,pop_with_sons)
                  daughters_ordered_tbl<-bind_rows(tibble(Population_ID=0,
                                                          daughters=list(root)),daughters_ordered_tbl)
                  root<-0
                }
                
                Clones_df<-obs_tumor_tibble_clones%>%
                  filter(clone%in%root)%>%
                  mutate(y_lower=-Ncells_clone/2,
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
                      filter(clone==p)%>%
                      mutate(center=(y_upper+y_lower)/2)%>%
                      dplyr::select(time,center)
                    
                    Ncells_parent<-obs_tumor_tibble[obs_tumor_tibble$Population_ID==p,c("time","Ncells")]%>%
                      rename("Ncells_p"="Ncells")
                    
                    Ncells_parent_clone<-obs_tumor_tibble_clones[obs_tumor_tibble_clones$clone==p,c("time","Ncells_clone")]%>%
                      rename("Ncells_p_clone"="Ncells_clone")
                    
                    Ncells_older_siblings<-obs_tumor_tibble_clones[obs_tumor_tibble_clones$clone%in%older_siblings,]%>%
                      group_by(time)%>%
                      summarize(Ncells_s=sum(Ncells_clone))
                    
                    d_df<-full_join(Ncells_parent,Ncells_parent_clone,by="time")%>%
                      full_join(Ncells_older_siblings,by="time")%>%
                      full_join(obs_tumor_tibble_clones%>%
                                  filter(clone==d),by="time")%>%
                      full_join(center_parent,by="time")%>%
                      filter(!is.na(clone))%>%
                      rowwise()%>%
                      mutate(y_lower=sum(center,prop_y_start*Ncells_p,-Ncells_p_clone/2,Ncells_s,na.rm=TRUE),
                             y_upper=y_lower+Ncells_clone,
                             time_appearance=time_appearance_d)%>%
                      dplyr::select(clone,time,Ncells_clone,"fun_eff",y_lower,y_upper,time_appearance)
                    
                    Clones_df<-rbind(Clones_df,d_df)
                    
                  }
                }
                Clones_df<-Clones_df%>%
                  filter(clone!=0)
                
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
                  filter(clone!=0)
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


setGeneric("get_muller_plot_show",function(Clones_df,freq,palette) standardGeneric("get_muller_plot_show"))
setMethod("get_muller_plot_show",
          signature(Clones_df="tbl",
                    freq="logical",
                    palette="vector"),
          function(Clones_df,freq,palette){
            
            tryCatch(expr = {
              
              fun_eff<-Clones_df$fun_eff
              
              palette<-palette[sort(unique(fun_eff))]
              
              palette_light<-lighten(palette,amount = 0.1)
              names(palette_light)<-names(palette)
              palette_dark<-darken(palette,amount = 0.3)
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

setGeneric("get_muller_plot_download",function(Clones_df,freq,palette) standardGeneric("get_muller_plot_download"))
setMethod("get_muller_plot_download",
          signature(Clones_df="tbl",
                    freq="logical",
                    palette="vector"),
          function(Clones_df,freq,palette){
            
            tryCatch(expr = {
              
              fun_eff<-Clones_df$fun_eff
              
              palette<-palette[sort(unique(fun_eff))]
              
              palette_light<-lighten(palette,amount = 0.1)
              names(palette_light)<-names(palette)
              palette_dark<-darken(palette,amount = 0.3)
              names(palette_dark)<-names(palette)
              
              if(freq==FALSE){

                latex_label_k<-scales::scientific(max(Clones_df$Ncells_clone), digits = 1)
                latex_label_k<-str_split(latex_label_k,"e",simplify = TRUE)
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
                             label = TeX(latex_label_k),
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

elbowed_link_right_down<-function(x1,y1,x2,y2,r){
  theta <- seq(pi/2, 0, length.out = 20)  # 0 to 90 degrees
  arc <- data.frame(
    x = x2 -r + r * cos(theta),
    y = y1 -r  + r * sin(theta)
  )
  
  path <- bind_rows(
    data.frame(x = c(x1, x2-r), y = c(y1, y1)),  # horizontal
    arc,
    data.frame(x = c(x2, x2), y = c(y1 -r, y2))   # vertical
  )
  return(path)
}

elbowed_link_up_right<-function(x1,y1,x2,y2,r){
  theta <- seq(pi, pi/2, length.out = 20)  # 0 to 90 degrees
  arc <- data.frame(
    x = x1 +r + r * cos(theta),
    y = y2 -r  + r * sin(theta)
  )
  
  path <- bind_rows(
    data.frame(x = c(x1, x1), y = c(y1, y2-r)),  # horizontal
    arc,
    data.frame(x = c(x1+r, x2), y = c(y2, y2))   # vertical
  )
  return(path)
}

elbowed_link_down_right<-function(x1,y1,x2,y2,r){
  theta <- seq(pi, 3*pi/2, length.out = 20)  # 0 to 90 degrees
  arc <- data.frame(
    x = x1 + r + r * cos(theta),
    y = y2 + r  + r * sin(theta)
  )
  
  path <- bind_rows(
    data.frame(x = c(x1, x1), y = c(y1, y2+r)),  # horizontal
    arc,
    data.frame(x = c(x1+r, x2), y = c(y2, y2))   # vertical
  )
  return(path)
}

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
        #label=paste("mut",row_number(),sep=""))%>%
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
    
    x_lim<-range(nodes_coord$x)+c(-0.1,0.1)
    y_lim<-range(nodes_coord$y)+c(-0.1,0.1)
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


setGeneric("get_ordered_clones_sequencing",function(Zprovv) standardGeneric("get_ordered_clones_sequencing"))
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
            
            Clones_df<-tibble(clone=Pop_ID,Ncells_clone=ncells_clone,Ncells=ncells)%>%
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
              Pop_ID<-Pop_ID[-1]
            }
            
            trasl<-min(Clones_df$y_lower)
            Clones_df$y_lower<-Clones_df$y_lower-trasl
            Clones_df$y_upper<-Clones_df$y_upper-trasl
            Clones_df$clone<-factor(Clones_df$clone,levels=Pop_ID)
            
            Clones_df<-distinct(Clones_df)
            return(Clones_df)
          })

setGeneric("get_tree_plot_download",function(path_sim,depth,parameters,palette) standardGeneric("get_tree_plot_download"))
setMethod("get_tree_plot_download",
          signature(path_sim="character",depth="numeric",palette="character",parameters="Parameters"),
          function(path_sim,depth,parameters,palette){
  Zprovvs<-list.files(path_sim)
  num_seq<-which.max(as.numeric(str_remove(str_remove(Zprovvs,"Zprovv"),".RData")))
  load(paste(path_sim,Zprovvs[num_seq],sep="/"))
  time_provv<-parameters@print_time[which.min(abs(time_provv-parameters@print_time))]
  
  pop<-lapply(Zprovv,Population)
  gen<-lapply(pop,genotype)
  fun_eff_num<-lapply(pop,functional_effect)
  fun_eff<-lapply(fun_eff_num,function(n){names(parameters@functional_effects)[n]})
  ncells<-sapply(Zprovv,Ncells)
  tot_ncells<-sum(ncells)
  unique_mut_id<-lapply(gen,function(g){
    unique_id_mut<-vector()
    for(j in 1:length(g)){
      unique_id_mut<-c(unique_id_mut,paste(g[1:j],collapse="_"))
    }
    return(unique_id_mut)
  })
  
  all_mut<-unique(unlist(unique_mut_id))
  used_nums<-as.numeric(gsub(pattern = "Mut",x = names(mut_names),replacement = ""))
  mut_nums<-1:length(all_mut)
  not_used_mut_nums<-mut_nums[!mut_nums%in%used_nums]
  mut_names_tbl<-tibble(mut=c(mut_names,all_mut[!all_mut%in%mut_names]),
                        names=c(names(mut_names),paste0("Mut",not_used_mut_nums,sep="",recycle0 = TRUE)))
  
  parents<-lapply(unique_mut_id,lag)
  mut_generation<-lapply(gen, function(g){1:length(g)})
  
  composition<-tibble(mut=unique_mut_id,parents,fun_eff,ncells,mut_generation)%>%
    unnest(c(mut,parents,fun_eff,mut_generation))%>%
    group_by(mut,parents,fun_eff,mut_generation)%>%
    summarise(ncells=sum(ncells),
              frequency=ncells/tot_ncells)%>%
    filter(frequency>depth)%>%
    ungroup()%>%
    arrange(desc(frequency))%>%
    merge(mut_names_tbl)
  
  
  roots<-composition$mut[is.na(composition$parents)]
  
  
  if(length(roots)>1){
    composition<-composition%>%
      mutate(parents=ifelse(is.na(parents),"0",parents))%>%
      bind_rows(tibble(mut="0"))
  }
  
  if(all(is.na(composition$parents))){
    wanted_mut<-composition$mut
    composition$y<-1:length(wanted_mut)
    
    plot<-ggplot(composition)+
      geom_label(aes(x=0,y=y,label=names,
                     fill=fun_eff),
                 label=composition$names,
                 color="white",
                 size=6,
                 label.r =unit(0.5,"lines"),
                 label.size = 0)+
      xlim(-0.1,0.1)+
      ylim(0.9,length(composition$mut)+0.1)+
      coord_fixed()+
      scale_fill_manual(values=palette,na.value = "white")+
      theme_void()+
      theme(legend.position = "none")
  }else{
    if(length(composition$mut)<10){
      wanted_mut<-composition$mut
    }else{
      wanted_mut<-unname(
        unlist(
          lapply(split(composition$mut,
                       composition$parents),
                 function(v){
                   v_filt<-v[(v%in%unique(composition$parents))]
                   v_filt<-c(v_filt,v[!v%in%v_filt][1])
                   if(length(v_filt)==length(v)-1){v_filt<-v}
                   return(v_filt)
                 }
          )
        )
      )
    }
    wanted_mut<-unique(c(wanted_mut,unique(composition$mut[is.na(composition$parents)])))
    
    mut_info<-composition%>%
      mutate(wanted=ifelse(mut%in%wanted_mut,TRUE,FALSE))%>%
      group_by(parents)%>%
      summarise(ndaught_drawn=sum(wanted),
                ndaught_hidden=n()-ndaught_drawn)%>%
      filter(!is.na(parents)&ndaught_hidden>0)%>%
      rowwise()%>%
      mutate(mut=paste(parents,"n",sep="_"),
             label=paste(ndaught_hidden,"more"),
             fun_eff=as.character(NA),
             mut_generation=unique(composition$mut_generation[composition$mut==parents])+1)%>%
      dplyr::select(mut,parents,label,mut_generation,fun_eff)%>%
      ungroup()
    
    mut_info<-composition%>%
      filter(mut%in%wanted_mut)%>%
      mutate(label=names)%>%
      #label=paste("mut",row_number(),sep=""))%>%
      dplyr::select(mut,parents,label,fun_eff,mut_generation)%>%
      bind_rows(mut_info)
    
    Zprovvs_ordered<-Zprovvs
    Zprovvs_ordered[as.numeric(str_remove(str_remove(Zprovvs,"Zprovv"),".RData"))+1]<-Zprovvs
    
    mut_info$time_appearance<-sapply(1:nrow(mut_info),function(j){
      mut_id<-mut_info$mut[j]
      if(grepl(x = mut_info$label[j],pattern = "more")){
        return(NA)
      }
      t<-time_provv
      for(k in 1:length(Zprovvs_ordered)){
        Z<-rev(Zprovvs_ordered)[k]
        load(paste(path_sim,Z,sep="/"))
        pop<-lapply(Zprovv,Population)
        gen<-lapply(pop,genotype)
        unique_mut_id<-lapply(gen,function(g){
          unique_id_mut<-vector()
          for(j in 1:length(g)){
            unique_id_mut<-c(unique_id_mut,paste(g[1:j],collapse="_"))
          }
          return(unique_id_mut)
        })
        
        if(!mut_id%in%unique(unlist(unique_mut_id))){break}else{
          t<-time_provv
        }
      }
      return(t)
    })
    NA_rows<-which(is.na(mut_info$time_appearance))
    for(j in NA_rows){
      parent<-mut_info$parents[j]
      mut_info$time_appearance[j]<-min(mut_info$time_appearance[-NA_rows][mut_info$parents[-NA_rows]==parent],na.rm = TRUE)
    }
    
    g <- graph_from_data_frame(mut_info%>%dplyr::select(parents,mut)%>%filter(!is.na(parents)), directed = TRUE)
    layout_composition <- create_layout(g, layout = "dendrogram", circular = FALSE)
    
    x_grid<-length(unique(layout_composition$y))
    y_grid<-sum(layout_composition$leaf)
    x_range<-c(0,1.5)
    y_range<-c(0,1)
    
    layout_composition<-layout_composition%>%
      filter(name!=0)%>%
      rowwise()%>%
      mutate(y=x_grid-unique(mut_info$mut_generation[mut_info$mut==name]))%>%
      ungroup()
    
    nodes_coord <- tibble(
      mut = layout_composition$name,
      y = rescale(-layout_composition$x, to = y_range)
    ) %>%
      left_join(
        mut_info %>% dplyr::select(mut, time_appearance),
        by = "mut"
      ) %>%
      mutate(
        x = rescale(time_appearance, to = x_range)
      ) %>%
      left_join(
        mut_info %>% dplyr::select(mut, label, fun_eff, mut_generation),
        by = "mut"
      ) %>%
      group_by(mut_generation) %>%
      mutate(n_mut_layer = n()) %>%
      ungroup()
    
    
    size_label<-min(40/max(nodes_coord$n_mut_layer),6)
    x_lim<-range(nodes_coord$x)+c(-0.1,0.1)
    y_lim<-range(nodes_coord$y)+c(-0.1,0.1)
    
    if(sum(!is.na(unique(layout_composition$y)))==1){
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
    }else{
      
      
      edges_coord<-tibble(
        x1=sapply(mut_info$parents[!is.na(mut_info$parents)],
                  function(parent){nodes_coord$x[nodes_coord$mut==parent]}),
        y1=sapply(mut_info$parents[!is.na(mut_info$parents)],
                  function(parent){nodes_coord$y[nodes_coord$mut==parent]}),
        x2=sapply(mut_info$mut[!is.na(mut_info$parents)],
                  function(mut){nodes_coord$x[nodes_coord$mut==mut]}),
        y2=sapply(mut_info$mut[!is.na(mut_info$parents)],
                  function(mut){nodes_coord$y[nodes_coord$mut==mut]}),
        x_mid=sapply(mut_info$parents[!is.na(mut_info$parents)],
                     function(parent){
                       daught<-mut_info$mut[!is.na(mut_info$parents)&mut_info$parents==parent]
                       min_dist_daught<-min(nodes_coord$x[nodes_coord$mut%in%daught])
                       return((min_dist_daught+nodes_coord$x[nodes_coord$mut==parent])/2)})
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
      
      x_lim<-range(nodes_coord$x)+c(-0.1,0.1)
      y_lim<-range(nodes_coord$y)+c(-0.1,0.1)
      plot<-ggplot()+
        geom_path(data=link,
                  aes(x=x,y=y,group=mut))+
        geom_point(data=nodes_coord%>%filter(!is.na(fun_eff)),
                   aes(x=x,y=y,fill=fun_eff),color="black",shape=21,size=5)+
        geom_label(data=nodes_coord%>%filter(is.na(fun_eff)),
                   aes(x=x,y=y,label=label),color="black",border.colour = "white",fill="white")+
        
        xlim(x_lim)+
        ylim(y_lim)+
        coord_fixed()+
        scale_fill_manual(values=palette,na.value = "white")+
        theme_void()+
        theme(legend.position = "none")
      
    }}
  
  return(plot)
})
