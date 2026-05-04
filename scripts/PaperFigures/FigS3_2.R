if(!require(purrr)){
  install.packages("purrr")
  library(purrr) 
}
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(latex2exp)

load("Data/Experiment1/ExpansivePhase/VAF_FunEff.RData")
VAFs_pre<-VAFs
load("Data/Experiment1/StablePhase/VAF_FunEff.RData")
VAFs_post<-VAFs
rm(VAFs)


VAFs_pre[2:4]<-mapply(FUN = function(tibble_pre,tibble_post){
  post_par<-tibble_post%>%select(-VAF)%>%distinct()
  return(tibble_pre%>%right_join(post_par))
},VAFs_pre[-1],VAFs_post[-1])
VAFs_post[2:4]<-mapply(FUN = function(tibble_pre,tibble_post){
  pre_par<-tibble_pre%>%select(-VAF)%>%distinct()
  return(tibble_post%>%right_join(pre_par))
},VAFs_pre[-1],VAFs_post[-1])

prob_surv<-mapply(function(tibble_pre,tibble_post){
  tibble_pre%>%
    group_by(across(c(-VAF)))%>%
    summarise(prek=sum(VAF>0,na.rm = TRUE),
              n_sim_prek=length(VAF),
              prob_surv_prek=prek/n_sim_prek)%>%
    merge(
      tibble_post%>%
        group_by(across(c(-VAF)))%>%
        summarise(postk=sum(VAF>0,na.rm = TRUE),
                  n_sim_postk=length(VAF),
                  prob_surv_postk=postk/n_sim_postk)
    )%>%
    group_by(across(c(-prek,postk,n_sim_postk,n_sim_prek,prob_surv_postk,prob_surv_prek)))%>%
    mutate(prop_test=(prop.test(x=c(prek,postk),n=c(n_sim_prek,n_sim_postk)))$p.value)%>%
    ungroup()%>%
    dplyr::select(-c(n_sim_prek,n_sim_postk))%>%
    pivot_longer(c(prob_surv_prek,prob_surv_postk),
                 names_to = "acquisition_prob",values_to = "prob_surv")%>%
    pivot_longer(c(prek,postk),
                 names_to = "acquisition_count",values_to = "count_surv")%>%
    mutate(acquisition_prob=str_remove(acquisition_prob,"prob_surv_"),
           acquisition_count=str_remove(acquisition_count,"count_surv_"))%>%
    filter(acquisition_prob==acquisition_count)%>%
    select(-acquisition_count)%>%
    rename("acquisition"="acquisition_prob")%>%
    distinct()
  
},VAFs_pre,VAFs_post)

count_surv_passenger<-prob_surv[[1]]$count_surv

# GROWTH vs passenger
count_surv_growth<-split(prob_surv[[2]]$count_surv,f = rep(1:5,each=2))
names(count_surv_growth)<-unique(prob_surv[[2]]$sel_adv)
lapply(count_surv_growth,function(counts){
  names(counts)<-c("preK","postK")
})
growth_prob_surv<-prob_surv[[2]]

p_values_passvsgrowth<-lapply(X = count_surv_growth,
                              FUN = function(c_pass){
                                sapply(FUN = function(x,n){prop.test(x,n)$p.value},
                                       X=list(c(count_surv_passenger[1],c_pass[1]),
                                              c(count_surv_passenger[2],c_pass[2])),
                                       n=c(1000,1000))
                              })

tiles_1<-Insite:::make_split_tiles(growth_prob_surv,
                                   x_col = sel_adv,
                                   value_col = prob_surv,
                                   half_col = acquisition,
                                   half_levels = c("prek","postk"),
                                   direction = "bl_tr")%>%
  full_join(as.data.frame(do.call(rbind,p_values_passvsgrowth))%>%
              rownames_to_column()%>%
              mutate(rowname=as.numeric(rowname))%>%
              rename(sel_adv=rowname,
                     prek=V1,
                     postk=V2)%>%
              pivot_longer(c(prek,postk),
                           names_to = "acquisition",
                           values_to = "p_value_cfrpass"))

plot<-ggplot() +
  geom_polygon(data=tiles_1, aes(x = px,
                                 y = py,
                                 group = interaction(tile_idx, acquisition),
                                 fill = prob_surv),color = "white", linewidth = 0.3) +
  geom_text(
    data=tiles_1%>%
      filter(acquisition=="prek")%>%
      group_by(sel_adv,p_value_cfrpass)%>%
      summarise(x=max(px)-0.1,y=min(py))%>%
      mutate(significance=ifelse(is.na(p_value_cfrpass)|(p_value_cfrpass>0.01),
                                 "",
                                 ifelse(p_value_cfrpass>0.001,
                                        "*",
                                        ifelse(p_value_cfrpass>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 1,
    color="white",size = 11
  ) +
  geom_text(
    data=tiles_1%>%
      filter(acquisition=="postk")%>%
      group_by(sel_adv,p_value_cfrpass)%>%
      summarise(x=min(px)+0.1,y=max(py)-0.25)%>%
      mutate(significance=ifelse(is.na(p_value_cfrpass)|(p_value_cfrpass>0.01),
                                 "",
                                 ifelse(p_value_cfrpass>0.001,
                                        "*",
                                        ifelse(p_value_cfrpass>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 0,
    color="white",size = 11
  ) +
  scale_x_continuous(breaks = 1:length(unique(tiles_1$x_label)),
                     labels = c(TeX("$\\frac{1}{10}$ WT"),
                                TeX("$\\frac{1}{2}$ WT"),
                                "WT",
                                TeX("$2$ WT"),
                                TeX("$10$ WT"))) +
  coord_equal() +
  scale_fill_gradient(limits=c(0.001,0.5),transform="log10",
                      low = "#DDE4D9" ,na.value = "#EAEEE7",
                      high = "#5D4E60")+
  theme_void()+
  theme(
    legend.position = "none",
    axis.text.x = element_text()
  )
plot

# COMPETITION vs passenger

count_surv_comp<-split(prob_surv[[3]]$count_surv,f = rep(1:222,each=2))
names(count_surv_comp)<-unique(paste(prob_surv[[3]]$sus,prob_surv[[3]]$off,sep="_"))
lapply(count_surv_comp,function(counts){
  names(counts)<-c("preK","postK")
})

comp_prob_surv<-prob_surv[[3]]

p_values_passvscomp<-lapply(X = count_surv_comp,
                            FUN = function(c_pass){
                              sapply(FUN = function(x,n){prop.test(x,n)$p.value},
                                     X=list(c(count_surv_passenger[1],c_pass[1]),
                                            c(count_surv_passenger[2],c_pass[2])),
                                     n=c(1000,1000))
                            })


tiles_1<-Insite:::make_split_tiles(comp_prob_surv,
                                   x_col = off,
                                   y_col = sus,
                                   value_col = prob_surv,
                                   half_col = acquisition,
                                   half_levels = c("prek","postk"),
                                   direction = "bl_tr")%>%
  full_join(as.data.frame(do.call(rbind,p_values_passvscomp))%>%
              rownames_to_column()%>%
              separate_wider_delim(cols =rowname ,delim="_",names = c("sus","off"))%>%
              mutate(sus=as.numeric(sus),
                     off=as.numeric(off))%>%
              rename(prek=V1,
                     postk=V2)%>%
              pivot_longer(c(prek,postk),names_to = "acquisition",values_to = "p_value_cfrpass"))

plot<-ggplot() +
  geom_polygon(data=tiles_1, aes(x = px,
                                 y = py,
                                 group = interaction(tile_idx, acquisition),
                                 fill = prob_surv),color = "white", linewidth = 0.3) +
  geom_text(
    data=tiles_1%>%
      filter(acquisition=="prek")%>%
      group_by(sus,off,p_value_cfrpass)%>%
      summarise(x=max(px)-0.1,y=min(py))%>%
      mutate(significance=ifelse(is.na(p_value_cfrpass)|(p_value_cfrpass>0.01),
                                 "",
                                 ifelse(p_value_cfrpass>0.001,
                                        "*",
                                        ifelse(p_value_cfrpass>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 1,
    color="white",size = 6
  ) +
  geom_text(
    data=tiles_1%>%
      filter(acquisition=="postk")%>%
      group_by(sus,off,p_value_cfrpass)%>%
      summarise(x=min(px)+0.1,y=max(py)-0.25)%>%
      mutate(significance=ifelse(is.na(p_value_cfrpass)|(p_value_cfrpass>0.01),
                                 "",
                                 ifelse(p_value_cfrpass>0.001,
                                        "*",
                                        ifelse(p_value_cfrpass>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 0,
    color="white",size = 6
  ) +
  scale_x_continuous(breaks = 1:length(unique(tiles_1$x_label)),
                     labels = unique(tiles_1$x_label)[unique(tiles_1$x_pos)]) +
  scale_y_continuous(breaks = 1:length(unique(tiles_1$y_label)),
                     labels = unique(tiles_1$y_label)[unique(tiles_1$y_pos)]) +
  coord_equal() +
  scale_fill_gradient(limits=c(0.001,0.5),transform="log10",
                      low = "#DDE4D9" ,na.value = "#EAEEE7",
                      high = "#5D4E60")+
  theme_void()+
  theme(
    legend.position = "none",
    axis.text = element_text()
  )
plot


# SPACE vs passenger

count_surv_space<-split(prob_surv[[4]]$count_surv,f = rep(1:5,each=2))
names(count_surv_space)<-unique(prob_surv[[4]]$add_space)
lapply(count_surv_space,function(counts){
  names(counts)<-c("preK","postK")
})

space_prob_surv<-prob_surv[[4]]

p_values_passvsspace<-lapply(X = count_surv_space,
                             FUN = function(c_pass){
                               sapply(FUN = function(x,n){prop.test(x,n)$p.value},
                                      X=list(c(count_surv_passenger[1],c_pass[1]),
                                             c(count_surv_passenger[2],c_pass[2])),
                                      n=c(1000,1000))
                             })

tiles_1<-Insite:::make_split_tiles(space_prob_surv,
                                   x_col = add_space,
                                   value_col = prob_surv,
                                   half_col = acquisition,
                                   half_levels = c("prek","postk"),
                                   direction = "bl_tr")%>%
  full_join(as.data.frame(do.call(rbind,p_values_passvsspace))%>%
              rownames_to_column()%>%
              mutate(rowname=as.numeric(rowname))%>%
              rename(add_space=rowname,
                     prek=V1,
                     postk=V2)%>%
              pivot_longer(c(prek,postk),
                           names_to = "acquisition",
                           values_to = "p_value_cfrpass"))

plot<-ggplot() +
  geom_polygon(data=tiles_1, aes(x = px,
                                 y = py,
                                 group = interaction(tile_idx, acquisition),
                                 fill = prob_surv),color = "white", linewidth = 0.3) +
  geom_text(
    data=tiles_1%>%
      filter(acquisition=="prek")%>%
      group_by(add_space,p_value_cfrpass)%>%
      summarise(x=max(px)-0.1,y=min(py))%>%
      mutate(significance=ifelse(is.na(p_value_cfrpass)|(p_value_cfrpass>0.01),
                                 "",
                                 ifelse(p_value_cfrpass>0.001,
                                        "*",
                                        ifelse(p_value_cfrpass>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 1,
    color="white",size = 11
  ) +
  geom_text(
    data=tiles_1%>%
      filter(acquisition=="postk")%>%
      group_by(add_space,p_value_cfrpass)%>%
      summarise(x=min(px)+0.1,y=max(py)-0.25)%>%
      mutate(significance=ifelse(is.na(p_value_cfrpass)|(p_value_cfrpass>0.01),
                                 "",
                                 ifelse(p_value_cfrpass>0.001,
                                        "*",
                                        ifelse(p_value_cfrpass>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 0,
    color="white",size = 11
  ) +
  scale_x_continuous(breaks = 1:length(unique(tiles_1$x_label)),
                     labels = c(TeX("$\\frac{1}{10}$ K"),
                                TeX("$\\frac{1}{2}$ K"),
                                "K",
                                TeX("$2$ K"),
                                TeX("$10$ K"))) +
  coord_equal() +
  scale_fill_gradient(limits=c(0.001,0.5),transform="log10",
                      low = "#DDE4D9" ,na.value = "#EAEEE7",
                      high = "#5D4E60")+
  theme_void()+
  theme(
    legend.position = "none",
    axis.text = element_text(),
    axis.text.y = element_blank()
  )
plot

