if(!require(purrr)){
  install.packages(purrr)
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
    summarise(count_surv_prek=sum(VAF>0,na.rm = TRUE),
              n_sim_prek=length(VAF),
              prob_surv_prek=count_surv_prek/n_sim_prek)%>%
    merge(
      tibble_post%>%
        group_by(across(c(-VAF)))%>%
        summarise(count_surv_postk=sum(VAF>0,na.rm = TRUE),
                  n_sim_postk=length(VAF),
                  prob_surv_postk=count_surv_postk/n_sim_postk)
    )%>%
    group_by(across(c(-count_surv_prek,count_surv_postk,n_sim_postk,n_sim_prek,prob_surv_postk,prob_surv_prek)))%>%
    mutate(prop_test=(prop.test(x=c(count_surv_prek,count_surv_postk),n=c(n_sim_prek,n_sim_postk)))$p.value)%>%
    ungroup()%>%
    dplyr::select(-c(n_sim_prek,n_sim_postk))%>%
    pivot_longer(c(prob_surv_prek,prob_surv_postk),
                 names_to = "acquisition_prob",values_to = "prob_surv")%>%
    pivot_longer(c(count_surv_prek,count_surv_postk),
                 names_to = "acquisition_count",values_to = "count_surv")%>%
    mutate(acquisition_prob=str_remove(acquisition_prob,"prob_surv_"),
           acquisition_count=str_remove(acquisition_count,"count_surv_"))%>%
    filter(acquisition_prob==acquisition_count)%>%
    select(-acquisition_count)%>%
    rename("acquisition"="acquisition_prob")%>%
    distinct()
    
},VAFs_pre,VAFs_post)


max_prob<-max(sapply(prob_surv,function(tibble_prob){max(tibble_prob$prob_surv)}))

# PASSENGER
passenger_prob_surv<-prob_surv[[1]]

tiles <- Insite:::make_split_tiles(passenger_prob_surv%>%mutate(par=1),
                          x_col = par,
                          value_col = prob_surv,
                          half_col = acquisition,
                          half_levels = c("prek","postk"),
                          direction = "bl_tr")

plot<-ggplot() +
  geom_polygon(data=tiles, aes(x = px,
                               y = py,
                               group = interaction(tile_idx, acquisition),
                               fill = prob_surv),
               color = "white", linewidth = 0.3) +
  scale_x_continuous(breaks = 1:length(unique(tiles$x_label)),
                     labels = unique(tiles$x_label)
  ) +
  coord_equal() +
  scale_fill_gradient(limits=c(0.001,0.5),transform="log10",
                      low = "#DDE4D9" ,
                      high = "#5D4E60")+
  theme_void()+
  theme(
    legend.position = "none",
  )
plot

plot<-ggplot() +
  geom_polygon(data=tiles, aes(x = px,
                               y = py,
                               group = interaction(tile_idx, acquisition),
                               fill = prob_surv),
               color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(is.na(unique(tiles$prop_test))|(unique(tiles$prop_test)>0.01),
                       "",
                       ifelse(unique(tiles$prop_test)>0.001,
                              "*",
                              ifelse(unique(tiles$prop_test)>0.0001,
                                     "**","***"
                              ))),
        x=min(tiles$px)+0.1,
        y=max(tiles$py)-0.25),hjust = 0,
    color="white",size = 50
  ) +
  scale_x_continuous(breaks = 1:length(unique(tiles$x_label)),
                     labels = unique(tiles$x_label)
  ) +
  coord_equal() +
  scale_fill_gradient(limits=c(0.001,0.5),transform="log10",
                      low = "#DDE4D9" ,
                      high = "#5D4E60")+
  theme_void()+
  theme(
    legend.position = "none",
  )
plot

count_surv_passenger<-prob_surv[[1]]$count_surv


# PROLIFERATION

growth_prob_surv<-prob_surv[[2]]

tiles <- Insite:::make_split_tiles(growth_prob_surv,
                          x_col = sel_adv,
                          value_col = prob_surv,
                          half_col = acquisition,
                          half_levels = c("prek","postk"),
                          direction = "bl_tr")

plot<-ggplot() +
  geom_polygon(data=tiles, aes(x = px,
                               y = py,
                               group = interaction(tile_idx, acquisition),
                               fill = prob_surv),
               color = "white", linewidth = 0.3) +
  scale_x_continuous(breaks = 1:length(unique(tiles$x_label)),
                     labels = c(TeX("$\\frac{1}{10}$ WT"),
                                TeX("$\\frac{1}{2}$ WT"),
                                "WT",
                                TeX("$2$ WT"),
                                TeX("$10$ WT"))) +
  geom_text(
    data=tiles%>%
      group_by(sel_adv,prop_test)%>%
      summarise(x=min(px)+0.1,y=max(py)-0.25)%>%
      mutate(significance=ifelse(is.na(prop_test)|(prop_test>0.01),
                                 "",
                                 ifelse(prop_test>0.001,
                                        "*",
                                        ifelse(prop_test>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 0,size = 11,
    color="white"
  ) +
  coord_equal() +
  scale_fill_gradient(limits=c(0.001,0.5),
                      low = "#DDE4D9" ,
                      transform="log10",
                      high = "#5D4E60",
                      name = "Survival probability"
                      
  )+
  theme_void()+
  theme(
    legend.ticks = element_blank(),
    #legend.position = "none",
    axis.text.x = element_text()
  )
plot


# COMPETITION 

comp_prob_surv<-prob_surv[[3]]

tiles <- Insite:::make_split_tiles(comp_prob_surv,
                          x_col = off,
                          y_col = sus,
                          value_col = prob_surv,
                          half_col = acquisition,
                          half_levels = c("prek","postk"),
                          direction = "bl_tr")

plot<-ggplot() +
  geom_polygon(data=tiles, aes(x = px,
                               y = py,
                               group = interaction(tile_idx, acquisition),
                               fill = prob_surv),color = "white", linewidth = 0.3) +
  geom_text(
    data=tiles%>%
      group_by(sus,off,prop_test)%>%
      summarise(x=min(px)+0.1,y=max(py)-0.25)%>%
      mutate(significance=ifelse(is.na(prop_test)|(prop_test>0.01),
                                 "",
                                 ifelse(prop_test>0.001,
                                        "*",
                                        ifelse(prop_test>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 0,
    color="white",size = 6
  ) +
  scale_x_continuous(breaks = 1:length(unique(tiles$x_label)),
                     labels = unique(tiles$x_label)[unique(tiles$x_pos)]) +
  scale_y_continuous(breaks = 1:length(unique(tiles$y_label)),
                     labels = unique(tiles$y_label)[unique(tiles$y_pos)]) +
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




# SPACE

space_prob_surv<-prob_surv[[4]]

tiles <- Insite:::make_split_tiles(space_prob_surv,
                          x_col = add_space,
                          value_col = prob_surv,
                          half_col = acquisition,
                          half_levels = c("prek","postk"),
                          direction = "bl_tr")

plot<-ggplot() +
  geom_polygon(data=tiles, aes(x = px,
                               y = py,
                               group = interaction(tile_idx, acquisition),
                               fill = prob_surv),color = "white", linewidth = 0.3) +
  geom_text(
    data=tiles%>%
      group_by(add_space,prop_test)%>%
      summarise(x=min(px)+0.1,y=max(py)-0.25)%>%
      mutate(significance=ifelse(is.na(prop_test)|(prop_test>0.01),
                                 "",
                                 ifelse(prop_test>0.001,
                                        "*",
                                        ifelse(prop_test>0.0001,
                                               "**","***"
                                        )))),
    aes(label = significance,
        x=x,
        y=y),hjust = 0,
    color="white",size = 11
  ) +
  scale_x_continuous(breaks = 1:length(unique(tiles$x_label)),
                     labels = c(TeX("$\\frac{1}{10}$ K"),
                                TeX("$\\frac{1}{2}$ K"),
                                "K",
                                TeX("$2$ K"),
                                TeX("$10$ K"))) +
  coord_equal() +
  scale_fill_gradient(limits=c(0.001,0.5),
                      transform="log10",
                      low = "#DDE4D9" ,
                      high = "#5D4E60")+
  theme_void()+
  theme(
    #legend.position = "none",
    axis.text.x = element_text()
  )
plot
