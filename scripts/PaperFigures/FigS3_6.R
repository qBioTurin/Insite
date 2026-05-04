if(!require(purrr)){
  install.packages("purrr")
  library(purrr) 
}
if(!require(ggridges)){
  install.packages("ggridges")
  library(ggridges) 
}
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


VAFs_passenger<-bind_rows(VAFs_pre[[1]],VAFs_post[[1]],.id = "acquisition")
VAFs_passenger$acquisition<-factor(VAFs_passenger$acquisition,labels = c("Expansive","Stable"))
VAFs_growth<-bind_rows(VAFs_pre[[2]],VAFs_post[[2]],.id = "acquisition")
VAFs_growth$acquisition<-factor(VAFs_growth$acquisition,labels = c("Expansive","Stable"))
VAFs_growth$sel_adv<-factor(VAFs_growth$sel_adv,
                            levels=sort(unique(VAFs_growth$sel_adv)))
VAFs_comp<-bind_rows(VAFs_pre[[3]],VAFs_post[[3]],.id = "acquisition")
VAFs_comp$sus<-factor(VAFs_comp$sus,
                      levels=sort(unique(VAFs_comp$sus),decreasing = TRUE))
VAFs_comp$off<-factor(VAFs_comp$off,
                      levels=sort(unique(VAFs_comp$off)))
VAFs_comp$acquisition<-factor(VAFs_comp$acquisition,labels = c("Expansive","Stable"))
VAFs_space<-bind_rows(VAFs_pre[[4]],VAFs_post[[4]],.id = "acquisition")
VAFs_space$acquisition<-factor(VAFs_space$acquisition,labels=c("Expansive","Stable"))
VAFs_space$add_space<-factor(VAFs_space$add_space,
                             levels=sort(unique(VAFs_space$add_space)))
min_VAF<-min(VAFs_passenger$VAF[VAFs_passenger$VAF>0],
             VAFs_growth$VAF[VAFs_growth$VAF>0],
             VAFs_comp$VAF[VAFs_comp$VAF>0],
             VAFs_space$VAF[VAFs_space$VAF>0])
max_VAF<-max(VAFs_passenger$VAF[VAFs_passenger$VAF>0],
             VAFs_growth$VAF[VAFs_growth$VAF>0],
             VAFs_comp$VAF[VAFs_comp$VAF>0],
             VAFs_space$VAF[VAFs_space$VAF>0])

x<-VAFs_passenger$VAF[VAFs_passenger$acquisition=="Stable"]

tbl_post<-VAFs_comp%>%
  filter(acquisition=="Stable",
         VAF>0)%>%
  mutate(par=paste(sus,off,sep="_"))
ks_pvalue_pass_comp_post<-sapply(X=split(tbl_post$VAF,tbl_post$par),
                                 FUN=function(y){
                                   if(length(x[x>0])>1&length(y[y>0])>1){
                                     return(ks.test(x[x>0],y[y>0])$p.value)
                                   }else{return(NA)}
                                 })%>%
  as_tibble(rownames="par")%>%
  separate_wider_delim(delim = "_",cols = par,names = c("sus","off"))%>%
  rename("pvalue"="value")%>%
  mutate(
    significance=ifelse(is.na(pvalue)|(pvalue>0.01),
                        "",
                        ifelse(pvalue>0.001,
                               "*",
                               ifelse(pvalue>0.0001,
                                      "**","***"
                               ))),
    sus=factor(sus,levels =c(2,1.5,1.1,1,0.9,0.5,0.2,0.1,0.01,0,-0.01,-0.1,-0.5,-1,-2)),
    off=factor(off,levels =rev(c(2,1.5,1.1,1,0.9,0.5,0.2,0.1,0.01,0,-0.01,-0.1,-0.5,-1,-2))))

x<-VAFs_passenger$VAF[VAFs_passenger$acquisition=="Expansive"]

tbl_pre<-VAFs_comp%>%
  filter(acquisition=="Expansive",
         VAF>0)%>%
  mutate(par=paste(sus,off,sep="_"))
ks_pvalue_pass_comp_pre<-sapply(X=split(tbl_pre$VAF,tbl_pre$par),
                                FUN=function(y){
                                  if(length(x[x>0])>1&length(y[y>0])>1){
                                    return(ks.test(x[x>0],y[y>0])$p.value)
                                  }else{return(NA)}
                                })%>%
  as_tibble(rownames="par")%>%
  separate_wider_delim(delim = "_",cols = par,names = c("sus","off"))%>%
  rename("pvalue"="value")%>%
  mutate(
    significance=ifelse(is.na(pvalue)|(pvalue>0.01),
                        "",
                        ifelse(pvalue>0.001,
                               "*",
                               ifelse(pvalue>0.0001,
                                      "**","***"
                               ))),
    sus=factor(sus,levels =c(2,1.5,1.1,1,0.9,0.5,0.2,0.1,0.01,0,-0.01,-0.1,-0.5,-1,-2)),
    off=factor(off,levels =rev(c(2,1.5,1.1,1,0.9,0.5,0.2,0.1,0.01,0,-0.01,-0.1,-0.5,-1,-2))))

df<-VAFs_comp%>%
  filter(VAF>0,acquisition=="Expansive")%>%
  select(-acquisition)%>%
  bind_rows(VAFs_passenger%>%
              filter(VAF>0,acquisition=="Expansive")%>%
              pull(VAF)%>%
              expand_grid(VAFs_comp%>%select(sus,off))%>%
              rename("VAF"="."),.id="FE")%>%
  mutate(FE=factor(FE,labels = c("Resource Control","Passenger")))

palette<-c(rbind(colorRampPalette(c("#C89FA3","#89696E"))(length(unique(VAFs_comp$sus))),
                 rep("#D0E1D4",length(unique(VAFs_comp$sus)))))

plot<-ggplot(df,
             aes(x=VAF,
                 group=FE,
                 fill=interaction(FE,sus),
                 color=interaction(FE,sus)
             ))+
  geom_rect(data=ks_pvalue_pass_comp_pre,
            aes(xmin =0.0000001,
                xmax = 0.5,
                ymin = -Inf,
                ymax=Inf,
                linetype = significance),linewidth = 1,
            fill="transparent",
            color="#6C534E",
            inherit.aes = FALSE)+
  coord_cartesian(clip = "off")+
  scale_linetype_manual(values = c(0,5,3,1),breaks = c("","*","**","***"))+
  geom_density(alpha=0.5)+
  scale_fill_manual(values = palette)+
  scale_color_manual(values = palette)+
  facet_grid(sus~off,scales="free",switch="both")+
  xlab("Offensiveness")+
  ylab("Susceptibility")+
  scale_x_continuous(transform = "log2",
                     limits = c(0.0000001,0.5),
                     breaks = c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.5),
                     labels = c("1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.1","0.5"))+
  theme_minimal()+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size=5,angle=90,hjust = 1,vjust = 0.5)
  )
plot

df<-VAFs_comp%>%
  filter(VAF>0,acquisition=="Stable")%>%
  select(-acquisition)%>%
  bind_rows(VAFs_passenger%>%
              filter(VAF>0,acquisition=="Stable")%>%
              pull(VAF)%>%
              expand_grid(VAFs_comp%>%select(sus,off))%>%
              rename("VAF"="."),.id="FE")%>%
  mutate(FE=factor(FE,labels = c("Resource Control","Passenger")))

palette<-c(rbind(rep("#D0E1D4",length(unique(VAFs_comp$sus))),
                 colorRampPalette(c("#C89FA3","#89696E"))(length(unique(VAFs_comp$sus)))
))

plot<-ggplot(df,
             aes(x=VAF,
                 group=FE,
                 fill=interaction(FE,sus),
                 color=interaction(FE,sus)
             ))+
  geom_rect(data=ks_pvalue_pass_comp_post,
            aes(xmin =0.0000001,
                xmax = 0.5,
                ymin = -Inf,
                ymax=Inf,
                linetype = significance),linewidth = 1,
            fill="transparent",
            color="#6C534E",
            inherit.aes = FALSE)+
  coord_cartesian(clip = "off")+
  scale_linetype_manual(values = c(0,5,3,1),breaks = c("","*","**","***"))+
  geom_density(alpha=0.5)+
  facet_grid(sus~off,scales="free",switch="both")+
  xlab("Offensiveness")+
  ylab("Susceptibility")+
  scale_fill_manual(values = palette)+
  scale_color_manual(values = palette)+
  scale_x_continuous(transform = "log2",
                     limits = c(0.0000001,0.5),
                     breaks = c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.5),
                     labels = c("1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.1","0.5"))+
  theme_minimal()+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size=5,angle=90,hjust = 1,vjust = 0.5)
  )
plot
