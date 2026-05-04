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

tbl_pre<-VAFs_pre[[3]]%>%
  mutate(par=paste(sus,off,sep="_"))%>%
  select(par,VAF)
tbl_post<-VAFs_post[[3]]%>%
  mutate(par=paste(sus,off,sep="_"))%>%
  select(par,VAF)
x<-split(tbl_pre$VAF,tbl_pre$par)
y<-split(tbl_post$VAF,tbl_post$par)
ks_pvalue_prepost_comp<-as_tibble(
  mapply(function(x,y){
    if(length(x[x>0])>1&length(y[y>0])>1){
      return(ks.test(x[x>0],y[y>0])$p.value)
    }else{return(NA)}
  },
  x,y),
  rownames="par")%>%
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

min_VAF<-min(VAFs_passenger$VAF[VAFs_passenger$VAF>0],
             VAFs_growth$VAF[VAFs_growth$VAF>0],
             VAFs_comp$VAF[VAFs_comp$VAF>0],
             VAFs_space$VAF[VAFs_space$VAF>0])
max_VAF<-max(VAFs_passenger$VAF[VAFs_passenger$VAF>0],
             VAFs_growth$VAF[VAFs_growth$VAF>0],
             VAFs_comp$VAF[VAFs_comp$VAF>0],
             VAFs_space$VAF[VAFs_space$VAF>0])

plot<-ggplot(data=VAFs_comp%>%
               filter(VAF>0),
             aes(x=VAF,y=acquisition,fill=after_stat(x)))+
  geom_rect(data=ks_pvalue_prepost_comp,
            aes(xmin =0.0000001,
                xmax = 0.5,
                ymin = -Inf,
                ymax=Inf,linetype = significance),
            fill="transparent",
            color="#CBCA83",
            inherit.aes = FALSE)+
  coord_cartesian(clip = "off")+
  scale_linetype_manual(values = c(0,3,5,1),breaks = c("","*","**","***"))+
  geom_density_ridges_gradient()+
  theme_minimal()+
  scale_fill_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                       midpoint = 0.25,limits=c(min_VAF,max_VAF))+
  scale_y_discrete(position = "right")+
  scale_x_continuous(limits = c(0,0.5),
                     breaks = c(0,0.1,0.2,0.3,0.4,0.5),
                     #labels = c("1e-6","1e-4","0.01","0.1","0.5")
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    legend.position = "none",
    axis.text.x = element_text(size=5,angle=90,hjust = 1,vjust = 0.5)
  )+
  facet_grid(sus~off,switch="both")+
  xlab("Offensiveness")+
  ylab("Susceptibility")

plot
