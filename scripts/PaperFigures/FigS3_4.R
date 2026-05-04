if(!require(purrr)){
  install.packages(purrr)
  library(purrr) 
}
if(!require(ggridges)){
  install.packages(ggridges)
  library(ggridges) 
}
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(latex2exp)
library(patchwork)

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

x<-VAFs_passenger$VAF[VAFs_passenger$acquisition=="Expansive"]

tbl_pre<-VAFs_growth%>%
  filter(acquisition=="Expansive",
         VAF>0)
ks_pvalue_pass_growth_pre<-sapply(X=split(tbl_pre$VAF,tbl_pre$sel_adv),
                                  FUN=function(y){
                                    if(length(x[x>0])>1&length(y[y>0])>1){
                                      return(ks.test(x[x>0],y[y>0])$p.value)
                                    }else{return(NA)}
                                  })%>%
  as_tibble(rownames="sel_adv")%>%
  rename("pvalue"="value")%>%
  mutate(
    significance=ifelse(is.na(pvalue)|(pvalue>0.01),
                        "",
                        ifelse(pvalue>0.001,
                               "*",
                               ifelse(pvalue>0.0001,
                                      "**","***"
                               ))),
    sel_adv=factor(sel_adv,levels = sort(as.numeric(sel_adv)))
    
  )


tbl_pre<-VAFs_space%>%
  filter(acquisition=="Expansive",
         VAF>0)
ks_pvalue_pass_space_pre<-sapply(X=split(tbl_pre$VAF,tbl_pre$add_space),
                                 FUN=function(y){
                                   if(length(x[x>0])>1&length(y[y>0])>1){
                                     return(ks.test(x[x>0],y[y>0])$p.value)
                                   }else{return(NA)}
                                 })%>%
  as_tibble(rownames="add_space")%>%
  rename("pvalue"="value")%>%
  mutate(
    significance=ifelse(is.na(pvalue)|(pvalue>0.01),
                        "",
                        ifelse(pvalue>0.001,
                               "*",
                               ifelse(pvalue>0.0001,
                                      "**","***"
                               )))
  )
ks_pvalue_pass_space_pre$add_space=factor(ks_pvalue_pass_space_pre$add_space,levels = sort(as.integer(ks_pvalue_pass_space_pre$add_space)))

x<-VAFs_passenger$VAF[VAFs_passenger$acquisition=="Stable"]

tbl_post<-VAFs_growth%>%
  filter(acquisition=="Stable",
         VAF>0)
ks_pvalue_pass_growth_post<-sapply(X=split(tbl_post$VAF,tbl_post$sel_adv),
                                   FUN=function(y){
                                     if(length(x[x>0])>1&length(y[y>0])>1){
                                       return(ks.test(x[x>0],y[y>0])$p.value)
                                     }else{return(NA)}
                                   })%>%
  as_tibble(rownames="sel_adv")%>%
  rename("pvalue"="value")%>%
  mutate(
    significance=ifelse(is.na(pvalue)|(pvalue>0.01),
                        "",
                        ifelse(pvalue>0.001,
                               "*",
                               ifelse(pvalue>0.0001,
                                      "**","***"
                               )))
  )

tbl_post<-VAFs_space%>%
  filter(acquisition=="Stable",
         VAF>0)
ks_pvalue_pass_space_post<-sapply(X=split(tbl_post$VAF,tbl_post$add_space),
                                  FUN=function(y){
                                    if(length(x[x>0])>1&length(y[y>0])>1){
                                      return(ks.test(x[x>0],y[y>0])$p.value)
                                    }else{return(NA)}
                                  })%>%
  as_tibble(rownames="add_space")%>%
  rename("pvalue"="value")%>%
  mutate(
    significance=ifelse(is.na(pvalue)|(pvalue>0.01),
                        "",
                        ifelse(pvalue>0.001,
                               "*",
                               ifelse(pvalue>0.0001,
                                      "**","***"
                               )))
  )
ks_pvalue_pass_space_post$add_space=factor(ks_pvalue_pass_space_post$add_space,levels = sort(as.integer(ks_pvalue_pass_space_post$add_space)))


df<-VAFs_growth%>%
  filter(VAF>0,acquisition=="Expansive")%>%
  select(-acquisition)%>%
  bind_rows(VAFs_passenger%>%
              filter(VAF>0,acquisition=="Expansive")%>%
              pull(VAF)%>%
              expand.grid(VAFs_growth$sel_adv)%>%
              rename("VAF"=Var1,"sel_adv"=Var2),.id="FE")%>%
  mutate(FE=factor(FE,labels = c("Proliferative","Passenger")))

palette<-c(rbind(colorRampPalette(c("#C89FA3","#89696E"))(length(unique(VAFs_growth$sel_adv))),
                 rep("#D0E1D4",length(unique(VAFs_growth$sel_adv)))))

plot_growth_expansive<-ggplot(df,
             aes(x=VAF,group=FE,fill=interaction(FE,sel_adv),color=interaction(FE,sel_adv)))+
  geom_rect(data=ks_pvalue_pass_growth_pre,
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
  geom_density(alpha=0.7)+
  scale_fill_manual(values = palette)+
  scale_color_manual(values = palette)+
  scale_x_continuous(transform = "log2",
                     limits = c(0.0000001,0.5),
                     breaks = c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.5),
                     labels = c("1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.1","0.5"))+
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=5,angle=90,hjust = 1,vjust = 0.5)
  )+
  facet_grid(.~sel_adv,scales="free_x")


df<-VAFs_space%>%
  filter(VAF>0,acquisition=="Expansive")%>%
  select(-acquisition)%>%
  bind_rows(VAFs_passenger%>%
              filter(VAF>0,acquisition=="Expansive")%>%
              pull(VAF)%>%
              expand.grid(VAFs_space$add_space)%>%
              rename("VAF"=Var1,"add_space"=Var2),.id="FE")%>%
  mutate(FE=factor(FE,labels = c("Proliferative","Passenger")))

palette<-c(rbind(colorRampPalette(c("#C89FA3","#89696E"))(length(unique(VAFs_space$add_space))),
                 rep("#D0E1D4",length(unique(VAFs_space$add_space)))))

plot_space_expansive<-ggplot(df,
             aes(x=VAF,group=FE,fill=interaction(FE,add_space),color=interaction(FE,add_space)))+
  geom_rect(data=ks_pvalue_pass_space_pre,
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
  geom_density(alpha=0.7)+
  scale_fill_manual(values = palette)+
  scale_color_manual(values = palette)+
  scale_x_continuous(transform = "log2",
                     limits = c(0.0000001,0.5),
                     breaks = c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.5),
                     labels = c("1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.1","0.5"))+
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=5,angle=90,hjust = 1,vjust = 0.5)
  )+
  facet_grid(.~add_space,scales="free_x")



# POST

df<-VAFs_growth%>%
  filter(VAF>0,acquisition=="Stable")%>%
  select(-acquisition)%>%
  bind_rows(VAFs_passenger%>%
              filter(VAF>0,acquisition=="Stable")%>%
              pull(VAF)%>%
              expand.grid(VAFs_growth$sel_adv)%>%
              rename("VAF"=Var1,"sel_adv"=Var2),.id="FE")%>%
  mutate(FE=factor(FE,labels = c("Proliferative","Passenger")))

palette<-c(rbind(colorRampPalette(c("#C89FA3","#89696E"))(length(unique(VAFs_growth$sel_adv))),
                 rep("#D0E1D4",length(unique(VAFs_growth$sel_adv)))))
Growth_labels <- c(
  "0.01" = "\\frac{1}{10} WT",
  "0.05" = "\\frac{1}{2} WT",
  "0.1"  = "WT",
  "0.2"  = "2 WT",
  "1"    = "10 WT"
)
growth_plotmath <- sapply(Growth_labels, function(x) latex2exp::TeX(x, output = "character"))

plot_growth_stable<-ggplot(df,
             aes(x=VAF,group=FE,fill=interaction(FE,sel_adv),color=interaction(FE,sel_adv)))+
  geom_rect(data=ks_pvalue_pass_growth_post,
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
  geom_density(alpha=0.7)+
  scale_fill_manual(values = palette)+
  scale_color_manual(values = palette)+
  scale_x_continuous(transform = "log2",
                     limits = c(0.0000001,0.5),
                     breaks = c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.5),
                     labels = c("1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.1","0.5"))+
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank()
  )+
  facet_grid(.~sel_adv,scales="free_x",
             labeller = labeller(
               sel_adv = as_labeller(growth_plotmath,default = label_parsed)
             ))


df<-VAFs_space%>%
  filter(VAF>0,acquisition=="Stable")%>%
  select(-acquisition)%>%
  bind_rows(VAFs_passenger%>%
              filter(VAF>0,acquisition=="Stable")%>%
              pull(VAF)%>%
              expand.grid(VAFs_space$add_space)%>%
              rename("VAF"=Var1,"add_space"=Var2),.id="FE")%>%
  mutate(FE=factor(FE,labels = c("Proliferative","Passenger")))

palette<-c(rbind(colorRampPalette(c("#C89FA3","#89696E"))(length(unique(VAFs_space$add_space))),
                 rep("#D0E1D4",length(unique(VAFs_space$add_space)))))

space_labels <- c(
  "100000" = "\\frac{1}{10} K",
  "500000" = "\\frac{1}{2} K",
  "1000000"  = "K",
  "2000000"  = "2 K",
  "10000000"    = "10 K"
)
space_plotmath <- sapply(space_labels, function(x) latex2exp::TeX(x, output = "character"))

plot_space_stable<-ggplot(df,
             aes(x=VAF,group=FE,fill=interaction(FE,add_space),color=interaction(FE,add_space)))+
  geom_rect(data=ks_pvalue_pass_space_post,
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
  geom_density(alpha=0.7)+
  scale_fill_manual(values = palette)+
  scale_color_manual(values = palette)+
  scale_x_continuous(transform = "log2",
                     limits = c(0.0000001,0.5),
                     breaks = c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.5),
                     labels = c("1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.1","0.5"))+
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank()
  )+
  facet_grid(.~add_space,scales="free_x",
             labeller = labeller(
               add_space = as_labeller(space_plotmath,default = label_parsed)
             ))


(plot_growth_stable/plot_growth_expansive)+
  plot_annotation(title="Proliferative Deregulation",
                  theme = theme(plot.title = element_text(face="bold")))

(plot_space_stable/plot_space_expansive)+
  plot_annotation(title="Limit Evasion",
                  theme = theme(plot.title = element_text(face="bold")))
