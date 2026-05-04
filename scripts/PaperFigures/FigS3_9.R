library(ggplot2)
library(dplyr)
library(ggridges)
if(!require(MoMAColors)){
  devtools::install_github("BlakeRMills/MoMAColors")
  library(MoMAColors)
}
library(scales)

load("Data/Experiment2/comp_vcf_synth.RData")

regression_lines<-tibble()
for(sus_i in unique(comp_synth_vcf_sampled$sus)){
  data<-comp_synth_vcf_sampled%>%
    group_by(sus,off)%>%
    dplyr::summarise(n_mut=dplyr::n())%>%
    filter(sus==sus_i,n_mut>5)
  
  if(sum(!is.na(data$n_mut))>5){
    lin_mod<-lm(formula = n_mut~off,data=data)
    regression_lines<-bind_rows(regression_lines,
                                tibble(sus=sus_i,
                                       off=sort(data$off[!is.na(data$n_mut)]),
                                       n_mut_fit=predict(lin_mod),
                                       adj_Rsq=summary(lin_mod)$adj.r.squared,
                                       Rsq=summary(lin_mod)$r.squared,
                                       intercept=summary(lin_mod)$coefficients["off",1]
                                ))
  }
}

plot<-comp_synth_vcf_sampled%>%
  group_by(sus,off)%>%
  dplyr::summarise(n_mut=dplyr::n())%>%
  filter(n_mut>5,sus<0.9)%>%
  ggplot()+
  geom_point(aes(x=off,y=n_mut,color=as.factor(sus)))+
  geom_line(data=regression_lines%>%
              filter(sus<0.9),
            aes(x=off,y=n_mut_fit,color=as.factor(sus)))+
  geom_text(data=regression_lines%>%
              filter(sus<0.9)%>%
              select(sus,Rsq)%>%
              distinct(),
            aes(x=0,
                y=0.32,
                label=paste("Rsquared ",round(Rsq,2)),
                color=as.factor(sus)),
            size=2)+
  ylab("Number of mutations detected")+
  xlab("Suceptibility")+
  guides(color = "none", 
         alpha = guide_legend("Number of data to compute the sd"))+
  scale_color_moma_d("Levine2")+
  facet_grid(.~sus,
             switch="both",
             labeller = labeller(
               sus = function(x) {paste("Offensiveness \n \n", x)}
             ))+
  theme_minimal()+
  theme(
    panel.grid.major = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.minor.y = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.minor.x = element_blank(),
    strip.placement = "outside",
    axis.title = element_text(size=10,),
    legend.position = "bottom",
    axis.text.x = element_text(size=8)
  )
plot

regression_lines<-tibble()
for(off_i in unique(comp_synth_vcf_sampled$off)){
  data<-comp_synth_vcf_sampled%>%
    group_by(sus,off)%>%
    dplyr::summarise(n_mut=dplyr::n())%>%
    filter(off==off_i,n_mut>5,sus<0.9)
  
  if(sum(!is.na(data$n_mut))>5){
    lin_mod<-lm(formula = n_mut~sus,data=data)
    regression_lines<-bind_rows(regression_lines,
                                tibble(off=off_i,
                                       sus=sort(data$sus[!is.na(data$n_mut)]),
                                       n_mut_fit=predict(lin_mod),
                                       adj_Rsq=summary(lin_mod)$adj.r.squared,
                                       Rsq=summary(lin_mod)$r.squared,
                                       
                                       intercept=summary(lin_mod)$coefficients["sus",1]
                                ))
  }
}


plot<-comp_synth_vcf_sampled%>%
  group_by(sus,off)%>%
  dplyr::summarise(n_mut=dplyr::n())%>%
  filter(n_mut>5,sus<0.9)%>%
  ggplot()+
  geom_point(aes(x=sus,y=n_mut,color=as.factor(off)))+
  geom_line(data=regression_lines%>%
              filter(sus<0.9),
            aes(x=sus,y=n_mut_fit,color=as.factor(off)))+
  geom_text(data=regression_lines%>%
              filter(sus<0.9)%>%
              select(off,Rsq)%>%distinct(),
            aes(x=-1,y=0.32,label=paste("Rsquared ",round(Rsq,2)),color=as.factor(off)),
            size=1.5)+
  ylab("Number of mutations detected")+
  xlab("Offensiveness")+
  guides(color = "none", 
         alpha = guide_legend("Number of data to compute the sd"))+
  scale_color_moma_d("Levine1")+
  facet_grid(.~off,
             switch="both",
             labeller = labeller(
               off = function(x) {paste("Susceptibility \n \n", x)}
             ))+
  theme_minimal()+
  theme(
    panel.grid.major = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.minor.y = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.minor.x = element_blank(),
    strip.placement = "outside",
    axis.title = element_text(size=10),
    legend.position = "bottom",
    axis.text.x = element_text(size=4)
  )
plot
