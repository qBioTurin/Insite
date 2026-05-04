library(ggplot2)
library(dplyr)
library(ggridges)
if(!require(MoMAColors)){
  devtools::install_github("BlakeRMills/MoMAColors")
  library(MoMAColors)
}
library(scales)

load("Data/Experiment2/comp_vcf_synth.RData")

plot<-comp_synth_vcf_sampled%>%
  filter(sus<0.9)%>%
  ungroup()%>%
  mutate(thrsh=(VAF>0.1))%>%
  group_by(thrsh,sus,off)%>%
  dplyr::summarize(n_mut=dplyr::n())%>%
  ggplot(aes(x=sus,y=n_mut,col=thrsh))+
  scale_color_manual(values=c("#364156","#D66853"),breaks = c(FALSE,TRUE),labels=c("VAF<0.1","VAF>0.1"))+
  guides(color=guide_legend(title=""))+
  geom_point()+
  ylab("Number of mutations detected")+
  xlab("Offensivenes")+
  facet_grid(
    . ~ off,
    switch = "both",
    labeller = labeller(
      off = function(x){paste("Susceptibility \n \n", x)}
    )
  )+
  theme_minimal()+
  theme(panel.grid.major = element_line(colour = "#F5F5F5",linewidth = 0.1),
        panel.grid.minor.y = element_line(colour = "#F5F5F5",linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size=10),
        legend.position = "bottom",
        axis.text.x = element_text(size=4))
plot

plot<-comp_synth_vcf_sampled%>%
  filter(sus<0.9)%>%
  ungroup()%>%
  mutate(thrsh=(VAF>0.1))%>%
  group_by(thrsh,sus,off)%>%
  dplyr::summarize(n_mut=dplyr::n())%>%
  ggplot(aes(x=off,y=n_mut,col=thrsh))+
  scale_color_manual(values=c("#364156","#D66853"),breaks = c(FALSE,TRUE),labels=c("VAF<0.1","VAF>0.1"))+
  guides(color=guide_legend(title=""))+
  geom_point()+
  ylab("Number of mutations detected")+
  xlab("Susceptibility")+
  facet_grid(
    . ~ sus,
    switch = "both",
    labeller = labeller(
      sus = function(x){paste("Offensivenes \n \n", x)}
    )
  )+
  theme_minimal()+
  theme(panel.grid.major = element_line(colour = "#F5F5F5",linewidth = 0.1),
        panel.grid.minor.y = element_line(colour = "#F5F5F5",linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size=10),
        legend.position = "bottom",
        axis.text.x = element_text(size=4))
plot
