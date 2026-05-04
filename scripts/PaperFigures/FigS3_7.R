library(ggplot2)
library(dplyr)
library(ggridges)
if(!require(MoMAColors)){
  devtools::install_github("BlakeRMills/MoMAColors")
  library(MoMAColors)
}
library(scales)
library(latex2exp)

load("Data/Experiment2/comp_vcf_synth.RData")
load("Data/Experiment2/growth_vcf_synth.RData")
load("Data/Experiment2/space_vcf_synth.RData")
load("Data/Experiment2/passenger_vcf_synth.RData")

plot<-passenger_synth_vcf_sampled%>%
  ggplot()+
  geom_dotplot(aes(x=VAF,fill=after_stat(x),color=after_stat(x)),
               dotsize=0.8,
               binwidth = 0.02)+
  xlim(0,1)+
  scale_fill_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                       midpoint = 0.5,limits=c(0,1))+
  scale_color_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                        midpoint = 0.5,limits=c(0,1))+
  theme_minimal()+
  ggtitle("Passenger")+
  theme(
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    title = element_text(face="bold",size=14),
    panel.background = element_rect(fill = "transparent",color="black",linewidth = 0.2),
    strip.placement = "outside",
    axis.title = element_text(size=10,face = "plain"),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=8)
  )
plot

Growth_labels <- c(
  "0.01" = "\\frac{1}{10} WT",
  "0.05" = "\\frac{1}{2} WT",
  "0.1"  = "WT",
  "0.2"  = "2 WT",
  "1"    = "10 WT"
)
growth_plotmath <- sapply(Growth_labels, function(x) latex2exp::TeX(x, output = "character"))
plot<-growth_synth_vcf_sampled%>%
  ggplot()+
  geom_dotplot(aes(x=VAF,
                   group=sel_adv,
                   fill=after_stat(x),
                   color=after_stat(x)),
               dotsize=0.8,
               binwidth = 0.02)+
  xlim(0,1)+
  ylim(0,1)+
  ylab("Proliferative advantage")+
  facet_grid(sel_adv~.,
             scales="free",
             switch="both",
             labeller = labeller(
               sel_adv = as_labeller(growth_plotmath,default = label_parsed)
             ))+
  scale_fill_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                       midpoint = 0.5,limits=c(0,1))+
  scale_color_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                        midpoint = 0.5,limits=c(0,1))+
  ggtitle("Proliferative Deregulation")+
  theme_minimal()+
  theme(
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    title = element_text(face="bold",size=14),
    panel.background = element_rect(fill = "transparent",color="black",linewidth = 0.2),
    axis.title = element_text(size=10,face = "plain"),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=8),
    strip.placement = "outside",
    strip.text.y.left = element_text(size=8,angle = 0)
  )
plot

plot<-comp_synth_vcf_sampled%>%
  mutate(sus=factor(sus,levels = sort(unique(comp_synth_vcf_sampled$sus))),
         off=factor(off,levels = sort(unique(comp_synth_vcf_sampled$off),decreasing = TRUE)))%>%
  group_by(sus,off)%>%
  dplyr::mutate(n_mut=dplyr::n())%>%
  ggplot()+
  geom_dotplot(aes(x=VAF,fill=after_stat(x),
                   color=after_stat(x)),
               dotsize=0.8,
               binwidth = 0.02)+
  scale_fill_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                       midpoint = 0.5,limits=c(0,1))+
  scale_color_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                        midpoint = 0.5,limits=c(0,1))+
  xlab("Susceptibility")+
  ylab("Offensiveness")+
  facet_grid(off~sus,switch="both")+
  ggtitle("Resource Control")+
  theme_minimal()+
  theme(
    title = element_text(face="bold",size=14),
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",color="black",linewidth = 0.2),
    axis.title = element_text(size=10,face = "plain"),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=8),
    strip.placement = "outside",
    strip.text.y.left = element_text(size=8,angle = 0)
  )
plot

space_labels <- c(
  "1e+05" = "\\frac{1}{10} K",
  "5e+05" = "\\frac{1}{2} K",
  "1e+06"  = "K",
  "2e+06"  = "2 K",
  "1e+07"    = "10 K"
)
space_plotmath <- sapply(space_labels, function(x) latex2exp::TeX(x, output = "character"))

plot<-space_synth_vcf_sampled%>%
  ggplot()+
  geom_dotplot(aes(x=VAF,group=add_space,
                   fill=after_stat(x),
                   color=after_stat(x)),
               dotsize=0.8,
               binwidth = 0.019)+
  ggtitle("Limit Evasion")+
  ylab("additional space")+
  facet_grid(add_space~.,
             scales="free",
             switch="both",
             labeller = labeller(
               add_space = as_labeller(space_plotmath,default = label_parsed)
             ))+
  xlim(0,1)+
  scale_fill_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                       midpoint = 0.5,limits=c(0,1))+
  scale_color_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                        midpoint = 0.5,limits=c(0,1))+
  theme_minimal()+
  theme(
    title = element_text(face="bold",size=14),
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",color="black",linewidth = 0.2),
    axis.title = element_text(size=10,face = "plain"),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=8),
    strip.placement = "outside",
    strip.text.y.left = element_text(size=8,angle = 0)
  )
plot
