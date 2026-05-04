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
  mutate(sus=factor(sus,levels = sort(unique(comp_synth_vcf_sampled$sus))),
         off=factor(off,levels = sort(unique(comp_synth_vcf_sampled$off))))%>%
  group_by(sus,off)%>%
  dplyr::mutate(n_mut=dplyr::n())%>%
  ggplot()+
  geom_density_ridges_gradient(aes(x=VAF,
                                   y = off,
                                   fill=after_stat(x)
  ))+
  scale_fill_gradient2(low = "#364156",mid = "#7D4E57",high = "#D66853",
                       midpoint = 0.5,limits=c(0,1))+
  scale_x_continuous(limits = c(0,1),
                     breaks = c(0,0.25,0.5,0.75,1),name = "Suceptibility")+
  ylab("Offensiveness")+
  facet_grid(.~sus,switch="x")+
  guides(fill = "none", 
         alpha = guide_legend("Number of data to compute the density"))+
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
