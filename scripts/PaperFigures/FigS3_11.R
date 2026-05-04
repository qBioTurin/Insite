library(ggplot2)
library(dplyr)
library(ggridges)
if(!require(MoMAColors)){
  devtools::install_github("BlakeRMills/MoMAColors")
  library(MoMAColors)
}
library(scales)

load("Data/Experiment2/passenger_vcf_synth_nmut.RData")

mean_and_sd<-passenger_synth_vcf_sampled%>%
  group_by(nmut)%>%
  summarise(mean=mean(VAF),
            sd=sd(VAF))%>%
  bind_rows(passenger_synth_vcf%>%
              group_by(nmut)%>%
              summarise(mean=mean(VAF),
                        sd=sd(VAF)),.id="VAF_sampling")


plot<-
  ggplot()+
  geom_violin(data=passenger_synth_vcf_sampled%>%
                bind_rows(passenger_synth_vcf,.id="VAF_sampling")%>%
                mutate(VAF_sampling=factor(VAF_sampling,labels = c("Sequence-alike procedure (sampled VAF)",
                                                                   "Half cell prevalence (theoretical VAF)"))),
              aes(x = VAF,y = VAF_sampling, fill=VAF_sampling, color=VAF_sampling),alpha=0.8)+
  facet_grid(nmut~.,scales="free",switch="both")+
  scale_x_continuous(transform = "log",
                     breaks = c(0,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1))+
  scale_fill_manual(values=c("#F5B841","#6A041D"),
                    labels=label_wrap(30))+
  scale_color_manual(values=c("#F5B841","#6A041D"),
                     labels=label_wrap(30))+
  ylab("Number of base pairs considered")+
  geom_text(data = mean_and_sd,
            aes(x=mean,
                y=VAF_sampling,
                label = paste("Mean ",signif(mean,2),"\n","sd ",signif(sd,2))),size=2,nudge_y = -1)+
  guides(fill=guide_legend(title="VAF computing method"),
         color=guide_legend(title="VAF computing method"))+
  theme_minimal()+
  theme(
    panel.grid.major.x = element_line(colour = "#F5F5F5",linewidth = 0.1),
    panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",color="black",linewidth = 0.2),
    axis.title = element_text(size=10),
    legend.position = "right",
    legend.key.spacing.y = unit(10,"pt"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=8),
    strip.placement = "outside",
    strip.text.y.left = element_text(size=8,angle = 0)
  )
plot
