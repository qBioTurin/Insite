library(dplyr)
library(tidyr)
library(scales)
if(!require(patchwork)){
  install.packages("patchwork")
  library(patchwork)
}
library(stringr)
library(ggplot2)


time_table<-read.table("Data/nD_exp/time_report.txt",header = TRUE)%>%
  mutate(EndSize_FE=str_remove(par_set,"params_"))%>%
  separate_wider_delim(EndSize_FE,"_",names = c("EndSize","FE","param"),too_few = "align_start")%>%
  mutate(param=ifelse(FE=="growth",paste("sel_adv",param,sep="_"),
                      ifelse(FE=="space",paste("additional_space",param,sep="_"),
                             ifelse(FE=="competition",paste("sus",param,"off_1",sep="_"),"None"))))


time_table_means<-time_table%>%
  group_by(par_set,EndSize,FE,param)%>%
  summarise(time=mean(time))

time_table_means$param<-factor(time_table_means$param,levels = c("None",
                                                                 "sel_adv_0d1",
                                                                 "sel_adv_0d2",
                                                                 "sel_adv_1",
                                                                 "sus_0d7_off_1",
                                                                 "sus_0d5_off_1",
                                                                 "sus_0d2_off_1",
                                                                 "sus_0d1_off_1",
                                                                 "additional_space_1e+06",
                                                                 "additional_space_2e+06",
                                                                 "additional_space_1e+07"),
                               labels=c("Passenger",
                                        "+WT",
                                        "+2WT",
                                        "+10WT",
                                        "0.7",
                                        "0.5",
                                        "0.2",
                                        "0.1",
                                        "+K",
                                        "+2K",
                                        "+10K"))

colors_palette1<-c("+WT"="#E5D5C9",
                  "+2WT"="#DDA77B",
                  "+10WT"="#DE8339",
                  "0.7"="#C6A5A5",
                  "0.5"="#B28686",
                  "0.2"="#926767",
                  "0.1"="#714747",
                  "Passenger"="#A8BCDF",
                  "+K"="#89CAC2",
                  "+2K"="#1B998B",
                  "+10K"="#225C56"
)

p<-ggplot(time_table_means%>%
            mutate(EndSize=ifelse(EndSize=="EndSize0d1Carring",
                                  "M=0.1K",
                                  ifelse(EndSize=="EndSize2Carring",
                                         "M=2K",
                                         "M=10K"))))+
  geom_bar(aes(x = time,y=param,fill=param),position = "dodge",stat = "identity")+
  scale_fill_manual(values = colors_palette1)+
  scale_x_continuous(breaks = c(10/60,30/60,1,2,3),labels = c("10sec","30sec","1min","2min","3min"))+
  facet_grid(.~EndSize)+
  xlab("Average time required for each simulation run")+
  theme_void()+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size=7,hjust = 1),
        strip.text = element_text(face="bold"),
        legend.position = "none",
        axis.title.x = element_text(size = 10),
        plot.title = element_text(face="bold"),
        panel.background = element_rect(fill="#FAFAFA"),
        panel.grid.major.x = element_line(colour = "white",linewidth = 1))
p

