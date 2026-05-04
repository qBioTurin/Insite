library(Insite)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
if(!require(patchwork)){
  install.packages(patchwork)
  library(patchwork)
}

load("Data/nD_exp/nD_pts.RData")

nD_pts<-nD_pts%>%
  separate_wider_delim(cols = "par_set",delim = "_",
                       names = c("Intro","EndSize","FE","par_value"),
                       too_few = "align_start")%>%
  dplyr::select(-Intro)%>%
  mutate(EndSize=as.numeric(str_replace(string = str_remove(
    str_remove(EndSize,
               "EndSize"),
    "Carring"),
    pattern = "d",
    ".")),
    par_value=as.numeric(str_replace(par_value,pattern = "d",".")),
    color_par=ifelse(FE=="growth",
                     paste("sel_adv",par_value,sep="_"),
                     ifelse(FE=="competition",
                            paste("sus",par_value,"off_1",sep="_"),
                            ifelse(FE=="space",
                                   paste("additional_space",par_value,sep="_"),
                                   "None"
                            ))
    ),
    det_size_reached=ifelse(det_size_reached,"Stopping size reached","Stopping size NOT reached"))


colors_palette<-c("sel_adv_0.1"="#E5D5C9",
                  "sel_adv_0.2"="#DDA77B",
                  "sel_adv_1"="#DE8339",
                  "sus_0.7_off_1"="#C6A5A5",
                  "sus_0.5_off_1"="#B28686",
                  "sus_0.2_off_1"="#926767",
                  "sus_0.1_off_1"="#714747",
                  "None"="#A8BCDF",
                  "additional_space_1e+06"="#89CAC2",
                  "additional_space_2e+06"="#1B998B",
                  "additional_space_1e+07"="#225C56"
)

area_grey<-tibble(
  n=seq(1,2,by=0.01)
)%>%
  mutate(D=(2-n)^(-2))%>%
  filter(D<ceiling(max(nD_pts$D)))

mean_branch<-tibble(
  n=seq(1,ceiling(max(nD_pts$n)),by=0.01)
)%>%
  mutate(D=9*(2*n-1)/8)

lin_tree<-tibble(
  n=seq(1,ceiling(max(nD_pts$n)),by=0.01)
)%>%
  mutate(
    frac_n=n-floor(n),
    D=((1-frac_n)^2+frac_n^2)^(-1))

jitter_width<-sd(nD_pts$n)/50
jitter_height<-sd(nD_pts$D)/50

EndSize_labels <- c(
  "0.1" = "0.1 K",
  "2" = "2 K",
  "5" = "5 K",
  "10" = "10 K"
)

plot<-ggplot(nD_pts)+
  geom_ribbon(data=area_grey,aes(x=n,ymin=D,ymax = ceiling(max(nD_pts$D))),fill="#DFE2E0")+
  geom_line(data=lin_tree,aes(x=n,y=D),color="#DFE2E0")+
  geom_jitter(aes(x=n,y=D,fill=color_par,shape=FE),
              width = jitter_width,
              height = jitter_height,
              color="white",stroke = 0.2,size=4,alpha=0.8)+
  scale_shape_manual(values = 21:24)+
  scale_fill_manual(values = colors_palette)+
  scale_y_continuous(limits = c(0,ceiling(max(nD_pts$D))))+
  scale_x_continuous(limits = c(0.5,ceiling(max(nD_pts$n))),breaks = 1:ceiling(max(nD_pts$n)))+
  facet_wrap(~EndSize, labeller = labeller(EndSize = EndSize_labels))+
  theme_minimal()+
  theme(legend.position = "none",
        panel.spacing.x = unit(30,"pt"),
        panel.grid.major = element_line(color="#DFE2E0",linewidth = 0.2),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color="#DFE2E0",linewidth = 0.2),
        strip.text = element_text(face="bold"))
plot

plot_det_size<-ggplot(nD_pts)+
  geom_ribbon(data=area_grey,aes(x=n,ymin=D,ymax = ceiling(max(nD_pts$D))),fill="#DFE2E0")+
  geom_line(data=lin_tree,aes(x=n,y=D),color="#DFE2E0")+
  geom_jitter(aes(x=n,y=D,fill=color_par,shape=FE),
              width = jitter_width,
              height = jitter_height,
              color="white",stroke = 0.2,size=4,alpha=0.8)+
  scale_shape_manual(values = 21:24)+
  scale_fill_manual(values = colors_palette)+
  scale_y_continuous(limits = c(0,ceiling(max(nD_pts$D))))+
  scale_x_continuous(limits = c(0.5,ceiling(max(nD_pts$n))),breaks = 1:ceiling(max(nD_pts$n)))+
  facet_grid(det_size_reached~EndSize, labeller = labeller(EndSize = EndSize_labels))+
  theme_minimal()+
  theme(legend.position = "none",
        panel.spacing.x = unit(30,"pt"),
        panel.grid.major = element_line(color="#DFE2E0",linewidth = 0.2),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color="#DFE2E0",linewidth = 0.2),
        strip.text = element_text(face="bold"))
plot_det_size



nD_pts_EndSize<-nD_pts%>%
  nest_by(EndSize)%>%
  pull(data)
names(nD_pts_EndSize)<-unique(nD_pts$EndSize)
  
miniplot_list<-lapply(nD_pts_EndSize,
                      function(nD_pts_tbl){
                        jitter_width<-sd(nD_pts_tbl$n)/10
                        jitter_height<-sd(nD_pts_tbl$D)/10
                        
                        if(ceiling(max(nD_pts_tbl$n))<=2){
                          x_breaks<-round(seq(from=1,
                                              to=ceiling(max(nD_pts_tbl$n)),
                                              length.out=3),digits = 2)
                        }else{
                          x_breaks<-1:ceiling(max(nD_pts_tbl$n))
                        }
                        
                        
                        nD_pts_tbl_comp<-nD_pts_tbl%>%
                          filter(FE=="competition")
                        nD_pts_tbl_NOcomp<-nD_pts_tbl%>%
                          filter(FE!="competition")%>%
                          dplyr::select(-par_value)
                        colors_palette_comp<-colors_palette[grepl(x=names(colors_palette),pattern = "sus")]
                        nD_pts_tbl_comp$par_value<-factor(nD_pts_tbl_comp$par_value,levels = c(0.7,0.5,0.2,0.1))
                        plot_comp_enlighted<-ggplot(nD_pts_tbl_comp)+
                          geom_ribbon(data=area_grey,aes(x=n,ymin=D,ymax = ceiling(max(nD_pts$D))),fill="#DFE2E0")+
                          geom_line(data=lin_tree,aes(x=n,y=D),color="#DFE2E0")+
                          geom_jitter(data=nD_pts_tbl_NOcomp,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,shape=FE),
                                      fill="#DFE2E0",
                                      alpha=0.1,
                                      color="white",
                                      stroke = 0.2,
                                      size=2)+
                          geom_jitter(data=nD_pts_tbl_comp,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,fill=color_par,shape=FE),
                                      alpha=1,color="white",
                                      stroke = 0.2,
                                      size=3)+
                          scale_shape_manual(values = 21:24)+
                          facet_grid(.~par_value)+
                          scale_fill_manual(values = colors_palette_comp)+
                          scale_y_continuous(limits = c(min(nD_pts_tbl$D),ceiling(max(nD_pts_tbl$D))))+
                          scale_x_continuous(limits = c(min(nD_pts_tbl$n),ceiling(max(nD_pts_tbl$n))),breaks = x_breaks)+
                          theme_minimal()+
                          theme(legend.position = "none",
                                panel.grid = element_line(color="#DFE2E0",linewidth = 0.2),
                                panel.grid.minor = if(length(x_breaks)>3){element_blank()}else{element_line(color="#DFE2E0",linewidth = 0.2)},
                                axis.text = element_text(size=7),
                                axis.title = element_blank(),
                                strip.text = element_text(size=7,face="bold"))
                        plot_comp_enlighted
                        
                        nD_pts_tbl_growth<-nD_pts_tbl%>%
                          filter(FE=="growth")
                        nD_pts_tbl_NOgrowth<-nD_pts_tbl%>%
                          filter(FE!="growth")%>%
                          dplyr::select(-par_value)
                        colors_palette_growth<-colors_palette[grepl(x=names(colors_palette),pattern = "sel_adv")]
                        nD_pts_tbl_growth$par_value<-factor(nD_pts_tbl_growth$par_value,
                                                            levels = c(0.1,0.2,1),
                                                            labels = c("+MUT","+2MUT","+10MUT"))
                        plot_growth_enlighted<-ggplot(nD_pts_tbl_growth)+
                          geom_ribbon(data=area_grey,aes(x=n,ymin=D,ymax = ceiling(max(nD_pts$D))),fill="#DFE2E0")+
                          geom_line(data=lin_tree,aes(x=n,y=D),color="#DFE2E0")+
                          geom_jitter(data=nD_pts_tbl_NOgrowth,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,shape=FE),
                                      fill="#DFE2E0",
                                      alpha=0.1,
                                      color="white",
                                      stroke = 0.2,
                                      size=2)+
                          geom_jitter(data=nD_pts_tbl_growth,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,fill=color_par,shape=FE),
                                      alpha=1,color="white",
                                      stroke = 0.2,
                                      size=3)+
                          scale_shape_manual(values = 21:24)+
                          facet_grid(.~par_value)+
                          scale_fill_manual(values = colors_palette_growth)+
                          scale_y_continuous(limits = c(min(nD_pts_tbl$D),ceiling(max(nD_pts_tbl$D))))+
                          scale_x_continuous(limits = c(min(nD_pts_tbl$n),ceiling(max(nD_pts_tbl$n))),breaks =x_breaks)+
                          theme_minimal()+
                          theme(legend.position = "none",
                                panel.grid = element_line(color="#DFE2E0",linewidth = 0.2),
                                panel.grid.minor = if(length(x_breaks)>3){element_blank()}else{element_line(color="#DFE2E0",linewidth = 0.2)},
                                axis.text = element_text(size=7),
                                axis.title = element_blank(),
                                strip.text = element_text(size=7,face="bold"))
                        plot_growth_enlighted
                        
                        nD_pts_tbl_space<-nD_pts_tbl%>%
                          filter(FE=="space")
                        nD_pts_tbl_NOspace<-nD_pts_tbl%>%
                          filter(FE!="space")%>%
                          dplyr::select(-par_value)
                        colors_palette_space<-colors_palette[grepl(x=names(colors_palette),pattern = "space")]
                        nD_pts_tbl_space$par_value<-factor(nD_pts_tbl_space$par_value,
                                                      levels=sort(unique(nD_pts_tbl_space$par_value)),
                                                      labels=c("+K","+2K","+10K"))
                        plot_space_enlighted<-ggplot(nD_pts_tbl_space)+
                          geom_ribbon(data=area_grey,aes(x=n,ymin=D,ymax = ceiling(max(nD_pts$D))),fill="#DFE2E0")+
                          geom_line(data=lin_tree,aes(x=n,y=D),color="#DFE2E0")+
                          geom_jitter(data=nD_pts_tbl_NOspace,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,shape=FE),
                                      fill="#DFE2E0",
                                      alpha=0.1,
                                      color="white",
                                      stroke = 0.2,
                                      size=2)+
                          geom_jitter(data=nD_pts_tbl_space,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,fill=color_par,shape=FE),
                                      alpha=1,color="white",
                                      stroke = 0.2,
                                      size=3)+
                          scale_shape_manual(values = 21:24)+
                          facet_grid(.~par_value)+
                          scale_fill_manual(values = colors_palette_space)+
                          scale_y_continuous(limits = c(min(nD_pts_tbl$D),ceiling(max(nD_pts_tbl$D))))+
                          scale_x_continuous(limits = c(min(nD_pts_tbl$n),ceiling(max(nD_pts_tbl$n))),breaks = x_breaks)+
                          theme_minimal()+
                          theme(legend.position = "none",
                                panel.grid = element_line(color="#DFE2E0",linewidth = 0.2),
                                panel.grid.minor = if(length(x_breaks)>3){element_blank()}else{element_line(color="#DFE2E0",linewidth = 0.2)},
                                axis.text = element_text(size=7),
                                axis.title = element_blank(),
                                strip.text = element_text(size=7,face="bold"))
                        plot_space_enlighted
                        
                        nD_pts_tbl_passenger<-nD_pts_tbl%>%
                          filter(FE=="passenger")
                        nD_pts_tbl_NOpassenger<-nD_pts_tbl%>%
                          filter(FE!="passenger")%>%
                          dplyr::select(-par_value)
                        colors_palette_passenger<-colors_palette[grepl(x=names(colors_palette),pattern = "None")]
                        plot_passenger_enlighted<-ggplot(nD_pts_tbl_passenger)+
                          geom_ribbon(data=area_grey,aes(x=n,ymin=D,ymax = ceiling(max(nD_pts$D))),fill="#DFE2E0")+
                          geom_line(data=lin_tree,aes(x=n,y=D),color="#DFE2E0")+
                          geom_jitter(data=nD_pts_tbl_NOpassenger,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,shape=FE),
                                      fill="#DFE2E0",
                                      alpha=0.1,
                                      color="white",
                                      stroke = 0.2,
                                      size=2)+
                          geom_jitter(data=nD_pts_tbl_passenger,
                                      width = jitter_width,
                                      height = jitter_height,
                                      aes(x=n,y=D,fill=color_par,shape=FE),
                                      alpha=1,color="white",
                                      stroke = 0.2,
                                      size=3)+
                          scale_shape_manual(values = 21:24)+
                          scale_fill_manual(values = colors_palette_passenger)+
                          scale_y_continuous(limits = c(min(nD_pts_tbl$D),ceiling(max(nD_pts_tbl$D))))+
                          scale_x_continuous(limits = c(min(nD_pts_tbl$n),ceiling(max(nD_pts_tbl$n))),breaks = x_breaks)+
                          theme_minimal()+
                          theme(
                            legend.position = "none",
                            panel.grid.major = element_line(color="#DFE2E0",linewidth = 0.2),
                            panel.grid.minor = if(length(x_breaks)>3){element_blank()}else{element_line(color="#DFE2E0",linewidth = 0.2)},
                            axis.text = element_text(size=7),
                            axis.title = element_blank(),
                            strip.text = element_text(size=7))
                        plot_passenger_enlighted
                        
                        (plot_growth_enlighted/
                            plot_space_enlighted/
                            plot_comp_enlighted/
                            plot_passenger_enlighted)+
                          plot_layout(design=c("AAA#
                                               BBB#
                                               CCCC
                                               D###"))
                      })

plot_allFE_div<-wrap_plots(miniplot_list)
plot_allFE_div

ggsave(plot_allFE_div,
       device = "pdf",
       height = 5,
       width = 15,path="Data/nD_exp",
       file = "nD_FunEff_divided.pdf")

plot_allFE<-plot/wrap_plots(miniplot_list)/plot_det_size+
  plot_layout(heights = c(1,3,2))
plot_allFE

ggsave(plot_allFE,
        device = "pdf",
        height = 12,
        width = 15,
        path="Data/nD_exp",
        file = "nD_FunEff.pdf")
