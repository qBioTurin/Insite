library(Insite)
library(stringr)
library(colorspace)
library(tibble)

tum_types<-c("AML","breast","kidney","lung","mesothelioma","uveal")
load("Data/best_fit/nD_points_bestfit.RData")

sim_division<-list()
for(tum_type in c("breast","kidney","lung")){
  clus<-sim_data%>%
    filter(dataset==tum_type)%>%
    pull(n)%>%
    kmeans(centers = 2)
  
  names_clus<-vector()
  names_clus[which.max(clus$centers)]<-"high n"
  names_clus[which.min(clus$centers)]<-"low n"
  
  sim_division[[tum_type]]<-sim_data%>%
    filter(dataset==tum_type)%>%
    bind_cols(cluster=clus$cluster)%>%
    mutate(name=names_clus[cluster])
  
}
d<-bind_rows(sim_division)%>%
  bind_rows(sim_data%>%
              filter(dataset %in% c("AML","mesothelioma","uveal"))%>%
              mutate(name=""))%>%
  mutate(time_end=time_end/365)%>%
  dplyr::select(dataset,time_end,name)



times<-d%>%
  ungroup()%>%
  nest_by(dataset,name)%>%
  pull(data)%>%
  lapply(function(t){pull(t,time_end)})

names(times)<-d%>%
  ungroup()%>%
  nest_by(dataset,name)%>%
  dplyr::select(dataset,name)%>%
  mutate(name=paste(dataset,name,sep="_"))%>%
  pull(name)

pvals<-sapply(
  times,
  function(v1){
    sapply(times,function(v2){
      t_test<-t.test(v1,v2)
      return(t_test$p.value)
    })
  }
)



d1<-d%>%
  group_by(name,dataset)%>%
  summarise(bar_y=max(time_end)*1.1)%>%
  filter(name!="")%>%
  pivot_wider(values_from = bar_y,names_from = name)%>%
  rename("high_n"="high n",
         "low_n"="low n")%>%
  rowwise()%>%
  mutate(hbar_y= max(high_n,low_n)*1.05)


pvals<-pvals[!rownames(pvals)%in%c("AML_","uveal_","mesothelioma_"),!colnames(pvals)%in%c("AML_","uveal_","mesothelioma_")] %>%
  as.data.frame() %>%
  rownames_to_column(var = "row_name") %>%
  pivot_longer(
    cols = -row_name, 
    names_to = "col_name", 
    values_to = "pval"
  )%>%
  filter(as.character(row_name) < as.character(col_name))%>%
  separate_wider_delim(cols = row_name,delim = "_",names = c("dataset","cluster"))%>%
  separate_wider_delim(cols = col_name,delim = "_",names = c("dataset1","cluster1"))%>%
  filter(dataset==dataset1)%>%
  merge(d1)%>%
  dplyr::select(dataset,pval,hbar_y)%>%
  mutate(pval_stars=ifelse(pval<10^(-4),"***",
                           ifelse(pval<10^(-3),"**",
                                  ifelse(pval<10^(-2),"*",round(pval,2)))))

p<-ggplot()+
  geom_boxplot(data = d%>%
                 mutate(name=factor(name,levels=c("low n","high n",""))),
               aes(x=as.factor(name),y=time_end,fill=dataset,color=dataset),alpha=0.8)+
  geom_segment(data = d1,aes(x="low n", xend="high n",y = hbar_y,yend=hbar_y),
               linejoin = "round",linewidth = 0.2)+
  geom_segment(data = d1,aes(x="low n", xend="low n",y = hbar_y,yend=low_n),
               linejoin = "round",linewidth = 0.3)+
  geom_segment(data = d1,aes(x="high n", xend="high n",y = hbar_y,yend=high_n),
               linejoin = "round",linewidth = 0.3)+
  geom_label(data = pvals,aes(x=1.5,y=hbar_y,label=pval_stars),vjust=0.7,fill = "white",border.color = "white")+
  facet_wrap("dataset",scales="free")+
  scale_x_discrete(breaks=c("","low n","high n"))+
  scale_fill_manual(values=c("#e76f51",
                             "#f4a261",
                             "#e9c46a",
                             "#8ab17d",
                             "#2a9d8f",
                             "#264653"))+
  scale_color_manual(values=c("#e76f51",
                              "#f4a261",
                              "#e9c46a",
                              "#8ab17d",
                              "#2a9d8f",
                              "#264653"))+
  theme_void()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7,colour = "darkgrey"),
        strip.text = element_text(size = 14),
        panel.grid.major.y = element_line(colour = "lightgrey",linewidth = 0.1),
        panel.grid.minor.y = element_line(colour = "lightgrey",linewidth = 0.1))
p
