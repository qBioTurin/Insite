library(Insite)
if (!require("dplyr")) install.packages("dplyr",repos = "https://cloud.r-project.org")
if (!require("tidyr")) install.packages("tidyr",repos = "https://cloud.r-project.org")
if (!require("ggplot2")) install.packages("ggplot2",repos = "https://cloud.r-project.org")
if (!require("stringr")) install.packages("stringr",repos = "https://cloud.r-project.org")

load("Data/best_fit/Noble_tum_data.RData")
real_points$dataset[real_points$dataset=="breast_SC"]<-"breast"

load("Data/best_fit/nD_points_bestfit.RData")

sim_division<-list()
for(tum_type in c("breast","kidney","lung")){
  clus<-sim_data%>%
    filter(dataset==tum_type)%>%
    pull(n)%>%
    kmeans(centers = 2)
  
  names_clus<-vector()
  names_clus[which.max(clus$centers)]<-"Right Group"
  names_clus[which.min(clus$centers)]<-"Left Group"
  
  sim_division[[tum_type]]<-sim_data%>%
    filter(dataset==tum_type)%>%
    bind_cols(cluster=clus$cluster)%>%
    mutate(name=names_clus[cluster])
    
}

v_line_df<-lapply(sim_division,function(sim_data_clustered){
  vline_df<-sim_data_clustered%>%
    group_by(name)%>%
    summarise(max_n=max(n),
           min_n=min(n))
  vline_min<-vline_df$max_n[vline_df$name=="Left Group"]
  vline_max<-vline_df$min_n[vline_df$name=="Right Group"]
  vline<-(vline_min+vline_max)/2
  return(vline)
})%>%
  as_tibble()%>%
  pivot_longer(cols = c("breast","kidney","lung"),
               names_to = "dataset",
               values_to="x_intercept")

plot<-ggplot()+
  geom_point(
    data=sim_data%>%arrange(centroid),
    aes(x=n,y=D,group=sim,fill=dataset,alpha=centroid,color=interaction(centroid,dataset),size=centroid),
    shape=21
  )+
  geom_vline(data = v_line_df,aes(xintercept = x_intercept),
             color="grey",linetype="dashed")+
  geom_point(
    data=real_points,
    aes(x=n,y=D),colour = "black",shape=4,stroke=1
  )+
  facet_wrap("dataset",scales='free')+
  scale_fill_manual(values=c("#e76f51",
                             "#f4a261",
                             "#e9c46a",
                             "#8ab17d",
                             "#2a9d8f",
                             "#264653"))+
  scale_alpha_manual(values=c(0.5,1))+
  scale_size_manual(values = c(2,3))+
  scale_color_manual(values = c("#f7c8be","white",
                                "#fbdbc4","white",
                                "#f6e9c7","white",
                                "#d3e1cd","white",
                                "#b6dbd4","white",
                                "#adb9be","white"))+
  scale_x_continuous(limits=c(0,14),
                     breaks = seq(from=0,to=14,by=2)) + 
  scale_y_continuous(limits=c(0,20),
                     breaks = seq(from=0,to=20,by=5))+
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

plot

for(tum_type in unique(sim_data$dataset)){
  load(paste("Data/best_fit",tum_type,"vcf_multiregional_allsims.RData",sep = "/"))
  
  mut_charact<-vcf_multiregional_all%>%
    dplyr::group_by(mut,sim)%>%
    mutate(mean_VAF=mean(VAF))%>%
    ungroup()%>%
    dplyr::count(mut,sim,fun_eff,mean_VAF)%>%
    mutate(patterns=ifelse(n>=9,"Public",
                           ifelse(n>=4,"Public Regional",
                                  ifelse(n>1,"Private Regional","Private Unique"))))%>%
    dplyr::count(patterns,sim,fun_eff,mean_VAF)
  
  if(tum_type=="AML"){
    fe_lev<-sort(unique(mut_charact$fun_eff))
  fe_lev<-fe_lev[fe_lev!="WT"]
  fe_lev<-c("WT",fe_lev)
  
  plot<-mut_charact%>%
    group_by(patterns,fun_eff)%>%
    mutate(mean_VAF=mean(mean_VAF),
           n=sum(n),
           fun_eff=factor(fun_eff,levels = fe_lev),
           patterns=factor(patterns,levels = c("Public","Public Regional","Private Regional","Private Unique")))%>%
    ggplot()+
    geom_col(aes(x=fun_eff,color=patterns,fill=patterns,y = n),position="fill",width=0.5)+
    scale_fill_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
    scale_color_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
    
    theme_void()+
    coord_fixed(ratio = 4)+
    labs(title = tum_type)+
    theme(axis.text.x = element_text(angle=90,hjust = 1,size=9),
          title = element_text(face = "bold"),
          legend.position = "none")
  
  print(plot)
  
  }else if(tum_type=="mesothelioma"){
    plot<-mut_charact%>%
      mutate(fe_group=factor(ifelse(fun_eff=="WT","WT",
                                    ifelse(grepl("Strong",fun_eff),"Strong Res",
                                           ifelse(grepl("Res",fun_eff),"Mild Res",
                                                  ifelse(grepl("S",fun_eff),"Limit Evasion",
                                                         "Proliferation")))
      ),
      levels=c("WT","Proliferation","Limit Evasion","Mild Res","Strong Res"))
      )%>%
      group_by(patterns,fe_group)%>%
      mutate(mean_VAF=mean(mean_VAF),
             n=sum(n))%>%
      mutate(patterns=factor(patterns,levels = c("Public","Public Regional","Private Regional","Private Unique")))%>%
      ggplot()+
      geom_col(aes(x=fe_group,color=patterns,fill=patterns,y = n),position="fill",width=0.5)+
      scale_fill_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
      scale_color_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
      theme_void()+
      coord_fixed(ratio = 4)+
      labs(title = tum_type)+
      theme(axis.text.x = element_text(angle=90,hjust = 1,size=9),
            title = element_text(face = "bold"),
            legend.position = "none")
    
    print(plot)
  }else if(tum_type=="uveal"){
    plot<-mut_charact%>%
      mutate(fe_group=factor(ifelse(fun_eff=="WT","WT",
                                    ifelse(fun_eff=="Res1","Strong Res",
                                           ifelse(fun_eff=="Res2","Med Res",
                                                  ifelse(grepl("Res",fun_eff),"Mild Res",
                                                         "Limit Evasion")))
      ),
      levels=c("WT","Limit Evasion","Mild Res","Med Res","Strong Res"))
      )%>%
      group_by(patterns,fe_group)%>%
      mutate(mean_VAF=mean(mean_VAF),
             n=sum(n))%>%
      mutate(patterns=factor(patterns,levels = c("Public","Public Regional","Private Regional","Private Unique")))%>%
      ggplot()+
      geom_col(aes(x=fe_group,color=patterns,fill=patterns,y = n),position="fill",width=0.5)+
      scale_fill_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
      scale_color_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
      theme_void()+
      coord_fixed(ratio = 4)+
      labs(title = tum_type)+
      theme(axis.text.x = element_text(angle=90,hjust = 1,size=9),
            title = element_text(face = "bold"),
            legend.position = "none")
    
    print(plot)
  }else{
    plot<-mut_charact%>%
      left_join(sim_division[[tum_type]]%>%
                  ungroup()%>%
              mutate(sim=str_remove(sim,"sim"))%>%
                dplyr::select(sim,name),by = "sim")%>%
      mutate(fe_group=factor(ifelse(fun_eff=="WT","WT",
                                    ifelse(grepl("Strong",fun_eff),"Strong Res",
                                           ifelse(grepl("Res",fun_eff),"Mild Res",
                                                  ifelse(grepl("S",fun_eff),"Limit Evasion",
                                                         "Proliferation")))),
                             levels = c("WT","Proliferation","Limit Evasion","Mild Res","Strong Res"))
      )%>%
      dplyr::group_by(patterns,fe_group,name)%>%
      dplyr::mutate(mean_VAF=mean(mean_VAF),
             n=sum(n))%>%
      mutate(patterns=factor(patterns,levels = c("Public","Public Regional","Private Regional","Private Unique")))%>%
      ggplot()+
      geom_col(aes(x=name,color=patterns,fill=patterns,y = n),position="fill")+
      facet_wrap("fe_group",ncol=4)+
      scale_fill_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
      scale_color_manual(values = c("#5F021F","#9A5B62","#C1BEA3","#FAF8D4"))+
      theme_void()+
      coord_fixed(ratio = 7)+
      labs(title = tum_type)+
      theme(axis.text.x = element_text(angle=90,hjust = 1,size=9),
            title = element_text(face = "bold"),
            legend.position = "none")
    
    print(plot)
  }
}


