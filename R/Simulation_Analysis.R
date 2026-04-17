#' Create a downloadable mutation tree plot
#'
#' Builds a mutation tree plot from simulation output, formatted for export or
#' download.
#'
#' @param path_sim Character string giving the path to the simulation output
#'   directory.
#' @param depth Detection threshold used to filter mutations by frequency (relative to tumor size).
#'  Only clones overcoming this threshold at least once in their lifetime will be plotted.
#' @param parameters A Parameters object corresponding to the parameters used in the simulation.
#' @param palette Named character vector of colors associated with phenotype.
#' @param freq if TRUE, the relative abundance is use as y-axis for the plot,
#'  if FALSE, absolute abundance is implemented.
#' @return A `ggplot2` object representing the mutation tree.
#'
#' @rdname get_muller_plot
#' @export
#' @import ggplot2
#' @import dplyr
#' @import colorspace
#' @importFrom tibble tibble
#' @import tidyr
#' @import colorspace
#' @importFrom latex2exp TeX
#' @importFrom scales scientific
#' @import stringr
#' @import patchwork
setGeneric("get_muller_plot", function(path_sim, depth, parameters,freq=FALSE, palette=NULL) standardGeneric("get_muller_plot"))

setMethod("get_muller_plot",
          signature(path_sim="character",
                    depth="numeric",
                    parameters="Parameters"),
          function(path_sim,
                   depth,
                   parameters,
                   freq=FALSE,
                   palette=NULL){
            obs_tumor<-get_obs_tum(path_sim, depth, parameters)
            Clones_df<-get_muller_plot_info(obs_Pop_ID = obs_tumor$obs_Pop_ID,
                                            obs_tumor_tibble = obs_tumor$obs_tumor_tibble,
                                            freq = freq,
                                            functional_effects = parameters@functional_effects)
            
            if(is.null(palette)){
              hues <- seq(0, 360, length.out = length(unique(Clones_df$fun_eff)) + 1)[-1]  
              palette <- hsv(h = hues / 360, s = 0.3, v = 0.7) 
              names(palette)<-unique(Clones_df$fun_eff)
            }
            
            p<-get_muller_plot_download(Clones_df = Clones_df,
                                        freq = relative,
                                        palette = palette)
            
            p <- p + plot_layout(guides = "collect") & 
              theme(legend.position = "bottom")
            
            return(p)
          })


#' Create a downloadable mutation tree plot
#'
#' Builds a mutation tree plot from simulation output, formatted for export or
#' download.
#'
#' @param path_sim Character string giving the path to the simulation output
#'   directory.
#' @param depth Detection threshold used to filter mutations by abundance or
#'   frequency.
#' @param parameters A Parameters object corresponding to the parameters used in the simulation.
#' @param palette Named character vector of colors associated with functional
#'   effects.
#' @return A `ggplot2` object representing the mutation tree.
#'
#' @rdname get_tree_plot
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr unnest
#' @importFrom stringr str_remove
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph create_layout
#' @importFrom tibble tibble
#' @importFrom scales rescale
setGeneric("get_tree_plot", function(path_sim, depth, parameters, palette=NULL) standardGeneric("get_tree_plot"))

setMethod("get_tree_plot",
          signature(path_sim="character",depth="numeric",parameters="Parameters"),
          function(path_sim,depth,parameters,palette=NULL){
            Zprovvs<-list.files(path_sim)
            num_seq<-which.max(as.numeric(str_remove(str_remove(Zprovvs,"Zprovv"),".RData")))
            load(paste(path_sim,Zprovvs[num_seq],sep="/"))
            time_provv<-parameters@print_time[which.min(abs(time_provv-parameters@print_time))]
            
            pop<-lapply(Zprovv,Population)
            gen<-lapply(pop,genotype)
            fun_eff_num<-lapply(pop,functional_effect)
            fun_eff<-lapply(fun_eff_num,function(n){names(parameters@functional_effects)[n]})
            ncells<-sapply(Zprovv,Ncells)
            tot_ncells<-sum(ncells)
            unique_mut_id<-lapply(gen,function(g){
              unique_id_mut<-vector()
              for(j in 1:length(g)){
                unique_id_mut<-c(unique_id_mut,paste(g[1:j],collapse="_"))
              }
              return(unique_id_mut)
            })
            
            all_mut<-unique(unlist(unique_mut_id))
            used_nums<-as.numeric(gsub(pattern = "Mut",x = names(mut_names),replacement = ""))
            mut_nums<-1:length(all_mut)
            not_used_mut_nums<-mut_nums[!mut_nums%in%used_nums]
            mut_names_tbl<-tibble::tibble(mut=c(mut_names,all_mut[!all_mut%in%mut_names]),
                                          names=c(names(mut_names),paste0("Mut",not_used_mut_nums,sep="",recycle0 = TRUE)))
            
            parents<-lapply(unique_mut_id,lag)
            mut_generation<-lapply(gen, function(g){1:length(g)})
            
            composition<-tibble::tibble(mut=unique_mut_id,parents,fun_eff,ncells,mut_generation)%>%
              unnest(c(mut,parents,fun_eff,mut_generation))%>%
              group_by(mut,parents,fun_eff,mut_generation)%>%
              summarise(ncells=sum(ncells),
                        frequency=ncells/tot_ncells)%>%
              filter(frequency>depth)%>%
              ungroup()%>%
              arrange(desc(frequency))%>%
              merge(mut_names_tbl)
            
            
            roots<-composition$mut[is.na(composition$parents)]
            
            
            if(length(roots)>1){
              composition<-composition%>%
                mutate(parents=ifelse(is.na(parents),"0",parents))%>%
                bind_rows(tibble(mut="0"))
            }
            
            if(all(is.na(composition$parents))){
              wanted_mut<-composition$mut
              composition$y<-1:length(wanted_mut)
              
              plot<-ggplot(composition)+
                geom_label(aes(x=0,y=y,label=names,
                               fill=fun_eff),
                           label=composition$names,
                           color="white",
                           size=6,
                           label.r =unit(0.5,"lines"),
                           label.size = 0)+
                xlim(-0.1,0.1)+
                ylim(0.9,length(composition$mut)+0.1)+
                coord_fixed()+
                scale_fill_manual(values=palette,na.value = "white")+
                theme_void()+
                theme(legend.position = "none")
            }else{
              if(length(composition$mut)<10){
                wanted_mut<-composition$mut
              }else{
                wanted_mut<-unname(
                  unlist(
                    lapply(split(composition$mut,
                                 composition$parents),
                           function(v){
                             v_filt<-v[(v%in%unique(composition$parents))]
                             v_filt<-c(v_filt,v[!v%in%v_filt][1])
                             if(length(v_filt)==length(v)-1){v_filt<-v}
                             return(v_filt)
                           }
                    )
                  )
                )
              }
              wanted_mut<-unique(c(wanted_mut,unique(composition$mut[is.na(composition$parents)])))
              
              mut_info<-composition%>%
                mutate(wanted=ifelse(mut%in%wanted_mut,TRUE,FALSE))%>%
                group_by(parents)%>%
                summarise(ndaught_drawn=sum(wanted),
                          ndaught_hidden=n()-ndaught_drawn)%>%
                filter(!is.na(parents)&ndaught_hidden>0)%>%
                rowwise()%>%
                mutate(mut=paste(parents,"n",sep="_"),
                       label=paste(ndaught_hidden,"more"),
                       fun_eff=as.character(NA),
                       mut_generation=unique(composition$mut_generation[composition$mut==parents])+1)%>%
                dplyr::select(mut,parents,label,mut_generation,fun_eff)%>%
                ungroup()
              
              mut_info<-composition%>%
                filter(mut%in%wanted_mut)%>%
                mutate(label=names)%>%
                dplyr::select(mut,parents,label,fun_eff,mut_generation)%>%
                bind_rows(mut_info)
              
              Zprovvs_ordered<-Zprovvs
              Zprovvs_ordered[as.numeric(str_remove(str_remove(Zprovvs,"Zprovv"),".RData"))+1]<-Zprovvs
              
              mut_info$time_appearance<-sapply(1:nrow(mut_info),function(j){
                mut_id<-mut_info$mut[j]
                if(grepl(x = mut_info$label[j],pattern = "more")){
                  return(NA)
                }
                t<-time_provv
                for(k in 1:length(Zprovvs_ordered)){
                  Z<-rev(Zprovvs_ordered)[k]
                  load(paste(path_sim,Z,sep="/"))
                  pop<-lapply(Zprovv,Population)
                  gen<-lapply(pop,genotype)
                  unique_mut_id<-lapply(gen,function(g){
                    unique_id_mut<-vector()
                    for(j in 1:length(g)){
                      unique_id_mut<-c(unique_id_mut,paste(g[1:j],collapse="_"))
                    }
                    return(unique_id_mut)
                  })
                  
                  if(!mut_id%in%unique(unlist(unique_mut_id))){break}else{
                    t<-time_provv
                  }
                }
                return(t)
              })
              NA_rows<-which(is.na(mut_info$time_appearance))
              for(j in NA_rows){
                parent<-mut_info$parents[j]
                mut_info$time_appearance[j]<-min(mut_info$time_appearance[-NA_rows][mut_info$parents[-NA_rows]==parent],na.rm = TRUE)
              }
              
              g <- graph_from_data_frame(mut_info%>%dplyr::select(parents,mut)%>%filter(!is.na(parents)), directed = TRUE)
              layout_composition <- create_layout(g, layout = "dendrogram", circular = FALSE)
              
              x_grid<-length(unique(layout_composition$y))
              y_grid<-sum(layout_composition$leaf)
              x_range<-c(0,1.5)
              y_range<-c(0,1)
              
              layout_composition<-layout_composition%>%
                filter(name!=0)%>%
                rowwise()%>%
                mutate(y=x_grid-unique(mut_info$mut_generation[mut_info$mut==name]))%>%
                ungroup()
              
              nodes_coord <- tibble(
                mut = layout_composition$name,
                y = rescale(-layout_composition$x, to = y_range)
              ) %>%
                left_join(
                  mut_info %>% dplyr::select(mut, time_appearance),
                  by = "mut"
                ) %>%
                mutate(
                  x = rescale(time_appearance, to = x_range)
                ) %>%
                left_join(
                  mut_info %>% dplyr::select(mut, label, fun_eff, mut_generation),
                  by = "mut"
                ) %>%
                group_by(mut_generation) %>%
                mutate(n_mut_layer = n()) %>%
                ungroup()
              
              
              size_label<-min(40/max(nodes_coord$n_mut_layer),6)
              x_lim<-range(nodes_coord$x)+c(-0.1,0.1)
              y_lim<-range(nodes_coord$y)+c(-0.1,0.1)
              
              if(sum(!is.na(unique(layout_composition$y)))==1){
                plot<-ggplot()+
                  geom_label(data=nodes_coord,
                             aes(x=x,
                                 y=y,
                                 label=label,
                                 fill=fun_eff,
                                 color=fun_eff),
                             size=size_label,
                             label.r =unit(0.5,"lines"),
                             label.size = 0)+
                  xlim(x_lim)+
                  ylim(y_lim)+
                  coord_fixed()+
                  scale_color_manual(values=rep("white",length(parameters@functional_effects)),na.value = "black")+
                  scale_fill_manual(values=palette,na.value = "white")+
                  theme_void()+
                  theme(legend.position = "none")
              }else{
                
                
                edges_coord<-tibble(
                  x1=sapply(mut_info$parents[!is.na(mut_info$parents)],
                            function(parent){nodes_coord$x[nodes_coord$mut==parent]}),
                  y1=sapply(mut_info$parents[!is.na(mut_info$parents)],
                            function(parent){nodes_coord$y[nodes_coord$mut==parent]}),
                  x2=sapply(mut_info$mut[!is.na(mut_info$parents)],
                            function(mut){nodes_coord$x[nodes_coord$mut==mut]}),
                  y2=sapply(mut_info$mut[!is.na(mut_info$parents)],
                            function(mut){nodes_coord$y[nodes_coord$mut==mut]}),
                  x_mid=sapply(mut_info$parents[!is.na(mut_info$parents)],
                               function(parent){
                                 daught<-mut_info$mut[!is.na(mut_info$parents)&mut_info$parents==parent]
                                 min_dist_daught<-min(nodes_coord$x[nodes_coord$mut%in%daught])
                                 return((min_dist_daught+nodes_coord$x[nodes_coord$mut==parent])/2)})
                )
                
                radius<-edges_coord%>%
                  rowwise()%>%
                  mutate(r=min(abs(x1-x2),abs(y1-y2))/2)%>%
                  pull(r)
                
                r<-min(radius[radius>0])
                
                link<-edges_coord%>%
                  rowwise()%>%
                  mutate(link=list(as_tibble(elbowed_link(x1,y1,x2,y2,r,x_mid))))%>%
                  ungroup()%>%
                  pull(link)%>%bind_rows(.id="mut")
                
                x_lim<-range(nodes_coord$x)+c(-0.1,0.1)
                y_lim<-range(nodes_coord$y)+c(-0.1,0.1)
                
                if(is.null(palette)){
                  fun_eff<-parameters@functional_effects
                  base_colors <- c(
                    "growth"      = "#C2D4DD",
                    "mutation"    = "#E4DDC3",
                    "competition" = "#E6D1D1",
                    "space"       = "#AEAEC2",
                    "passenger"   = "#B6CBBC"
                  )
                  
                  jitter_color <- function(hex_color) {
                    rgb_obj <- hex2RGB(hex_color)
                    hls_obj <- as(rgb_obj, "HLS")
                    
                    coords <- hls_obj@coords
                    coords[1, "L"] <- max(0, min(1, coords[1, "L"] + runif(1, -0.06, 0.06)))
                    coords[1, "S"] <- max(0, min(1, coords[1, "S"] + runif(1, -0.04, 0.04)))
                    hls_mod <- HLS(coords[1, "H"], coords[1, "L"], coords[1, "S"])
                    return(hex(hls_mod))
                  }
                  
                  clamp <- function(x, min, max) pmin(pmax(x, min), max)
                  
                  palette <- sapply(fun_eff, function(x) jitter_color(base_colors[x]))
                  
                  names(palette) <- names(fun_eff)
                }
                
                plot<-ggplot()+
                  geom_path(data=link,
                            aes(x=x,y=y,group=mut))+
                  geom_point(data=nodes_coord%>%filter(!is.na(fun_eff)),
                             aes(x=x,y=y,fill=fun_eff),color="black",shape=21,size=5)+
                  geom_label(data=nodes_coord%>%filter(is.na(fun_eff)),
                             aes(x=x,y=y,label=label),color="black",border.colour = "white",fill="white")+
                  
                  xlim(x_lim)+
                  ylim(y_lim)+
                  coord_fixed()+
                  scale_fill_manual(values=palette,na.value = "white")+
                  theme_void()+
                  theme(legend.position = "none")
                
              }}
            
            return(plot)
          })

#' Simulate sequencing from tumor data
#'
#' This function simulates bulk or multi-regional sequencing starting from either
#' a tumor slice (`Zprovv`) or a precomputed ordered clone table (`Clones_df`).
#' It generates one or more VCF-like tibbles containing sampled variant allele
#' frequencies (VAF), read depths (DP), and alternative allele counts (AD).
#'
#' @param Zprovv A list representing a sequenced tumor slice. Required if
#'   `Clones_df` is not provided.
#' @param Clones_df A tibble of ordered clones as returned by
#'   \code{get_ordered_clones_sequencing()}. If provided, `Zprovv` is not required.
#' @param parameters An object of class \code{Parameters} containing model parameters.
#' @param seed Optional numeric seed for reproducibility. If \code{NULL}, a seed is generated.
#' @param n_regions Integer. Number of spatial regions to sequence (1 = bulk).
#' @param n_seq_cells Integer. Number of cells sequenced per region. Use \code{Inf}
#'   to sequence the entire tumor mass.
#' @param Nrep Integer. Number of independent sequencing repetitions.
#' @param dens An object with class "density" containing the density used to sample sequencing depth (DP).
#'   If NULL, a default density derived from TCGA data is used.
#'
#' @return A list of tibbles (length \code{Nrep}), each representing a simulated
#'   sequencing experiment. Each tibble contains:
#'   \itemize{
#'     \item \code{region}: region identifier
#'     \item \code{mut}: mutation identifier
#'     \item \code{fun_eff}: functional effect label
#'     \item \code{fun_eff_category}: functional effect category
#'     \item \code{sample_DP}: sampled read depth
#'     \item \code{sample_AD}: sampled alternative allele count
#'     \item \code{VAF}: variant allele frequency
#'   }
#'
#' @details
#' If \code{Clones_df} is not provided, it is computed internally from \code{Zprovv}.
#' The function assumes that clone genotypes and functional effects are consistent
#' with the provided \code{parameters}.
#'
#' @examples
#' \dontrun{
#' res <- sequencing(
#'   Zprovv = Zprovv,
#'   parameters = parameters,
#'   n_regions = 2,
#'   n_seq_cells = 500,
#'   Nrep = 3
#' )
#' }
#'
#' @export
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
#' @importFrom stats rnorm rbinom
#' @import methods
sequencing <- function(
    Zprovv = NULL,
    Clones_df = NULL,
    parameters,
    seed = NULL,
    n_regions = 1,
    n_seq_cells = 1000,
    Nrep = 1,
    dens = NULL
) {
  
  if (is.null(Zprovv) && is.null(Clones_df)) {
    stop("Either Zprovv or Clones_df must be provided.")
  }
  
  if (is.null(Clones_df)) {
    Clones_df <- get_ordered_clones_sequencing(Zprovv)
  }
  
  if (is.null(dens)) {
    dens <- default_dens
  }
  
  if (is.null(seed)) {
    seed_selected <- as.integer(Sys.time())
  } else {
    seed_selected <- seed
  }
  set.seed(seed_selected)
  
  pop <- lapply(Zprovv, Population)
  Pop_ID <- seq_along(pop)
  ncells <- sapply(Zprovv, Ncells)
  
  gen <- lapply(pop, genotype)
  fun_eff <- lapply(pop, functional_effect)
  fun_eff_label <- lapply(fun_eff, function(f) {
    names(parameters@functional_effects)[f]
  })
  
  unique_mut_id <- lapply(gen, function(g) {
    sapply(seq_along(g), function(i) paste(g[1:i], collapse = "_"))
  })
  
  results <- lapply(seq_len(Nrep), function(nrep) {
    
    if (n_regions * n_seq_cells >= sum(ncells)) {
      n_cells_eff <- sum(ncells)
      seq_min_y_list <- min(Clones_df$y_lower)
    } else {
      L <- (max(Clones_df$y_upper) - min(Clones_df$y_lower)) -
        (n_regions - 1) * n_seq_cells
      
      u <- sort(runif(n_regions, 0, L))
      seq_min_y_list <- min(Clones_df$y_lower) +
        u + (0:(n_regions - 1)) * n_seq_cells
      n_cells_eff <- n_seq_cells
    }
    
    vcf_multiregional <- dplyr::bind_rows(
      lapply(seq_min_y_list, function(seq_min_y) {
        
        seq_max_y <- seq_min_y + n_cells_eff
        
        Clones_df_seq <- Clones_df %>%
          dplyr::filter(y_lower < seq_max_y, y_upper > seq_min_y)
        
        mut_fe_to_add <- lapply(Clones_df_seq$clone, function(c) {
          g <- gen[[c]]
          lg <- length(g)
          
          if (lg == 1) return(tibble::tibble())
          
          f <- fun_eff[[c]]
          
          out <- purrr::map_dfr(seq_len(lg - 1), function(i) {
            newg <- head(g, -i)
            newf <- names(parameters@functional_effects)[f[lg - i]]
            
            if (any(sapply(gen[Clones_df_seq$clone], identical, newg))) {
              return(NULL)
            }
            
            tibble::tibble(
              mut = paste(newg, collapse = "_"),
              fun_eff = newf,
              y_lower = seq_min_y,
              y_upper = seq_max_y
            )
          })
          
          out
        }) %>% dplyr::bind_rows() %>% dplyr::distinct()
        
        Clones_df_seq <- merge(
          Clones_df_seq,
          tibble::tibble(
            clone = Pop_ID,
            mut = sapply(unique_mut_id, tail, 1),
            fun_eff = sapply(fun_eff_label, tail, 1)
          )
        ) %>%
          dplyr::bind_rows(mut_fe_to_add) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            fun_eff_category = parameters@functional_effects[fun_eff],
            y2 = min(y_upper, seq_max_y),
            y1 = max(y_lower, seq_min_y),
            Ncells_seq = round(y2 - y1),
            prob = Ncells_seq / (2 * (seq_max_y - seq_min_y))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(mut, fun_eff, fun_eff_category, Ncells_seq, prob)
        
        sample_DP <- round(
          sample(dens$x, nrow(Clones_df_seq),
                 prob = dens$y, replace = TRUE) +
            stats::rnorm(1, 0, dens$bw)
        )
        sample_DP[sample_DP < 0] <- 0
        
        sample_AD <- mapply(
          stats::rbinom,
          prob = Clones_df_seq$prob,
          size = sample_DP,
          MoreArgs = list(n = 1)
        )
        
        Clones_df_seq %>%
          dplyr::mutate(
            sample_DP = sample_DP,
            sample_AD = sample_AD
          ) %>%
          dplyr::filter(sample_AD > 0) %>%
          dplyr::mutate(VAF = sample_AD / sample_DP) %>%
          dplyr::select(-prob, -Ncells_seq)
        
      }),
      .id = "region"
    )
    
    vcf_multiregional
  })
  
  return(results)
}
