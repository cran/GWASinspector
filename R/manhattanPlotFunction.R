plot_manhattan_standalone<-function(study.sample, plot.title, plot.subtitle, sig.threshold.log, fileName){


  if(nrow(study.sample) < 1)
    stop('Manhattan plot skipped! not enough variables.', call. = FALSE)


  ## p < 10^-300 => 10^-300
  study.sample <- correct_extreme_pvalues(study.sample)

  study.sample <- study.sample[order(CHR,POSITION,decreasing = FALSE)]
  ###

  ### calculations
  y.scale.min = floor(-log10(max(study.sample$PVALUE, na.rm = TRUE)))  ## minimum value for Y axis

  # plot to -log10 == 10 even if no pvalue is as low as that
  p.Min = ceiling(-log10(min(study.sample$PVALUE, na.rm = TRUE)))
  p.Min = ifelse(p.Min < 10 , 10 , p.Min)

  if (p.Min %% 2 == 0){
    plot.scale = seq(2, p.Min , 2)
  }else{
    plot.scale = seq(2, p.Min+1 ,2)
  }
  ###

  #----
  # mathematical expression for y axis titles
  log10P <- expression(paste("-log"[10], "(",italic(p),")"))

  study.sample <- fill_missing_chr_positions(study.sample)

  man_plot<-ggplot(study.sample,aes(study.sample$POSITION,-log10(study.sample$PVALUE) ,colour = factor(study.sample$CHR))) +
    scale_y_continuous(log10P, limits=c(y.scale.min,p.Min)) + #, breaks = plot.scale) +
    labs(x="Chromosome",y=" Observed -log(P-value)",
         subtitle=plot.subtitle , title=plot.title) +
    facet_grid(.~study.sample$CHR, scales = "free",space = 'free', switch = "both", margins = TRUE) +
    geom_point() +
    theme_grey() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size=12, face = "bold")
          ,strip.text.y = element_text(size=14, face = "bold")
          ,axis.text.x = element_blank()
          ,axis.ticks.x = element_blank()
          ,legend.position="none"
          ,plot.title=element_text(size=30, face="bold",hjust = 0.5)
          ,plot.subtitle=element_text(size=15, face="bold",hjust = 0.5)
          ,axis.text.y=element_text(size=15)
          ,axis.title.x=element_text(size=25)
          ,axis.title.y=element_text(size=25)) +
    geom_hline(data = study.sample, aes(yintercept = sig.threshold.log),color="red")


  ggsave(plot=man_plot,
         device = .QC$graphic.device,
         filename = fileName,
         width = 400,
         height = 200,
         units = 'mm',
         dpi = 200)


  #### write log ####
  message("Manhattan plot saved!")
  #### remove variables from RAM

  rm(study.sample,
     man_plot)

  #invisible(gc())

}



plot_manhattan<-function(study.sample, plot.title, plot.subtitle, sig.threshold.log, fileName){


  if(nrow(study.sample) < 1)
  {
    print_and_log('Manhattan plot skipped! not enough variables.','warning',display=.QC$config$debug$verbose)
    return(NULL)
  }

  ## p < 10^-300 => 10^-300
  study.sample <- correct_extreme_pvalues(study.sample)

  study.sample <- study.sample[order(CHR,POSITION,decreasing = FALSE)]
  ###

  ### calculations
  y.scale.min = floor(-log10(max(study.sample$PVALUE, na.rm = TRUE)))  ## minimum value for Y axis

  # plot to -log10 == 10 even if no pvalue is as low as that
  p.Min = ceiling(-log10(min(study.sample$PVALUE, na.rm = TRUE)))
  p.Min = ifelse(p.Min < 10 , 10 , p.Min)

  if (p.Min %% 2 == 0){
    plot.scale = seq(2, p.Min , 2)
  }else{
    plot.scale = seq(2, p.Min+1 ,2)
  }
  ###

  #----
  # mathematical expression for y axis titles
  log10P <- expression(paste("-log"[10], "(",italic(p),")"))

  study.sample <- fill_missing_chr_positions(study.sample)

  man_plot<-ggplot(study.sample,aes(study.sample$POSITION,-log10(study.sample$PVALUE) ,colour = factor(study.sample$CHR))) +
    scale_y_continuous(log10P, limits=c(y.scale.min,p.Min)) + #, breaks = plot.scale) +
    labs(x="Chromosome",y=" Observed -log(P-value)",
         subtitle=plot.subtitle , title=plot.title) +
    facet_grid(.~study.sample$CHR, scales = "free",space = 'free', switch = "both", margins = TRUE) +
    geom_point() +
   # theme_grey() +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size=12, face = "bold"),
          strip.text.y = element_text(size=14, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          plot.title=element_text(size=30, face="bold",hjust = 0.5),
          plot.subtitle=element_text(size=15, face="bold",hjust = 0.5),
          axis.text.y=element_text(size=15),
          axis.title.x=element_text(size=25),
          axis.title.y=element_text(size=25),
          panel.spacing = unit(0, "points"),
          panel.border = element_rect(fill=NA,color="black", size=0.1, linetype="solid")
          ) +
    geom_hline(data = study.sample, aes(yintercept = sig.threshold.log),color="red")


  ggsave(plot=man_plot,
         device = .QC$graphic.device,
         filename = fileName,
         width = 400,
         height = 200,
         units = 'mm',
         dpi = 200)


  #### write log ####
  print_and_log("Manhattan plot saved!",
                'info')
  #### remove variables from RAM
  rm(study.sample,
     man_plot)

  #invisible(gc())

}

fill_missing_chr_positions <- function(plot.data) {

  # check if map file exists in package folder
  chr_pos_map.file <- system.file("extdata", "chr_pos_map.rds", package = "GWASinspector")

  if(!file.exists(chr_pos_map.file))
    return(plot.data)

  # load the map file
  chr_pos_map <- readRDS(chr_pos_map.file)

  # this column might be missing , if the function is run independently from the package
  if(!is.element('PVALUE.calculated',names(plot.data)))
    plot.data[,PVALUE.calculated := NA]


  for(i in 1:nrow(chr_pos_map))
  {

    plot.data <- rbind(plot.data,data.table::data.table(CHR = chr_pos_map[i]$CHR,
                                                        PVALUE = NA ,
                                                        POSITION = chr_pos_map[i]$`min(POS)` - 2 ,
                                                        PVALUE.calculated = NA))

    plot.data <- rbind(plot.data,data.table::data.table(CHR = chr_pos_map[i]$CHR,
                                                        PVALUE = NA ,
                                                        POSITION = chr_pos_map[i]$`max(POS)` + 2 ,
                                                        PVALUE.calculated = NA))
  }



  return(plot.data)
}
