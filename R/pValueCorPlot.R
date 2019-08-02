plot.observedP.vs.expectedP<-function(input.data,plot_cutoff_p,PVcor,pvalCorPlotPath,plot.subtitle){


  study.sample <- input.data[PVALUE < plot_cutoff_p]

  if(nrow(study.sample) < 1)
  {
    print.and.log('No valid P-value below threshold exist. P-value correlation plot skipped!',
                  'warning',display=.QC$config$debug$verbose)
    return(NULL)
  }

  ## p < 10^-300 => 10^-300
  study.sample <- correct.extreme.pvalues(study.sample)
  study.sample <- correct.extreme.calculated.pvalues(study.sample)



  boolScale <- scale_colour_manual(values=c("0"="darkred","1"="darkblue"),
                                   name="Variants",
                                   labels=c("0"="LQ","1"= "HQ"))

  #----
  # mathematical expression for x , y axis titles
  log10Pe <- expression(paste("Expected -log"[10],  "(",italic(p),")"))
  log10Po <- expression(paste("Observed -log"[10],  "(",italic(p),")"))

  plot <- ggplot(data = NULL) +
    geom_abline(slope=1, intercept=0,colour='grey')+
    {
      if(nrow(study.sample[HQ == FALSE]) > 0)
        geom_point(data = study.sample[HQ == FALSE] ,aes(y=-log10(PVALUE),x=-log10(PVALUE.calculated),color = '0'))
    } +{
      if(nrow(study.sample[HQ == TRUE]) > 0)
        geom_point(data = study.sample[HQ == TRUE] ,aes(y=-log10(PVALUE),x=-log10(PVALUE.calculated),color = '1'))
    } +
    labs(title="P-value correlation plot", subtitle=plot.subtitle,
         x=log10Pe,
         y=log10Po) +
    theme_classic() +
    theme(strip.background = element_blank()
          ,axis.title.y=element_text(size=10)
          ,axis.title.x=element_text(size=10)
          ,legend.position = "top",
          legend.justification='left',
          legend.direction='horizontal'
          ,legend.text = element_text(colour="black", size=8)
          ,legend.key.size = unit(.2, "cm")
          ,legend.key.width = unit(.5, "cm"),
          legend.key.height=unit(0.5, "cm"),
          legend.title = element_text(colour="black", size=10,face="bold")
          ,legend.background = element_rect(fill="gray90", size=0.2, linetype="solid",
                                            colour ="darkblue")
          ,plot.title=element_text(size=14, face="bold",hjust = 0.5)
          ,plot.subtitle=element_text(size=10, face="bold",hjust = 0.5)
          )+
    annotate("text",x=1,y=max(-log10(input.data$PVALUE))
             ,hjust =0.1,
             label=sprintf('italic(r) == %.3f', PVcor),parse= TRUE) +
    boolScale


  ggsave(plot=plot,
		 device = .QC$graphic.device,
         filename = pvalCorPlotPath,
         units = c('mm'),
         width = 180,
         height =120,
         dpi = 300)




  print.and.log("P-value correlation plot saved! ",
                'info')

  #### remove variables from RAM
  rm(study.sample, plot ,boolScale)
  invisible(gc())
}
