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
    annotate("text",x=-log10(plot_cutoff_p),y=max(-log10(study.sample$PVALUE))
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


plot.observedP.vs.expectedP_dual<-function(input.data,plot_cutoff_p,PVcor,pvalCorPlotPath,plot.subtitle){


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
    annotate("text",x=-log10(plot_cutoff_p),y=max(-log10(study.sample$PVALUE))
             ,hjust =0.1,
             label=sprintf('italic(r) == %.3f', PVcor),parse= TRUE) +
    boolScale


  plot_HQ <- ggplot(data = NULL) +
    geom_abline(slope=1, intercept=0,colour='grey')+
    {
      if(nrow(study.sample[HQ == TRUE]) > 0)
        geom_point(data = study.sample[HQ == TRUE] ,aes(y=-log10(PVALUE),x=-log10(PVALUE.calculated),color = '1'))
    } +
    labs(title="P-value correlation plot (HQ only)", subtitle=plot.subtitle,
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
    # annotate("text",x=1,y=max(-log10(study.sample[HQ == TRUE]$PVALUE))
    #          ,hjust =0.1,
    #          label=sprintf('italic(r) == %.3f', PVcor),parse= TRUE) +
    boolScale

  ggsave(plot=arrangeGrob(plot,
                          plot_HQ,
                          ncol = 2,
                          #top = textGrob("plot title", gp=gpar(fontsize=12, fontface='bold')) ,
                          padding = unit(1,'cm')),
         device = .QC$graphic.device,
         filename = pvalCorPlotPath,
         units = 'cm',
         width = 24,
         height =10,
         dpi = 600)



  print.and.log("P-value correlation plot saved! ",
                'info')

  #### remove variables from RAM
  rm(study.sample, plot,plot_HQ ,boolScale)
  invisible(gc())
}


plotScatterSmooth.observedP.vs.expectedP<-function(input.data,plot_cutoff_p,PVcor,pvalCorSmPlotPath,plot.subtitle){


  study.sample <- input.data[PVALUE < plot_cutoff_p]

  if(nrow(study.sample) < 1)
  {
    print.and.log('No valid P-value below threshold exist. P-value correlation plot skipped!',
                  'warning',display=.QC$config$debug$verbose)
    return(NULL)
  }

  # mathematical expression for x , y axis titles
  log10Pe <- expression(paste("Expected -log"[10],  "(",italic(p),")"))
  log10Po <- expression(paste("Observed -log"[10],  "(",italic(p),")"))


  if(.QC$img.extension == 'png')
    png(pvalCorSmPlotPath , width = 3300,height = 1500,units = 'px',res=200)
  else
    jpeg(pvalCorSmPlotPath , width = 3300,height = 1500,units = 'px',res=200)

  par(mfrow = c(1, 2),mar=c(5, 6, 6, 2),cex.main=1.2, cex.lab=1.2, font.main=1)


  smoothScatter(x = -log10(study.sample[PVALUE < plot_cutoff_p & HQ == TRUE]$PVALUE),
                y = -log10(study.sample[PVALUE < plot_cutoff_p & HQ == TRUE]$PVALUE.calculated),
                pch=20,
                xlab = log10Po,
                ylab = log10Pe,
                nrpoints = 20000)
  text(sprintf('r = %.3f',PVcor) ,x = -log10(plot_cutoff_p)+0.2,y= -log10(min(study.sample$PVALUE)))
  title("HQ variants", adj = 0.5, line = 1)

  if(.QC$thisStudy$LQ.count > 0)
  {
    smoothScatter(x = -log10(study.sample[PVALUE < plot_cutoff_p & HQ == FALSE]$PVALUE),
                  y = -log10(study.sample[PVALUE < plot_cutoff_p & HQ == FALSE]$PVALUE.calculated),
                  pch=20,
                  xlab = log10Po,
                  ylab = log10Pe,
                  nrpoints = 20000)
    title("LQ variants", adj = 0.5, line = 1)
  }


  title("P-value correlation plot",line=-2, outer = TRUE,cex.main=2)
  title(plot.subtitle,line=-3.5, outer = TRUE,cex.main=1.6)

  dev.off()

  print.and.log("P-value correlation plot saved! ",
                'info')

}
