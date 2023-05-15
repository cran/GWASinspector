QQ_plots<-function(input.data,plot_cutoff_p = 0.01,plot.title.text)
{


  #
  # 1- missing P-values in the file are replaced by calculated P-values (STDERR, EFFECT calculation)
  # 2- extreme P-values are corrected
  # 3- dataset is sorted by P-values (from low to high values)
  # 4- -log10(pvalue) is added to dataset
  # 5- -lo10(expected pvalue) is added to dataset. i.e. -log10(ppoints(rowCount))
  # 6- upper/lower limits (confidence interval lines in plot) are calculated.
  #
  # Filtering:
  #   7- variants with NA values for each filter are removed from dataset. e.g. variants with missing HWE-P are first removed when plotting HWE-P qqplot.
  # 8- variants are separated to groups based on a filter value. e.g.
  #       all variants to group 1
  #       HWE-P >= 10^-6 to group 2
  #       HWE-P >= 10^-4 to group 3
  # 9- expected pvalues are again calculated in each of the subgroups (2,3)
  # 10 - percentage of each subgroup to total variant count in step 7 are calculated (NA values are not included in count).
  # 11 - variants with pvalues that are higher than a thresold are removed. e.g. -log10(pvalue) > 2
  # 12 - plot is created from each group.




  #replace missing pvalues with calculated pvalues
  input.data[is.na(PVALUE) , PVALUE := PVALUE.calculated]

  ## p < 10^-300 => 10^-300
  input.data <- correct_extreme_pvalues(input.data)

  pvalue.count <- nrow(input.data)

  # sort the dataset on pvalue
  input.data <- input.data[order(PVALUE,decreasing=FALSE)]
  # calculate log of P value
  input.data[,logP := -log10(PVALUE)]
  # calculate log of expected P value
  input.data[, expected.logP := -log10(ppoints(pvalue.count))]



  # calculate limit for confidence interval
  # used from https://gist.github.com/slowkow/9041570
  # ci.data <- data.table('cupper' = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:pvalue.count, shape2 = pvalue.count:1)),
  #                       'clower' = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:pvalue.count, shape2 = pvalue.count:1)))
  #
  input.data[, cupper := -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:pvalue.count, shape2 = pvalue.count:1))]
  input.data[, clower := -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:pvalue.count, shape2 = pvalue.count:1))]


  # pValue threshold for selecting variants (from config file , default= 0.01)
  plotting.threshold <- -log10(plot_cutoff_p)


  qq.af.plot <- plot_qqplot_af(input.data,plotting.threshold)


  ## set x and y labels for consistency between plots
  if(!is.null(qq.af.plot$name)) # empty plot (rectgrob) and not ggplot
  {
    xlab <- seq(0,1,0.25)
    ylab <- seq(0,1,0.25)

  }else{

    # this works ggplot2 v 2.2
    xlab <- as.numeric(ggplot_build(qq.af.plot)$layout$panel_ranges[[1]]$x.labels)
    ylab <- as.numeric(ggplot_build(qq.af.plot)$layout$panel_ranges[[1]]$y.labels)

    if(length(xlab) == 0)
    {
      # this workd gor ggplot2 v 3
      xlab <- as.numeric(ggplot_build(qq.af.plot)$layout$panel_params[[1]]$x.major_source)
      ylab <- as.numeric(ggplot_build(qq.af.plot)$layout$panel_params[[1]]$y.major_source)
    }

  }

  # xlim <- ggplot_build(qq.af.plot)$layout$panel_ranges[[1]]$x.range
  # ylim <- ggplot_build(qq.af.plot)$layout$panel_ranges[[1]]$y.range
  ##


  qq.hwe.plot <-  plot_qqplot_hwe(input.data,plotting.threshold ,xlab,ylab)
  qq.call.plot <- plot_qqplot_call(input.data,plotting.threshold,xlab,ylab)
  qq.impq.plot <- plot_qqplot_impq(input.data,plotting.threshold,xlab,ylab)




  ggsave(filename = .QC$thisStudy$QQPlotPath,
		 device = .QC$graphic.device,
         plot = arrangeGrob(qq.af.plot,
                            qq.hwe.plot,
                            qq.call.plot,
                            qq.impq.plot,
                            ncol = 2,
                            top = textGrob(plot.title.text, gp=gpar(fontsize=16, fontface='bold')) ,
                            padding = unit(1,'cm')),
         width = 25,
         height = 25,
         units = 'cm',
         dpi = 600)



  ### write log #####
  print_and_log("QQ plot is saved!",'info')
  #### remove plots from RAM ####
  rm(qq.af.plot,
     qq.hwe.plot,
     qq.call.plot,
     qq.impq.plot
  )

  invisible(gc())
}



plot_qqplot_af <- function(input.data,plotting.threshold) {


  #create an emoty plot , return if column not found in data or no valid variant
  plot <- textGrob('Insufficient data for "AF" variable QQ-plot',gp = gpar(fontsize=12,col='red', fontface='bold'))

  # <-ggdraw() +
  #    draw_label('Insufficient data for "allele frequency" variable', x = .5, y = .5,
  #               vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')



  if('EFF_ALL_FREQ' %notin% colnames(input.data))
    return(plot)


  input.data <-subset(input.data[!is.na(EFF_ALL_FREQ)],select=c('logP','expected.logP','EFF_ALL_FREQ','cupper','clower'))


  input.data[,frq := 1]

  input.data[EFF_ALL_FREQ < 0.05 | EFF_ALL_FREQ > 0.95, frq := 2]
  input.data[EFF_ALL_FREQ < 0.03 | EFF_ALL_FREQ > 0.97, frq := 3]
  input.data[EFF_ALL_FREQ < 0.01 | EFF_ALL_FREQ > 0.99 , frq := 4]


  input.data[frq < 2 , expected.logP2 := -log10(ppoints(nrow(input.data[frq < 2])))]
  input.data[frq < 3 , expected.logP3 := -log10(ppoints(nrow(input.data[frq < 3])))]
  input.data[frq < 4 , expected.logP4 := -log10(ppoints(nrow(input.data[frq < 4])))]


  boolScale <- scale_colour_manual(values=c("1"="black","2"="orange","3"="green","4"="yellow"),
                                   name="Variants",
                                   labels=c('1'= 'All',
                                            '2'=sprintf('> 0.01 (%s)',calculatePercent(nrow(input.data[ frq < 4 ]),
                                                                                       nrow(input.data))),
                                            '3'=sprintf('> 0.03 (%s)',calculatePercent(nrow(input.data[ frq < 3 ]),
                                                                                       nrow(input.data))),
                                            '4'= sprintf('> 0.05 (%s)',calculatePercent(nrow(input.data[ frq < 2 ]),
                                                                                        nrow(input.data)))
                                   )
  )


  input.data <- input.data[logP > plotting.threshold]
  ## at least one valid variant should exist
  if(nrow(input.data) < 1)
    return(plot)


  #----
  # mathematical expression for x , y axis titles
  log10Pe <- expression(paste("Expected -log"[10],  "(",italic(p),")"))
  log10Po <- expression(paste("Observed -log"[10],  "(",italic(p),")"))

  .QC$QQ.min <- min(input.data$logP, na.rm = TRUE)
  .QC$QQ.max <- max(input.data$logP, na.rm = TRUE)

  .QC$QQ_exp.min <- min(input.data$expected.logP, na.rm = TRUE)
  .QC$QQ_exp.max <- max(input.data$expected.logP, na.rm = TRUE)


  plot <- ggplot(data=input.data) +
    geom_line(aes(expected.logP, cupper), linetype = 2) +
    geom_line(aes(expected.logP, clower), linetype = 2) +
    geom_point(data=input.data, aes(x=expected.logP,y=logP,color='1'),shape=1)+
    geom_point(data=input.data[ frq < 4 ], aes(x=expected.logP4,y=logP,color='2'),shape=1)+
    geom_point(data=input.data[ frq < 3 ], aes(x=expected.logP3 ,y=logP,color='3'),shape=1)+
    geom_point(data=input.data[ frq < 2 ], aes(x=expected.logP2,y=logP,color='4'),shape=1)+
    geom_abline(slope=1, intercept=0,colour='grey') +
    labs(title="QQ plot - allele frequency",
         x=log10Pe,
         y=log10Po) +
    theme_classic() +
    theme(axis.title.y=element_text(size=10)
          ,axis.title.x=element_text(size=10)
          ,axis.text.x=element_text(size=8)
          ,axis.text.y=element_text(size=8),legend.position="top"
          ,legend.text = element_text(colour="black", size=8)
          ,legend.key.size = unit(.2, "cm")
          ,legend.key.width = unit(.5, "cm")
          ,legend.title = element_text(colour="black", size=10,face="bold")
          ,legend.background = element_rect(fill="gray90")
          ,plot.title=element_text(size=12, face="bold",hjust = 0.5)
          ,plot.margin=unit(c(1,1,1,1),"cm"))+
    annotate("text",
             x = plotting.threshold + 0.5,
             y = max(input.data$logP),
             label=sprintf('lambda == %s', .QC$thisStudy$lambda),
             parse= TRUE)+
    boolScale + coord_cartesian(xlim = c(.QC$QQ_exp.min,.QC$QQ_exp.max), ylim = c(.QC$QQ.min,.QC$QQ.max))

  return(plot)

}

plot_qqplot_hwe <- function(input.data,plotting.threshold,xlab,ylab) {

  #create an emoty plot , return if column not found in data or no valid variant
  plot <- textGrob('Insufficient data for "HWE P-value" variable QQ-plot',gp = gpar(fontsize=12,col='red', fontface='bold'))


  if('HWE_PVAL' %notin%  colnames(input.data))
    return(plot)



  input.data <-subset(input.data[!is.na(HWE_PVAL)],select=c('logP','expected.logP','HWE_PVAL','cupper','clower'))
  input.data <-subset(input.data,select=c('logP','expected.logP','HWE_PVAL','cupper','clower'))


  input.data[,hwp := 1]
  input.data[HWE_PVAL >= 10^-6, hwp := 2 ]
  input.data[HWE_PVAL >= 10^-4, hwp := 3 ]

  input.data[hwp > 1 , expected.logP2 := -log10(ppoints(nrow(input.data[hwp > 1])))]
  input.data[hwp > 2 , expected.logP3 := -log10(ppoints(nrow(input.data[hwp > 2])))]



  boolScale <- scale_colour_manual(values=c("1"="black","2"="blue","3"="orange"),
                                   name="Variants",
                                   breaks=c('1','2','3'),
                                   labels=c('1'= 'All',
                                            '2'= sprintf('> 1e-06 (%s)',calculatePercent(nrow(input.data[hwp > 1]),
                                                                                         nrow(input.data))),
                                            '3'=sprintf('> 1e-04 (%s)',calculatePercent(nrow(input.data[hwp >2]),
                                                                                        nrow(input.data)))))

  ## only variants with low p values are selected
  input.data <- input.data[logP > plotting.threshold]

  ## at least one valid variant should exist
  if(nrow(input.data) < 1)
    return(plot)


  #----
  # mathematical expression for x , y axis titles
  log10Pe <- expression(paste("Expected -log"[10],  "(",italic(p),")"))
  log10Po <- expression(paste("Observed -log"[10],  "(",italic(p),")"))


  plot <- ggplot(data=input.data) +
    geom_line(aes(expected.logP, cupper), linetype = 2) +
    geom_line(aes(expected.logP, clower), linetype = 2) +
    geom_point(data = input.data, aes(x= expected.logP, y =logP , color = '1'),shape=1) +
    {
      if(nrow(input.data[hwp >1]) > 0)
        geom_point(data = input.data[hwp >1], aes(x= expected.logP2, y =logP , color = '2'),shape=1)
    }+{
      if(nrow(input.data[hwp >2]) > 0)
        geom_point(data = input.data[hwp >2], aes(x= expected.logP3, y =logP , color = '3'),shape=1)
    }+
    geom_abline(slope=1, intercept=0,colour='grey') +
    labs(title="QQ plot - HWE P-value",
         x=log10Pe,
         y=log10Po) +
    theme_classic() +
    theme(axis.title.y=element_text(size=10)
          ,axis.title.x=element_text(size=10)
          ,axis.text.x=element_text(size=8)
          ,axis.text.y=element_text(size=8),legend.position="top"
          ,legend.text = element_text(colour="black", size=8)
          ,legend.key.size = unit(.2, "cm")
          ,legend.key.width = unit(.5, "cm")
          ,legend.title = element_text(colour="black", size=10,face="bold")
          ,legend.background = element_rect(fill="gray90")
          ,plot.title=element_text(size=12, face="bold",hjust = 0.5)
          ,plot.margin=unit(c(1,2,1,1),"cm")
          ,axis.line.x = element_line(color="black", size = 0.5)
          ,axis.line.y = element_line(color="black", size = 0.5)) +
    boolScale +
    coord_cartesian (xlim = c(.QC$QQ_exp.min,.QC$QQ_exp.max), ylim = c(.QC$QQ.min,.QC$QQ.max)) +
    { if(length(xlab) > 0)
      scale_x_continuous(breaks = xlab, labels = xlab)
    }+
    { if(length(ylab) > 0)
      scale_y_continuous(breaks = ylab, labels = ylab )
    }


  return(plot)

}

plot_qqplot_call <- function(input.data,plotting.threshold,xlab,ylab) {

  #create an emoty plot , return if column not found in data or no valid variant
  plot <- textGrob('Insufficient data for "call rate" variable QQ-plot',gp = gpar(fontsize=12,col='red', fontface='bold'))


  if('CALLRATE' %notin%  colnames(input.data))
    return(plot)


  input.data <-subset(input.data[!is.na(CALLRATE)],select=c('logP','expected.logP','CALLRATE','cupper','clower'))

  input.data[, call := 1]
  input.data[CALLRATE >= 0.95, call := 2 ]
  input.data[CALLRATE >= 0.98, call := 3 ]
  input.data[CALLRATE >= 0.99, call := 4 ]

  input.data[call >1 , expected.logP2 := -log10(ppoints(nrow(input.data[call >1])))]
  input.data[call >2 , expected.logP3 := -log10(ppoints(nrow(input.data[call >2])))]
  input.data[call >3 , expected.logP4 := -log10(ppoints(nrow(input.data[call >3])))]


  boolScale <- scale_colour_manual(values=c("1"="black","2"="blue","3"="green","4"="red"),
                                   name="Variants",
                                   breaks=c('1','2','3','4'),
                                   labels=c('1'='All',
                                            '2'=sprintf('> 0.95 (%s)',calculatePercent(nrow(input.data[call > 1]),
                                                                                       nrow(input.data))),
                                            '3'=sprintf('> 0.98 (%s)',calculatePercent(nrow(input.data[call >2]),
                                                                                       nrow(input.data))),
                                            '4'=sprintf('> 0.99 (%s)',calculatePercent(nrow(input.data[call >3]),
                                                                                       nrow(input.data)))))

  input.data <- input.data[logP > plotting.threshold]
  ## at least one valid variant should exist
  if(nrow(input.data) < 1)
    return(plot)


  #----
  # mathematical expression for x , y axis titles
  log10Pe <- expression(paste("Expected -log"[10],  "(",italic(p),")"))
  log10Po <- expression(paste("Observed -log"[10],  "(",italic(p),")"))


  plot <- ggplot(data=input.data) +
    geom_line(aes(expected.logP, cupper), linetype = 2) +
    geom_line(aes(expected.logP, clower), linetype = 2) +
    geom_point(data = input.data, aes(x= expected.logP, y =logP , color = '1'),shape=1) +
    {
      if(nrow(input.data[call > 1]) > 0)
        geom_point(data = input.data[call >1], aes(x= expected.logP2, y =logP , color = '2'),shape=1)
    } +
    {
      if(nrow(input.data[call >2]) > 0)
        geom_point(data = input.data[call >2], aes(x= expected.logP3, y =logP , color = '3'),shape=1)
    } +
    {
      if(nrow(input.data[call >3]) > 0)
        geom_point(data = input.data[call >3], aes(x= expected.logP4, y =logP , color = '4'),shape=1)
    }+
    geom_abline(slope=1, intercept=0,colour='grey') +
    labs(title="QQ plot - call rates",
         x=log10Pe,
         y=log10Po) +
    theme_classic() +
    theme(axis.title.y=element_text(size=10)
          ,axis.title.x=element_text(size=10)
          ,axis.text.x=element_text(size=8)
          ,axis.text.y=element_text(size=8),legend.position="top"
          ,legend.text = element_text(colour="black", size=8)
          ,legend.key.size = unit(.2, "cm")
          ,legend.key.width = unit(.5, "cm")
          ,legend.title = element_text(colour="black", size=10,face="bold")
          ,legend.background = element_rect(fill="gray90")
          ,plot.title=element_text(size=12, face="bold",hjust = 0.5)
          ,plot.margin=unit(c(1,1,1,1),"cm")
          ,axis.line = element_line(color="black", size = 0.5)) +
    boolScale +
    coord_cartesian (xlim = c(.QC$QQ_exp.min,.QC$QQ_exp.max), ylim = c(.QC$QQ.min,.QC$QQ.max)) +
    { if(length(xlab) > 0)
      scale_x_continuous(breaks = xlab, labels = xlab)
    }+
    { if(length(ylab) > 0)
      scale_y_continuous(breaks = ylab, labels = ylab )
    }


  return(plot)

}

plot_qqplot_impq <- function(input.data,plotting.threshold,xlab,ylab) {

  #create an emoty plot , return if column not found in data or no valid variant
  plot <- textGrob('Insufficient data for "imputation quality QQ-plot', gp = gpar(fontsize=12,col='red', fontface='bold'))

  if('IMP_QUALITY' %notin%  colnames(input.data))
    return(plot)

  input.data <-subset(input.data[!is.na(IMP_QUALITY)],select=c('logP','expected.logP','IMP_QUALITY','cupper','clower'))

  input.data[, impq := 1]
  input.data[IMP_QUALITY >= 0.3, impq := 2 ]
  input.data[IMP_QUALITY >= 0.5, impq := 3 ]
  input.data[IMP_QUALITY >= 0.7, impq := 4 ]
  input.data[IMP_QUALITY >= 0.9, impq := 5 ]

  input.data[impq >1 , expected.logP2 := -log10(ppoints(nrow(input.data[impq >1])))]
  input.data[impq >2 , expected.logP3 := -log10(ppoints(nrow(input.data[impq >2])))]
  input.data[impq >3 , expected.logP4 := -log10(ppoints(nrow(input.data[impq >3])))]
  input.data[impq >4 , expected.logP5 := -log10(ppoints(nrow(input.data[impq >4])))]


  boolScale <- scale_colour_manual(values=c("1"="black","2"="blue","3"="orange","4"="green","5"="yellow"),
                                   name="Variants",
                                   breaks=c('1','2','3','4','5'),
                                   labels=c('1'='All',
                                            '2'=sprintf('> 0.3 (%s)',calculatePercent(nrow(input.data[impq >1]),
                                                                                      nrow(input.data))),
                                            '3'=sprintf('> 0.5 (%s)',calculatePercent(nrow(input.data[impq >2]),
                                                                                      nrow(input.data))),
                                            '4'=sprintf('> 0.7 (%s)',calculatePercent(nrow(input.data[impq >3]),
                                                                                      nrow(input.data))),
                                            '5'=sprintf('> 0.9 (%s)',calculatePercent(nrow(input.data[impq >4]),
                                                                                      nrow(input.data)))))

  input.data <- input.data[logP > plotting.threshold]
  ## at least one valid variant should exist
  if(nrow(input.data) < 1)
    return(plot)


  #----
  # mathematical expression for x , y axis titles
  log10Pe <- expression(paste("Expected -log"[10],  "(",italic(p),")"))
  log10Po <- expression(paste("Observed -log"[10],  "(",italic(p),")"))


  plot <- ggplot(data=input.data) +
    geom_line(aes(expected.logP, cupper), linetype = 2) +
    geom_line(aes(expected.logP, clower), linetype = 2) +
    geom_abline(slope=1, intercept=0,colour='grey') +
    geom_point(data = input.data, aes(x= expected.logP, y =logP , color = '1'),shape=1) +
    {
      if(nrow(input.data[impq >1]) > 0)
        geom_point(data = input.data[impq >1], aes(x= expected.logP2, y =logP , color = '2'),shape=1)
    }+{
      if(nrow(input.data[impq >2]) > 0)
        geom_point(data = input.data[impq >2], aes(x= expected.logP3, y =logP , color = '3'),shape=1)
    }+{
      if(nrow(input.data[impq >3]) > 0)
        geom_point(data = input.data[impq >3], aes(x=expected.logP4, y =logP , color = '4'),shape=1)
    }+{
      if(nrow(input.data[impq >4]) > 0)
        geom_point(data = input.data[impq >4], aes(x= expected.logP5, y =logP , color = '5'),shape=1)
    }+
    labs(title="QQ plot - imputation quality",
         x=log10Pe,
         y=log10Po) +
    theme_classic() +
    theme(axis.title.y=element_text(size=10)
          ,axis.title.x=element_text(size=10)
          ,axis.text.x=element_text(size=8)
          ,axis.text.y=element_text(size=8),legend.position="top"
          ,legend.text = element_text(colour="black", size=8)
          ,legend.key.size = unit(.2, "cm")
          ,legend.key.width = unit(.5, "cm")
          ,legend.title = element_text(colour="black", size=10,face="bold")
          ,legend.background = element_rect(fill="gray90")
          ,plot.title=element_text(size=12, face="bold",hjust = 0.5)
          ,plot.margin=unit(c(1,2,1,1),"cm")
          ,axis.line.x = element_line(color="black", size = 0.5)) +
    boolScale +
    coord_cartesian (xlim = c(.QC$QQ_exp.min,.QC$QQ_exp.max), ylim = c(.QC$QQ.min,.QC$QQ.max)) +
    { if(length(xlab) > 0)
      scale_x_continuous(breaks = xlab, labels = xlab)
    }+
    { if(length(ylab) > 0)
      scale_y_continuous(breaks = ylab, labels = ylab )
    }

  return(plot)

}
