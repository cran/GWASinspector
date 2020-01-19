plot.variable.frequency.histograms<-function(input.data, study , plot.title.text)
{

  ###FIXME only draw plot for variable with more than a threshold
  # at least one varaible should exist for plotting
  plot.itemcount.threshold <- 1  #calculatePlotThreshold(input.data)
  data.colnames<-colnames(input.data)
  data.row.count <- nrow(input.data)

  ####
  tmp.data<-NULL   #temp value for keeping input.data

  #### declare plots inorder to check if they exist in data ######
  effect.plot <- textGrob('Insufficient data for "EFFECT" variable', gp = gpar(fontsize=12,col='red', fontface='bold'))

  # ggdraw() +
  #   draw_label('Insufficient data for "EFFECT" variable', x = .5, y = .5,
  #              vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')

  maf.plot <- textGrob('Insufficient data for "AF" variable', gp = gpar(fontsize=12,col='red', fontface='bold'))

  # <-ggdraw() +
  #   draw_label('Insufficient data for "AF" variable' , x = .5, y = .5,
  #              vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')

  err.plot <- textGrob('Insufficient data for "STDERR" variable', gp = gpar(fontsize=12,col='red', fontface='bold'))

  # <-ggdraw() +
  #   draw_label('Insufficient data for "STDERR" variable', x = .5, y = .5,
  #              vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')

  call.plot <- textGrob('Insufficient data for "CALLRATE" variable', gp = gpar(fontsize=12,col='red', fontface='bold'))

  # <-ggdraw() +
  #   draw_label('Insufficient data for "CALLRATE" variable', x = .5, y = .5,
  #              vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')

  hwe.plot <- textGrob('Insufficient data for "HWE P-value" variable', gp = gpar(fontsize=12,col='red', fontface='bold'))

  # <-ggdraw() +
  #   draw_label('Insufficient data for "HWE P-value" variable', x = .5, y = .5,
  #              vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')

  impq.plot <- textGrob('Insufficient data for "Imputaion Quality" variable', gp = gpar(fontsize=12,col='red', fontface='bold'))

  # <-ggdraw() +
  #   draw_label('Insufficient data for "Imputaion Quality" variable', x = .5, y = .5,
  #              vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')

  maf.good.items<- data.row.count - length(study$EFF_ALL_FREQ)
  callrate.good.items<- data.row.count - length(study$CALLRATE)
  hwe.good.items<- data.row.count - length(study$HWE_PVAL)
  impq.good.items<- data.row.count - length(study$IMP_QUALITY)

  # these two columns are crucial, so missing values have been removed
  effect.good.items<- data.row.count ##- length(study$column.NA.list$EFFECT)
  err.good.items<- data.row.count ##- length(study$STDERR)


  #####create plots #####

  if('EFFECT' %in% data.colnames && effect.good.items >= plot.itemcount.threshold )
  {
    tmp.data<-data.table(input.data[!is.na(EFFECT)]$EFFECT)
    colnames(tmp.data) <- 'EFFECT'


    # convert BETA to Beta
    # effect.title <- tools::toTitleCase(tolower(.QC$config$input_parameters$effect_type_string))
    effect.title <- chartr('ETA','eta', .QC$config$input_parameters$effect_type_string)


    effect.plot<-ggplot(data=tmp.data,aes(x=EFFECT)) +
      geom_histogram(bins=25,aes(y=..density..), color="darkblue", fill="lightblue") +
      labs(x= effect.title , y="Density") +
      ggtitle("Effect size") +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.text.x = element_text(size=8, face = "bold")
            ,strip.text.y = element_text(size=8, face = "bold")
            ,legend.position="none"
            ,plot.title=element_text(size=10, face="bold",hjust = 0.5))
  }



  if('EFF_ALL_FREQ' %in% data.colnames && maf.good.items >= plot.itemcount.threshold)
  {
    tmp.data<- data.table(input.data[!is.na(EFF_ALL_FREQ)]$EFF_ALL_FREQ)
    colnames(tmp.data) <- 'EFF_ALL_FREQ'

    maf.plot<-ggplot(tmp.data,aes(x=EFF_ALL_FREQ)) +
      geom_histogram(bins=25,aes(y=..density..), color="darkgreen", fill="green") +
      labs(x="Frequency",y="Density") +
      ggtitle("Allele Frequency") +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.text.x = element_text(size=8, face = "bold")
            ,strip.text.y = element_text(size=8, face = "bold")
            ,legend.position="none"
            ,plot.title=element_text(size=10, face="bold",hjust = 0.5))+
      scale_x_continuous(limits=c(min(tmp.data$EFF_ALL_FREQ, na.rm = TRUE),max(tmp.data$EFF_ALL_FREQ, na.rm = TRUE)))
  }



  if('IMP_QUALITY' %in% data.colnames && impq.good.items >= plot.itemcount.threshold)
  {
    tmp.data<-data.table(input.data[!is.na(IMP_QUALITY)]$IMP_QUALITY)
    # tmp.data<-tmp.data[tmp.data$IMPUTED == TRUE,]
    colnames(tmp.data) <- 'IMP_QUALITY'



    if(nrow(tmp.data) > plot.itemcount.threshold)
    {
      impq.plot<-ggplot(tmp.data,aes(x=IMP_QUALITY)) +
        geom_histogram(bins=25,aes(y=..density..), color="black", fill="gray") +
        labs(x="Imputation Quality",
             y="Density",
             # subtitle="Imputed SNPs only",
             title="Imputation Quality") +
        coord_cartesian(xlim = c(0, 1.2)) +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.text.x = element_text(size=8, face = "bold")
              ,strip.text.y = element_text(size=8, face = "bold")
              ,legend.position="none"
              ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
              ,plot.subtitle=element_text(size=8,vjust = 2, hjust=0.5, face="italic", color="darkblue")
        )
    }
  }



  if('CALLRATE' %in% data.colnames && callrate.good.items >= plot.itemcount.threshold)
  {
    tmp.data<-data.table(input.data[!is.na(CALLRATE)]$CALLRATE)
    # tmp.data<-tmp.data[tmp.data$IMPUTED == FALSE,]
    colnames(tmp.data) <- 'CALLRATE'



    if(nrow(tmp.data) > plot.itemcount.threshold)
    {
      call.plot<-ggplot(tmp.data,aes(x=CALLRATE)) +
        geom_histogram(bins=25,aes(y=..density..), color="navy", fill="slategray") +
        labs(x="Call Rate",
             y="Density",
             # subtitle="Genotyped SNPs only",
             title="Call Rate") +
        coord_cartesian(xlim = c(0, 1.2)) +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.text.x = element_text(size=8, face = "bold")
              ,strip.text.y = element_text(size=8, face = "bold")
              ,legend.position="none"
              ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
              ,plot.subtitle=element_text(size=8,vjust = 2, hjust=0.5, face="italic", color="darkblue")
        )
      # scale_x_continuous(limits=c(min(tmp.data$CALLRATE),max(tmp.data$CALLRATE)))
    }
  }



  if('HWE_PVAL' %in% data.colnames && hwe.good.items >= plot.itemcount.threshold)
  {
    tmp.data<-data.table(input.data[!is.na(HWE_PVAL)]$HWE_PVAL)
    # tmp.data<-tmp.data[tmp.data$IMPUTED == FALSE,]
    colnames(tmp.data) <- 'HWE_PVAL'




    if(nrow(tmp.data) > plot.itemcount.threshold)
    {
      hwe.plot<-ggplot(tmp.data,aes(x=HWE_PVAL)) +
        geom_histogram(bins=25,aes(y=..density..), color="yellow", fill="yellow4") +
        labs(x="HWE P-value",
             y="Density",
             # subtitle="Genotyped SNPs only",
             title="HWE P-value") +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.text.x = element_text(size=8, face = "bold")
              ,strip.text.y = element_text(size=8, face = "bold")
              ,legend.position="none"
              ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
              ,plot.subtitle=element_text(size=8,vjust = 2, hjust=0.5, face="italic", color="darkblue")
        )
    }
  }



  if('STDERR' %in% data.colnames && err.good.items >= plot.itemcount.threshold)
  {
    tmp.data<-data.table(input.data[!is.na(STDERR)]$STDERR)
    colnames(tmp.data) <- 'STDERR'



    err.plot<-ggplot(tmp.data,aes(x=STDERR)) +
      geom_histogram(bins=25,aes(y=..density..), color="darkred", fill="brown") +
      labs(x="Standard Error",y="Density") +
      ggtitle("Standard Error") +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.text.x = element_text(size=8, face = "bold")
            ,strip.text.y = element_text(size=8, face = "bold")
            ,legend.position="none"
            ,plot.title=element_text(size=10, face="bold",hjust = 0.5))
  }




  ggsave(filename   = study$histPlotPath,
		 device = .QC$graphic.device,
         plot = arrangeGrob(effect.plot,
                            maf.plot,
                            call.plot,
                            err.plot,
                            hwe.plot,
                            impq.plot,ncol=3,
                            top = textGrob(plot.title.text, gp=gpar(fontsize=16, fontface='bold')) ,
                            padding = unit(1,'cm')),
         width = 35,
         height = 16,
         units = 'cm',
         dpi = 600)


  #dev.off()

  ### write log #####
  print.and.log("histogram plot is saved!",'info')

  #### remove plots from RAM ####
  rm(err.plot,
     hwe.plot,
     call.plot,
     impq.plot,
     maf.plot,
     effect.plot,
     plot.title,
     final.plot,
     tmp.data)

  invisible(gc())

}
