plot.DataMAF.vs.RefMAF<-function(input.data,mafPlotPath,AFcor,AFcor.palindromic,AFcor.INDEL,plot.title.text)
{
  if(!is.numeric(AFcor) || !is.numeric(AFcor.palindromic)){
    print.and.log('Allele frequency plot is skipped!','warning',display=.QC$config$debug$verbose)
    return(NULL)
  }

  # split data to snp and indel
  snp.input.data <- NULL
  indel.input.data <- NULL
  indel.ggplot <- NULL

  split.input.data <- split(input.data,input.data$VT)

  snp.input.data <- split.input.data[[1]] #as.data.table(input.data[VT == 1 , ])

  if(length(split.input.data) == 2)
    indel.input.data <- split.input.data[[2]] # as.data.table(input.data[VT == 2 , ])


  # manual scale for HQ / LQ variants -- used for both plots
  boolScale_palindromic <- scale_colour_manual(name="Variants",
                                               values=c('0' = "darkred", '1' = "darkblue"),
                                               labels=c( "0"="LQ","1"= "HQ"))


  AF_dif_threshold <- ifelse(AFcor > 0.97 , 0.05 , 0.1)
  ######## AMbiguous alleles ####
  ## ============================

  #### select a sample of 50000 aalleles with maf-maf < 0.1
  # because plotting millions of data will take a long time
  # there is no need to plot alleles with similar AF

  # check if user wants all variants to be plotted
  # this is a bad idea and should not be used normally
  if(.QC$config$debug$reduced.AF.plot)
  {
    similar.palindromic.variants <- which(snp.input.data$palindromic == TRUE &
                                            abs(snp.input.data$EFF_ALL_FREQ - snp.input.data$AF) < AF_dif_threshold)
  }
  else
  {
    similar.palindromic.variants <- data.table()
  }

  ## if similar items are above 100,000 only 100,000 are selected randomly, if  not all ambiguous alleles are plotted
  if(length(similar.palindromic.variants) > 100000)
  {
    amb.alleles <- rbind(snp.input.data[palindromic == TRUE & abs(EFF_ALL_FREQ - AF) > AF_dif_threshold],
                         snp.input.data[sample(similar.palindromic.variants,100000)])
  }else
  {
    amb.alleles <- snp.input.data[palindromic == TRUE] 
  }


  # palindromic plot
  amb.ggplot<-ggplot(data = NULL) +{
    if(nrow(amb.alleles[HQ == FALSE]) > 0)
      geom_point(data = amb.alleles[HQ == FALSE] ,aes(x = EFF_ALL_FREQ,y = AF,colour = '0'), size = .8)
  }+{
    if(nrow(amb.alleles[HQ == TRUE]) > 0)
      geom_point(data = amb.alleles[HQ == TRUE] , aes(x = EFF_ALL_FREQ,y = AF,colour = '1'), size = .8)
  }+
    labs(title="Allele Frequency Correlation",
         x="Reported allele frequency",
         y="Reference allele frequency",
         subtitle="Palindromic SNPs") +
    # theme_classic() +
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size=8, face = "bold")
          ,strip.text.y = element_text(size=8, face = "bold")
          ,legend.position = "top",
          legend.justification='left',
          legend.direction='horizontal'
          ,legend.text = element_text(colour="black", size=8)
          ,legend.key.size = unit(.2, "cm")
          ,legend.key.width = unit(.5, "cm"),
          legend.title = element_text(colour="black", size=10,face="bold")
          ,legend.background = element_rect(fill="gray90", size=0.2, linetype="solid",
                                            colour ="darkblue")
          ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
          ,plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="darkblue"))+
    annotate("text",
             x=0.1,
             y=1,
             label=sprintf('italic(r) == %.3f' ,AFcor.palindromic),
             parse= TRUE) +
    boolScale_palindromic



  ############### All alleles #######

  # manual scale for HQ / LQ variants -- used for both plots
  boolScale<- scale_colour_manual(name="Variants",
                                  values=c('0' = "darkred", '1' = "darkblue" ),
                                  labels=c( "0"="LQ","1"= "HQ" ))



  # check if user wants all variants to be plotted
  # this is a bad idea and should not be used normally
  if(.QC$config$debug$reduced.AF.plot)
  {
    similar.all.variants <- which(abs(snp.input.data$EFF_ALL_FREQ - snp.input.data$AF) < AF_dif_threshold)
  }
  else
  {
    similar.all.variants <- data.table()
  }


  ## if similar items are above 100,000 only 100,000 are selected randomly, if not, all alleles are plotted
  if(length(similar.all.variants) > 100000)
  {
    all.alleles <- rbind(snp.input.data[abs(EFF_ALL_FREQ - AF) > AF_dif_threshold],
                         snp.input.data[sample(similar.all.variants,100000)])
  }else
  {
    all.alleles <- snp.input.data
  }



  #all varaints plot
  all.ggplot<-ggplot(data = NULL) +{
    if(nrow(all.alleles[HQ == FALSE]) > 0)
      geom_point(data = all.alleles[HQ == FALSE] ,aes(x = EFF_ALL_FREQ,y = AF,colour = '0'), size = .8)
  }+{
    if(nrow(all.alleles[HQ == TRUE]) > 0)
      geom_point(data = all.alleles[HQ == TRUE] , aes(x = EFF_ALL_FREQ,y = AF,colour = '1'), size = .8)
  }+
    labs(title="Allele Frequency Correlation",
         x="Reported allele frequency",
         y="Reference allele frequency",
         subtitle="All SNPs") +
    # theme_classic() +
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size=6, face = "bold")
          ,strip.text.y = element_text(size=6, face = "bold")
          ,legend.position = "top",
          legend.justification='left',
          legend.direction='horizontal'
          ,legend.text = element_text(colour="black", size=8)
          ,legend.key.size = unit(.2, "cm")
          ,legend.key.width = unit(.5, "cm")
          ,legend.title = element_text(colour="black", size=10,face="bold")
          ,legend.background = element_rect(fill="gray90", size=0.2, linetype="solid",
                                            colour ="darkblue")
          ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
          ,plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="darkblue")) +
    annotate("text",
             x=0.1,
             y=1,
             label=sprintf('italic(r) == %.3f' ,AFcor),
             parse= TRUE) +
    boolScale


  #### 3 INDEL plot ####
  if(!is.null(indel.input.data))
  {

    # manual scale for HQ / LQ variants -- used for both plots
    boolScale<- scale_colour_manual(name="Variants",
                                    values=c('0' = "darkred", '1' = "darkblue" ),
                                    labels=c( "0" = "LQ", "1" = "HQ"))


    # check if user wants all variants to be plotted
    # this is a bad idea and should not be used normally
    if(.QC$config$debug$reduced.AF.plot)
    {
      similar.indel.variants <- which(abs(indel.input.data$EFF_ALL_FREQ - indel.input.data$AF) < AF_dif_threshold)
    }
    else
    {
      similar.indel.variants <- data.table()
    }



    ## if similar items are above 100,000 only 100,000 are selected randomly, if not, all alleles are plotted
    if(length(similar.indel.variants) > 100000)
    {
      indel.alleles <- rbind(indel.input.data[abs(EFF_ALL_FREQ - AF) > AF_dif_threshold],
                             indel.input.data[sample(similar.all.variants,100000)])
    }else
    {
      indel.alleles <- indel.input.data
    }


    #INDEL varaints plot
    indel.ggplot <- ggplot(data = NULL) +{
      if(nrow(indel.alleles[HQ == FALSE]) > 0)
        geom_point(data = indel.alleles[HQ == FALSE] ,aes(x = EFF_ALL_FREQ,y = AF,colour = '0'), size = .8)
    }+{
      if(nrow(indel.alleles[HQ == TRUE]) > 0)
        geom_point(data = indel.alleles[HQ == TRUE] , aes(x = EFF_ALL_FREQ,y = AF,colour = '1'), size = .8)
    }+
      labs(title="Allele Frequency Correlation",
           x="Reported allele frequency",
           y="Reference allele frequency",
           subtitle="All non-SNPs") +
      # theme_classic() +
      scale_x_continuous(limits=c(0,1)) +
      scale_y_continuous(limits=c(0,1)) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.text.x = element_text(size=6, face = "bold")
            ,strip.text.y = element_text(size=6, face = "bold")
            ,legend.position = "top",
            legend.justification='left',
            legend.direction='horizontal'
            ,legend.text = element_text(colour="black", size=8)
            ,legend.key.size = unit(.2, "cm")
            ,legend.key.width = unit(.5, "cm")
            ,legend.title = element_text(colour="black", size=10,face="bold")
            ,legend.background = element_rect(fill="gray90", size=0.2, linetype="solid",
                                              colour ="darkblue")
            ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
            ,plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="darkblue")) +
      annotate("text",
               x=0.1,
               y=1,
               label=sprintf('italic(r) == %.3f' ,AFcor.INDEL),
               parse= TRUE) +
      boolScale


    rm(similar.indel.variants,
       indel.alleles)

  }


  if(is.null(indel.ggplot)) # 2 subplots with no indel
    ggsave(filename = mafPlotPath,
           device = .QC$graphic.device,
           plot = arrangeGrob(all.ggplot,
                              amb.ggplot,
                              ncol = 2,
                              top = textGrob(plot.title.text, gp=gpar(fontsize=12, fontface='bold')) ,
                              padding = unit(1,'cm')),
           width = 24,
           height = 10,
           units = 'cm',
           dpi = 600)
  else                      # 3 sub plots with indel
    ggsave(filename = mafPlotPath,
           device = .QC$graphic.device,
           plot = arrangeGrob(all.ggplot,
                              amb.ggplot,
                              indel.ggplot,
                              ncol = 3,
                              top = textGrob(plot.title.text, gp=gpar(fontsize=12, fontface='bold')) ,
                              padding = unit(1,'cm')),
           width = 36,
           height = 10,
           units = 'cm',
           dpi = 600)


  #### write log ####
  print.and.log(paste('AF correlation plot nodes (all):', nrow(all.alleles)),'info',display=.QC$config$debug$verbose)
  print.and.log(paste('AF correlation plot nodes (ambiguous) :', nrow(amb.alleles)),'info',display=.QC$config$debug$verbose)

  print.and.log(sprintf("Allele frequency correlation plot saved!"),
                'info')


  ##### remove plots from RAM ####
  rm(amb.ggplot,
     all.ggplot,
     indel.ggplot,
     similar.palindromic.variants,
     similar.all.variants,
     all.alleles,
     amb.alleles)

  invisible(gc())

}
