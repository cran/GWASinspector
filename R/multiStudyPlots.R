multi.study.precision.plot <- function(study.list, graphic.device , figure.path){

  study.names <- sapply(study.list, get.study.name) # studyFunctions.R

  # find the maximum length of file names for plot aspect and width
  # it should be at least 160
  p.width <-  max(sapply(study.names,nchar)) * 2
  p.width <- ifelse(p.width < 160 , 160 , p.width)

  precision.table <- as.data.table(t(sapply(study.list, get.precision.plot.values))) #multiStudyFunctions.R
  precision.table <- cbind(precision.table,study.names)

  precision.table$SE.mean.HQ <- unlist(precision.table$SE.mean.HQ)
  precision.table$sqrt.n <- unlist(precision.table$sqrt.n)
  precision.table$study.names <- unlist(precision.table$study.names)

  # check if it has missing or NA values - No HQ variants
  precision.table[is.na(SE.mean.HQ) | is.na(sqrt.n) , missing := TRUE]

  if(any(!is.na(precision.table$missing))){
    print.and.log(sprintf("studies removed from 'Precision Plot' due to NA value (possibly no HQ variants or missing N_Total column)! (%s)",
                          paste(precision.table[missing == TRUE]$study.names,collapse = ' | ')),'warning',display=.QC$config$debug$verbose)
    precision.table <- precision.table[is.na(missing)] ## if 'missing' column is not set to TRUE it is null and not missing
  }

  if(nrow(precision.table) > 1)
  {
    #add number to display in plot beside each point
    # precision.table$order <- 1:nrow(precision.table)
    precision.table$order <- sapply(.QC$qc.study.list, function(x) return(x$number))
    precision.table[,ordered.names := sprintf('%s- %s',order,study.names)]

    # minimum and maximum values for x axis
    xmin <- min(precision.table$sqrt.n) - 5
    xmax <- max(precision.table$sqrt.n) + 5


    plot <- ggplot(precision.table,aes(y=SE.mean.HQ,x=sqrt.n,label=order)) +
      geom_point(aes(color = ordered.names))+
      labs(title="Precision by Sample Size",
           x=expression(sqrt('sample size')),
           y=expression('1/mean(Standard Error)'),
           subtitle="High quality SNPs only!") +
      geom_text(hjust = -1.2, nudge_x = 0.1,  size=3) +
      coord_cartesian(xlim = c(xmin,xmax))+
      theme_classic() +
      theme(strip.background = element_blank(),
            axis.title.x = element_text(size=8, face = "bold"),
            axis.text.x = element_text(size=6, face = "bold"),
            axis.title.y = element_text(size=8, face = "bold"),
            axis.text.y = element_text(size=6, face = "bold"),legend.key.height=unit(.5,"line")
            ,legend.position = 'right'
            ,legend.text = element_text(colour="black", size=5)
            ,legend.title.align = 0.5
            ,legend.title = element_text(colour="black", size=5,
                                         face="bold")
            ,legend.background = element_rect(fill="gray90",
                                              size=0.5, linetype="solid",
                                              colour ="darkblue")
            ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
            ,plot.subtitle=element_text(size=7, hjust=0.5, face="italic", color="darkblue"))+
      guides(col = guide_legend(ncol = 1, byrow = TRUE)) +
      labs(shape="study names", colour="study names")


    ggsave(plot=plot,
           device = graphic.device ,
           filename = figure.path,
           units = c('mm'),
           width = p.width,
           height =100,
           dpi = 300)




    rm(plot)
    invisible(gc())

    print.and.log("Precision plot is saved!",'info')

  }else
  {
    print.and.log("Precision plot not saved due to insufficient data!",'warning',display=.QC$config$debug$verbose)
  }

}






multi.study.skew.kurt.plot <- function(study.list, graphic.device , figure.path){

  study.names <- sapply(study.list, get.study.name) # studyFunctions.R

  # find the maximum length of file names for plot aspect and width
  # it should be at least 160
  p.width <-  max(sapply(study.names,nchar)) * 2
  p.width <- ifelse(p.width < 160 , 160 , p.width)


  skew.kurt.table <- as.data.table(t(sapply(study.list, get.skewness.kurtosis))) #multiStudyFunctions.R
  skew.kurt.table <- cbind(skew.kurt.table,study.names)

  skew.kurt.table$kurtosis <- unlist(skew.kurt.table$kurtosis)
  skew.kurt.table$skewness <- unlist(skew.kurt.table$skewness)
  skew.kurt.table$study.names <- unlist(skew.kurt.table$study.names)

  # check if it has missing or NA values - No HQ variants
  skew.kurt.table[is.na(skewness) | is.na(kurtosis) , missing := TRUE]
  if(nrow(skew.kurt.table[missing == TRUE]) > 0){
    print.and.log(sprintf("studies removed from 'Skewness-Kurtosis Plot' due to NA value (possibly no HQ variants)! (%s)",
                          paste(skew.kurt.table[missing==TRUE]$study.names,collapse = ' | ')),'warning',display=.QC$config$debug$verbose)
    skew.kurt.table <- skew.kurt.table[is.na(missing)]
  }


  if(nrow(skew.kurt.table) > 1)
  {
    ## add order values to table from 1 to row count
    skew.kurt.table$order <- 1:nrow(skew.kurt.table)



    if(!is.null(.QC$qc.study.list))  # if function is called from inspect()
      skew.kurt.table$order <- sapply(.QC$qc.study.list, function(x) return(x$number))
    else if(!is.null(study.list)) # if function is called independently
      skew.kurt.table$order <- sapply(study.list, function(x) return(x$number))



        skew.kurt.table[,ordered.names := sprintf('%s- %s',order,study.names)]

    xmin <- min(skew.kurt.table$skewness) - 0.05
    xmax <- max(skew.kurt.table$skewness) + 0.05

    plot <- ggplot(skew.kurt.table,aes(y=kurtosis,x=skewness,label=order)) +
      geom_point(aes(color=ordered.names))+
      labs(title="Skewness vs Kurtosis plot",
           x="Skewness",
           y="Kurtosis",
           subtitle="High quality SNPs only!") +
      geom_text(hjust = 0.005, nudge_x = 0.003, size =2) +
      coord_cartesian(xlim = c(xmin,xmax))+
      theme_classic() +
      theme(strip.background = element_blank(),
            axis.title.x = element_text(size=8, face = "bold"),
            axis.text.x = element_text(size=6, face = "bold"),
            axis.title.y = element_text(size=8, face = "bold"),
            axis.text.y = element_text(size=6, face = "bold"),legend.key.height=unit(.5,"line")
            ,legend.position = 'right'
            ,legend.text = element_text(colour="black", size=5)
            ,legend.title.align = 0.5
            ,legend.title = element_text(colour="black", size=5,
                                         face="bold")
            ,legend.background = element_rect(fill="gray90",
                                              size=0.5, linetype="solid",
                                              colour ="darkblue")
            ,plot.title=element_text(size=10, face="bold",hjust = 0.5)
            ,plot.subtitle=element_text(size=7, hjust=0.5, face="italic", color="darkblue")) +
      guides(col = guide_legend(ncol = 1, byrow = TRUE)) +
      labs(shape="study names", colour="study names")


    ggsave(plot=plot,
           filename = figure.path,
           device = graphic.device ,
           units = c('mm'),
           width = p.width,
           height =100,
           dpi = 300)


    rm(plot)
    invisible(gc())

    print.and.log("Skewness vs kurtosis plot is saved!",'info')
  }else
  {
    print.and.log("Skewness vs kurtosis plot not saved due to insufficient data!",'warning',display=.QC$config$debug$verbose)
  }
}






multi.study.eff.plot <- function(study.list, graphic.device , figure.path)
{

  # order files on max sample size
  study.list = tryCatch(
    study.list[order(sapply(study.list,'[[','MAX_N_TOTAL'))],
    error = function(err) {
      print.and.log(paste('Error in ordering list of files:',err$message),'warning',display=.QC$config$debug$verbose)
      return(study.list)
    }
  )

  effect.plot.list <- lapply(study.list, function(x) return(x$effect.plot))


  ggsave(plot = arrangeGrob(grobs = effect.plot.list,
                            ncol = length(effect.plot.list)),
         filename = figure.path ,
         device = graphic.device ,
         units = c('mm'),
         # width = length(effect.plot.list) * 120,
         width = length(effect.plot.list) * 40,
         height =50,
         dpi=150,
         limitsize = FALSE )


  invisible(gc())


  print.and.log("Effect-size box plot is saved!",'info')
}
