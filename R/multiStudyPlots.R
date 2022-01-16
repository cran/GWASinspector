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

  if(nrow(precision.table) >= 1)
  {
    #add number to display in plot beside each point
    #precision.table$order <- 1:nrow(precision.table)
    # if(!is.null(study.list)) # if function is called independently
    #   precision.table$order <- sapply(study.list, function(x) return(x$number))
    # else if(!is.null(.QC$qc.study.list))  # if function is called from inspect()
    #   precision.table$order <- sapply(.QC$qc.study.list, function(x) return(x$number))


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


  if(nrow(skew.kurt.table) >= 1)
  {
    ## add order values to table from 1 to row count
    #skew.kurt.table$order <- 1:nrow(skew.kurt.table)

    # if(!is.null(study.list)) # if function is called independently
    #   skew.kurt.table$order <- sapply(study.list, function(x) return(x$number))
    # else if(!is.null(.QC$qc.study.list))  # if function is called from inspect()
    #   skew.kurt.table$order <- sapply(.QC$qc.study.list, function(x) return(x$number))



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
    {
      ## If all studies have N_CASES column, use this column for sorting
      if(all(sapply(study.list, function(x) ifelse("N_CASES" %in% x$renamed.File.Columns, TRUE , FALSE))))
      {
        study.list[order(sapply(study.list,'[[','MAX_N_CASES'))]
        print.and.log("N_CASES column was used for ordering the plots.")
      }
      else
      {
        study.list[order(sapply(study.list,'[[','MAX_N_TOTAL'))]
      }

    },
    error = function(err) {
      print.and.log(paste('Error in ordering list of files:',err$message),'warning',display=.QC$config$debug$verbose)
      return(study.list)
    }
  )

  ## Generate effect-size plot instead of loading it
  ##effect.plot.list <- lapply(study.list, function(x) return(x$effect.plot))
  effect.plot.list <- lapply(study.list, function(x) return(generateEffectSizePlot(x)))

  # remove null studies, with no HQ variants
  effect.plot.list2 <- effect.plot.list[!sapply(effect.plot.list,is.null)]

  #if(any(sapply(effect.plot.list,is.null)))
  if(length(effect.plot.list2) != length(effect.plot.list))
  {
    #print.and.log("Effect-size box plot NOT saved due to insufficient data!",'warning',display=.QC$config$debug$verbose)
    print.and.log(sprintf("%s studies are removed from effect-size box plot!",length(effect.plot.list) - length(effect.plot.list2)),
                  'warning',
                  display=.QC$config$debug$verbose)
  }

  if(length(effect.plot.list2) >= 1)
  {
    ggsave(plot = arrangeGrob(grobs = effect.plot.list2,
                              ncol = length(effect.plot.list2)),
           filename = figure.path ,
           device = graphic.device ,
           units = c('mm'),
           # width = length(effect.plot.list) * 120,
           width = length(effect.plot.list2) * 40,
           height =50,
           dpi=150,
           limitsize = FALSE )


    invisible(gc())


    print.and.log("Effect-size box plot is saved!",'info')
  }else
  {
    print.and.log("Effect-size box plot not saved due to insufficient data!",'warning',display=.QC$config$debug$verbose)
  }
}


generateEffectSizePlot <- function(study)
{

  if(is.null(study$effect.plot.df) | any(sapply(study$effect.plot.df,is.na)))
  {
    return(NULL)
  }else
  {
    df <- study$effect.plot.df

    file.number = study$number

    if("N_CASES" %in% study$renamed.File.Columns)
    {
      file.N.max = study$MAX_N_CASES
      print.and.log("N_CASES will be used for MAX_N value.")
    }
    else
      file.N.max = study$MAX_N_TOTAL



    plot = ggplot(df, aes(x)) +
      geom_boxplot(
        aes(ymin = y0,
            lower = y25,
            middle = y50,
            upper = y75,
            ymax = y100),
        stat = "identity") +
      labs(x= file.number ,
           y="effect size",
           subtitle = sprintf('N = %s', file.N.max)) +
      theme_classic(base_size = 8)+
      coord_cartesian(ylim = c(-0.6,0.6)) +
      geom_hline(yintercept = 0.1,
                 linetype = 2,
                 color='red') +
      geom_hline(yintercept = -0.1,
                 linetype = 2,
                 color='red') +
      theme(axis.text.x=element_blank()) +
      if(!is.null(study$effect.plot.df_y_upper) && !is.null(study$effect.plot.df_y_lower))
      {
        geom_errorbar(aes(ymin=study$effect.plot.df_y_lower,
                          ymax=study$effect.plot.df_y_upper),
                      width=0.6,
                      position=position_dodge(.9))
      }


    if(is.null(study$effect.plot.df_y_upper) || is.null(study$effect.plot.df_y_lower))
      print.and.log("Error-bars missing for effect-size plot. Horizontal lines are at -0.1 and 0.1.","warning")

    return(plot)
  }
}
