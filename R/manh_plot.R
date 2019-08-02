#' Creates the Manhattan plot
#'
#' A function to generate Manhattan plots.
#'
#' @param dataset Data frame or data table containing the below columns
#' @param chr Name of chromosome column
#' @param pvalue Name of P-value column
#' @param position Name of position column
#' @param fileName Full name and path of file to be saved (file extension should be 'png'). e.g. “c:/users/researcher/study/man_plot.png”
#' @param plot.title Title of the plot, default value is 'Manhattan plot'
#' @param plot.subtitle Subtitle of the plot
#' @param sig.threshold.log The -log10 transformed significance threshold, used for plotting a threshold line (e.g. 8 = 10^-8)
#' @param p.threshold Threshold for plotting variants (i.e. p-values > 0.01 will not be plotted). Setting a higher threshold will significantly increase plotting time
#' @param beta (optional) Name of the effect-size column
#' @param std.error (optional) Name of the standard error column
#' @param check.columns Whether to check input columns for invalid values
#' @return Generates and saves a Manhattan plot for the provided data.
#' @examples
#' input.data = read.table(gzfile(system.file("extdata", "sample.txt.gz", package = "GWASinspector")),
#'                         header = TRUE,
#'                         stringsAsFactors = FALSE,
#'                         fill = TRUE)
#' tmpPlotFile = paste(tempfile(),'png',sep = '.')
#' man.plot(dataset = input.data, chr = 'CHR', pvalue = 'PVALUE', position = 'POSITION',
#'          plot.title = 'Manhattan plot', plot.subtitle = 'This data is fabricated!',
#'          fileName = tmpPlotFile , p.threshold = '0.5')
#'
man.plot <- function(dataset,
                     chr,
                     pvalue,
                     position,
                     fileName,
                     plot.title = 'Manhattan Plot',
                     plot.subtitle = '',
                     p.threshold = 0.01,
                     sig.threshold.log = -log10(5*10^-8),
                     beta=NULL,
                     std.error = NULL,
                     check.columns = TRUE) {


  if(!is.data.table(dataset))
    dataset <- as.data.table(dataset)

  # check if columns are present
  if(!is.element(chr,names(dataset))){
    print.and.log('skipping Manhattan plot! \'CHR\' column not found.','warning')
    return(NULL)
  }

  if(!is.element(pvalue,names(dataset))){
    print.and.log('skipping Manhattan plot! \'P-value\' column not found.','warning')
    return(NULL)
  }

  if(!is.element(position,names(dataset))){
    print.and.log('skipping Manhattan plot! \'Position\' column not found.','warning')
    return(NULL)
  }

  ## get data column names
  col.names <- colnames(dataset)

  selected.col.names <- c(chr,pvalue,position,beta,std.error)

  if('PVALUE.calculated' %in% col.names) # if function is called from package and contains calculated P-values
    selected.col.names <- c(selected.col.names, 'PVALUE.calculated')

  ## 1- check if function parameters are correctly set
  missing.beta <- FALSE
  missing.std.err <- FALSE


  if(is.null(beta))
    missing.beta <- TRUE

  if(is.null(std.error))
    missing.std.err <- TRUE




  wrong.col.index <- which(selected.col.names %notin% col.names)
  if(length(wrong.col.index) > 0)
  {
    wrong.col.names <- selected.col.names[wrong.col.index]

    stop(sprintf('Manhattan plot error: selected column not found in dataset (%s)!',
                 paste(wrong.col.names,collapse = ' | ')))

  }


  ## check thresholds
  if(p.threshold > 1 | p.threshold <= 0 )
    stop(sprintf('Manhattan plot error: Invalid P-value threshold: %s!',p.threshold))

  if(sig.threshold.log > 300 | sig.threshold.log <= 1 )
    stop(sprintf('Manhattan plot error: Invalid significant threshold: %s!',sig.threshold.log))


  ##2- subset input data and calculate missing pvalues

  if('HQ' %in% col.names) # filter data if function is called from package and contains HQ column
    plot.data <- subset(dataset[HQ==TRUE], select = selected.col.names)
  else
    plot.data <- subset(dataset, select = selected.col.names)

  ## 3- rename columns
  names(plot.data)[names(plot.data) == pvalue] <- 'PVALUE'
  names(plot.data)[names(plot.data) == beta] <- 'EFFECT'
  names(plot.data)[names(plot.data) == std.error] <- 'STDERR'
  names(plot.data)[names(plot.data) == chr] <- 'CHR'
  names(plot.data)[names(plot.data) == position] <- 'POSITION'


  ## 4- check columns for invalid values and chromosome numbers
  # this step is not required if function is called from package and as part of GWASinspector
  # check.columns parameter is false from GWASinspector
  if(check.columns){
    plot.data <- process.column.CHR(plot.data)
    plot.data <- process.column.POSITION(plot.data)
    plot.data <- process.column.PVALUE(plot.data)

    if(!missing.beta)
      plot.data <- process.column.EFFECT(plot.data)
    if(!missing.std.err)
      plot.data <- process.column.STDERR(plot.data)
  }
  # 5- calcualte missing p
  missing.p.count <- nrow(plot.data[is.na(PVALUE)])

  if( missing.p.count > 0){ # if there are missing pvaluse, try to replace them from calculated  values or calculate now
    print.and.log(sprintf('missing P-values are found in dataset (%s from %s variants)!',
                          missing.p.count,
                          nrow(plot.data)),
                  'warning')

    if('PVALUE.calculated' %in% selected.col.names) # if function is called from package
    {
      plot.data[is.na(PVALUE), PVALUE := PVALUE.calculated]
      print.and.log('missing P-values are calculated and replaced!','info')
    }
    else{
      if(!missing.std.err & !missing.beta)
      {
        plot.data[is.na(PVALUE), PVALUE := pchisq((EFFECT/STDERR)^2, 1, lower.tail=FALSE)]
        print.and.log('missing P-values are calculated from STDERR and Beta values!'
                      ,'warning')
      }
      else
      {
        print.and.log('Define STDERR and Beta columns for calculating missing items.'
                      ,'warning')
      }
    }
  }

  # 6- plot data
  # only select variant with existing Chromosome and low P values
  plot.data <- plot.data[!is.na(CHR) & PVALUE < p.threshold]

  ## manhattanPlotFunction.R
  plot.manhattan(plot.data, plot.title,plot.subtitle, sig.threshold.log, fileName)




}


