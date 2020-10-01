#' Compare summary statistics of GWAS result files
#'
#' This function compares the key metrics of previously inspected files.
#' This allows the user to check if the results of these studies are comparable (important when running a meta-analysis) and that there are no significant anomalies.
#'
#' @param input.file.list list, full path of the Inspector result files. This file is in RDS format and will be generated for each GWAS result file during the inspection algorithm.
#' @param output.path character, full path to the folder where output files should be saved.
#' @return Key metrics report and plots of previously inspected files are generated and saved in the specified output folder.
#'
compare.GWASs <- function(input.file.list, output.path)
{
  if(missing(input.file.list) || missing(output.path))
    stop('Function arguments are not set.', call. = FALSE)


  on.exit({

  })


  ##==========================

  output.path <- check.output.dir(output.path)


  #
  .QC <- new.env()

  .QC$package.name <- 'GWASinspector'

  .QC$package.description <- 'Comprehensive, efficient and easy to use quality control of genome-wide association study results'

  .QC$script.version <- '1.4.8.2'

  .QC$verbose <- FALSE
  #

  .QC$txt.report.path  <-  paste(output.path,'report.txt',sep = '/')
  .QC$precisionPlotPath <- paste(output.path,'precision_plot.png',sep = '/')
  .QC$skew_kurtPlotPath <- paste(output.path,'skew_kurt.plot.png',sep = '/')
  .QC$effsizePlotPath <- paste(output.path,'effect_size.png',sep = '/')



  .QC$input.file.list  <- check.if.files.exist.and.more.than.one(input.file.list)

  .QC$input.file.list <- lapply(.QC$input.file.list,
                            readRDS)

  for(i in 1:length(.QC$input.file.list))
  {
    .QC$input.file.list[[i]]$number = i
    #input.file.list[[i]]$effect.plot$plot_env$file.number = i
    .QC$input.file.list[[i]]$effect.plot$labels$x = i


    ## compatibility with older versions
    if(is.null(.QC$input.file.list[[i]]$effect_type_string))
      .QC$input.file.list[[i]]$effect_type_string = 'BETA'

    if(is.null(.QC$input.file.list[[i]]$tables$variable.summary.HQ))
      .QC$input.file.list[[i]]$tables$variable.summary.HQ = data.table()
  }


  create.comparison.plots(.QC$input.file.list,
                          'png',
                          .QC$precisionPlotPath,
                          .QC$skew_kurtPlotPath,
                          .QC$effsizePlotPath)

  multi.file.txt.report.file(.QC$input.file.list,
                             .QC$txt.report.path)


  multi.file.html.report.file(.QC,
                              output.path ,
                              'report.html')




}



############################################################

create.comparison.plots <- function(input.file.list,
                                    graphic.device ,
                                    precisionPlotPath ,
                                    skew_kurt ,
                                    effsizePlotPath)
{


  tryCatch( suppressWarnings(multi.study.precision.plot(input.file.list ,
                                                        graphic.device ,
                                                        precisionPlotPath)),
            error = function(err){
              message(paste('Warning in precision plot:',err$message))
            }
  )

  # skew-kurt plot
  tryCatch(multi.study.skew.kurt.plot(input.file.list, graphic.device , skew_kurt),
           error = function(err){
             message(paste('Warning in skewness-kurtosis plot:',err$message))
           }
  )

  # boxplot effects
  tryCatch(multi.study.eff.plot(input.file.list , graphic.device , effsizePlotPath),
           error = function(err){
             message(paste('Warning in effect-size comparison plot:',err$message))
           }
  )
}


check.if.files.exist.and.more.than.one <- function(input.file.list)
{

  if(length(input.file.list) < 2)
    stop('At least 2 study files should be selected.')



  input.file.exists = sapply(input.file.list, check.each.GWAS.file)
  input.file.list = input.file.list[which(input.file.exists)]

  if(length(input.file.list) < 2)
    stop('Not enough file for comparison',call. = FALSE)
  else
    return(input.file.list)
}


check.each.GWAS.file <- function(GWAS.file.path)
{
  if(!file.exists(GWAS.file.path))
  {
    message(paste('GWAS file not found:' , GWAS.file.path))
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}


check.output.dir <- function(output.path)
{
  new.dir <- NULL

  if(!dir.exists(output.path))
  {
    stop(sprintf("%s: %s", "Output directory does not exist.",output.path))
  }
  else
  {
    # remove the trailing backslash from path
    output.path <- gsub('/+$', '' , output.path)

    new.dir <- paste(output.path,'comparing_GWASs',sep = '/')

    tryCatch({
      if (!dir.exists(new.dir))
      {
         dir.create(new.dir)
      }
      else
      {
        new.dir <- sprintf("%s_%s", new.dir, as.numeric(Sys.time()))
        dir.create(new.dir)
      }
    }, error = function(err)
    {
      stop(sprintf("%s: %s", 'Error creating output directory.', err$message))
    })
  }

  message(paste("Report directory:",new.dir))
  return(new.dir)
}

check.output.report.files <- function(output.path, filename)
{
  report.path = paste(output.path,filename,sep = '/')

  # if(file.exists(report.path))
  #   stop(paste('output file exists:',report.path),call. = FALSE)
  # else
    return(report.path)
}


write.to.report.file <- function(message,txt.report.path){

  write(message,
        file= txt.report.path,
        append=TRUE)

}


multi.file.txt.report.file <- function(study.list, txt.report.path)
{

  write.to.report.file('==================================================',txt.report.path)
  write.to.report.file('============== Quality Check Report ==============',txt.report.path)
  write.to.report.file('==================================================',txt.report.path)
  write.to.report.file(' ',txt.report.path)


  ## ======================================
  report.table <- data.table(sapply(study.list, function(x) return(x$file.name)))
  report.table <- cbind(seq(1:nrow(report.table)) , report.table)
  colnames(report.table) <- c('Number', 'Study Name')
  row.names(report.table) <- seq(1:nrow(report.table))

  write.to.report.file(kable(report.table,format = 'rst'),txt.report.path)


  write.to.report.file(' ',txt.report.path)
  write.to.report.file(' ',txt.report.path)
  write.to.report.file(' ',txt.report.path)

  ## ======================================
  report.table <- t(data.table(
    ## 'File Names' = sapply(study.list, function(x) return(x$file.name)),
    'Sample Size (Max)' = sapply(study.list, function(x) return(x$MAX_N_TOTAL)),
    'Missing Columns' = sapply(study.list, function(x) return(paste(x$missing.Columns, collapse = ' | '))),
    'SNPs in input file' = sapply(study.list, function(x) return(x$input.data.rowcount)),
    'Variant count after step 1 *' = sapply(study.list, function(x) return(x$rowcount.step1)),
    'Variant count after step 2 **' = sapply(study.list, function(x) return(x$rowcount.step2)),
    'Variant count after step 3 ***' = sapply(study.list, function(x) return(x$rowcount.step3)),
    'SNP variants' = sapply(study.list, function(x) return(calculatePercent(sum(as.numeric(gsub(x=x$tables$multi_allele_count_preProcess[1:3],
                                                                                                pattern = ",",
                                                                                                replacement = ""))),x$rowcount.step3,pretty=TRUE))),
    'Non-SNP variants' = sapply(study.list, function(x) return(calculatePercent(sum(as.numeric(gsub(x=x$tables$multi_allele_count_preProcess[4:6],
                                                                                                    pattern = ",",
                                                                                                    replacement = "")))
                                                                                ,x$rowcount.step3,
                                                                                pretty=TRUE))),
    'Monomorphic' = sapply(study.list, function(x) return(calculatePercent(x$monomorphic.count,x$rowcount.step3,pretty=TRUE))),
#/	'Duplicates' = sapply(study.list, function(x) return(calculatePercent(x$duplicate.count,x$rowcount.step3,pretty=TRUE))),
    'Palindromics' = sapply(study.list, function(x) return(calculatePercent(x$palindromic.rows,x$rowcount.step3,pretty=TRUE))),
    'Genotyped variants' = sapply(study.list, function(x) return(ifelse(x$tables$imputed.tbl[IMPUTED=='genotyped',.N] == 0,
                                                                        "NA",
                                                                        calculatePercent(as.numeric(x$tables$imputed.tbl[IMPUTED=="genotyped",N]),x$rowcount.step3,pretty=TRUE)))),
    'Imputed variants' = sapply(study.list, function(x) return(ifelse(x$tables$imputed.tbl[IMPUTED=='imputed',.N] == 0,
                                                                      "NA",
                                                                      calculatePercent(as.numeric(x$tables$imputed.tbl[IMPUTED=="imputed",N]),x$rowcount.step3,pretty=TRUE)))),
    'Negative-strand SNPs' = sapply(study.list, function(x) return(calculatePercent(x$neg.strand.count,x$rowcount.step3,pretty=TRUE))),
    'Allele Frequency Correlation (Standard Ref)' = sapply(study.list, function(x) return(x$AFcor.std_ref)),
    'Palindromic Allele Frequency correlation (Standard Ref)' = sapply(study.list, function(x) return(x$AFcor.palindromic.std_ref)),
    'Allele Frequency Correlation (Alternative Ref)' = sapply(study.list, function(x) return(x$AFcor.alt_ref)),
    'Palindromic Allele Frequency correlation (Alternative Ref)' = sapply(study.list, function(x) return(x$AFcor.palindromic.alt_ref)),
    'Lambda - Total' = sapply(study.list, function(x) return(x$lambda)),
    'Lambda - Genotyped' = sapply(study.list, function(x) return(x$lambda.gen)),
    'Lambda - Imputed' = sapply(study.list, function(x) return(x$lambda.imp)),
    'P-value Correlation' = sapply(study.list, function(x) return(x$PVcor)),
    "Visscher's Statistic (HQ variants)" = sapply(study.list, function(x) return(x$Visschers.stat.HQ)),
    "Fixed HWE P-value" = sapply(study.list, function(x) return(x$fixed.hwep)),
    "Fixed Imputation Quality" = sapply(study.list, function(x) return(x$fixed.impq)),
    "Fixed Sample Size" = sapply(study.list, function(x) return(x$fixed.n_total)),
    "Fixed Call Rate" = sapply(study.list, function(x) return(x$fixed.callrate)),
    "Effect size (All variants)" = " ",
    "-     Min." = sapply(study.list, function(x) return(x$tables$variable.summary['Min.', x$effect_type_string])),
    "-     1st Qu." = sapply(study.list, function(x) return(x$tables$variable.summary['1st Qu.',x$effect_type_string])),
    "-     Median" = sapply(study.list, function(x) return(x$tables$variable.summary['Median',x$effect_type_string])),
    "-     Mean" = sapply(study.list, function(x) return(x$tables$variable.summary['Mean',x$effect_type_string])),
    "-     3rd Qu." = sapply(study.list, function(x) return(x$tables$variable.summary['3rd Qu.',x$effect_type_string])),
    "-     Max." = sapply(study.list, function(x) return(x$tables$variable.summary['Max.',x$effect_type_string])),
    "-     Min.(HQ variants)" = sapply(study.list, function(x) return(ifelse(nrow(x$tables$variable.summary.HQ) > 0,
                                                                x$tables$variable.summary.HQ['Min.', x$effect_type_string],
                                                                NA))),
    "-     Max.(HQ variants)" = sapply(study.list, function(x) return(ifelse(nrow(x$tables$variable.summary.HQ) > 0,
                                                                x$tables$variable.summary.HQ['Max.',x$effect_type_string],
                                                                NA))),
    "Standard Error (median)" = sapply(study.list, function(x) return(x$tables$variable.summary['Median','STDERR'])),
    "Standard Error (median) (HQ variants)" = sapply(study.list, function(x) return(ifelse(nrow(x$tables$variable.summary.HQ) > 0,
                                                                                  x$tables$variable.summary.HQ['Median','STDERR'],
                                                                                  NA)))
  )
)

  colnames(report.table) <- seq(1:length(study.list))

  write.to.report.file(knitr::kable(report.table, align = "c" ,format = 'rst'),txt.report.path)


  write.to.report.file(' ',txt.report.path)
  write.to.report.file('* step1: removing variants with missing crucial values.',txt.report.path)
  write.to.report.file('** step2: removing monomorphic or duplicated variants, and specified chromosomes.',txt.report.path)
  write.to.report.file('*** step3: removing mismatched, ambiguous and multi-allelic variants that could not be verified.',txt.report.path)

  write.to.report.file('\n==============================================',txt.report.path)
  write.to.report.file(sprintf('Generated by GWASinspector package - v.%s', utils::packageVersion('GWASinspector')),txt.report.path)
  write.to.report.file(as.character(Sys.time()),txt.report.path)

  message(sprintf("Output file saved from %s input files." ,length(study.list)))
}

multi.file.html.report.file <- function(.QC, html.report.folder , html.report.file)
{
  input.file.list <- .QC$input.file.list

  if(!rmarkdown::pandoc_available())
  {
    message('pandoc module is required for converting report to Html format! check the manual on how to install.')
    return(NULL)
  }


  if(!is.element('kableExtra',installed.packages()))
  {
    message('kableExtra package is suggested for pretty Html format! check the manual on how to install.')
    report.template.file <- system.file("rmd", "multiFileReport.rmd", package = "GWASinspector")
  }
  else
  {
   # `%>%` <- NULL
   # `%>%` <<- kableExtra::`%>%`
    report.template.file <- system.file("rmd", "multiFileReport_alone_extra.rmd", package = "GWASinspector")
  }




    render.success <- tryCatch({
      # clear cache and RAM
      knitr::knit_meta(class=NULL, clean = TRUE)
      invisible(gc())

      render(report.template.file,
             output_dir = html.report.folder,
             output_file = html.report.file,
             quiet = TRUE)

      return(TRUE)

    },
    error=function(err){

      print.and.log(paste('---[ERROR saving main html file!---]\nThe result is also saved as txt and is in the output folder.',err$message),
                    'warning')

      return(FALSE)
    } )

    if(render.success)
      message(paste('HTML report file saved as ',report.output.path))


  }
