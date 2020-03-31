multi.file.comparison <- function() {



  # return if only one file exists
  if(length(.QC$qc.study.list) < 2)
    return(NULL)

  print.and.log('============== Comparing Input Files ==============',
                'info')


  # order files on max sample size
  # .QC$qc.study.list <- tryCatch(
  #   .QC$qc.study.list[order(sapply(.QC$qc.study.list,'[[','MAX_N_TOTAL'))],
  #   error = function(err) {
  #     print.and.log(paste('Error in ordering list of files:',err$message),'warning')
  #     return(.QC$qc.study.list)
  #   }
  # )



  ## text report
  # ========================================
  finalReport.to.txt.file(.QC$config,
                          .QC$qc.study.list)


  ## plots
  # ========================================

  # precision plot
  if(.QC$config$plot_specs$make_plots)
  {

    graphic.device = .QC$graphic.device

    tryCatch( multi.study.precision.plot(.QC$qc.study.list , graphic.device , .QC$config$paths$precisionPlotPath),
              error = function(err){
                print.and.log(paste('error in plotting precision plot:',err$message),'warning',display=.QC$config$debug$verbose)
              }
    )

    # skew-kurt plot
    tryCatch(multi.study.skew.kurt.plot(.QC$qc.study.list, graphic.device , .QC$config$paths$skew_kurt),
             error = function(err){
               print.and.log(paste('error in plotting skewness-kurtosis plot:',err$message),'warning',display=.QC$config$debug$verbose)
             }
    )

    # boxplot effects
    tryCatch(multi.study.eff.plot(.QC$qc.study.list , graphic.device , .QC$config$paths$effsizePlotPath),
             error = function(err){
               print.and.log(paste('error in plotting effect-size comparison plot:',err$message),'warning',display=.QC$config$debug$verbose)
             }
    )

  }
  else
  {
    print.and.log('Plots are skipped!','warning',display=.QC$config$debug$verbose)
  }

  ## free RAM
  invisible(gc())
}

finalReport.to.txt.file <- function(config,study.list)
{

  #Save the final report for multiple file comparisons
  # remove old report file if exists
  if(file.exists(config$paths$txt.report))
    file.remove(config$paths$txt.report)


  writeFileComparisonTXTreport('==================================================')
  writeFileComparisonTXTreport('============= Quality check Report ===============')
  writeFileComparisonTXTreport('==================================================')
  writeFileComparisonTXTreport(' ')
  writeFileComparisonTXTreport(sprintf('Start time: %s', format(config$new_items$starttime, "%b %d %Y - %X")))
  writeFileComparisonTXTreport(sprintf('End time: %s', format(config$new_items$endtime, "%b %d %Y - %X")))

  writeFileComparisonTXTreport(sprintf('Input directory: \'%s\'',config$paths$dir_data))

  writeFileComparisonTXTreport(sprintf('Output directory: \'%s\'',config$paths$dir_output))

  writeFileComparisonTXTreport(sprintf('References directory: \'%s\'',config$paths$dir_references))

  writeFileComparisonTXTreport(sprintf('Alternate header file: \'%s\'',config$supplementaryFiles$header_translations))


  writeFileComparisonTXTreport(sprintf('Allele frequency Reference set: %s', basename(config$supplementaryFiles$allele_ref_std)))

  if(!is.na(config$supplementaryFiles$allele_ref_alt))
    writeFileComparisonTXTreport(sprintf('Alternate Allele frequency Reference set: %s', basename(config$supplementaryFiles$allele_ref_alt)))

  if(!is.na(config$supplementaryFiles$beta_ref_std))
    writeFileComparisonTXTreport(sprintf('Effect size Reference set: %s', basename(config$supplementaryFiles$beta_ref_std)))


  writeFileComparisonTXTreport(' ')
  writeFileComparisonTXTreport(' ')
  writeFileComparisonTXTreport(' ')

  ## ======================================
  report.table <- data.table(sapply(study.list, function(x) return(x$file.name)))
 # report.table <- cbind(seq(1:nrow(report.table)) , report.table)
  report.table <- cbind(sapply(.QC$qc.study.list, function(x) return(x$number)) , report.table)
  colnames(report.table) <- c('Number', 'Study Name')
  row.names(report.table) <- seq(1:nrow(report.table))

  writeFileComparisonTXTreport(kable(report.table,format = 'rst'))


  writeFileComparisonTXTreport(' ')
  writeFileComparisonTXTreport(' ')
  writeFileComparisonTXTreport(' ')

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
    "Effect size" = " ",
    "-     Min." = sapply(study.list, function(x) return(x$tables$variable.summary['Min.', .QC$config$input_parameters$effect_type_string])),
    "-     1st Qu." = sapply(study.list, function(x) return(x$tables$variable.summary['1st Qu.',.QC$config$input_parameters$effect_type_string])),
    "-     Median" = sapply(study.list, function(x) return(x$tables$variable.summary['Median',.QC$config$input_parameters$effect_type_string])),
    "-     Mean" = sapply(study.list, function(x) return(x$tables$variable.summary['Mean',.QC$config$input_parameters$effect_type_string])),
    "-     3rd Qu." = sapply(study.list, function(x) return(x$tables$variable.summary['3rd Qu.',.QC$config$input_parameters$effect_type_string])),
    "-     Max." = sapply(study.list, function(x) return(x$tables$variable.summary['Max.',.QC$config$input_parameters$effect_type_string])),
    "Standard Error (median)" = sapply(study.list, function(x) return(x$tables$variable.summary['Median','STDERR'])),
    "Fixed HWE P-value" = sapply(study.list, function(x) return(x$fixed.hwep)),
    "Fixed Imputation Quality" = sapply(study.list, function(x) return(x$fixed.impq)),
    "Fixed Sample Size" = sapply(study.list, function(x) return(x$fixed.n_total)),
    "Fixed Call Rate" = sapply(study.list, function(x) return(x$fixed.callrate))
  ))

  # colnames(report.table) <- seq(1:length(study.list))
  colnames(report.table) <- sapply(.QC$qc.study.list, function(x) return(x$number))

  writeFileComparisonTXTreport(kable(report.table, align = "c" ,format = 'rst'))


  writeFileComparisonTXTreport(' ')
  writeFileComparisonTXTreport('* step1: removing variants with missing crucial values.')
  writeFileComparisonTXTreport('** step2: removing monomorphic or duplicated variants, and specified chromosomes.')
  writeFileComparisonTXTreport('*** step3: removing mismatched, ambiguous and multi-allelic variants that could not be verified.')



  print.and.log(sprintf('Report file saved as \'%s\'', .QC$config$paths$txt.report),
                'info')

}


