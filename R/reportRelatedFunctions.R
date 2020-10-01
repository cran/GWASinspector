create.report.files <- function() {


  # print.and.log('\n','info')
  print.and.log('============== Creating Report Files ==============',
                'info')

  if(!.QC$pandoc.exists)
    print.and.log('pandoc module is required for converting report to Html format! check the manual on how to install.','warning',display=.QC$config$debug$verbose)

  if(!.QC$kableExtra.package.exists)
    print.and.log('kableExtra package is suggested for pretty Html format! check the manual on how to install.','warning',display=.QC$config$debug$verbose)



  if(.QC$pandoc.exists){

    tryCatch(
      # multi file comparison report - html
      writeMainReportFile(), #reportRelatedFunctions.R
      error = function(err) print.and.log(paste0('Error in converting main report to html format. %s ',err$message))
    )

    tryCatch(
      # file specific report - html
      writeStudyReportFile(), #reportRelatedFunctions.R
      error = function(err) print.and.log(paste0('Error in converting input file report to html format. %s ',err$message))
    )
  }
  else
    print.and.log('Writing Html report is skipped! required packages not found.','warning',display=.QC$config$debug$verbose)





  # EXCEL
  writeExcelReportFile()


}




writeMainReportFile <- function() {

  # FIXME do nothing and return if only one file is selected
  # create report of only one file exists !!!??
  if(length(.QC$qc.study.list) == 1)
    return(NULL)


  # check if template file exists and get the path
  multi.file.report.template <- get.multi.file.report.template()

  # if user wants the report file and template file exists
  if(.QC$config$output_parameters$html_report & !is.null(multi.file.report.template))
  {

    # path of html file
    report.output.path <- .QC$config$paths$html.report

    render.success <- tryCatch({
      # clear cache and RAM
      knitr::knit_meta(class=NULL, clean = TRUE)
      invisible(gc())

      render(multi.file.report.template,
             output_dir = .QC$config$paths$dir_output,
             output_file = report.output.path,
             quiet = TRUE)

      print.and.log(sprintf('HTML report file saved as %s!',report.output.path),
                    'info')
      return(TRUE)

    },
    error=function(err){

      print.and.log(paste('---[ERROR saving main html file!---]\nThe result is also saved as txt and is in the output folder.',err$message),
                    'warning',display=.QC$config$debug$verbose)

      return(FALSE)
    }
    )

    if(render.success)
      print.and.log(sprintf('HTML report file saved as %s!',report.output.path),
                    'info')

  }else
  {
    print.and.log('Writing the report file is skipped!','info')
  }
}


writeStudyReportFile <- function(){

  # check if template file exists and get the path
  report.template <- get.study.specific.report.template()



  # if user wants the report file and template file exists
  if(.QC$config$output_parameters$html_report & !is.null(report.template))
  {


    sapply(.QC$qc.study.list, function(study){
      tryCatch({
        .QC$thisStudy <- study

        # path of html file
        report.output.path <- study$html.report.path

        # clear cache and RAM
        knitr::knit_meta(class=NULL, clean = TRUE)
        invisible(gc())

        render(report.template,
               output_dir = .QC$config$paths$dir_output,
               output_file = report.output.path,
               quiet = TRUE)


        print.and.log(sprintf('HTML report file saved as %s!',report.output.path),
                      'info')
      }
      ,error=function(err){
        print.and.log(paste('---[ERROR saving study html file!---]\nThe result is also saved as txt and is in the output folder.',err$message),
                      'warning',display=.QC$config$debug$verbose)
      }
      )
    })
  }else
  {
    print.and.log('Writing the report file is skipped!','info')
  }
}


get.study.specific.report.template <- function() {


  # check if package default report template file is present and accessible. report is skipped if template file not found
  if(.QC$kableExtra.package.exists)
  {
    report.template.file <- system.file("rmd", "mainReport_extra.rmd", package = "GWASinspector")
  }
  else
  {
    report.template.file <- system.file("rmd", "mainReport.rmd", package = "GWASinspector")
  }

  if(file.exists(report.template.file))
    return(report.template.file)
  else
  {
    print.and.log('Report template file is not found in package! try re-installing GWASinspector package.','warning',display=.QC$config$debug$verbose)
    print.and.log('Writing the report file is skipped!','info')
    return(NULL)
  }
}


get.multi.file.report.template <- function() {
  # check if package default report template file is present and accessible. report is skipped if template file not found
  if(.QC$kableExtra.package.exists)
    report.template.file <- system.file("rmd", "multiFileReport_extra.rmd", package = "GWASinspector")
  else
    report.template.file <- system.file("rmd", "multiFileReport.rmd", package = "GWASinspector")


  if(file.exists(report.template.file))
    return(report.template.file)
  else
  {
    print.and.log('Main-Report template file is not found in package! try re-installing GWASinspector package.','warning',display=.QC$config$debug$verbose)
    print.and.log('Writing the main-report file is skipped!','info')
    return(NULL)
  }
}


# display a report table to user for each input file
report.to.txt.file <- function(study) {

  # remove old report file if exists
  if(file.exists(study$txt.report.path))
    file.remove(study$txt.report.path)

  # report intro
  writeTXTreport('============================================================')
  writeTXTreport(sprintf('================= %s v.%s ==================',
                         .QC$package.name,
                         .QC$script.version))
  writeTXTreport('============================================================')

  writeTXTreport(' ')
 # writeTXTreport(paste('Script version:', .QC$script.version))
  writeTXTreport(paste('System Information:', .QC$r.version))
  writeTXTreport(sprintf('Start Time: %s', format(study$starttime, "%b %d %Y - %X")))
  writeTXTreport(sprintf('End Time: %s', format(study$endtime, "%b %d %Y - %X")))


  writeTXTreport(' ')


  ### ==================================
  writeTXTreport(' ')
  writeTXTreport('==========================================================')
  writeTXTreport('==================== User preferences ====================')
  writeTXTreport('==========================================================')
  writeTXTreport(' ')
  writeTXTreport(sprintf('Alterative header file: %s', basename(.QC$config$supplementaryFiles$header_translations)))
  writeTXTreport(sprintf('Allele frequency standard reference dataset: %s', basename(.QC$config$supplementaryFiles$allele_ref_std)))


  if(!is.na(.QC$config$supplementaryFiles$allele_ref_alt))
    writeTXTreport(sprintf('Allele frequency alternate reference dataset: %s', basename(.QC$config$supplementaryFiles$allele_ref_alt)))

  if(!is.na(.QC$config$supplementaryFiles$beta_ref_std))
    writeTXTreport(sprintf('Effect size reference dataset: %s', basename(.QC$config$supplementaryFiles$beta_ref_std)))

  writeTXTreport(' ')

  # ===================================
  writeTXTreport('==========================================================')
  writeTXTreport('= Filter values for selecting High-Quality (HQ) variants =')

  filter.table <- data.table(
    "Allele frequency" =  format(.QC$config$filters$HQfilter_FRQ,scientific = FALSE))

  if("HWE_PVAL" %in% study$renamed.File.Columns.sorted)
    filter.table <- cbind(filter.table, "HWE p-value" = format(.QC$config$filters$HQfilter_HWE,scientific = FALSE))
  else
    filter.table <- cbind(filter.table, "HWE p-value" = "Not included")


  if("CALLRATE" %in% study$renamed.File.Columns.sorted)
    filter.table <- cbind(filter.table, "Call-rate" = format(.QC$config$filters$HQfilter_cal,scientific = FALSE))
  else
    filter.table <- cbind(filter.table, "Call-rate" = "Not included")


  if("IMP_QUALITY" %in% study$renamed.File.Columns.sorted)
    filter.table <- cbind(filter.table, "Imputation quality" = format(.QC$config$filters$HQfilter_imp,scientific = FALSE))
  else
    filter.table <- cbind(filter.table, "Imputation quality" = "Not included")

  # filter.table <- t(data.table(
  #   "Allele frequency" =  format(.QC$config$filters$HQfilter_FRQ,scientific = FALSE),
  #   {
  #
  #   },
  #   {
  #     if("CALLRATE" %in% study$renamed.File.Columns.sorted)
  #        "Call-rate" =  format(.QC$config$filters$HQfilter_cal,scientific = FALSE)
  #   },
  #   {
  #     if("IMP_QUALITY" %in% study$renamed.File.Columns.sorted)
  #        "Imputation quality" = format(.QC$config$filters$HQfilter_imp, scientific = FALSE)
  #   }
  #   ))
  filter.table <- t(filter.table)

  colnames(filter.table) <- 'Value'
  writeTXTreport(kable(filter.table,format = "rst"))


  writeTXTreport(' ')

  writeTXTreport(paste('Effect type:', .QC$config$input_parameters$effect_type))

  ### ==================================
  writeTXTreport(' ')
  writeTXTreport(' ')
  writeTXTreport('==========================================================')
  writeTXTreport('================= Input file description =================')
  writeTXTreport('==========================================================')
  writeTXTreport(' ')

  # input file spec
  writeTXTreport(sprintf('Input file name: %s', basename( study$file.path)))
  writeTXTreport(sprintf('Input file line count (including header): %s', study$file.line.count))
  writeTXTreport(sprintf('Input file ends with a new line: %s', study$file.endsWithNewLine))
  # writeTXTreport(sprintf('Duplicated lines: %s', format(.QC$thisStudy$dup_lines_count,big.mark = ',',scientific = FALSE)))

  # it is mentioned in log file
  # if(study$hanNoneBaseAlleles)
  #   writeTXTreport('WARNING: Input file has unknown character for INDEL variants!')

  writeTXTreport(' ')


  ## column names and translations
  writeTXTreport(' ')
  writeTXTreport('========== Column names and translations ================')
  writeTXTreport(' ')
  column.tbl <- rbind(.QC$thisStudy$original.File.Columns.sorted,
                      .QC$thisStudy$renamed.File.Columns.sorted)
  rownames(column.tbl) <- c('Original', 'Renamed')

  writeTXTreport(kable(t(column.tbl),format = "rst"))

  writeTXTreport(' ')
  writeTXTreport(' ')

  writeTXTreport('================== Column report  ======================')
  writeTXTreport(' ')

  ### invalid items
  b <- t(data.frame('CHR' = c(abs(study$column.NA.list$CHR - study$column.INVALID.list$CHR) ,
                              study$column.INVALID.list$CHR,
                              ' '),

                    'POSITION' = c(abs(study$column.NA.list$POSITION - study$column.INVALID.list$POSITION) ,
                                   study$column.INVALID.list$POSITION,
                                   ' '),

                    'EFFECT_ALL' = c(abs(study$column.NA.list$EFFECT_ALL - study$column.INVALID.list$EFFECT_ALL) ,
                                     study$column.INVALID.list$EFFECT_ALL,
                                     ' '),

                    'OTHER_ALL' = c(abs(study$column.NA.list$OTHER_ALL - study$column.INVALID.list$OTHER_ALL) ,
                                    study$column.INVALID.list$OTHER_ALL,
                                    ' '),

                    'EFFECT' = c(abs(study$column.NA.list$EFFECT - study$column.INVALID.list$EFFECT) ,
                                # study$column.INVALID.list$EFFECT,
                                 ' ',
                                 ' '),

                    'STDERR' = c(abs(study$column.NA.list$STDERR - study$column.INVALID.list$STDERR - study$column.INVALID.list$zero.STDERR) ,
                                 study$column.INVALID.list$STDERR,
                                 study$column.INVALID.list$zero.STDERR),

                    'EFF_ALL_FREQ' = c(abs(study$column.NA.list$EFF_ALL_FREQ - study$column.INVALID.list$EFF_ALL_FREQ - study$column.INVALID.list$minusone.EFF_ALL_FREQ),
                                       study$column.INVALID.list$EFF_ALL_FREQ,
                                       study$column.INVALID.list$minusone.EFF_ALL_FREQ),

                    'HWE_PVAL' = c(abs(study$column.NA.list$HWE_PVAL - study$column.INVALID.list$HWE_PVAL - study$column.INVALID.list$minusone.HWE_PVAL) ,
                                   study$column.INVALID.list$HWE_PVAL,
                                   study$column.INVALID.list$minusone.HWE_PVAL),

                    'PVALUE' = c(abs(study$column.NA.list$PVALUE - study$column.INVALID.list$PVALUE - study$column.INVALID.list$minusone.PVALUE) ,
                                 study$column.INVALID.list$PVALUE,
                                 study$column.INVALID.list$minusone.PVALUE),

                    'IMPUTED' = c(abs(study$column.NA.list$IMPUTED - study$column.INVALID.list$IMPUTED),
                                  study$column.INVALID.list$IMPUTED,
                                  ' '),

                    'IMP_QUALITY' = c(abs(study$column.NA.list$IMP_QUALITY - study$column.INVALID.list$IMP_QUALITY) ,
                                      study$column.INVALID.list$IMP_QUALITY,
                                      ' '),

                    'MARKER' = c(abs(study$column.NA.list$MARKER - study$column.INVALID.list$MARKER) ,
                                 ' ',
                                 ' '),

                    'N_TOTAL' = c(abs(study$column.NA.list$N_TOTAL - study$column.INVALID.list$N_TOTAL) ,
                                  study$column.INVALID.list$N_TOTAL,
                                  ' '),

                    'STRAND' = c(abs(study$column.NA.list$STRAND - study$column.INVALID.list$STRAND) ,
                                 study$column.INVALID.list$STRAND,
                                 ' '),

                    'CALLRATE' = c(abs(study$column.NA.list$CALLRATE - study$column.INVALID.list$CALLRATE - study$column.INVALID.list$minusone.CALLRATE),
                                   study$column.INVALID.list$CALLRATE ,
                                   study$column.INVALID.list$minusone.CALLRATE)


  ))


  colnames(b) <- c('NA values','Invalid values','Uncertain values')
  writeTXTreport(kable(b,format = "rst"))


  ## ===================================
  writeTXTreport(' ')
  writeTXTreport('=======================================================')
  writeTXTreport('================= Variant processing ==================')
  writeTXTreport('=======================================================')

  writeTXTreport('* step1: removing variants with missing crucial values and duplicated lines.')
  writeTXTreport('** step2: removing monomorphic variants and specified chromosomes.')
  writeTXTreport('*** step3: removing mismatched, ambiguous and multi-allelic variants that could not be verified.')
  writeTXTreport(' ')


  count.table <- t(data.table(
    "input variant count" = format(study$input.data.rowcount, big.mark="," , scientific = FALSE),
    'Missing crucial variable' = calculatePercent(study$missing.crucial.rowcount,
                                                  study$input.data.rowcount,
                                                  pretty = TRUE),
    'Duplicated variants' = calculatePercent(study$duplicate.count,
                                             study$input.data.rowcount,
                                             pretty = TRUE),
    "variant count after step 1 *"= calculatePercent(study$rowcount.step1,
                                                     study$input.data.rowcount,
                                                     decimal.place=3,
                                                     pretty = TRUE),
    'Monomorphic variants' = calculatePercent(study$monomorphic.count,
                                              study$input.data.rowcount,
                                              pretty = TRUE),
    "variant count after step 2 **"= calculatePercent(study$rowcount.step2,
                                                      study$input.data.rowcount,
                                                      decimal.place=3,
                                                      pretty = TRUE),
    "variant count after step 3 ***"= calculatePercent(study$rowcount.step3,
                                                       study$input.data.rowcount,
                                                       decimal.place=3,
                                                       pretty = TRUE)))


  colnames(count.table) <- 'count'



  writeTXTreport(kable(count.table,format = "rst"))



  writeTXTreport(' ')
  writeTXTreport('NOTE: All further reports are based on variants after step3 (which will be saved as output file).')
  writeTXTreport(' ')
  writeTXTreport(' ')

  ##==============================================
  writeTXTreport('==================================================')
  writeTXTreport('============ Description of variants =============')
  count.table <- t(data.table(
    'High Quality variants' = calculatePercent(study$HQ.count,
                                               study$rowcount.step3,
                                               pretty = TRUE),
    'Low Quality variants' = calculatePercent(study$LQ.count,
                                              study$rowcount.step3,
                                              pretty = TRUE),
    'Palindromic variants' = calculatePercent(study$palindromic.rows,
                                              study$rowcount.step3,
                                              pretty = TRUE),
    'Non-Palindromic variants' = calculatePercent(study$non.palindromic.rows,
                                                  study$rowcount.step3,
                                                  pretty = TRUE),
    'variants +' = calculatePercent(study$palindormicHighDiffEAF,
                                    study$palindromic.rows,
                                    pretty = TRUE),
    'variants ++' = calculatePercent(study$nonpalindormicHighDiffEAF ,
                                      study$non.palindromic.rows,
                                      pretty = TRUE),
    'variants +++' =  calculatePercent(study$palindormicExtremeDiffEAF ,
                                      study$palindromic.rows,
                                      pretty = TRUE)))

  colnames(count.table) <- 'count'



  writeTXTreport(kable(count.table,format = "rst"))

  writeTXTreport(sprintf('+ Palindromic variants with high allele frequency difference (> %s)',
                         .QC$config$filters$threshold_diffEAF))

 writeTXTreport(sprintf('++ Non-palindromic variants with high allele frequency difference (> %s)',
                         .QC$config$filters$threshold_diffEAF))

  writeTXTreport('+++ Palindromic variants with opposite allele frequency "compared to the reference" (> 0.65 for the input file and < 0.35 for the reference, or vice versa)')

  writeTXTreport(' ')

  ###
  writeTXTreport(' ')
  writeTXTreport(paste('Negative strand variants:',study$neg.strand.count))
  writeTXTreport(' ')
  writeTXTreport(paste('Allele frequency = 0 :',study$column.INVALID.list$zero.EFF_ALL_FREQ))
  writeTXTreport(' ')
  writeTXTreport(paste('Allele frequency = 1 :',study$column.INVALID.list$one.EFF_ALL_FREQ))
  writeTXTreport(' ')



  ### imputation table
  writeTXTreport('Imputation status')

  tbl = study$tables$imputed.tbl
  colnames(tbl) <- c('','Count')


  writeTXTreport(kable(tbl, align = "l",format = "rst"))

  writeTXTreport(' ')


  writeTXTreport(' ')
  writeTXTreport('========================================================')
  writeTXTreport('= Result from matching with standard reference dataset =')
  writeTXTreport('========================================================')

  ## not helpful anymore
  # match.table1 <- study$tables$match.ref.table
  # colnames(match.table1)[colnames(match.table1) == 'Std_ref'] <- 'Standard Reference'
  #
  #
  #
  # match.table <- data.table(apply(match.table1,2, function(x)
  #   return(calculatePercent(x,
  #                           study$rowcount.step2,
  #                           pretty = TRUE,
  #                           decimal.place = 3)
  #   )
  # ))
  #
  # match.table <- cbind(colnames(match.table1),match.table)
  # colnames(match.table) <- c('Reference' ,'Count')
  #
  #
  # writeTXTreport(kable(match.table,format = "rst"))

  writeTXTreport(' ')


  #writeTXTreport('Variant types after matching with reference datasets\n')
  writeTXTreport(kable(study$tables$multi_allele_count_preProcess,format = "rst"))

  writeTXTreport(' ')

  ##========================================
  # print.and.log('--------[Result from matching with standard reference file!]--------','info', cat = FALSE)

  writeTXTreport(' ')
  # writeTXTreport('========================================================')
  # writeTXTreport('= Result from matching with standard reference dataset =')
  # writeTXTreport('========================================================')

  count.table <- t(data.table(
    'Verified variants' = calculatePercent(study$found.rows.std,
                                        study$rowcount.step2,
                                        decimal.place=3,
                                        pretty=TRUE),
    'Not-found variants' = calculatePercent(study$not.found.rows.std,
                                            study$rowcount.step2,
                                            decimal.place=3,
                                            pretty=TRUE),
    # 'Mismatch variants' = calculatePercent(study$mismatched.rows.std,
    #                                        study$found.rows.std,
    #                                        decimal.place=3,
    #                                        pretty=TRUE),
    # 'Non-verified multiallelic variants' = calculatePercent(study$multiAlleleVariants.rowcount,
    #                                                         study$found.rows.std,
    #                                                         decimal.place=3,
    #                                                         pretty=TRUE),
    # 'Ambiguous variants' = calculatePercent(study$ambiguos.rows,
    #                                         study$found.rows.std,
    #                                         pretty=TRUE),
    'Flipped variants' = calculatePercent(study$flipped.rows.std,
                                          study$found.rows.std,
                                          pretty=TRUE),
    'Switched variants' = calculatePercent(study$switched.rows.std,
                                           study$found.rows.std,
                                           pretty=TRUE),
    '============================' ='==============',
    'Allele frequency correlation' = '',
    '   r (all variants)' = study$AFcor.std_ref,
    '   r (palindromic)' = study$AFcor.palindromic.std_ref,
    '   r (non-palindromic)' = study$AFcor.non.palindromic.std_ref,
    '   r (INDEL)' = study$AFcor.std_ref.indel))

  colnames(count.table) <- 'count'


  writeTXTreport(kable(count.table,format = "rst"))
  writeTXTreport(' ')

  ##=========================================

  if(!is.na(.QC$config$supplementaryFiles$allele_ref_alt))
  {
    # print.and.log('-------[Result from matching with alternate reference file!]-------','info', cat = FALSE)
    writeTXTreport(' ')
    writeTXTreport('=========================================================')
    writeTXTreport('= Result from matching with alternate reference dataset =')
    writeTXTreport('=========================================================')

    count.table <- t(data.table(
      'Verified variants' = calculatePercent(study$found.rows.alt ,
                                          study$not.found.rows.std,
                                          decimal.place=3,
                                          pretty=TRUE),
      'Not-found variants' = calculatePercent(study$not.found.rows.alt ,
                                              study$not.found.rows.std,
                                              decimal.place=3,
                                              pretty=TRUE),
      # 'Mismatch variants' = calculatePercent(study$mismatched.rows.alt ,
      #                                        study$found.rows.alt,
      #                                        decimal.place=3,
      #                                        pretty=TRUE),
      'Flipped variants' = calculatePercent(study$flipped.rows.alt ,
                                            study$found.rows.alt,
                                            pretty=TRUE),
      'Switched variants' = calculatePercent(study$switched.rows.alt ,
                                             study$found.rows.alt,
                                             pretty=TRUE),
      '============================' ='==============',
      'Allele frequency correlation' = '',
      '   r (all variants)' = study$AFcor.alt_ref,
      '   r (palindromic)' = study$AFcor.palindromic.alt_ref,
      '   r (non-palindromic)' = study$AFcor.non.palindromic.alt_ref))

    colnames(count.table) <- 'count'


    writeTXTreport(kable(count.table,format = "rst"))
    writeTXTreport(' ')

  }

  ##========================================
  writeTXTreport(' ')
  writeTXTreport('AF correlation for each chromosome')
  writeTXTreport(kable(study$AFcor.std_ref.CHR ,format = "rst",align = "c"))


  ##=========================================
  # print.and.log('-------[Calculated variables]-------','info', cat = FALSE)
  writeTXTreport(' ')
  writeTXTreport('==============================================')
  writeTXTreport('============ QC summary statistics ===========')
  writeTXTreport('==============================================')
  writeTXTreport(' ')

  writeTXTreport('Pvalue correlation (observed vs expected)')
  writeTXTreport('Note: Only variants with a valid P-value are used for P-value correlation calculation.')
  count.table <- t(data.table(
    'included variants' = calculatePercent(study$rownum.PVcor,
                                           study$rowcount.step3,
                                           pretty = TRUE),
    '   r' = study$PVcor
  ))



  colnames(count.table) <- 'value'
  writeTXTreport(kable(count.table,format = "rst"))
  writeTXTreport(' ')


  writeTXTreport(' ')

  count.table <- t(data.table(
    'Skewness' = study$skewness,
    'Skewness (HQ)' = study$skewness.HQ,
    'Kurtosis' = study$kurtosis,
    'Kurtosis (HQ)'= study$kurtosis.HQ,
    "Visscher's stat" = study$Visschers.stat ,
    "Visscher's stat (HQ)" = study$Visschers.stat.HQ,
    "Lambda - total" = study$lambda ,
    'Lambda - genotyped' = study$lambda.gen,
    'Lambda - imputed' = study$lambda.imp,
    '============================' = '==============',
    'Sample Size (Max)' = study$MAX_N_TOTAL,
    "Fixed HWE P-value" = study$fixed.hwep,
    "Fixed Imputation Quality" = study$fixed.impq,
    "Fixed Call Rate" = study$fixed.callrate,
    "Fixed Sample Size" = study$fixed.n_total
  ))

  colnames(count.table) <- 'value'

  writeTXTreport(kable(count.table,format = "rst"))
  writeTXTreport(' ')


  ##=========================================
  # print.and.log('-------[Calculated variables]-------','info', cat = FALSE)
  writeTXTreport(' ')
  writeTXTreport('==============================================')
  writeTXTreport('========== Distribution statistics  ==========')
  writeTXTreport('==============================================')
  writeTXTreport(' ')

  writeTXTreport(sprintf('All variants (%s)' , prettyNum(.QC$thisStudy$rowcount.step3,big.mark = ",")))
  writeTXTreport(kable(t(study$tables$variable.summary), format = "rst"))
  writeTXTreport(' ')

  if(nrow(study$tables$variable.summary.HQ ) > 0 & study$HQ.count != study$rowcount.step3)
  {
    writeTXTreport(sprintf('HQ variants (%s)' , prettyNum(.QC$thisStudy$HQ.count,big.mark = ",")))
    writeTXTreport(kable(t(study$tables$variable.summary.HQ), format = "rst"))
    writeTXTreport(' ')
  }


  ##========================================
  # writeTXTreport(' ')
  # writeTXTreport('==============================================')
  # writeTXTreport('============= Column statistics  =============')
  # writeTXTreport(' ')


  ### chromosome table


  if(!all(is.na(study$tables$CHR.tbl)))
  {
    writeTXTreport(' ')
    writeTXTreport('Variant count for each chromosome')
    tbl = study$tables$CHR.tbl
    colnames(tbl) <- c('Chromosome','Variant count')


    writeTXTreport(kable(tbl, align = "c",format = "rst"))

  }

  if(length(study$missing_chromosomes) >0 )
  {
    writeTXTreport(' ')
    writeTXTreport(sprintf("%s %s","Missing chromosome(s) number",paste(.QC$thisStudy$missing_chromosomes,collapse = ", ")))
  }

  writeTXTreport(' ')

  ### alleles
  writeTXTreport(' ')
  writeTXTreport('Effect allele distribution in SNP variants')

  tbl = merge(study$tables$EFFECT_ALL.tbl,
              study$tables$EFFECT_ALL.post.matching.tbl,
              by="EFFECT_ALL",
              all = TRUE)
  tbl = t(tbl)

  rownames(tbl) <- c('Allele','Count (input file)','Count (post-matching)')
  colnames(tbl) <- tbl[1,]
  writeTXTreport(kable(tbl[-1,], align = "c",format = "rst"))

  writeTXTreport(' ')


  ###
  writeTXTreport(' ')
  writeTXTreport('Other allele distribution in SNP variants')
  tbl = merge(study$tables$OTHER_ALL.tbl,
              study$tables$OTHER_ALL.post.matching.tbl,
              by="OTHER_ALL",
              all = TRUE)
  tbl = t(tbl)

  rownames(tbl) <- c('Allele','Count (input file)','Count (post-matching)')
  colnames(tbl) <- tbl[1,]
  writeTXTreport(kable(tbl[-1,], align = "c",format = "rst"))

  ##

  ## END OF REPORT
  # =============
  print.and.log(sprintf('Report file saved as \'%s\'',study$txt.report.path),
                'info')
}



# save each study object as rdata file
# to compare different files after each is run separately
save.rds.file <- function(study) {


  # rm(list=setdiff(ls(envir =  study$effect.plot$plot_env),
  #                 c('y_lower','y_upper','df','file.N.max','file.number')),
  #    envir =  study$effect.plot$plot_env)
  # #
  # rm('.QC' , envir =  study$effect.plot$plot_env)
  # rm('study' , envir =  study$effect.plot$plot_env)
  # rm('input.data' , envir =  study$effect.plot$plot_env)


  tryCatch(
    {
      if(.QC$config$output_parameters$object_file)
        saveRDS(object = study, file =  study$rds.study.rds.path, version = '2')
    },
    error = function(err)
    {
      print.and.log(paste('Could not save study RDS object file:',err$message),'warning',display=.QC$config$debug$verbose)
    }
  )
}
