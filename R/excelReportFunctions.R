writeExcelReportFile <- function()
{


  report.success <- tryCatch(create_xlsx_report(.QC$config,.QC$qc.study.list),
                             error = function(err)  {
                               print_and_log(paste0('Excel report could not be saved. ', err$message),'warning')
                               return(FALSE)
                             }
  )

  if(report.success)
    print_and_log(sprintf('Excel report file saved at: \'%s\'', .QC$config$paths$xlsx.report),'info')
}


create_xlsx_report <- function(config,study.list){


  ##### styles ####
  style1 <- openxlsx::createStyle(
    fontSize = 16,
    textDecoration = "bold",
    halign = "left",
    fontColour = "darkblue"
  )

  style2 <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center"
  )

  style2_left <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "left"
  )

  style3 <- openxlsx::createStyle(
    textDecoration = "bold",
    fontColour="darkgreen"
  )

  style4 <- openxlsx::createStyle(
    fontSize = 14,
    textDecoration = c("bold", "italic"),
    fontColour="darkred"
  )

  style5 <- openxlsx::createStyle(
    textDecoration = c("bold"),
    border = c('top','bottom'),
    borderStyle = c('thin','medium'),
    halign = "center"
  )


  styleList <- list(
    'style1'=style1,
    'style2'=style2,
    'style2_left'=style2_left,
    'style3'=style3,
    'style4'=style4,
    'style5'=style5
  )
  #######

  ## initialize excel file
  wb <- openxlsx::createWorkbook(
    creator = 'GWASinspector_package',
    title = 'QC output',
    subject = 'Quality Check'
  )


  ## first sheet
  first_excel_sheet(wb,styleList,config,study.list)

  ## second sheet
  second_excel_sheet(wb,styleList,config,study.list)

  ## third sheet
  third_excel_sheet(wb,styleList,config,study.list)

  ## sheet per file
  for(i in 1:length(study.list))
  {
    file_excel_sheet(wb,styleList,i,config,study.list)
  }

  ## save output file
  openxlsx::saveWorkbook(
    wb = wb,
    file = config$paths$xlsx.report,
    overwrite = TRUE
  )

  return(TRUE)

}

first_excel_sheet <- function(wb,styleList,config,study.list) {

  sheet <- openxlsx::addWorksheet(wb, sheetName = "Report variables")

  #### title
  openxlsx::writeData(wb ,sheet ,x = paste0(.QC$package.name ," Report v.", .QC$script.version),startRow = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style1']],rows = 1,cols = 1)
  openxlsx::setColWidths(wb,sheet,cols = c(1,2),c(50,100))

  #### intro table
  tbl <- t(data.table(
    format(config$new_items$starttime, "%b %d %Y - %X"),
    format(config$new_items$endtime, "%b %d %Y - %X"),
    .QC$script.version,
    .QC$r.version,
    config$paths$dir_data,
    config$paths$dir_output,
    config$paths$dir_references,
    basename(config$supplementaryFiles$allele_ref_std),
    basename(config$supplementaryFiles$header_translations)))

  tbl<- cbind(c('Start time',
                'End time',
                'Script version',
                'System Information',
                'Input directory',
                'Output directory',
                'References directory',
                'Allele frequency Reference dataset' ,
                'Alternate header file'), tbl)

  if(!is.na(config$supplementaryFiles$allele_ref_alt)){
    tbl <- rbind(tbl ,c( 'Alternate Allele frequency Reference dataset ',
                         data.table(basename(config$supplementaryFiles$allele_ref_alt))))
  }

  if(!is.na(config$supplementaryFiles$beta_ref_std)){
    tbl <- rbind(tbl ,c( 'Effect size Reference dataset',
                         data.table(basename(config$supplementaryFiles$beta_ref_std))))
  }

  openxlsx::writeData(wb ,sheet ,x = tbl,startRow = 3,colNames = FALSE)

}

second_excel_sheet <- function(wb,styleList,config,study.list) {


  sheet <- openxlsx::addWorksheet(wb, sheetName = "File Names")

  ## title
  openxlsx::writeData(wb ,sheet ,x = "Input File Names",startRow = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style1']],rows = 1,cols = 1)

  ## table
  report.table <- data.table('file_number'=sapply(study.list, function(x) return(x$number)))
  report.table$file_name <-data.table(sapply(study.list, function(x) return(basename(x$file.path))))
  colnames(report.table) <- c('File number','File Name')


  openxlsx::writeData(
    wb,
    sheet,
    x = report.table,
    startRow = 3,
    rowNames = FALSE,
    colNames = FALSE
  )
  openxlsx::addStyle(wb,sheet,style = styleList[['style2']],rows = seq(3,3+nrow(report.table)),cols = 1)

  # col width
  openxlsx::setColWidths(wb,sheet,cols = c(1,2),c(30,60))

}

third_excel_sheet <- function(wb,styleList,config,study.list) {


  sheet <- openxlsx::addWorksheet(wb, sheetName = "Multi File Report")

  ## header
  openxlsx::writeData(wb ,sheet ,x = "Comparing input files",xy = c(1,1))
  openxlsx::addStyle(wb,sheet,style = styleList[['style1']],rows = 1,cols = 1)


  ## table
  report.table <- t(data.table(
    sapply(study.list, function(x) return(x$MAX_N_TOTAL)),
    sapply(study.list, function(x) return(paste(x$missing.Columns, collapse = ' | '))),
    sapply(study.list, function(x) return(format(x$input.data.rowcount, big.mark="," , scientific = FALSE))),
    sapply(study.list, function(x) return(format(x$rowcount.step1, big.mark="," , scientific = FALSE))),
    sapply(study.list, function(x) return(format(x$rowcount.step2, big.mark="," , scientific = FALSE))),
    sapply(study.list, function(x) return(format(x$rowcount.step3, big.mark="," , scientific = FALSE))),
    sapply(study.list, function(x) return(calculatePercent(sum(as.numeric(gsub(x=x$tables$multi_allele_count_preProcess[1:3],
                                                                               pattern = ",",
                                                                               replacement = ""))),x$rowcount.step3,pretty=TRUE))),
    sapply(study.list, function(x) return(calculatePercent(sum(as.numeric(gsub(x=x$tables$multi_allele_count_preProcess[4:6],
                                                                               pattern = ",",
                                                                               replacement = "")))
                                                           ,x$rowcount.step3,
                                                           pretty=TRUE))),
    sapply(study.list, function(x) return(calculatePercent(x$monomorphic.count,x$rowcount.step3,pretty=TRUE))),
    #sapply(study.list, function(x) return(calculatePercent(x$duplicate.count,x$rowcount.step3,pretty=TRUE))),
    sapply(study.list, function(x) return(calculatePercent(x$palindromic.rows,x$rowcount.step3,pretty=TRUE))),
    sapply(study.list, function(x) return(ifelse(x$tables$imputed.tbl[IMPUTED=='genotyped',.N] == 0,
                                                 "NA",
                                                 calculatePercent(as.numeric(x$tables$imputed.tbl[IMPUTED=="genotyped",N]),x$rowcount.step3,pretty=TRUE)))),
    sapply(study.list, function(x) return(ifelse(x$tables$imputed.tbl[IMPUTED=='imputed',.N] == 0,
                                                 "NA",
                                                 calculatePercent(as.numeric(x$tables$imputed.tbl[IMPUTED=="imputed",N]),x$rowcount.step3,pretty=TRUE)))),
    sapply(study.list, function(x) return(calculatePercent(x$neg.strand.count,x$rowcount.step3,pretty=TRUE))),
    sapply(study.list, function(x) return(x$AFcor.std_ref)),
    sapply(study.list, function(x) return(x$AFcor.palindromic.std_ref)),
    sapply(study.list, function(x) return(x$AFcor.alt_ref)),
    sapply(study.list, function(x) return(x$AFcor.palindromic.alt_ref)),
    sapply(study.list, function(x) return(x$lambda)),
    sapply(study.list, function(x) return(x$lambda.gen)),
    sapply(study.list, function(x) return(x$lambda.imp)),
    sapply(study.list, function(x) return(x$PVcor)),
    sapply(study.list, function(x) return(x$Visschers.stat.HQ)),
    sapply(study.list, function(x) return(x$fixed.hwep)),
    sapply(study.list, function(x) return(x$fixed.impq)),
    sapply(study.list, function(x) return(x$fixed.n_total)),
    sapply(study.list, function(x) return(x$fixed.callrate)),
    rep(' ', length(study.list)),
    sapply(study.list, function(x) return(x$tables$variable.summary['Min.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['1st Qu.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['Median', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['Mean', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['3rd Qu.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['Max.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(ifelse(nrow(x$tables$variable.summary.HQ) > 0,
                                                 x$tables$variable.summary.HQ['Min.', .QC$config$input_parameters$effect_type_string],
                                                 NA))),
    sapply(study.list, function(x) return(ifelse(nrow(x$tables$variable.summary.HQ) > 0,
                                                 x$tables$variable.summary.HQ['Max.', .QC$config$input_parameters$effect_type_string],
                                                 NA))),
    sapply(study.list, function(x) return(x$tables$variable.summary['Median','STDERR'])),
    sapply(study.list, function(x) return(ifelse(nrow(x$tables$variable.summary.HQ) > 0,
                                                 x$tables$variable.summary.HQ['Median','STDERR'],
                                                 NA)))

  ))



  #report.table <- rbind(sapply(.QC$qc.study.list, function(x) return(x$number)),report.table)

  row.names(report.table) <- c(#"File number",
    "Sample Size (Max)",
    "Missing Columns",
    "SNPs in input file",
    "Variant count after step 1 *",
    "Variant count after step 2 **",
    "Variant count after step 3 ***",
    "SNP variants",
    "Non-SNP variants",
    "Monomorphic",
    #/	"Duplicates",
    "Palindromics",
    "Genotyped variants",
    "Imputed variants",
    "Negative-strand SNPs",
    "Allele Frequency Correlation (Standard Ref)",
    "Palindromic Allele Frequency correlation (Standard Ref)",
    "Allele Frequency Correlation (Alternative Ref)",
    "Palindromic Allele Frequency correlation (Alternative Ref)",
    "Lambda - Total",
    "Lambda - Genotyped",
    "Lambda - Imputed",
    "P-value Correlation",
    "Visscher's Statistic (HQ variants)",
    "Fixed HWE P-value",
    "Fixed Imputation Quality",
    "Fixed Sample Size",
    "Fixed Call Rate",
    .QC$config$input_parameters$effect_type_string,
    "-      Min.",
    "-      1st Qu.",
    "-      Median",
    "-      Mean",
    "-      3rd Qu.",
    "-      Max.",
    "-      Min. (HQ variants)",
    "-      Max. (HQ variants)",
    "Standard Error (median)",
    "Standard Error (median) (HQ variants)")



  colnames(report.table) <- sapply(.QC$qc.study.list, function(x) return(x$number))


  ## write table
  openxlsx::writeData(
    wb,
    sheet,
    x = report.table,
    startRow = 3,
    rowNames = TRUE,
    colNames = TRUE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = seq(3,3+nrow(report.table)),cols = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style2']],rows = 3,cols = seq(2,2+nrow(report.table)))


  ###


  openxlsx::writeData(
    wb,
    sheet,
    x = "step1: removing variants with missing crucial values and duplicated lines.",
    startRow = 42
  )

  openxlsx::writeData(
    wb,
    sheet,
    x = "step2: removing monomorphic variants and specified chromosomes.",
    startRow = 43
  )

  openxlsx::writeData(
    wb,
    sheet,
    x = "step3: removing mismatched, ambiguous and multi-allelic variants that could not be verified.",
    startRow = 44
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style3']], cols = 1,rows = 42:44)

  # col width
  openxlsx::setColWidths(wb,sheet,cols = seq(1,(length(study.list) + 1)),c(70,rep(30,times=length(study.list))))

}

file_excel_sheet <- function(wb,styleList,i,config,study.list) {

  sheet <- openxlsx::addWorksheet(wb, sheetName = paste0('File ',i))
  study <- study.list[[i]]
  row.index <- 3

  ## title
  openxlsx::writeData(wb ,sheet ,x = basename(study$file.path),startRow = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style1']],rows = 1,cols = 1)

  ## table

  # 1 introduction
  tbl <- t(data.table(
    format(study$starttime, "%b %d %Y - %X"),
    format(study$endtime, "%b %d %Y - %X"),
    basename(study$file.path),
    study$file.line.count,
    study$file.endsWithNewLine))

  tbl<- cbind(c('Start time',
                'End time',
                'Input File Name',
                'Input File Line Count (including header)',
                'Input File ends with a new line'), tbl)


  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = FALSE
  )
  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(3:8),cols = 1)

  # col width
  openxlsx::setColWidths(wb,sheet,cols = 1:8,widths = c(40,rep(25,times=7)))

  # 2- filters
  filter.table <- data.table(
    "Allele frequency" =  format(.QC$config$filters$HQfilter_FRQ,scientific = FALSE))

  if("HWE_PVAL" %in% study$renamed.File.Columns.sorted)
    filter.table <- cbind(filter.table, "HWE p-value" = format(.QC$config$filters$HQfilter_HWE,scientific = FALSE)) else
      filter.table <- cbind(filter.table, "HWE p-value" = "Not included")

  if("CALLRATE" %in% study$renamed.File.Columns.sorted)
    filter.table <- cbind(filter.table, "Call-rate" = format(.QC$config$filters$HQfilter_cal,scientific = FALSE)) else
      filter.table <- cbind(filter.table, "Call-rate" = "Not included")


  if("IMP_QUALITY" %in% study$renamed.File.Columns.sorted)
    filter.table <- cbind(filter.table, "Imputation quality" = format(.QC$config$filters$HQfilter_imp,scientific = FALSE)) else
      filter.table <- cbind(filter.table, "Imputation quality" = "Not included")

  filter.table <- t(filter.table)

  colnames(filter.table) <- 'Value'


  row.index <- row.index + 6 #
  openxlsx::writeData(wb ,sheet ,x = "Filter values for selecting High-Quality (HQ) variants" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  row.index <- row.index + 1 #
  openxlsx::writeData(
    wb,
    sheet,
    x = filter.table,
    startRow = row.index,
    rowNames = TRUE,
    colNames = FALSE
  )
  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(filter.table))),cols = 1)


  ## 3 column names

  row.index <- row.index + 6 #
  openxlsx::writeData(wb ,sheet ,x = "Column names and translation" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)


  column.tbl <- rbind(study$original.File.Columns.sorted,
                      study$renamed.File.Columns.sorted)
  rownames(column.tbl) <- c('Original', 'Renamed')
  column.tbl <- t(column.tbl)
  row.index <- row.index + 1 #

  openxlsx::writeData(
    wb,
    sheet,
    x = column.tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = TRUE
  )
  openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 1:2)


  ## column report
  row.index <- row.index + nrow(column.tbl) + 2
  openxlsx::writeData(wb ,sheet ,x = "Column report" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  tbl <- t(data.frame('CHR' = c(abs(study$column.NA.list$CHR - study$column.INVALID.list$CHR) ,
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
                                   # study$column.INVALID.list$MARKER,
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



  colnames(tbl) <- c('NA values','Invalid values','Uncertain values')

  row.index <- row.index + 1 #

  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = TRUE,
    colNames = TRUE
  )
  openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 2:4)
  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)


  ## 4 variant processing
  row.index <- row.index + 19 #
  openxlsx::writeData(wb ,sheet ,x = "Variant processing" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  row.index <- row.index + 1 #
  openxlsx::writeData(wb ,sheet ,x = "step1: removing variants with missing crucial values and duplicated lines." ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)

  row.index <- row.index + 1 #
  openxlsx::writeData(wb ,sheet ,x = "step2: removing monomorphic variants and specified chromosomes." ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)

  row.index <- row.index + 1 #
  openxlsx::writeData(wb ,sheet ,x = "step3: removing mismatched, ambiguous and multi-allelic variants that could not be verified." ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)



  tbl <- t(data.table(
    format(study$input.data.rowcount, big.mark="," , scientific = FALSE),
    calculatePercent(study$duplicate.count,
                     study$input.data.rowcount,
                     pretty = TRUE),
    calculatePercent(study$missing.crucial.rowcount,
                     study$input.data.rowcount,
                     pretty = TRUE),
    calculatePercent(study$rowcount.step1,
                     study$input.data.rowcount,
                     decimal.place=3,
                     pretty = TRUE),
    calculatePercent(study$monomorphic.count,
                     study$input.data.rowcount,
                     pretty = TRUE),
    calculatePercent(study$rowcount.step2,
                     study$input.data.rowcount,
                     decimal.place=3,
                     pretty = TRUE),
    calculatePercent(study$rowcount.step3,
                     study$input.data.rowcount,
                     decimal.place=3,
                     pretty = TRUE)
  ))


  tbl<- cbind(c(   "input variant count",
                   'Duplicated variants',
                   'Missing crucial variable',
                   "variant count after step 1 *",
                   'Monomorphic variants',
                   "variant count after step 2 **",
                   "variant count after step 3 ***"),tbl)

  row.index <- row.index + 2 # 13
  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = FALSE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)

  ## description of variant
  row.index <- row.index + 10 # 8

  openxlsx::writeData(wb ,sheet ,x = "Description of variants" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)


  tbl <- t(data.table(
    calculatePercent(study$HQ.count,
                     study$rowcount.step3,
                     pretty = TRUE),
    calculatePercent(study$LQ.count,
                     study$rowcount.step3,
                     pretty = TRUE),
    calculatePercent(study$palindromic.rows,
                     study$rowcount.step3,
                     pretty = TRUE),
    calculatePercent(study$non.palindromic.rows,
                     study$rowcount.step3,
                     pretty = TRUE),
    calculatePercent(study$palindormicHighDiffEAF,
                     study$palindromic.rows,
                     pretty = TRUE),
    calculatePercent(study$nonpalindormicHighDiffEAF ,
                     study$non.palindromic.rows,
                     pretty = TRUE),
    calculatePercent(study$palindormicExtremeDiffEAF ,
                     study$palindromic.rows,
                     pretty = TRUE)
  ))


  tbl<- cbind(c( 'High Quality variants',
                 'Low Quality variants',
                 'Palindromic variants',
                 'Non-Palindromic variants',
                 'variants +',
                 'variants ++',
                 'variants +++'),tbl)

  row.index <- row.index + 1 # 23

  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = FALSE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)



  row.index <- row.index + 8 # 31
  openxlsx::writeData(wb ,sheet ,x = sprintf('+ palindromic variants with high allele frequency difference (> %s)',
                                             config$filters$threshold_diffEAF) ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)


  row.index <- row.index + 1 # 32
  openxlsx::writeData(wb ,sheet ,x = sprintf('++ Non-palindromic variants with high allele frequency difference (> %s)',
                                             config$filters$threshold_diffEAF) ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)


  row.index <- row.index + 1 # 33
  openxlsx::writeData(wb ,sheet ,x = '+++ palindromic variants with opposite allele frequency "compared to the reference" (> 0.65 for the input file and < 0.35 for the reference, or vice versa)' ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)

  ## indels
  if(nrow(study$tables$multi_allele_count_preProcess) > 1)
  {
    row.index <- row.index + 3 # 36
    openxlsx::writeData(wb ,sheet ,x = "Variant types" ,startRow = row.index)
    openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)


    row.index <- row.index + 1 # 37

    openxlsx::writeData(
      wb,
      sheet,
      x = study$tables$multi_allele_count_preProcess,
      startRow = row.index,
      rowNames = TRUE,
      colNames = TRUE
    )

    openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)

    ## matching with reference panels

    row.index <- row.index + 8 # 36

    openxlsx::writeData(wb ,sheet ,x = "Result from matching with standard reference file" ,startRow = row.index)
    openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)


    tbl <- t(data.table(
      calculatePercent(study$found.rows.std,
                       study$rowcount.step2,
                       decimal.place=3,
                       pretty=TRUE),
      calculatePercent(study$not.found.rows.std,
                       study$rowcount.step2,
                       decimal.place=3,
                       pretty=TRUE),
      # calculatePercent(study$mismatched.rows.std,
      #                  study$found.rows.std,
      #                  decimal.place=3,
      #                  pretty=TRUE),
      # calculatePercent(study$multiAlleleVariants.rowcount,
      #                  study$found.rows.std,
      #                  decimal.place=3,
      #                  pretty=TRUE),
      # calculatePercent(study$ambiguos.rows,
      #                  study$found.rows.std,
      #                  pretty=TRUE),
      calculatePercent(study$flipped.rows.std,
                       study$found.rows.std,
                       pretty=TRUE),
      calculatePercent(study$switched.rows.std,
                       study$found.rows.std,
                       pretty=TRUE),
      study$AFcor.std_ref,
      study$AFcor.palindromic.std_ref,
      study$AFcor.non.palindromic.std_ref,
      study$AFcor.std_ref.indel
    ))

    tbl <- cbind(c(
      'Verified variants',
      'Not-found variants',
      # 'Mismatch variants',
      # 'Non-verified multiallelic variants',
      # 'Ambiguous variants',
      'Flipped variants',
      'Switched variants',
      'Allele frequency correlation (all variants)',
      'Allele frequency correlation (palindromic variants)',
      'Allele frequency correlation (non-palindromic variants)',
      'Allele frequency correlation (INDEL)'
    ) , tbl)


    row.index <- row.index + 1 # 37
    openxlsx::writeData(
      wb,
      sheet,
      x = tbl,
      startRow = row.index,
      rowNames = FALSE,
      colNames = FALSE
    )

    openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)

    #

    if(!is.na(config$supplementaryFiles$allele_ref_alt))
    {
      row.index <- row.index + 10 # 47
      openxlsx::writeData(wb ,sheet ,x = "Result from matching with alternative reference file" ,startRow = row.index)
      openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)


      tbl <- t(data.table(
        calculatePercent(study$found.rows.alt,
                         study$not.found.rows.std,
                         decimal.place=3,
                         pretty=TRUE),
        calculatePercent(study$not.found.rows.alt,
                         study$not.found.rows.std,
                         decimal.place=3,
                         pretty=TRUE),
        # calculatePercent(study$mismatched.rows.alt,
        #                  study$found.rows.alt,
        #                  decimal.place=3,
        #                  pretty=TRUE),
        calculatePercent(study$flipped.rows.alt,
                         study$found.rows.alt,
                         pretty=TRUE),
        calculatePercent(study$switched.rows.alt,
                         study$found.rows.alt,
                         pretty=TRUE),
        study$AFcor.alt_ref,
        study$AFcor.palindromic.alt_ref,
        study$AFcor.non.palindromic.alt_ref
      ))

      tbl <- cbind(c(
        'Verified variants',
        'Not-found variants',
        # 'Mismatch variants',
        'Flipped variants',
        'Switched variants',
        'Allele frequency correlation (all variants)',
        'Allele frequency correlation (palindromic variants)',
        'Allele frequency correlation (non-palindromic variants)'
      ) , tbl)


      row.index <- row.index + 1 # 48
      openxlsx::writeData(
        wb,
        sheet,
        x = tbl,
        startRow = row.index,
        rowNames = FALSE,
        colNames = FALSE
      )

      openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)

    }
    #

  }

  ## AF correlation per chromosme
  row.index <- row.index + 11 # 58

  if( any(is.na(study$AFcor.std_ref.CHR)) && nrow(study$AFcor.std_ref.CHR) > 0)
  {
    openxlsx::writeData(wb ,sheet ,x = "AF correlation for each chromosome" ,startRow = row.index)
    openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

    row.index <- row.index + 1 # 48

    openxlsx::writeData(
      wb,
      sheet,
      x = study$AFcor.std_ref.CHR,
      startRow = row.index,
      rowNames = FALSE,
      colNames = TRUE
    )

    openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(study$AFcor.std_ref.CHR))),cols = 1)
    openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 1:2)

    row.index <- row.index + nrow(study$AFcor.std_ref.CHR) + 2
  }

  ## summary statistics

  openxlsx::writeData(wb ,sheet ,x = "QC summary statistics" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  tbl <- t(data.table(
    calculatePercent(study$rownum.PVcor,
                     study$rowcount.step3,
                     pretty = TRUE),
    study$PVcor,
    study$skewness,
    study$skewness.HQ,
    study$kurtosis,
    study$kurtosis.HQ,
    study$Visschers.stat,
    study$Visschers.stat.HQ,
    study$lambda,
    study$lambda.gen,
    study$lambda.imp,
    study$MAX_N_TOTAL,
    study$fixed.hwep,
    study$fixed.impq,
    study$fixed.callrate,
    study$fixed.n_total
  ))

  tbl <- cbind(c('Number of variants for P-value correlation',
                 'P-value correlation (all)',
                 'Skewness',
                 'Skewness (HQ)',
                 'Kurtosis',
                 'Kurtosis (HQ)',
                 "Visscher's stat",
                 "Visscher's stat (HQ)",
                 "Lambda - total",
                 'Lambda - genotyped',
                 'Lambda - imputed',
                 'Sample Size (Max)',
                 "Fixed HWE P-value",
                 "Fixed Imputation Quality",
                 "Fixed Call Rate",
                 "Fixed Sample Size"),tbl)


  row.index <- row.index + 1 # 59

  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = FALSE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)


  # variable summary statistics

  row.index <- row.index + 19 # 78
  openxlsx::writeData(wb ,sheet ,x = "Distribution statistics" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)


  row.index <- row.index + 2 # 80

  openxlsx::writeData(wb ,sheet ,x = "All variants" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)

  row.index <- row.index + 1

  row.names(study$tables$variable.summary) <- c("min.","first_quartile","median","mean","third_quartile","max." )

  openxlsx::writeData(
    wb,
    sheet,
    x = t(study$tables$variable.summary),
    startRow = row.index,
    rowNames = TRUE,
    colNames = TRUE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 2:7)

  ##
  if(nrow(study$tables$variable.summary.HQ ) > 0 & study$HQ.count != study$rowcount.step3)
  {
    row.index <- row.index + 10

    openxlsx::writeData(wb ,sheet ,x = "HQ variants only" ,startRow = row.index)
    openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)

    row.index <- row.index + 1
    row.names(study$tables$variable.summary.HQ) <- c("min.","first_quartile","median","mean","third_quartile","max." )
    openxlsx::writeData(
      wb,
      sheet,
      x = t(study$tables$variable.summary.HQ),
      startRow = row.index,
      rowNames = FALSE,
      colNames = FALSE
    )

    openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)
    openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 2:7)
  }
  # chromosome summary statistics

  row.index <- row.index + 10
  openxlsx::writeData(wb ,sheet ,x = "Variant count for each chromosome" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)


  ####
  row.index <- row.index + 1 # 89
  chr.tbl.length <- 0
  if(!is.null(study$tables$CHR.tbl))
  {
    if(ncol(study$tables$CHR.tbl) == 2)
    {

      chr.tbl <-  study$tables$CHR.tbl
      colnames(chr.tbl) <- c('Chromosome Number','Variant Count')
      chr.tbl.length <- nrow(chr.tbl)

      openxlsx::writeData(
        wb,
        sheet,
        x = chr.tbl,
        startRow = row.index,
        rowNames = FALSE,
        colNames = TRUE
      )

      openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+ chr.tbl.length)),cols = 1)
      openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 1:2)
    }
  }

  row.index <- row.index + chr.tbl.length + 2

  if(length(study$missing_chromosomes) >0 )
  {

    openxlsx::writeData(wb, sheet,
                        x = sprintf("%s %s", "Missing chromosome(s) number", paste(study$missing_chromosomes,collapse = ", ")),
                        startRow = row.index)

    openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

    row.index <- row.index + 2
  }


  ###############
  ###
  openxlsx::writeData(wb ,sheet ,x = "Effect allele distribution in SNP variants" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  row.index <- row.index + 1 #

  tbl = merge(study$tables$EFFECT_ALL.tbl,
              study$tables$EFFECT_ALL.post.matching.tbl,
              by="EFFECT_ALL",
              all = TRUE)

  colnames(tbl) <- c('Allele','Count (input file)','Count (post-matching)')

  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = TRUE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 1:3)


  ###
  row.index <- row.index + 6 #
  openxlsx::writeData(wb ,sheet ,x = "Other allele distribution in SNP variants" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  row.index <- row.index + 1 #

  tbl = merge(study$tables$OTHER_ALL.tbl,
              study$tables$OTHER_ALL.post.matching.tbl,
              by="OTHER_ALL",
              all = TRUE)

  colnames(tbl) <- c('Allele','Count (input file)','Count (post-matching)')

  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = TRUE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 1:3)

  ##


  row.index <- row.index + 6 #
  openxlsx::writeData(wb ,sheet ,x = "Imputation status" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  row.index <- row.index + 1 #
  tbl = study$tables$imputed.tbl
  # tbl$IMPUTED <- c('Genotyped','Imputed')
  colnames(tbl) <- c('Status','Count')

  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = TRUE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(tbl))),cols = 1)
  openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 1:2)

  ##


  row.index <- row.index + 5

  openxlsx::writeData(wb ,sheet ,x = "Negative strand variants" ,startRow = row.index)
  openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

  tbl <- t(data.frame('Negative strand variants' = study$neg.strand.count))
  colnames(tbl) <- c('Count')

  row.index <- row.index + 1
  openxlsx::writeData(
    wb,
    sheet,
    x = tbl,
    startRow = row.index,
    rowNames = FALSE,
    colNames = TRUE
  )

  openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = row.index:(row.index+1) ,cols = 1)



  # effect size comparison

  if(!is.null(study$tables$betaCor.tbl)  & !is.na(config$supplementaryFiles$beta_ref_std)){

    row.index <- row.index + 3 # 88

    openxlsx::writeData(wb ,sheet ,x = "Effect-size comparison" ,startRow = row.index)
    openxlsx::addStyle(wb,sheet,style = styleList[['style4']],rows = row.index,cols = 1)

    row.index <- row.index + 1
    openxlsx::writeData(
      wb,
      sheet,
      x = study$tables$betaCor.tbl ,
      startRow = row.index,
      rowNames = TRUE,
      colNames = TRUE
    )

    openxlsx::addStyle(wb,sheet,style = styleList[['style2_left']],rows = c(row.index:(row.index+nrow(study$tables$betaCor.tbl))),cols = 1)
    openxlsx::addStyle(wb,sheet,style = styleList[['style5']],rows = row.index,cols = 1:3)



    row.index <- row.index + 6 # 9
    openxlsx::writeData(wb ,sheet ,x = "* Data is presented as r(N). Variants were filtered on reference data P-values." ,startRow = row.index)
    openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)


    row.index <- row.index + 1
    openxlsx::writeData(wb ,sheet ,x = "** Data is presented as r(N). Variants were filtered on input result file P-values." ,startRow = row.index)
    openxlsx::addStyle(wb,sheet,style = styleList[['style3']],rows = row.index,cols = 1)

  }

}
