xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
  rows <-xlsx::createRow(sheet,rowIndex=rowIndex)
  sheetTitle <-xlsx::createCell(rows, colIndex=1)
  xlsx::setCellValue(sheetTitle[[1,1]], title)
  xlsx::setCellStyle(sheetTitle[[1,1]], titleStyle)
}


writeExcelReportFile <- function()
{

  # check if xlsx function is installed.
  # this package is removed from dependencies, because iot requires JAVA to be installed priorly
  # which is not always available
  if (!.QC$xlsx.package.exists) {
    print.and.log("Package \"xlsx\" is needed for creating excel reports. Please install it.",'warning',display=.QC$config$debug$verbose)
    return(NULL)
  }


  report.success <- tryCatch(create.xlsx.report(.QC$config,.QC$qc.study.list),
                             error = function(err)  {
                               print.and.log(paste0('Excel report could not be saved. ', err$message),'warning')
                               return(FALSE)
                             }
  )

  if(report.success)
    print.and.log(sprintf('Excel report file saved at: \'%s\'', .QC$config$paths$xlsx.report),'info')
}


create.xlsx.report <- function(config,study.list){


  # if(file.exists(config$paths$xlsx.report))
  #   file.remove(config$paths$xlsx.report)



  wb<-xlsx::createWorkbook(type="xlsx")


  # STYLES
  #=====================================
  # TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16,color="darkblue", isBold=TRUE, underline=1)
  TITLE_STYLE <- xlsx::CellStyle(wb)+ xlsx::Font(wb,  heightInPoints=16,
                                                 color="darkblue", isBold=TRUE)

  SUB_TITLE_STYLE <- xlsx::CellStyle(wb) +
    xlsx::Font(wb,  heightInPoints=14, color="darkred",
               isItalic=TRUE, isBold=FALSE)

  NOTE_TITLE_STYLE <- xlsx::CellStyle(wb) + xlsx::Font(wb,isBold=TRUE)

  NOTE_TITLE_STYLE2 <- xlsx::CellStyle(wb) + xlsx::Font(wb,isBold=TRUE,color="darkgreen")



  # Styles for the data table row/column names
  TABLE_ROWNAMES_STYLE <- xlsx::CellStyle(wb,alignment = xlsx::Alignment(horizontal = "ALIGN_CENTER")) +
    xlsx::Font(wb, isBold=TRUE)
  TABLE_ROWNAMES_STYLE_Left <- xlsx::CellStyle(wb,alignment = xlsx::Alignment(horizontal = "ALIGN_LEFT")) +
    xlsx::Font(wb, isBold=TRUE)

  TABLE_COLNAMES_STYLE <- xlsx::CellStyle(wb) + xlsx::Font(wb, isBold=TRUE) +
    xlsx::Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
    xlsx::Border(color="black", position=c("TOP", "BOTTOM"),
                 pen=c("BORDER_THIN", "BORDER_THICK"))


  xlsx::CellStyle(wb, dataFormat=NULL, alignment=NULL,
                  border=NULL, fill=NULL, font=NULL)

  #===========================================


  # SHEET 1
  # =====================
  sheet <- xlsx::createSheet(wb, sheetName = "Report variables")

  xlsx.addTitle(sheet, rowIndex=1, title=paste0(.QC$package.name ," Report v.", .QC$script.version),
                titleStyle = TITLE_STYLE)
  # xlsx.addTitle(sheet, rowIndex=8,
  #               title=sprintf('Alternate Allele frequency Reference set: %s', basename(config$supplementaryFiles$allele_ref_alt)),
  #               titleStyle = SUB_TITLE_STYLE)

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


  xlsx::addDataFrame(tbl, sheet, startRow=3, startColumn=1,
                     colnamesStyle = TABLE_COLNAMES_STYLE,row.names = FALSE,
                     rownamesStyle = TABLE_ROWNAMES_STYLE,col.names = FALSE)

  #

########removed to each file sheet
#
#   tbl <- t(data.table(
#     format(config$filters$HQfilter_FRQ,
#            scientific = FALSE),
#     format(config$filters$HQfilter_HWE,
#            scientific = FALSE),
#     format(config$filters$HQfilter_cal,
#            scientific = FALSE),
#     format(config$filters$HQfilter_imp,
#            scientific = FALSE)
#   ))
#
#
#   tbl <- cbind(c(
#     'Allele frequency',
#     'HWE p-value',
#     'Call-rate',
#     'Imputation quality'),tbl)
#
#   xlsx.addTitle(sheet, rowIndex=18, title="High Quality variant filter parameters",
#                 titleStyle = SUB_TITLE_STYLE)
#
#   xlsx::addDataFrame(tbl, sheet, startRow=19, startColumn=1,
#                      colnamesStyle = TABLE_COLNAMES_STYLE,row.names = FALSE,
#                      rownamesStyle = TABLE_ROWNAMES_STYLE,col.names = FALSE)



  xlsx::setColumnWidth(sheet, colIndex=1, colWidth=50)
  xlsx::setColumnWidth(sheet, colIndex=2, colWidth=100)
  #==================

  #========sheet 2 (file names)==========

  sheet2 <- xlsx::createSheet(wb, sheetName = "File Names")


  report.table <- data.table(sapply(study.list, function(x) return(basename(x$file.path))))
  colnames(report.table) <- c('Input File Name')
  # row.names(report.table) <- seq(1:nrow(report.table))
  row.names(report.table) <-sapply(study.list, function(x) return(x$number))


  xlsx.addTitle(sheet2, rowIndex=1, title="Input File Names",
                titleStyle = TITLE_STYLE)

  xlsx::addDataFrame(report.table, sheet2, startRow=3, startColumn=1,
                     colnamesStyle = TABLE_COLNAMES_STYLE,
                     rownamesStyle = TABLE_ROWNAMES_STYLE,col.names = FALSE,row.names = TRUE)



  xlsx::setColumnWidth(sheet2, colIndex=c(1), colWidth=30)
  xlsx::setColumnWidth(sheet2, colIndex=c(2), colWidth=100)

  #========sheet 3 (multi file report)==========
  sheet3 <- xlsx::createSheet(wb, sheetName = "Multi File Report")


  report.table <- t(data.table(
    sapply(study.list, function(x) return(x$MAX_N_TOTAL)),
    sapply(study.list, function(x) return(paste(x$missing.Columns, collapse = ' | '))),
    sapply(study.list, function(x) return(format(x$input.data.rowcount, big.mark="," , scientific = FALSE))),
    sapply(study.list, function(x) return(format(x$rowcount.step1, big.mark="," , scientific = FALSE))),
    sapply(study.list, function(x) return(format(x$rowcount.step2, big.mark="," , scientific = FALSE))),
    sapply(study.list, function(x) return(format(x$rowcount.step3, big.mark="," , scientific = FALSE))),
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
    rep(' ', length(study.list)),
    sapply(study.list, function(x) return(x$tables$variable.summary['Min.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['1st Qu.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['Median', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['Mean', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['3rd Qu.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['Max.', .QC$config$input_parameters$effect_type_string])),
    sapply(study.list, function(x) return(x$tables$variable.summary['Median','STDERR'])),
    sapply(study.list, function(x) return(x$fixed.hwep)),
    sapply(study.list, function(x) return(x$fixed.impq)),
    sapply(study.list, function(x) return(x$fixed.n_total)),
    sapply(study.list, function(x) return(x$fixed.callrate))))

  report.table <- rbind(seq(1:length(study.list)),report.table)

  row.names(report.table)<- c('File number',
                              'Sample Size (Max)',
                              'Missing Columns',
                              'SNPs in input file',
                              'Variant count after step 1 *',
                              'Variant count after step 2 **',
                              'Variant count after step 3 ***',
                              'Monomorphic',
#/	'Duplicates',
                              'Palindromics',
                              'Genotyped variants',
                              'Imputed variants',
                              'Negative-strand SNPs',
                              'Allele Frequency Correlation (Standard Ref)',
                              'Palindromic Allele Frequency correlation (Standard Ref)',
                              'Allele Frequency Correlation (Alternative Ref)',
                              'Palindromic Allele Frequency correlation (Alternative Ref)',
                              'Lambda - Total',
                              'Lambda - Genotyped',
                              'Lambda - Imputed',
                              'P-value Correlation',
                              "Visscher's Statistic (HQ variants)",
                              .QC$config$input_parameters$effect_type_string,
                              "      Min.",
                              "      1st Qu.",
                              "      Median",
                              "      Mean",
                              "      3rd Qu.",
                              "      Max.",
                              "Standard Error (median)",
                              "Fixed HWE P-value",
                              "Fixed Imputation Quality",
                              "Fixed Sample Size",
                              "Fixed Call Rate")

  # colnames(report.table) <- c(seq(1:length(study.list)))
  colnames(report.table) <- sapply(.QC$qc.study.list, function(x) return(x$number))


  xlsx.addTitle(sheet3, rowIndex=1, title="Comparing input files",
                titleStyle = TITLE_STYLE)

  xlsx::addDataFrame(report.table, sheet3, startRow=3, startColumn=1,
                     col.names = FALSE,row.names = TRUE,
                     colnamesStyle = TABLE_COLNAMES_STYLE,
                     rownamesStyle = TABLE_ROWNAMES_STYLE_Left)


  xlsx.addTitle(sheet3, rowIndex= 40, title="step1: removing variants with missing crucial values and duplicated lines.",
                titleStyle = NOTE_TITLE_STYLE2)

  xlsx.addTitle(sheet3, rowIndex= 41, title="step2: removing monomorphic variants and specified chromosomes.",
                titleStyle = NOTE_TITLE_STYLE2)

  xlsx.addTitle(sheet3, rowIndex= 42, title="step3: removing mismatched, ambiguous and multi-allelic variants that could not be verified.",
                titleStyle = NOTE_TITLE_STYLE2)

  xlsx::setColumnWidth(sheet3, colIndex= 1, colWidth = 70)
  xlsx::setColumnWidth(sheet3, colIndex = c(2:(length(study.list) + 1)), colWidth=30)


  #========sheets per input file ==========

  for(i in 1:length(study.list))
  {
    fileSheet <- xlsx::createSheet(wb, sheetName = sprintf('File %s',i))

    study <- study.list[[i]]
    row.index <- 1

    xlsx.addTitle(fileSheet, rowIndex = row.index, title=basename(study$file.path),
                  titleStyle = TITLE_STYLE)


    # introduction
    tbl <- t(data.table(
      format(study$starttime, "%b %d %Y - %X"),
      format(study$endtime, "%b %d %Y - %X"),
      basename(study$file.path),
      study$file.line.count))

    tbl<- cbind(c('Start time',
                  'End time',
                  'Input File Name',
                  'Input File Line Count (including header)'), tbl)


    row.index <- row.index + 2 # 3
    xlsx::addDataFrame(tbl, fileSheet, startRow = row.index , startColumn=1,
                       colnamesStyle = TABLE_COLNAMES_STYLE,row.names = FALSE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE,col.names = FALSE)


    #
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

    filter.table <- t(filter.table)

    colnames(filter.table) <- 'Value'

    row.index <- row.index + 6 #

    xlsx.addTitle(fileSheet, rowIndex=row.index, title="Filter values for selecting High-Quality (HQ) variants",
                  titleStyle = SUB_TITLE_STYLE)

    row.index <- row.index + 1 #
    xlsx::addDataFrame(filter.table, fileSheet, startRow=row.index, startColumn=1,
                       colnamesStyle = TABLE_COLNAMES_STYLE,row.names = TRUE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left,col.names = FALSE)


    #

    row.index <- row.index + 6 #
    xlsx.addTitle(fileSheet, rowIndex= row.index, title='Column names and translation',
                  titleStyle = SUB_TITLE_STYLE)


    column.tbl <- rbind(.QC$thisStudy$original.File.Columns.sorted,
                        .QC$thisStudy$renamed.File.Columns.sorted)
    rownames(column.tbl) <- c('Original', 'Renamed')

    row.index <- row.index + 2 # 3
    xlsx::addDataFrame(t(column.tbl), fileSheet, startRow = row.index , startColumn=1,
                       colnamesStyle = TABLE_COLNAMES_STYLE,row.names = FALSE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE,col.names = TRUE)



   #
    row.index <- row.index + ncol(column.tbl) + 2
    xlsx.addTitle(fileSheet, rowIndex= row.index, title='Column report',
                  titleStyle = SUB_TITLE_STYLE)

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



    colnames(b) <- c('NA values','Invalid values','Uncertain values')

    row.index <- row.index + 1 #
    xlsx::addDataFrame(b , fileSheet, startRow= row.index, startColumn=1,
                       col.names = TRUE ,row.names = TRUE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)

    row.index <- row.index + 19 # 8

    xlsx.addTitle(fileSheet, rowIndex= row.index, title="Variant processing",
                  titleStyle = SUB_TITLE_STYLE)


    row.index <- row.index + 1 # 9
    xlsx.addTitle(fileSheet, rowIndex=row.index, title="step1: removing variants with missing crucial values and duplicated lines.",
                  titleStyle = NOTE_TITLE_STYLE2)

    row.index <- row.index + 1 # 10
    xlsx.addTitle(fileSheet, rowIndex= row.index, title="step2: removing monomorphic variants and specified chromosomes.",
                  titleStyle = NOTE_TITLE_STYLE2)

    row.index <- row.index + 1 # 11
    xlsx.addTitle(fileSheet, rowIndex= row.index, title="step3: removing mismatched, ambiguous and multi-allelic variants that could not be verified.",
                  titleStyle = NOTE_TITLE_STYLE2)


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
    xlsx::addDataFrame(tbl, fileSheet, startRow=row.index, startColumn=1,
                       col.names = FALSE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)




    #



    row.index <- row.index + 10 # 8

    xlsx.addTitle(fileSheet, rowIndex= row.index, title="Description of variants",
                  titleStyle = SUB_TITLE_STYLE)


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
    xlsx::addDataFrame(tbl, fileSheet, startRow= row.index, startColumn=1,
                       col.names = FALSE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)


    row.index <- row.index + 8 # 31
    xlsx.addTitle(fileSheet, rowIndex= row.index, title=sprintf('+ palindromic variants with high allele frequency difference (> %s)',
                                                                config$filters$threshold_diffEAF),
                  titleStyle = NOTE_TITLE_STYLE)

    row.index <- row.index + 1 # 32
    xlsx.addTitle(fileSheet, rowIndex= row.index, title=sprintf('++ Non-palindromic variants with high allele frequency difference (> %s)',
                                                                config$filters$threshold_diffEAF),
                  titleStyle = NOTE_TITLE_STYLE)

	row.index <- row.index + 1 # 33
    xlsx.addTitle(fileSheet, rowIndex= row.index, title='+++ palindromic variants with opposite allele frequency "compared to the reference" (> 0.65 for the input file and < 0.35 for the reference, or vice versa)',
                  titleStyle = NOTE_TITLE_STYLE)



    ### indel-snp type
    if(nrow(study$tables$multi_allele_count_preProcess) > 1)
    {
      row.index <- row.index + 3 # 36
      xlsx.addTitle(fileSheet, rowIndex= row.index, title="Variant types",
                    titleStyle = SUB_TITLE_STYLE)

      # row.index <- row.index + 1 # 37
      # xlsx::addDataFrame(study$tables$VT.tbl, fileSheet, startRow= row.index, startColumn=1,
      #                    col.names = FALSE ,row.names = FALSE,
      #                    colnamesStyle = TABLE_COLNAMES_STYLE,
      #                    rownamesStyle = TABLE_ROWNAMES_STYLE_Left)

      row.index <- row.index + 1 # 37


      xlsx::addDataFrame(study$tables$multi_allele_count_preProcess, fileSheet, startRow= row.index, startColumn=1,
                         col.names = TRUE ,row.names = TRUE,
                         colnamesStyle = TABLE_COLNAMES_STYLE,
                         rownamesStyle = TABLE_ROWNAMES_STYLE_Left)

    }
    #
    row.index <- row.index + 8 # 36
    xlsx.addTitle(fileSheet, rowIndex= row.index, title="Result from matching with standard reference file",
                  titleStyle = SUB_TITLE_STYLE)

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
    xlsx::addDataFrame(tbl, fileSheet, startRow= row.index, startColumn=1,
                       col.names = FALSE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)

    #

    if(!is.na(config$supplementaryFiles$allele_ref_alt))
    {
      row.index <- row.index + 10 # 47
      xlsx.addTitle(fileSheet, rowIndex= row.index, title="Result from matching with alternative reference file",
                    titleStyle = SUB_TITLE_STYLE)

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
      xlsx::addDataFrame(tbl, fileSheet, startRow= row.index, startColumn=1,
                         col.names = FALSE ,row.names = FALSE,
                         colnamesStyle = TABLE_COLNAMES_STYLE,
                         rownamesStyle = TABLE_ROWNAMES_STYLE_Left)
    }
    #

    row.index <- row.index + 11 # 58
    xlsx.addTitle(fileSheet, rowIndex = row.index, title="AF correlation for each chromosome",
                  titleStyle = SUB_TITLE_STYLE)
    row.index <- row.index + 1 # 48



    xlsx::addDataFrame(study$AFcor.std_ref.CHR , fileSheet, startRow= row.index, startColumn=1,
                       col.names = TRUE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)


    #
    row.index <- row.index + nrow(study$AFcor.std_ref.CHR) + 2
    xlsx.addTitle(fileSheet, rowIndex = row.index, title="QC summary statistics",
                  titleStyle = SUB_TITLE_STYLE)

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
    xlsx::addDataFrame(tbl, fileSheet, startRow= row.index, startColumn=1,
                       col.names = FALSE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)



    # variable summary statistics

    row.index <- row.index + 19 # 78
    xlsx.addTitle(fileSheet, rowIndex= row.index, title="Distribution statistics",
                  titleStyle = SUB_TITLE_STYLE)

    row.index <- row.index + 2 # 80

    row.names(study$tables$variable.summary) <- c("min.","first_quartile","median","mean","third_quartile","max." )
    xlsx::addDataFrame(t(study$tables$variable.summary), fileSheet, startRow= row.index, startColumn=1,
                       col.names = TRUE ,row.names = TRUE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)




    # column summary statistics

    row.index <- row.index + 10
    xlsx.addTitle(fileSheet, rowIndex= row.index, title="Variant count for each chromosome",
                  titleStyle = SUB_TITLE_STYLE)


    ####
    row.index <- row.index + 1 # 89
    chr.tbl.length <- 0
    if(!is.na(study$tables$CHR.tbl))
    {
      chr.tbl <-  study$tables$CHR.tbl
      colnames(chr.tbl) <- c('Chromosome Number','Variant Count')
      chr.tbl.length <- nrow(chr.tbl)

      xlsx::addDataFrame(chr.tbl, fileSheet, startRow= row.index, startColumn=1,
                         col.names = TRUE ,row.names = FALSE,
                         colnamesStyle = TABLE_COLNAMES_STYLE,
                         rownamesStyle = TABLE_ROWNAMES_STYLE_Left)

    }

    if(length(.QC$thisStudy$missing_chromosomes) >0 )
    {

      row.index <- row.index + nrow(chr.tbl) + 2
      xlsx.addTitle(fileSheet, rowIndex= row.index,
                    title=sprintf("%s %s","Missing chromosome(s) number",paste(.QC$thisStudy$missing_chromosomes,collapse = ", ")),
                    titleStyle = SUB_TITLE_STYLE)

      row.index <- row.index + 2
    }
    else
    {
      row.index <- row.index + nrow(chr.tbl) + 2 #
    }

    ###############
    ###
    xlsx.addTitle(fileSheet, rowIndex= row.index, title="Effect allele distribution in SNP variants",
                  titleStyle = SUB_TITLE_STYLE)



    row.index <- row.index + 1 #

    tbl = merge(study$tables$EFFECT_ALL.tbl,
                study$tables$EFFECT_ALL.post.matching.tbl,
                by="EFFECT_ALL",
                all = TRUE)

    colnames(tbl) <- c('Allele','Count (input file)','Count (post-matching)')

    xlsx::addDataFrame(tbl , fileSheet, startRow= row.index, startColumn=1,
                       col.names = TRUE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)


    ###
    row.index <- row.index + 6 #
    xlsx.addTitle(fileSheet, rowIndex= row.index, title="Other allele distribution in SNP variants",
                  titleStyle = SUB_TITLE_STYLE)



    row.index <- row.index + 1 #

    tbl = merge(study$tables$OTHER_ALL.tbl,
                study$tables$OTHER_ALL.post.matching.tbl,
                by="OTHER_ALL",
                all = TRUE)

    colnames(tbl) <- c('Allele','Count (input file)','Count (post-matching)')

    xlsx::addDataFrame(tbl , fileSheet, startRow= row.index, startColumn=1,
                       col.names = TRUE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)



    ##


    row.index <- row.index + 6 #
    xlsx.addTitle(fileSheet, rowIndex= row.index, title='Imputation status',
                  titleStyle = SUB_TITLE_STYLE)



    row.index <- row.index + 1 #
    tbl = study$tables$imputed.tbl
    # tbl$IMPUTED <- c('Genotyped','Imputed')
    colnames(tbl) <- c('Status','Count')

    xlsx::addDataFrame(tbl , fileSheet, startRow= row.index, startColumn=1,
                       col.names = TRUE ,row.names = FALSE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)


    ##





    #
    b <- t(data.frame('Negative strand variants' = study$neg.strand.count))
    colnames(b) <- c('Count')

    row.index <- row.index + 5
    xlsx::addDataFrame(b , fileSheet, startRow= row.index, startColumn=1,
                       col.names = TRUE ,row.names = TRUE,
                       colnamesStyle = TABLE_COLNAMES_STYLE,
                       rownamesStyle = TABLE_ROWNAMES_STYLE_Left)


    # effect size comparison

    if(!is.na(config$supplementaryFiles$beta_ref_std)){

      row.index <- row.index + 3 # 88
      xlsx.addTitle(fileSheet, rowIndex= row.index, title="Effect-size comparison",
                    titleStyle = SUB_TITLE_STYLE)


      b <- t(data.frame('r' = study$effect.rho_4))
      colnames(b) <- c('Value')

      row.index <- row.index + 1
      xlsx::addDataFrame(b , fileSheet, startRow= row.index, startColumn=1,
                         col.names = TRUE ,row.names = TRUE,
                         colnamesStyle = TABLE_COLNAMES_STYLE,
                         rownamesStyle = TABLE_ROWNAMES_STYLE_Left)
    }

    ## END OF EXCEL REPORT
    xlsx::setColumnWidth(fileSheet, colIndex= 1, colWidth = 30)
    xlsx::setColumnWidth(fileSheet, colIndex = c(2:10), colWidth=20)

  }



  # ========= write to file ===========

  xlsx::saveWorkbook(wb, config$paths$xlsx.report)

  return(TRUE)

}
