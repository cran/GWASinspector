#### added columns to dataset
# found	, whether found in any reference files
# palindromic
# match	, no nedd to flip or switch
# flip
# switch
# wrong	, wrong variant based on alleles and matching with references
# highDiffEAF
# HQ
# PVALUE.calculated

process.each.file <- function(study){


  print.and.log(sprintf('============== [ File %s from %s ] ==============',
                        .QC$file.counter,
                        length(.QC$qc.study.list)),
                'info',
                cat = TRUE)





  ## processing crucial columns and creating report variables
  config <- .QC$config

  # study variable contains all variables and paths
  .QC$thisStudy <- study

  .QC$thisStudy$starttime <- Sys.time()

  .QC$thisStudy$number <- .QC$file.counter
  .QC$thisStudy$effect_type_string <- .QC$config$input_parameters$effect_type_string

  .QC$file.counter <- .QC$file.counter + 1 ## counter that displays how many files are processed

  #print.and.log(mem_used(),'info',cat= FALSE)

  ### uplolading and processing file
  ## ==============================================
  input.data <- uploadInputFile()

  #print.and.log(mem_used(),'info',cat= FALSE)


  # return null if file could not be read
  #====================================
  if(is.null(input.data)){
    print.and.log('File removed from QC analysis due to error in loading!','warning')
    addEmptyStudy_studyInstance(.QC$thisStudy)
    return(NULL)
  }

  invisible(gc())


  # remove duplicate lines of data
  # study$dup_lines_count is set
  #TODO ignored because it is very time consuming
  #input.data <- removeDuplicatedLines(input.data)



  # IMPORTANT: variants may be removed in 2 steps
  # 1- missing crucial variables
  # 2- monomorphic, duplicate, chromosomes
  ## ==============================================

  input.data <- processInputFile(input.data)
  #print.and.log(mem_used(),'info',cat= FALSE)

  # return null if error encountered
  #====================================
  if(is.null(input.data)){
    print.and.log('File removed from QC analysis due to error in processing!','warning')
    addEmptyStudy_studyInstance(.QC$thisStudy)
    return(NULL)
  }

  # variable.statistics.pre.matching(input.data) moved to process.matched.data()

  invisible(gc())



  ### matching with reference dataset
  ## ==============================================
  print.and.log('Comparing input file with reference file (this might take long) ...','info')
  input.data<-tryCatch(compareInputfileWithReferenceData(input.data),
                       error = function(err) {
                         print.and.log(paste('Error in comparing process:' , err$message), 'warning')
                         return(NULL)
                       }
  )

  #====================================
  if(is.null(input.data)){
    print.and.log('File removed from QC analysis due to error in comparing process!','warning')
    addEmptyStudy_studyInstance(.QC$thisStudy)
    return(NULL)
  }

  invisible(gc())
  #print.and.log(mem_used(),'info',cat= FALSE)


  ## allele matching
  ## ==============================================
  print.and.log('Allele matching vs Allele Freq reference dataset ...','info')

  input.data <- tryCatch(input.data[,c("match_result","palindromic") :=
                                      variant.match(EFFECT_ALL,OTHER_ALL,ALT,REF,VT),
                                    by = list(EFFECT_ALL,OTHER_ALL,ALT,REF,VT)],
                         error = function(err) {
                           print.and.log(paste('Error in allele matching:' , err$message), 'warning')
                           return(NULL)
                         }
  )

  #====================================
  if(is.null(input.data)){
    print.and.log('File removed from QC analysis due to error in allele matching!','warning')
    addEmptyStudy_studyInstance(.QC$thisStudy)
    return(NULL)
  }

  invisible(gc())
  #print.and.log(mem_used(),'info',cat= FALSE)



  ## flip and switch and analyzing
  # IMPORTANT: variants may be removed in 1 step
  # 3- mismatched variants
  ## ==============================================
  print.and.log('Data processing ...','info')
  input.data <- tryCatch(process.matched.data(input.data),
                         error = function(err) {
                           print.and.log(paste('Error in variant processing:' , err$message), 'warning')
                           return(NULL)
                         }
  )

  #====================================
  if(is.null(input.data)){
    print.and.log('File removed from QC analysis due to error in variant processing!','warning')
    addEmptyStudy_studyInstance(.QC$thisStudy)
    return(NULL)
  }


  invisible(gc())



  ## Calculations
  ## ==============================================
  ## calculate Pvalue from stderr and effect
  input.data <- calculate.PVALUE(input.data) #pValueFunctions.R

  ## calculate correlation between the 2 Pvalues
  calculate.pvalue.correlation(input.data)

  # calcualte AF correlation between study and std ref
  calculate.af.correlation.std_ref(input.data) #calculationFunctions.R

  # calcualte AF correlation between study and alt ref
  if(!is.na(config$supplementaryFiles$allele_ref_alt))
    calculate.af.correlation.alt_ref(input.data) #calculationFunctions.R

  # selcet high quality variants based on 4 filters
  input.data <- applyHQfilter(input.data) #calculationFunctions.R


  # generate distribution statistcis for HQ variants
  # summary of important columns
  if(.QC$thisStudy$HQ.count == .QC$thisStudy$rowcount.step3)
  {
    .QC$thisStudy$tables$variable.summary.HQ <- .QC$thisStudy$tables$variable.summary
  }
  else if(.QC$thisStudy$HQ.count > 0)
  {
    .QC$thisStudy$tables$variable.summary.HQ <- variable.statistics.post.matching_HQ(input.data)
  }

  ##TODO  outliers are defined as skewness > 0.1 or < -0.1, or kurtosis > 10
  calculateSkewness(input.data)
  calculateSkewness.HQ(input.data)
  calculateKurtosis(input.data)
  calculateKurtosis.HQ(input.data)
  calculateVischerStats(input.data)
  calculateVischerStats.HQ(input.data)
  calculateLambda(input.data)




  ### ===================== variables for multifile comparison ==============
  ## --------- effect plot for box plot -------------
  eff.col <- as.data.table(input.data[HQ==TRUE]$EFFECT)
  ##------ variables for precision plot ---------
  ## ============================================
  .QC$thisStudy$STDERR.mean.HQ <- mean(input.data[HQ == TRUE]$STDERR,na.rm = TRUE)


  df <- data.frame(
    x = 1,
    y0 = min(eff.col$V1),
    y25 = stats::quantile(eff.col$V1, 0.25),
    y50 = stats::median(eff.col$V1),
    y75 = stats::quantile(eff.col$V1, 0.75),
    y100 = max(eff.col$V1)
  )

  # upper whisker = largest observation less than or equal to upper hinge + 1.5 * IQR
  y_upper = df$y75 + 1.5 * stats::IQR(eff.col$V1)

  #lower whisker = smallest observation greater than or equal to lower hinge - 1.5 * IQR
  y_lower = df$y25 - 1.5 * stats::IQR(eff.col$V1)

  #TODO maybe show values inside plot
  .QC$thisStudy$effect.plot.df <- df
  .QC$thisStudy$effect.plot.df_y_upper <- y_upper
  .QC$thisStudy$effect.plot.df_y_lower <- y_lower

  # generate the plot usin generateEffectSizePlot() function instead of saving the plot
  # file.number = .QC$thisStudy$number
  #
  # if("N_CASES" %in% .QC$thisStudy$renamed.File.Columns)
  # {
  #   file.N.max = .QC$thisStudy$MAX_N_CASES
  #   print.and.log("N_CASES will be used for MAX_N value.")
  # }
  # else
  #   file.N.max = .QC$thisStudy$MAX_N_TOTAL
  #
  #
  #
  # .QC$thisStudy$effect.plot = ggplot(df, aes(x)) +
  #   geom_boxplot(
  #     aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
  #     stat = "identity") +
  #   geom_errorbar(aes(ymin=y_lower, ymax=y_upper), width=0.6,
  #                 position=position_dodge(.9)) +
  #   labs(x= file.number ,y="effect size",subtitle = sprintf('N = %s', file.N.max)) +
  #   theme_classic(base_size = 8)+
  #   coord_cartesian(ylim = c(-0.6,0.6)) +
  #   geom_hline(yintercept = 0.1,linetype = 2,color='red') +
  #   geom_hline(yintercept = -0.1,linetype = 2,color='red') +
  #   theme(axis.text.x=element_blank())


  ### ===================== END of multi study variables =====================


  ## updating alternate refrence file
  ## ==============================================
  #cat('\n--- [evaluating alternate reference ...] ---',fill=TRUE)
  if(!is.na(config$supplementaryFiles$allele_ref_alt))
    update.alternate.reference(input.data)



  ## data cleaning 1
  ## ==============
  rm(eff.col) # was used for eff.plot

  # added while checking with reference sets
  input.data[, CHR.y := NULL]
  input.data[, POS := NULL]
  input.data[, ALT := NULL]
  input.data[, REF := NULL]
  # input.data[, AF := NULL] required for alele freq plot and correlation calculation
  input.data[, DATE_ADDED := NULL]
  # input.data[, SOURCE := NULL] required for alele freq plot and correlation calculation

  # added while matching with refrence sets
  input.data[, match := NULL]
  input.data[, flip := NULL]
  input.data[, switch := NULL]
  input.data[, wrong := NULL]

  invisible(gc())
  #print.and.log(mem_used(),'info',cat= FALSE)





  ## 6- plots
  ## ==============================================
  print.and.log('Plotting ...','info')
  drawPlots(input.data)

  invisible(gc())


  ## ==============================================
  print.and.log('Saving data set ...','info')
  ## missing pvalues set from calculated pvalues
  input.data <-fill.missing.pvalues.from.calculated.pvalues(input.data)


  saveDataSet.final(input.data) #'saveFilesFunctions.R'



  invisible(gc())
  #print.and.log(mem_used(),'info',cat= FALSE)


  ## data cleaning 2
  ## ==============
  # these items contain long list of missing and invalid variant indexes
  # keeps huge memory
  # not required after creating the report
  .QC$thisStudy$column.NA.list <- lapply(.QC$thisStudy$column.NA.list, length)
  .QC$thisStudy$column.INVALID.list <- lapply(.QC$thisStudy$column.INVALID.list, length)


  .QC$thisStudy$endtime <- Sys.time()

  # save txt report file
  report.to.txt.file(.QC$thisStudy)




  ## Compare Beta(Effect) values with refrence set and draw plots
  # =============================================================


  if(nrow(.QC$reference.data.effect) > 0 ){

    print.and.log('Comparing input file with Effect-Size reference file ...','info')
    input.data<-tryCatch(compareInputfileWithBetaReferenceFile(input.data),
                         error = function(err) {
                           print.and.log(paste('Error in comparing with Effect-Size reference:' , err$message), 'warning')
                           return(data.table())
                         }
    )

    # only variants that could be matched are returned.

    # gc not done in previous function
    invisible(gc())

    ## allele matching for effect size plot
    ## ==============================================
    if(nrow(input.data) > 0){



      # THIS STEP IS NOT REQUIRED, BECAUSE DATA IS ALREADY MATCHED WITH REFERENCE DATASET
      # print.and.log('Allele matching vs Effect-Size reference dataset ...','info')
      # input.data <- tryCatch(allele.match.effectPlot(input.data),
      #                        error = function(err) {
      #                          print.and.log(paste('Error in allele matching:' , err$message), 'warning')
      #                          return(NULL)
      #                        }
      # )




      if(config$plot_specs$make_plots & !is.null(input.data))
        tryCatch(plot.DataEFFECT.vs.RefEFFECT(input.data,
                                              .QC$thisStudy$effPlotPath,
                                              .QC$thisStudy$plot.title),
                 error = function(err) {
                   print.and.log(paste('Error plotting effect-size correlation plot:' , err$message), 'warning')
                   return(NULL)
                 }
        )


    }else
    {
      print.and.log('No variants were found in Effect-size reference dataset!','warning')
      print.and.log('Effect-size comparison plot is skipped!','warning')
      .QC$thisStudy$effect.rho_4 <- 'NA (no variants were found)'
    }


    # write(x = '\n\n[Comparing Input File with Effect-size reference dataset (HQ variants)]',
    #     file = .QC$thisStudy$txt.report.path ,
    #     sep =  '\n',
    #     append = TRUE)

    writeTXTreport('\n\nComparing Input File with Effect-size reference dataset (HQ variants)')
    writeTXTreport(kable(.QC$thisStudy$tables$betaCor.tbl, align = "l",format = "rst"))

    writeTXTreport('\n* Data is presented as r(N). Variants were filtered on reference data P-values. ')
    writeTXTreport('** Data is presented as r(N). Variants were filtered on input result file P-values. ')

#     write(x = paste('r (P-value < 0.001) =' ,  .QC$thisStudy$effect.rho_3),
#         file = .QC$thisStudy$txt.report.path ,
#         sep =  '\n',
#         append = TRUE)
#
#     write(x = paste('r (P-value < 0.0001) =' ,  .QC$thisStudy$effect.rho_4),
#         file = .QC$thisStudy$txt.report.path ,
#         sep =  '\n',
#         append = TRUE)
#
#     write(x = paste('r (P-value < 0.00001) =' ,  .QC$thisStudy$effect.rho_5),
# 	    file = .QC$thisStudy$txt.report.path ,
# 	    sep =  '\n',
# 	    append = TRUE)
#
#     write(x = paste('r (P-value < 0.000001) =' ,  .QC$thisStudy$effect.rho_6),
# 	    file = .QC$thisStudy$txt.report.path ,
# 	    sep =  '\n',
# 	    append = TRUE)

  }

  ## save significant variants ( p-value  < 1e-8  & HQ & rsID) ==> CHR-POS-MARKER-PVALUE
  #save.significant.variants(input.data)


  # save rds file
  studyClass <- create.Study(.QC$thisStudy)
  # .QC$StudyList <- append( .QC$StudyList , studyClass)
  .QC$StudyList@studyList <- append(.QC$StudyList@studyList , studyClass)
  .QC$StudyList@studyCount <- length(.QC$StudyList@studyList)

  save.rds.file(.QC$thisStudy)

  rm(input.data)
  invisible(gc())

  ## ended
  ##==============
  print.and.log('\n','info')

  ## add info to text report
  writeTXTreport("\n==============================================")
  writeTXTreport(sprintf('Generated by %s package - v.%s',.QC$package.name, .QC$script.version))
  writeTXTreport(as.character(Sys.time()))

  return(.QC$thisStudy)
}






create.file.specific.config <- function(file.name){

  ### these varaibles are used for each file that is QC'ed ###

  # study$starttime
  # study$endtime
  #
  #
  # study$AFcor.std_ref
  # study$AFcor.non.palindromic.std_ref
  # study$AFcor.palindromic.std_ref

  # study$AFcor.alt_ref
  # study$AFcor.non.palindromic.alt_ref
  # study$AFcor.palindromic.alt_ref


  # study$duplicate.count
  # study$HQ.count
  # study$input.data.rowcount
  # study$kurtosis
  # study$kurtosis.HQ
  # study$lambda
  # study$lambda.gen
  # study$lambda.imp
  # study$LQ.count
  # study$missing.crucial.rowcount
  # study$multiAlleleVariants.rowcount
  # study$missing.alleles.rowcount
  # study$monomorphic.count
  # study$neg.strand.count
  # study$PVcor
  # study$PVcor.palindromic
  # study$skewness
  # study$skewness.HQ
  # study$Visschers.stat
  # study$Visschers.stat.HQ
  # Study$MAX_N_TOTAL
  # study$hasINDEL
  # study$hID.added
  # check if there are characters in chromosome column.
  # they are converted to number and should be deconverted before saving final dataset
  # study$character.chromosome

  # if study has none base character as alleles => so they sholod not be regarded as mismatches
  # study$hanNoneBaseAlleles

  # study$column.INVALID.list$CALLRATE
  # study$column.INVALID.list$CHR
  # study$column.INVALID.list$EFF_ALL_FREQ
  # study$column.INVALID.list$EFFECT
  # study$column.INVALID.list$EFFECT_ALL
  # study$column.INVALID.list$HWE_PVAL
  # study$column.INVALID.list$IMP_QUALITY
  # study$column.INVALID.list$IMPUTED
  # study$column.INVALID.list$minusone.CALLRATE
  # study$column.INVALID.list$minusone.EFF_ALL_FREQ
  # study$column.INVALID.list$minusone.HWE_PVAL
  # study$column.INVALID.list$minusone.PVALUE
  # study$column.INVALID.list$N_TOTAL
  # study$column.INVALID.list$one.EFF_ALL_FREQ
  # study$column.INVALID.list$OTHER_ALL
  # study$column.INVALID.list$POSITION
  # study$column.INVALID.list$PVALUE
  # study$column.INVALID.list$STDERR
  # study$column.INVALID.list$zero.EFF_ALL_FREQ
  # study$column.INVALID.list$zero.STDERR
  #
  #
  #
  #
  # study$column.NA.list$CALLRATE
  # study$column.NA.list$CHR
  # study$column.NA.list$EFF_ALL_FREQ
  # study$column.NA.list$EFFECT
  # study$column.NA.list$EFFECT_ALL
  # study$column.NA.list$HWE_PVAL
  # study$column.NA.list$IMP_QUALITY
  # study$column.NA.list$IMPUTED
  # study$column.NA.list$MARKER
  # study$column.NA.list$N_TOTAL
  # study$column.NA.list$OTHER_ALL
  # study$column.NA.list$POSITION
  # study$column.NA.list$PVALUE
  # study$column.NA.list$STDERR
  # study$column.NA.list$STRAND
  #
  #
  #
  # study$tables$CHR.tbl
  # study$tables$VT.tbl
  # study$tables$EFFECT_ALL.tbl
  # study$tables$imputed.tbl
  # study$tables$OTHER_ALL.tbl
  # study$tables$match.ref.table  whcih ref is used for matching
  # study$tables$VT.ref.table     how many variants of each type and multi allelic or not


  ## multi file comparison
  # study$STDERR.mean.HQ
  # study$N.max
  # study$effect.plot

  # title if the plots . it is either the name of input file or set by user
  # it is automatrically set as file name if there are multiple files
  # study$plot.title


  study <- c()
  config <- .QC$config

  file.name <- as.character(file.name[[1]])

  study$file.path <- file.name

  study$file.extension <- file_ext(file.name)

  study$zipped.File <- FALSE

  if(study$file.extension %in% c('gz','zip','bz2'))
    study$zipped.File <- TRUE



  ##-- getting filename without extensions
  study$file.name <- tools::file_path_sans_ext(basename(file.name))

  # remove the extension before gz
  #if(study$file.extension == 'gz' && any(endsWith(x =  study$file.name,suffix = c('.txt','.csv')))) ## replaced to be compatible with older R
  if(study$file.extension %in% c('gz','bz2') && any( grepl("(txt|csv)$", x = study$file.name)))
  {
    study$file.name <- sub(x= study$file.name , pattern = '\\.(\\w{3})$', replacement = '')
  }


  ###------
  print.and.log(sprintf("Checking Study file : '%s'",file.name),
                'info')

  # get file line count if file is not zipped
  # if(config$test.run) # do not count the lines if it is a test run
  #   study$file.line.count <- 'NA (test run)'
  # else if(study$zipped.File)  # do not count the lines if it is zipped
  #   study$file.line.count <- 'NA (zipped file)'
  # else
  #   study$file.line.count <- get.file.line.count(file.name)

  # get file line count if file is not zipped
  if(study$file.extension != 'zip') # do not count the lines if it is a test run
  {
    fileInspection <- get.file.line.count_RUtils(file.name)
    study$file.line.count <- fileInspection[1]
    study$file.endsWithNewLine <- fileInspection[2]
  }


  ###### ==== file names ==== ##########
  # this prefix is put infront of all files. it does not have any extension
  # /home/user/  +  QC  +  _  +  studyfile
  files.prefix <- sprintf('%s/%s_%s',
                          config$paths$dir_output,
                          config$paths$filename_output_tag,
                          study$file.name)



  ## add plot paths
  study$manPlotPath<-paste0(files.prefix, '_graph_M', .QC$img.extension)


  study$histPlotPath<-paste0(files.prefix, '_graph_histogram' , .QC$img.extension)


  study$stdMafPlotPath<-paste0(files.prefix, '_graph_EAF_SR' , .QC$img.extension)

  study$stdMafSmPlotPath<-paste0(files.prefix, '_graph_EAF_SR_sm' , .QC$img.extension)

  study$altMafPlotPath<-paste0(files.prefix, '_graph_EAF_AR' , .QC$img.extension)


  study$pvalCorPlotPath<-paste0(files.prefix, '_graph_p_correlation' , .QC$img.extension)

  study$pvalCorSmPlotPath<-paste0(files.prefix, '_graph_p_correlation_sm' , .QC$img.extension)


  study$QQPlotPath<-paste0(files.prefix, '_graph_QQ' , .QC$img.extension)

  study$effPlotPath<-paste0(files.prefix, '_beta' , .QC$img.extension)

  study$plot.title <- ifelse(length(config$paths$input_files) > 1 | config$plot_specs$plot_title == 'none',
                             study$file.name,
                             config$plot_specs$plot_title)

  ## add dataset paths

  #variant with missing crucial values
  study$SNPs_invalid.path<-paste(files.prefix, 'vars_invalid_allele.txt' , sep='_')


  #variant with invalid allele values
  study$SNPs_removed.path<-paste(files.prefix, 'vars_removed.txt' , sep='_')

  #variant with invalid both allele values
  study$SNPs_invalid_both.path<-paste(files.prefix, 'vars_invalid_alleles.txt' , sep='_')

  #variant with missing non-crucial values
  study$SNPs_improbable_values.path<-paste(files.prefix, 'vars_improbable_values.txt' , sep='_')

  # duplicate variants
  study$SNPs_duplicates.path<-paste(files.prefix, 'vars_duplicates.txt' , sep='_')
  study$SNPs_duplicates_postMatch.path<-paste(files.prefix, 'vars_duplicates_post_match.txt' , sep='_')

  # mismatched bi-allelic variants
  study$SNPs_mismatches.path<-paste(files.prefix, 'vars_mismatches_BA.txt' , sep='_')

  # monomorphic variants
  study$SNPs_monomorphic.path<-paste(files.prefix, 'vars_monomorphic.txt' , sep='_')

  # multi_allelic variants
  # vairants that are found in reference but cannot be matched because of many variants on the same position
  # mismatched multi-allelic variants
  study$SNPs_multi_allelic.path<-paste(files.prefix, 'vars_mismatches_MA.txt' , sep='_')

  # duplicated match variants
  study$SNPs_ambiguous.path<-paste(files.prefix, 'vars_ambiguous.txt' , sep='_')

  # significant variants
  study$SNPs_significant.path<-paste(files.prefix, 'vars_significant.txt' , sep='_')


  # cleaned output file
  study$output.path<-paste(files.prefix, '.txt' , sep='')


  # effect-size reference dataset path
  study$effect_size_ref.output.path<-paste(files.prefix, 'effect-size_ref.rds' , sep='_')


  # html report file
  study$html.report.path<-paste(files.prefix, 'report.html' , sep='_')

  # txt report file
  study$txt.report.path<-paste(files.prefix, 'report.txt' , sep='_')

  # rdata object file
  study$rds.study.rds.path<-paste(files.prefix, 'object.rds' , sep='_')



  ## ====== column variables ======== ##

  study$original.File.Columns <- character(0)
  study$original.File.Columns.sorted <- character(0)
  study$renamed.File.Columns.sorted  <- character(0)
  study$wanted.columns.index <- character(0)
  study$renamed.File.Columns.classes <- character(0)
  study$missing.Columns <- 'none'

  study$missing.PVALUE.column <- FALSE



  study$hasINDEL <- FALSE
  study$hID.added <- FALSE
  study$hanNoneBaseAlleles <- FALSE
  study$HQ.count  <- 0
  study$LQ.count  <- 0

  study$input.data.rowcount <- 0 # originla file row number - before any processing
  study$rowcount.step1 <- 0      # after processing crucial columns
  study$rowcount.step2 <- 0      # aftere removing duplicates, monomorphics and chromosal variants
  study$rowcount.step3 <- 0      # aftere

  study$missing.crucial.rowcount<- 0
  study$missing.alleles.rowcount <- 0
  study$neg.strand.count<- 0

  study$monomorphic.count <- 0
  study$duplicate.count<- 0


  study$found.rows <- 0
  study$mismatched.rows<- 0
  study$ambiguos.rows<- 0
  study$not.found.rows <- 0
  study$switched.rows <- 0
  study$flipped.rows <- 0
  study$palindromic.rows <- 0
  study$non.palindromic.rows<- 0

  study$skewness <- 'NA'
  study$skewness.HQ <- 'NA'
  study$Visschers.stat <- 'NA'
  study$Visschers.stat.HQ <- 'NA'
  study$MAX_N_TOTAL <- 'NA'


  study$tables$imputed.tbl <-  as.data.table(matrix(NA,2,2))
  study$tables$multi_allele_count_preProcess <- data.table()
  study$tables$multi_allele_count_postProcess <- data.table()
  study$tables$variable.summary.HQ <- data.table()

  colnames(study$tables$imputed.tbl) <- c('IMPUTED','N')
  study$tables$imputed.tbl$IMPUTED[1] <- 'Imputed'
  study$tables$imputed.tbl$IMPUTED[2] <- 'Genotyped'

  study$lambda <- 0
  study$lambda.gen <- 'not available'
  study$lambda.imp <- 'not available'

  study$PVcor <- 'not available (missing P-value column)'
  study$PVcor.palindromic <- 'not available (missing P-value column)'
  study$rownum.PVcor <- 0


  study$fixed.hwep <- 'NA'
  study$fixed.callrate <- 'NA'
  study$fixed.n_total <- 'NA'
  study$fixed.impq <- 'NA'

  study$AFcor.std_ref <- 'NA'
  study$AFcor.non.palindromic.std_ref <- 'NA'
  study$AFcor.palindromic.std_ref <- 'NA'
  study$AFcor.std_ref.indel <- 'NA'
  study$AFcor.std_ref.CHR <- 'NA'


  study$AFcor.alt_ref <- 'NA'
  study$AFcor.non.palindromic.alt_ref <- 'NA'
  study$AFcor.palindromic.alt_ref <- 'NA'
  study$AFcor.alt_ref.indel <- 'NA'

  study$character.chromosome <- FALSE

  ## multi file comparison
  study$STDERR.mean.HQ <- 0
  study$effect.rho_4 <- 'NA (not calculated)'
  # study$column.INVALID.list$CALLRATE <- numeric(length = 0L)
  # study$column.INVALID.list$CHR <- numeric(length = 0L)
  # study$column.INVALID.list$EFF_ALL_FREQ <- numeric(length = 0L)
  # study$column.INVALID.list$EFFECT <- numeric(length = 0L)
  # study$column.INVALID.list$EFFECT_ALL <- numeric(length = 0L)
  # study$column.INVALID.list$HWE_PVAL <- numeric(length = 0L)
  # study$column.INVALID.list$IMP_QUALITY <- numeric(length = 0L)
  # study$column.INVALID.list$IMPUTED <- numeric(length = 0L)
  # study$column.INVALID.list$N_TOTAL <- numeric(length = 0L)
  # study$column.INVALID.list$OTHER_ALL <- numeric(length = 0L)
  # study$column.INVALID.list$POSITION <- numeric(length = 0L)
  # study$column.INVALID.list$PVALUE <- numeric(length = 0L)
  # study$column.INVALID.list$STDERR <- numeric(length = 0L)

  study$column.INVALID.list$minusone.CALLRATE <- numeric(length = 0L)
  study$column.INVALID.list$zero.EFF_ALL_FREQ <- numeric(length = 0L)
  study$column.INVALID.list$one.EFF_ALL_FREQ <- numeric(length = 0L)
  study$column.INVALID.list$minusone.EFF_ALL_FREQ <- numeric(length = 0L)
  study$column.INVALID.list$minusone.HWE_PVAL <- numeric(length = 0L)
  study$column.INVALID.list$minusone.PVALUE <- numeric(length = 0L)
  study$column.INVALID.list$zero.STDERR <- numeric(length = 0L)


  study$column.INVALID.list$CALLRATE = numeric(length = 0L)
  study$column.INVALID.list$CHR = numeric(length = 0L)
  study$column.INVALID.list$EFF_ALL_FREQ = numeric(length = 0L)
  study$column.INVALID.list$EFFECT = numeric(length = 0L)
  study$column.INVALID.list$EFFECT_ALL = numeric(length = 0L)
  study$column.INVALID.list$HWE_PVAL = numeric(length = 0L)
  study$column.INVALID.list$IMP_QUALITY = numeric(length = 0L)
  study$column.INVALID.list$IMPUTED = numeric(length = 0L)
  study$column.INVALID.list$MARKER = numeric(length = 0L)
  study$column.INVALID.list$N_TOTAL = numeric(length = 0L)
  study$column.INVALID.list$OTHER_ALL = numeric(length = 0L)
  study$column.INVALID.list$BOTH_ALL = numeric(length = 0L)
  study$column.INVALID.list$POSITION = numeric(length = 0L)
  study$column.INVALID.list$PVALUE = numeric(length = 0L)
  study$column.INVALID.list$STDERR = numeric(length = 0L)
  study$column.INVALID.list$STRAND = numeric(length = 0L)
  study$column.NA.list$CALLRATE = numeric(length = 0L)
  study$column.NA.list$CHR = numeric(length = 0L)
  study$column.NA.list$EFF_ALL_FREQ = numeric(length = 0L)
  study$column.NA.list$EFFECT = numeric(length = 0L)
  study$column.NA.list$EFFECT_ALL = numeric(length = 0L)
  study$column.NA.list$HWE_PVAL = numeric(length = 0L)
  study$column.NA.list$IMP_QUALITY = numeric(length = 0L)
  study$column.NA.list$IMPUTED = numeric(length = 0L)
  study$column.NA.list$MARKER = numeric(length = 0L)
  study$column.NA.list$N_TOTAL = numeric(length = 0L)
  study$column.NA.list$OTHER_ALL = numeric(length = 0L)
  study$column.NA.list$POSITION = numeric(length = 0L)
  study$column.NA.list$PVALUE = numeric(length = 0L)
  study$column.NA.list$STDERR = numeric(length = 0L)
  study$column.NA.list$STRAND = numeric(length = 0L)

  # study$column.NA.list$CALLRATE <- numeric(length = 0L)
  # study$column.NA.list$CHR <- numeric(length = 0L)
  # study$column.NA.list$EFF_ALL_FREQ <- numeric(length = 0L)
  # study$column.NA.list$EFFECT <- numeric(length = 0L)
  # study$column.NA.list$EFFECT_ALL <- numeric(length = 0L)
  # study$column.NA.list$HWE_PVAL <- numeric(length = 0L)
  # study$column.NA.list$IMP_QUALITY <- numeric(length = 0L)
  # study$column.NA.list$IMPUTED <- numeric(length = 0L)
  # study$column.NA.list$MARKER <- numeric(length = 0L)
  # study$column.NA.list$N_TOTAL <- numeric(length = 0L)
  # study$column.NA.list$OTHER_ALL <- numeric(length = 0L)
  # study$column.NA.list$POSITION <- numeric(length = 0L)
  # study$column.NA.list$PVALUE <- numeric(length = 0L)
  # study$column.NA.list$STDERR <- numeric(length = 0L)
  # study$column.NA.list$STRAND <- numeric(length = 0L)

  # specific chromosome count in file
  study$x.chr.count.removed <- 0
  study$y.chr.count.removed <- 0
  study$xy.chr.count.removed <- 0
  study$m.chr.count.removed <- 0


  study$dup_lines_count <- 0

  # get variable summary statistics from table
  # ===================
  # study$tables$variable.summary['Min.','EFFECT']
  # study$tables$variable.summary['1st Qu.','EFFECT']
  # study$tables$variable.summary['Median','EFFECT']
  # study$tables$variable.summary['Mean','EFFECT']
  # study$tables$variable.summary['3rd Qu.','EFFECT']
  # study$tables$variable.summary['Max.','EFFECT']


  study$effect.plot <- textGrob("Empty Effect Plot", gp = gpar(fontsize=12,col='red', fontface='bold'))

  # study$effect.plot <- ggdraw() +
  #   draw_label('Empty Effect Plot', x = .5, y = .5,
  #              vjust = .5, hjust = .5, size = 10, fontface = 'bold',colour = 'darkred')


  # check header of the input file
  study <-  tryCatch(checkRequiredColumnNames(file.name,study),
                     error = function(err)
                     {
                       print.and.log(paste('Could not read file header:',err$message),'warning')
                       return(NULL)
                     }
  )

  # print.and.log("\nColumn check is done!",'info')
  # cat('----------------------------------------------------',fill = TRUE)
  if(!is.null(study))
    return(study)

}

get.file.line.count <- function(file.path)
{

  ## get line count of input file
  # this values is used to compare count of file lines with count of vatriants in input datset
  # if they  do not match , it means some lines are not read
  # skip if input file is zipped (because line number is binary line count)
  if(.QC$wc.exists){
    line <- system(sprintf("wc -l %s",  file.path), intern = TRUE)

    lineCount.regex.value <- regexec('([0-9]*)([A-z ]*)',line)[[1]]
    start.char <- lineCount.regex.value[2]
    end.char <- start.char + attributes(lineCount.regex.value)$match.length[2] - 1


    # removed to exclude stringr package
    # file.line.count <- format(as.numeric(str_match(line, "([0-9]*)([A-z ]*)")[,2]),
    #                           big.mark="," ,
    #                           scientific = FALSE)

    file.line.count <- format(as.numeric(substr(line, start.char, end.char)),
                              big.mark="," ,
                              scientific = FALSE)

  }
  else
  {
    file.line.count <- 'NA (command not accessible)'
  }

  return(file.line.count)
}

get.file.line.count_RUtils <- function(file.path)
{

  ## get line count of input file
  # this values is used to compare count of file lines with count of vatriants in input datset
  # if they  do not match , it means some lines are not read
  # skip if input file is zipped (because line number is binary line count)

  file.line.count <- tryCatch({

    line <- R.utils::countLines(file.path)

    file_ends_withNewline <- as.character(attr(line,"lastLineHasNewline"))

    line <- format(line,
                   big.mark="," ,
                   scientific = FALSE)

    return(list(line,file_ends_withNewline))
  },
  error = function(err) return("NA (command not accessible)","NA")
  )

  return(file.line.count)
}

get.study.name <- function(study) {

  return(study$file.name)
}


get.skewness.kurtosis <- function(study) {
  vec <- data.table('kurtosis' = study$kurtosis.HQ,
                    'skewness' = study$skewness.HQ,
                    'order' = study$number)

  return(vec)
}

get.precision.plot.values <- function(study) {
  vec <- data.table('SE.mean.HQ' = 1 / study$STDERR.mean.HQ,
                    'sqrt.n' = sqrt(study$MAX_N_TOTAL),
                    'order' = study$number)

  return(vec)
}

verify.files.with.user <- function(qc.study.list, user.verification)
{

  study.table <- sapply(qc.study.list,function(study)
    return(c('File Name' = basename(study$file.path),
             'Lines (including header)'= study$file.line.count,
             'missing items (first 100 lines)'= study$file.header.na,
             'Missing Columns'= paste(study$missing.Columns,collapse = ' | ')
    )
    )
  )

  if(is.null(ncol(study.table)))
    print.and.log('Error in filenames. check if filename contains space or dash or another uknown character.','fatal')

  study.table <- cbind(seq(1:ncol(study.table)),t(study.table))
  colnames(study.table)[1] <- '#'

  print.and.log(kable(study.table,format = "rst"),
                'info',
                cat= FALSE)


  if(user.verification)
  {
    # warn user that existing files will be overwritten if output folder containts files that match the input file pattern
    if(.QC$config$new_items$non.empty.output.folder)
      warning('existing files in the output folder will be overwritten!!!')

    if(menu(c("Yes", "No"),
            title=sprintf('Do you want to continue (1 = \'Yes\' | 2 = \'No\' )?'))
       == 2)
      runStopCommand('QC ended by user!')
  }
}


variable.statistics.pre.matching <- function(input.data)
{

  # count alleles for report

  .QC$thisStudy$tables$EFFECT_ALL.tbl <- input.data[VT == 1,.(.N),keyby = EFFECT_ALL]  ## count alleles
  .QC$thisStudy$tables$OTHER_ALL.tbl <- input.data[VT == 1,.(.N),keyby = OTHER_ALL]  ## count alleles

}

variable.statistics.post.matching <- function(input.data)
{
  # pallindromics
  .QC$thisStudy$palindromic.rows <- input.data[palindromic ==TRUE,.N]
  .QC$thisStudy$non.palindromic.rows <- nrow(input.data) - .QC$thisStudy$palindromic.rows


  if(is.element('CHR' , names(input.data))){
    .QC$thisStudy$tables$CHR.tbl <- input.data[,.(.N),keyby = CHR]  ## count variants per chromosome


    ### looking for missing chromosomes
    chr_range = range(as.numeric(.QC$thisStudy$tables$CHR.tbl$CHR),na.rm = T)

    chr_range = if(chr_range[2] <= 23)
      seq(chr_range[1]:23)
    else
      seq(chr_range[1]:chr_range[2])

    .QC$thisStudy$missing_chromosomes <- which(chr_range %notin% .QC$thisStudy$tables$CHR.tbl$CHR)
    ###################################

    .QC$thisStudy$tables$CHR.tbl$CHR <- as.character(.QC$thisStudy$tables$CHR.tbl$CHR)
    .QC$thisStudy$tables$CHR.tbl[is.na(CHR), CHR:= 'missing']

    .QC$thisStudy$tables$CHR.tbl$N <- thousand.sep(.QC$thisStudy$tables$CHR.tbl$N)
  }else{
    .QC$thisStudy$tables$CHR.tbl <- NA
    .QC$thisStudy$tables$CHR.tbl$CHR <- ''
    .QC$thisStudy$tables$CHR.tbl$N <- '0'
  }
  # count alleles after flipping and switching and compare to pre-match allele count
  # for report
  .QC$thisStudy$tables$EFFECT_ALL.post.matching.tbl <- input.data[VT == 1,.(.N),keyby = EFFECT_ALL]
  .QC$thisStudy$tables$OTHER_ALL.post.matching.tbl <- input.data[VT == 1,.(.N),keyby = OTHER_ALL]  ## count alleles


  if('IMPUTED' %in% colnames(input.data))
  {
    .QC$thisStudy$tables$imputed.tbl <- input.data[,.(.N),keyby = IMPUTED]  ## count variants per chromosome

    .QC$thisStudy$tables$imputed.tbl$IMPUTED <- as.character(.QC$thisStudy$tables$imputed.tbl$IMPUTED)
    .QC$thisStudy$tables$imputed.tbl[is.na(IMPUTED) , IMPUTED := 'missing/invalid']
    .QC$thisStudy$tables$imputed.tbl[IMPUTED == 0 , IMPUTED := 'genotyped']
    .QC$thisStudy$tables$imputed.tbl[IMPUTED == 1 , IMPUTED := 'imputed']
  }

  ##
  ## =================
  existing.columns <- intersect(c('PVALUE','HWE_PVAL','CALLRATE','EFF_ALL_FREQ','IMP_QUALITY','EFFECT','STDERR') ,
                                names(input.data))

  .QC$thisStudy$tables$variable.summary <- cbind(sapply(input.data[,existing.columns,with= FALSE],
                                                        function(x) return(summary(x)[1:6])))

  .QC$thisStudy$tables$variable.summary <- formatC(.QC$thisStudy$tables$variable.summary, format = 'g')


  # change effect column name to BETA or LN(OR)
  colnames(.QC$thisStudy$tables$variable.summary)[colnames(.QC$thisStudy$tables$variable.summary) == 'EFFECT'] <- .QC$config$input_parameters$effect_type_string

  ## get frequency table for multi-allelic variants
  .QC$thisStudy$tables$multi_allele_count_preProcess <- getMultiAlleleCountTbl(input.data)






  # variables that are found in standard reference file
  .QC$thisStudy$found.rows.std <- input.data[match_result != 9L & SOURCE == 'Std_ref' , .N]
  .QC$thisStudy$switched.rows.std <- input.data[match_result == 3L & SOURCE == 'Std_ref' , .N]
  .QC$thisStudy$flipped.rows.std <- input.data[match_result == 2L & SOURCE == 'Std_ref' , .N]

  # variables that are found in alternate reference file
  .QC$thisStudy$found.rows.alt <- input.data[match_result != 9L & SOURCE != 'Std_ref' , .N]
  .QC$thisStudy$switched.rows.alt <- input.data[match_result == 3L & SOURCE != 'Std_ref' , .N]
  .QC$thisStudy$flipped.rows.alt <- input.data[match_result == 2L & SOURCE != 'Std_ref' , .N]

  # variables that are not found in standard reference file
  .QC$thisStudy$not.found.rows.std <- nrow(input.data) - .QC$thisStudy$found.rows.std

  # variables that are not found in either standard or alternate reference file
  .QC$thisStudy$not.found.rows.alt <- nrow(input.data) - .QC$thisStudy$found.rows.std - .QC$thisStudy$found.rows.alt


}

variable.statistics.post.matching_HQ <- function(input.data)
{

  existing.columns <- intersect(c('PVALUE','HWE_PVAL','CALLRATE','EFF_ALL_FREQ','IMP_QUALITY','EFFECT','STDERR') ,
                                names(input.data))

  tbl <- cbind(sapply(input.data[HQ == 1,existing.columns,with= FALSE],
                                                        function(x) return(summary(x)[1:6])))

  tbl <- formatC(tbl, format = 'g')


  # change effect column name to BETA or LN(OR)
  colnames(tbl)[colnames(tbl) == 'EFFECT'] <- .QC$config$input_parameters$effect_type_string

  return(tbl)
}


addEmptyStudy_studyInstance <- function(study)
{
  faultyStudy <- new("Study")
  faultyStudy@File$file.path <- study$file.path
  faultyStudy@Successful_run <- FALSE
  faultyStudy@Counts$input.data.rowcount <- study$input.data.rowcount
  faultyStudy@Counts$rowcount.step1 <- study$rowcount.step1
  faultyStudy@Counts$rowcount.step3 <- study$rowcount.step3
  faultyStudy@Counts$found.rows <- NA
  faultyStudy@Correlations$AFcor.std_ref <- NA
  faultyStudy@Correlations$PVcor <- NA
  faultyStudy@Statistics$lambda <- NA

  .QC$StudyList@studyList <- append(.QC$StudyList@studyList , faultyStudy)
  .QC$StudyList@studyCount <- length(.QC$StudyList@studyList)
}

addEmptyStudy_pathOnly <- function(study.path)
{
  faultyStudy <- new("Study")
  faultyStudy@File$file.path <- study.path
  faultyStudy@Successful_run <- FALSE
  faultyStudy@Counts$input.data.rowcount <- NA
  faultyStudy@Counts$rowcount.step1 <- NA
  faultyStudy@Counts$rowcount.step3 <- NA
  faultyStudy@Counts$found.rows <- NA
  faultyStudy@Correlations$AFcor.std_ref <- NA
  faultyStudy@Correlations$PVcor <- NA
  faultyStudy@Statistics$lambda <- NA

  .QC$StudyList@studyList <- append(.QC$StudyList@studyList , faultyStudy)
  .QC$StudyList@studyCount <- length(.QC$StudyList@studyList)
}
