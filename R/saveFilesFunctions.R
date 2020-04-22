#this file includes functions for saving data sets and logging the job
# 1- saveDataSet()  => general function
# 2- saveDataSet.final() => saving matched data at the end of algorithm IF TRUE in config file





# Save a dataset to output directory
saveDataSet <- function(dataset,
                        file.path,
                        columnSeparator="\t",
                        naValue = "NA",
                        decValue = '.',
                        zipped = FALSE,
                        ordered = FALSE)
{

  # ORDER THE  output rows on CHR-position
  # order on RS number if CHR is not present
  if(ordered == TRUE){
    if(is.element('CHR', colnames(dataset)) && is.element('POSITION', colnames(dataset)) )
      dataset <- dataset[order(CHR,POSITION)]
    else if(is.element('MARKER', colnames(dataset)))
      dataset <- dataset[order(MARKER)]
  }


  ## 1- save dataset
  # zip if selected by user and Rtools is installed
  if(zipped){
    # much slower tha fwrite
    tryCatch(
      # write.table(dataset,
      #             dec = decValue,
      #             na = naValue,
      #             quote = FALSE,
      #             sep = columnSeparator,
      #             row.names = FALSE,
      #             file = gzfile(file.path ,compression = 3)),
      fwrite(x = dataset,
             file = file.path,
             sep = columnSeparator,
             na = naValue,
             dec = decValue,
             quote = FALSE,
             compress = "auto"),
      error = function(err)
      {
        print.and.log(sprintf("Error saving file! check below message for more information:\n%s",
                              err),
                      'warning')
      }
    )
  }else {
    tryCatch(
      fwrite(x = dataset,
             file = file.path,
             sep = columnSeparator,
             na = naValue,
             dec = decValue,
             quote = FALSE),

      error = function(err)
      {
        print.and.log(sprintf("Error saving file! check below message for more information:\n%s",
                              err),
                      'warning')
      }
    )
  }


  ## 2- create log
  if(zipped)
    print.and.log(sprintf("A zipped data file is saved as '%s' with '%s' rows!",
                          file.path,
                          thousand.sep(nrow(dataset))),
                  'info')
  else
    print.and.log(sprintf("A data file is saved as '%s' with '%s' rows!",
                          file.path,
                          thousand.sep(nrow(dataset))),
                  'info')

}



## TODO use 'filename_output' for saved file
saveDataSet.final<-function(dataset)
{
  config <- .QC$config

  if(config$output_parameters$save_final_dataset){

    ## remove reference addded columns from final dataset
    ## MULTOI_ALLELIC column is needed if it is generated
    ## =========================================

    requiredColNames <- .QC$thisStudy$renamed.File.Columns.sorted

    if(is.element('MULTI_ALLELIC',names(dataset)) && config$output_parameters$add_column_multiallelic == TRUE)
      requiredColNames <- c(requiredColNames,'MULTI_ALLELIC')

    if(is.element('highDiffEAF',names(dataset)) && config$output_parameters$add_column_AFmismatch == TRUE)
      requiredColNames <- c(requiredColNames,'highDiffEAF')

    if(is.element('AF',names(dataset)) && config$output_parameters$add_column_AF == TRUE)
      requiredColNames <- c(requiredColNames,'AF')

    if(is.element('REF_RSID',names(dataset)) && config$output_parameters$add_column_rsid == TRUE)
      requiredColNames <- c(requiredColNames,'REF_RSID')



    dataset <- subset(dataset , select = requiredColNames)



    # check if there were characters in chromosme column that shold be deconverted to characters
    ## =========================================
    if(.QC$thisStudy$character.chromosome)
      dataset <- deconvert.column.CHR(dataset)

    # make sure position is not saved as scientific number ,e.g. 8e-6
    dataset$POSITION<-format(dataset$POSITION,scientific = FALSE,trim = TRUE)


    ## rename columns based on user preference
    #PLINK, GWAMA, GENABEL,...
    ## =========================================
    selected.header.format <- config$output_parameters$out_header

    if(selected.header.format == 'GENABEL'){

      if('build' %notin% colnames(data))
        dataset[, build := "unknown"]


      if("pgc" %notin% colnames(data)) {
        if('sebeta' %in% colnames(data)) # this is not in our standard columns
          dataset[, pgc := pchisq((EFFECT/(sebeta * sqrt(.QC$thisStudy$lambda) ))^2,
                                  1, lower.tail=FALSE)]
        else
          dataset[, pgc := NA]
      }

      if(!"lambda.estimate" %notin% colnames(data))
        dataset[, lambda.estimate := .QC$thisStudy$lambda]

      if(!"lambda.se" %notin% colnames(data))
        dataset[, lambda.se := NA]

    }

    new.col.names <- changeColumnNames(colnames(dataset), selected.header.format)
    colnames(dataset) <- new.col.names


    ## =========================================

    # add zip extension to output file path
    if(config$output_parameters$gzip_final_dataset)
      .QC$thisStudy$output.path <- paste0(.QC$thisStudy$output.path, '.gz')


    saveDataSet(dataset,
                .QC$thisStudy$output.path,
                columnSeparator = config$output_parameters$out_sep,
                naValue = config$output_parameters$out_na,
                decValue = config$output_parameters$out_dec,
                zipped = config$output_parameters$gzip_final_dataset,
				ordered = .QC$config$output_parameters$ordered)
  }else{
    print.and.log('Saving final dataset is skipped!','warning',display=.QC$config$debug$verbose)
  }
}


save.NA.Dataset <- function(input.data,input.data.backup) {

  #get the list of noncrucial columns that should be checked for NA values
  nonCrucialColumns <- getNonCrucialColumnNames()
  crucialColumns <- getCrucialColumnNames.onFileAnalysis()

  # get column names of input file
  current.col.names <- colnames(input.data)

  #
  row.count <- nrow(input.data)

  # we need to check the columns for NA values:
  # all column should not be NA
  # only required columns that are NOT all NA are checked for NA items
  wantedColumnsList <- sapply(nonCrucialColumns, function(x)
  {
    if(x %in% current.col.names){
      if(input.data[is.na(eval(parse(text = x))), .N] == row.count)
      {
        print.and.log(sprintf('Column \'%s\' removed from improbable_variant file because all where NA!',x)
                      ,'warning',display=.QC$config$debug$verbose)
        return(FALSE)
      }
      else
      {
        return(TRUE)
      }
    }else
    {
      return(FALSE)
    }
  })

  # get the list of required columns
  unwantedColumnsList <- names(wantedColumnsList[which(wantedColumnsList == FALSE)])
  wantedColumnsList <- names(wantedColumnsList[which(wantedColumnsList == TRUE)])



  na.rows <- which(!complete.cases(input.data[,wantedColumnsList,with=FALSE]))

  ## na.rows <- apply(input.data,1,anyNA) => very slow
  # na.rows <- which(!complete.cases(input.data))

  if(length(na.rows) > 0){

    # only require the first 100 rows of crucial columns and those that are not totally NA
    input.data <-  head(input.data[na.rows,union(crucialColumns,wantedColumnsList),with =FALSE],100)
    input.data.backup <-  head(input.data.backup[na.rows,union(crucialColumns,wantedColumnsList),with =FALSE],100)

    INVALID_COLUMN<-apply(input.data,
                          1,
                          function(x)
                          {
                            c<- which(is.na(x))
                            paste(names(c),collapse = '|')
                          })
    input.data.backup <- cbind(input.data.backup,INVALID_COLUMN)
    #	[filename_output]_SNPs_improbable_values.txt
    saveDataSet(input.data.backup,
                .QC$thisStudy$SNPs_improbable_values.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec,
				ordered = .QC$config$output_parameters$ordered)

  }
}


save.and.remove.unusable.variants <- function(input.data,input.data.backup) {

  #TODO break to two functions
  ## === 1 - first check allele columns => data is saved as SNPS_invalid_alleles.txt
  ## === 2 - check other crucial columns => data is saved as SNPS_removed.txt
  ## === 3 - remove from dataset

  ## 1
  ## save invalid allele variants as dataset

  invalid.allele.index.union <- union( .QC$thisStudy$column.INVALID.list$EFFECT_ALL,
                                       .QC$thisStudy$column.INVALID.list$OTHER_ALL)

  if(length(invalid.allele.index.union) > 0)
  {

    # [filename_output]_SNPs_invalid_allele.txt
    saveDataSet(input.data.backup[head(invalid.allele.index.union,100)],
                file.path = .QC$thisStudy$SNPs_invalid.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec,
				ordered = .QC$config$output_parameters$ordered)
  }

  ## find rows with NA in allele columns
  ## data is removed from dataset
  # invalid alleles are set to NA in processCOlumns() => process.column.EFFECT_ALL()
  # so, invalid.allele.index is a subset of missing.allele.index
  #


  ## 2
  ## find rows with NA in crucial columns (not alleles)
  ## first 100 is saved
  ## data is removed from dataset

  crucial.columns <- getCrucialColumnNames.onFileAnalysis()
  # allele columns are checked in previous phase (SNPs_invalid_allele) and not needed here
  crucial.columns <- crucial.columns[!crucial.columns %in% c('EFFECT_ALL','OTHER_ALL')]

  # add EFFECT to crucial columns if not present
  # this was removed from initial set because effect could be either beta , or
  if(!is.element('EFFECT',crucial.columns))
    crucial.columns <- append(crucial.columns,'EFFECT')


  crucial.col.index <- match(crucial.columns ,names(input.data))


  # get the list of NA values in each crucial column
  # column operaion is faster than row operation
  NA.list <- sapply(crucial.col.index, function(x)
    {
      which(is.na(input.data[,x,with=FALSE]))
    }
  )


  # check the union of missing rows indexes for each crucial column
  # this result set shows which rows have at least one missing values in a crucial column
  # NA.list.union <- unlist(NA.list[1])
  # for(i in 1:length(NA.list)){
  #   NA.list.union<-union(NA.list.union,
  #                        unlist(NA.list[[i]]))
  # }

  NA.list.union <- Reduce(union,NA.list) ## easier than above loop


  if(length(NA.list.union) > 0)
  {
    sample.data <- input.data.backup[head(NA.list.union,100)] ## 100 samples are saved from backup file, before invalid items are converted to NA

    # [filename_output]_SNPs_removed.txt
    saveDataSet(sample.data,
                file.path = .QC$thisStudy$SNPs_removed.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec,
				ordered = .QC$config$output_parameters$ordered)
  }

  ## 3
  # find union of missing values in crucial clumns with either allele columns for removing from dataset
  total.NA.list.union <- Reduce(union, list(NA.list.union,
                                            .QC$thisStudy$column.NA.list$OTHER_ALL,
                                            .QC$thisStudy$column.NA.list$EFFECT_ALL
  )
  )

  if(length(total.NA.list.union) > 0)
  {
    print.and.log(sprintf('%s variants that missed a crucial value were removed from dataset (step 1)!',thousand.sep(length(total.NA.list.union))),
                  'warning',display=.QC$config$debug$verbose)
    input.data <- input.data[!total.NA.list.union]
    .QC$thisStudy$missing.crucial.rowcount <- length(total.NA.list.union)
  }


  return(input.data)
}

save.significant.variants <- function(input.data) {

  sig.variants <- input.data[HQ == TRUE &
                               PVALUE < 5e-08 &
                               startsWith(MARKER,'rs'),
                             list(CHR,POSITION,MARKER,PVALUE)]

  saveDataSet(sig.variants,
              .QC$thisStudy$SNPs_significant.path,
              columnSeparator = .QC$config$output_parameters$out_sep,
				ordered = .QC$config$output_parameters$ordered)

}



# this is for debugging purposes only
save.pre_modification.file <- function(pre.modification.file)
{

  if(.QC$config$debug$save_pre_modification_file)
    data.table::fwrite(pre.modification.file ,
                       file = .QC$config$paths$save_pre_modification_file ,
                       sep = '\t')
}
