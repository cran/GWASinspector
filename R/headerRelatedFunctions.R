## this file includes functions for reading header file
# functions included are:
# 1- checkRequiredColumnNames => check if input file has the required columns and rename colnames
# 2- getCrucialColumnNames => colnames that are necessary and will cause the function to end
## TODO this function items (necessary colnames) might be added to config file

# 3- getRequiredColumnNames => colnames that are selected from input file (not all columns will be loaded)
## TODO this function items (colnames) might be added to config file

# 4- getFileHeaderKV => check if header file has two columns with no duplicates
# create a key-value variable for checking and renaming header columns

##below two functions get a hashset of keyvalues from header file for searching and replacing input colnames
# 5- replaceInputStringSetFromHashSet => compare all values of a list to hashset and return the keys
# 6- replaceInputStringFromHashSet => compare a single value to hashset and return the key



checkRequiredColumnNames <- function(inputFile, study){

  file.extension<-study$file.extension #file_ext(inputFile)

  ##set NA string set from config
  na.string <- .QC$config$input_parameters$na.string

  ##set column separator from config file
  sep.strings <- .QC$config$input_parameters$column_separator

  # if(file.extension %in% c('gz','zip')){
  #
  #   inputFilePath<-sprintf("gzip -dc %s", inputFile)
  #   data<-fread( inputFilePath,
  #                nrows=100,
  #                header = TRUE,
  #                fill=TRUE,
  #                sep='auto')

  if(file.extension == "zip") {

    # name if txt file inside zip file
    # FIXME how to know the exact name!?!?
    embeded.file <- paste0(tools::file_path_sans_ext(basename(inputFile)),'.txt')

    sep.strings <- ifelse(sep.strings == 'auto' , '' , sep.strings)

    data <- read.table(unz(description = inputFile,
                           filename = embeded.file),
                       sep = sep.strings ,
                       header = TRUE,
                       na.strings = na.string,
                       nrows=100,
                       fill = TRUE)

    close(unz(
      description = inputFile,
      filename = embeded.file)
    )

  } else if (file.extension == "gz"){

    sep.strings <- ifelse(sep.strings == 'auto' , '' , sep.strings)

    data <- read.table(gzfile(inputFile),
                       sep = sep.strings,
                       header = TRUE,
                       nrows = 100,
                       na.strings = na.string,
                       fill = TRUE)

    close(gzfile(inputFile))

  }else if(file.extension %in% c('txt','dat','csv'))
  {
    sep.strings <- ifelse(sep.strings == 'auto' , '' , sep.strings)

    data <- read.table(file = inputFile,
                       sep = sep.strings,
                       header = TRUE,
                       nrows = 100,
                       na.strings = na.string,
                       fill = TRUE)

    # data<-fread(inputFile,
    #             nrows=100,
    #             header = TRUE,
    #             fill = TRUE,
    #             sep='auto',
    #             na.strings = na.string)

  }

  ###
  # convert to data.table from read.table function
  if(!is.data.table(data))
    data <- as.data.table(data)


  wanted.columns<-getRequiredColumnNames()
  original.File.Columns.upper<-toupper(colnames(data)) ## upper is required for translating column names
  original.File.Columns <- colnames(data) ## for loading dataset and saving final data

  # headerKV<-getFileHeaderKV()

  ###
  renamed.File.Columns<-replaceInputStringSetFromHashSet(.QC$headerKV, original.File.Columns.upper)
  # check if file header is correctly translated, return NULL if there is a duplicate in header names
  if(is.null(renamed.File.Columns))
    return(NULL)




  #### check if all required comuns exist ####
  ## if not, script is stopped

  crucial.columns<-getCrucialColumnNames.onFileLoading()
  missing.crucial.column.indexes<-which(crucial.columns %notin% renamed.File.Columns)

  if(length(missing.crucial.column.indexes) > 0){
    print.and.log(sprintf("File will be ignored. columns \'%s\' were not found.",
                          paste(crucial.columns[missing.crucial.column.indexes], collapse = '|')),
                  'warning')
    addEmptyStudy(inputFile)
    return(NULL)
  }


  ## check if file already has an EFFECT column. this is problematic and can not be continued
  if(is.element('EFFECT',original.File.Columns.upper))
    print.and.log('Column EFFECT will be renamed to BETA!'
                          ,'warning')

  if(is.element('EFFECT',renamed.File.Columns))
    print.and.log('Input file can not have a column that is named EFFECT.'
                  ,'fatal')


  ## check if file has the selected effect column , either BETA or OR
  # if(!is.element(.QC$config$input_parameters$effect_type,renamed.File.Columns))
  #   print.and.log(sprintf('You have selected \"%s\" as the effect column, which is not found in the input file.' ,
  #                         .QC$config$input_parameters$effect_type ),'fatal')



  ### check other columns and notify that one is missing
  missing.columns<-wanted.columns[which(wanted.columns  %notin% renamed.File.Columns )]

  if(length(missing.columns) > 0 )
  {
    study$missing.Columns <- missing.columns # add to config fo report


    if('PVALUE' %in% missing.columns){
      study$missing.PVALUE.column <- TRUE

      print.and.log('PVALUE column is missing. calculated values will be used!',
                    'warning')
      print.and.log('Pvalue correlation plot will not be saved!',
                    'warning')
    }

    # CHR is a crucial column now
#     if('CHR' %in% missing.columns)
#       print.and.log('Manhattan plot will not be saved (CHR column is missing)!',
#                     'warning')

    if('EFF_ALL_FREQ' %in% missing.columns)
      print.and.log('Allele frequency plots will not be saved (EFF_ALL_FREQ column is missing)!',
                    'warning')


  }

  ##----
  ## find colmn classes of input file to use when loading the file
  renamed.File.Columns.classes <- sapply(renamed.File.Columns, function(x) setColumnClass(x))
  names(renamed.File.Columns.classes) <- colnames(data)
  # remove nulls that mean unmatched columns
  renamed.File.Columns.classes = renamed.File.Columns.classes[!sapply(renamed.File.Columns.classes,
                                                                      is.null)]
  #----

  ## saving variables in config file
  col.index.all <- c(1:ncol(data))
  wanted.columns.index<-na.omit(match(wanted.columns,renamed.File.Columns))
  ## add column names to config file
  study$original.File.Columns <- original.File.Columns
  study$renamed.File.Columns <- renamed.File.Columns

  ##used for loadinh the input file, columns are loaded as required order - unknown columns are pushed to the end
  study$original.File.Columns.sorted <- c(original.File.Columns[wanted.columns.index],
                                          original.File.Columns[col.index.all[-wanted.columns.index]])

  #used for renaming the original file header after loading it and saving datasets - it is sorted
  study$renamed.File.Columns.sorted <- c(renamed.File.Columns[wanted.columns.index],
                                         renamed.File.Columns[col.index.all[-wanted.columns.index]])

  ## index of columns , first the known and then the unknowm
  study$wanted.columns.index <- c(wanted.columns.index,
                                  col.index.all[-wanted.columns.index])

  # classes for loading the columns- for the 15 known classes, these items are set. for the rest it is based on fread function
  study$renamed.File.Columns.classes <- unlist(renamed.File.Columns.classes)


  # VERY IMPORTANT: number of NA values in the first 100 lines of file
  # this is so important because it shows if columns are correctly seperated or not
  # i.e. header row may be separated based on space character while other rows are separated by tab
  study$file.header.na <- count.NA(data)

  return(study)
}

count.NA <- function(input.data){

  na.count <- length(which(is.na(input.data)))
  d <- dim(input.data)
  na.percent <- calculatePercent(na.count,  (d[1] * d[2]))

  return(na.percent)
}

##algorithm stops if one of these columns are missing
getCrucialColumnNames.onFileLoading<-function()
{
  wanted.columns<-c(
    "CHR",
    "POSITION",
    "EFFECT_ALL",
    "OTHER_ALL",
   # "EFFECT",  added as a separate test (either BETA or OR)
    "STDERR"
  )

  if(.QC$config$input_parameters$effect_type == 'BETA')
    wanted.columns <- append(wanted.columns,'BETA')
  else  if(.QC$config$input_parameters$effect_type == 'OR')
    wanted.columns <- append(wanted.columns,'OR')

  return(wanted.columns)
}

getCrucialColumnNames.onFileAnalysis<-function()
{
  wanted.columns<-c(
    "CHR",
    "POSITION",
    "EFFECT_ALL",
    "OTHER_ALL",
    "EFFECT",
    "STDERR"
  )

  return(wanted.columns)
}

getNonCrucialColumnNames<-function()
{
  wanted.columns<-c("STRAND",
                    "PVALUE",
                    "EFF_ALL_FREQ",
                    "HWE_PVAL",
                    "IMP_QUALITY",
                    "IMPUTED",
                    "CALLRATE",
                    "N_TOTAL",
                    "MARKER")

  return(wanted.columns)
}

getRequiredColumnNames<-function()
{
  wanted.columns<-c("CHR",
                    "MARKER",
                    "POSITION",
                    "STRAND",
                    "EFFECT_ALL",
                    "OTHER_ALL",
                   # "EFFECT",
                    "STDERR",
                    "PVALUE",
                    "EFF_ALL_FREQ",
                    "HWE_PVAL",
                    "IMP_QUALITY",
                    "IMPUTED",
                    "CALLRATE",
                    "N_TOTAL")

  if(.QC$config$input_parameters$effect_type == 'BETA')
    wanted.columns <- append(wanted.columns,'BETA')
  else  if(.QC$config$input_parameters$effect_type == 'OR')
    wanted.columns <- append(wanted.columns,'OR')

  return(wanted.columns)
}


setColumnClass<-function(column){
  switch(column,
         "CHR"='character',
         "MARKER" ='character' ,
         "POSITION"='numeric',
         "STRAND"='character',
         "EFFECT_ALL"='character',
         "OTHER_ALL"='character',
         "EFFECT"='numeric',
         "STDERR"='numeric',
         "PVALUE"='numeric',
         "EFF_ALL_FREQ"='numeric',
         "HWE_PVAL"='numeric',
         "IMP_QUALITY"='numeric',
         "IMPUTED"='character',
         "CALLRATE"='numeric',
         "N_TOTAL"='numeric'
         ## TODO,'CHR' check this for unknown columns
  )
}

getFileHeaderKV<-function()
{
  headerTable <- read.table(file = .QC$config$supplementaryFiles$header_translations,
                            sep='',
                            header = FALSE,
                            stringsAsFactors = FALSE)
# TODO delete
#   ###checking header file
#   if(ncol(headerTable) != 2L) {
#     print.and.log(sprintf('\'headerTable\' should have two columns but has %s!',ncol(headerTable)),
#                   'fatal')
#   }
#
#   if(any(duplicated(headerTable[ ,2]))) {
#     print.and.log('Duplicated items found in header table!',
#                   'fatal')
#   }

  ##creating hashset of header elements
  headerKV <- hash::hash(keys = as.matrix(headerTable[, 2]),
                         values = as.matrix(headerTable[, 1]))
  # flog.info("header KV file loaded!", name= .QC$log.name.main)

  print.and.log('Header translation table created!',
                'info')

  return(headerKV)
}


### create header hashset
replaceInputStringSetFromHashSet <- function(hashset , inputStringSet)
{
  if (!hash::is.hash(hashset))
  {
    stop('hashset not loaded correctly')
  }

  returnSet<-NULL
  if (all(hash::has.key(inputStringSet, hashset)))
  {
    returnSet<-hash::values(hashset, inputStringSet)
  } else
  {
    returnSet<-sapply(inputStringSet, function(x)
      replaceInputStringFromHashSet(hashset, x))
  }

  returnSet <- as.vector(returnSet)

  tbl <- data.table('header' = inputStringSet , 'translated_header' = returnSet)

  print.and.log('Header translation in input file...','info',display=.QC$config$debug$verbose)
  print.and.log(kable(tbl,format = "rst"),
                'info',
                cat= FALSE,
                display= .QC$config$debug$verbose)


  if(any(duplicated(returnSet)))
  {

    print.and.log('duplicated values are found in translated header!','warning',display=.QC$config$debug$verbose)
    return(NULL)
  }
  else
  {
    return(as.vector(returnSet))
  }
}

replaceInputStringFromHashSet <- function(hashset, inputString)
{
  if (!hash::is.hash(hashset))
  {
    stop('hashset not loaded correctly')
  }

  if (hash::has.key(inputString, hashset))
  {
    hash::values(hashset, inputString)
  } else
  {
    inputString
  }
}
