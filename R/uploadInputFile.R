## upload only required columns of the input file for faster loading and lower memory usage
## input file must have an extension c('gz','zip','txt','dat') are included
## TODO fread function can not read files with comments
## read.table method with commnent.char='#' should be included too
## 1- rename colnames
## 2- convert alleles to upper format characters
## 3-set index on MARKER column for better searching

uploadInputFile<-function()
{

  ## FIXME if invalid values should be reported to users, setting column classes must be avoided.
  ## because by setting column classes , invalid vlaues will be converted to NA
  ## for example if 'inf' is in P-value coloumn , and this columns type is set as numeric, it will be converted to NA when loaded.
  ## other solution is to set them all as 'character' on loading

  config <- .QC$config
  study <- .QC$thisStudy
  file.path <- study$file.path
  file.extension <- study$file.extension

  ##set NA string set from config
  na.string<-config$input_parameters$na.string

  ##set column separator from config file
  sep.strings<-config$input_parameters$column_separator

  #check how the file was read
  read.method <- 'fread'

  ###
  print_and_log(sprintf("Loading Study file : '%s'",file.path),
                'info')

  if(file.extension == "zip")
  {


    input.data = tryCatch(
      {
        if(.QC$unzip.exists){

            input.data = fread(cmd = paste('unzip -cq',file.path),
                               na.strings = na.string,
                               sep = sep.strings,
                               colClasses = study$renamed.File.Columns.classes,
                               data.table = TRUE,
                               showProgress = FALSE,
                               nrows = {if(config$test.run) config$debug$test_row_count else -1},
                               fill = TRUE)

        }
        else # if RTools is not found
        {
          sep.strings <- ifelse(sep.strings == 'auto' , '' , sep.strings)
          # name if txt file inside zip file
          # FIXME how to know the exact name!?!?
          embeded.file <- paste0(tools::file_path_sans_ext(basename(file.path)),'.txt')

          read.method <- 'read.table'


            input.data =  read.table(unz(description = file.path,
                                         filename = embeded.file),
                                     sep = sep.strings,
                                     header = TRUE,
                                     na.strings = na.string,
                                     stringsAsFactors = FALSE,
                                     nrows = {if(config$test.run) config$debug$test_row_count else -1},
                                     fill = TRUE)



          close(unz(
            description = file.path,
            filename = embeded.file)
          )
        }



        # convert to data.table from read.table function
        if(!is.data.table(input.data))
          setDT(input.data)

        ####Rename columns from user-defined names to standard names
        colnames(input.data) = study$renamed.File.Columns

        input.data
      }
      ,error=function(err) {
        print_and_log(sprintf('Error found in input file - switching to method #2! %s',err$message),'warning',display=.QC$config$debug$verbose)
        return(NULL)
      }
    )
  }
  # else if (file.extension == "gz")
  # {
  #
  #   input.data = tryCatch(
  #     {
  #
  #       if(.QC$gzip.exists){
  #
  #           input.data = fread(cmd = sprintf('gzip -cd "%s"',file.path),
  #                              na.strings = na.string,
  #                              sep = sep.strings,
  #                              colClasses = study$renamed.File.Columns.classes,
  #                              data.table = TRUE,
  #                              showProgress = FALSE,
  #                              nrows = {if(config$test.run) config$debug$test_row_count else -1},
  #                              fill = TRUE)
  #
  #
  #       }else # if RTools is not found
  #       {
  #         sep.strings = ifelse(sep.strings == 'auto' , '' , sep.strings)
  #         read.method <- 'read.table'
  #
  #         if(config$test.run)
  #
  #
  #           input.data = read.table(gzfile(file.path),
  #                                   sep = sep.strings,
  #                                   header = TRUE,
  #                                   na.strings = na.string,
  #                                   stringsAsFactors = FALSE,
  #                                   fill = TRUE,
  #                                   nrows = {if(config$test.run) config$debug$test_row_count else -1})
  #
  #
  #         close(gzfile(file.path))
  #
  #       }
  #
  #
  #       #
  #
  #       # convert to data.table from read.table function
  #       if(!is.data.table(input.data))
  #         setDT(input.data)
  #
  #
  #       ####Rename columns from user-defined names to standard names
  #       if(ncol(input.data) == length(study$renamed.File.Columns))
  #         colnames(input.data) <- study$renamed.File.Columns
  #
  #       input.data
  #     }
  #     ,error=function(err) {
  #       print_and_log(sprintf('Error found in input file - switching to method #2! %s',err$message),'warning',display=.QC$config$debug$verbose)
  #       return(NULL)
  #     }
  #   )
  #
  # }
  else if(file.extension %in% c('txt','dat','csv','gz','bz2'))
  {
    # first try fread for loading the file
    input.data = tryCatch(
      {

        input.data = fread(file.path,
                            select = study$original.File.Columns.sorted,
                            na.strings = na.string,
                            sep = sep.strings,
                            colClasses = study$renamed.File.Columns.classes,
                            data.table = TRUE,
                            showProgress = FALSE,
                            fill = TRUE,
                            nrows = {if(config$test.run) config$debug$test_row_count else -1})

        colnames(input.data) <- study$renamed.File.Columns.sorted

        input.data
      }
      ,
      error=function(err)  {
        print_and_log(sprintf('Error found in input file - switching to method #2! %s',err$message),'warning',display=.QC$config$debug$verbose)

return(NULL)

      })
  }

  # read.table is used if fread creates an error
  if(is.null(input.data) || nrow(input.data) == 0)
  {
    # 'auto' is thhe default value for sep sttring in fread
    # '' should be the default value for read.table
    sep.strings <- ifelse(sep.strings == 'auto' , '' , sep.strings)
    read.method <- 'read.table'

    input.data <- tryCatch(
      {


        if(file.extension != "gz")
        {

          input.data = read.table(file.path,
                                  na.strings = na.string,
                                  sep = sep.strings,
                                  stringsAsFactors = FALSE,
                                  header = TRUE,
                                  blank.lines.skip = TRUE,
                                  nrows = {if(config$test.run) config$debug$test_row_count else -1},
                                  fill = TRUE)
        }
        else
        {
          input.data = read.table(gzfile(file.path),
                                  sep = sep.strings,
                                  header = TRUE,
                                  na.strings = na.string,
                                  stringsAsFactors = FALSE,
                                  nrows = {if(config$test.run) config$debug$test_row_count else -1},
                                  fill = TRUE)
        }

        colnames(input.data) <- study$renamed.File.Columns

        input.data
      },
      error=function(err) {
        print_and_log(sprintf('ignoring input file',err$message),'warning')
        return(NULL)
      }
    )

  }

  # return null if file could not be read
  #====================================
  if(is.null(input.data) || nrow(input.data) == 0)
    return(data.table(numeric(0)))


  # convert to data.table from read.table function if function could be read
  #====================================
  if(!is.data.table(input.data))
    setDT(input.data)
  #input.data <- as.data.table(input.data)

  ####Rename columns from user-defined names to standard names
  # colnames(input.data) <- config$new_items$final.selected.Columns ## this is for when only required columns are loaded

  if(ncol(input.data) != length(study$renamed.File.Columns.sorted))
  {
    print_and_log('Error in file name column names! column name can not be empty.','warning')
    return(NULL)
  }


  # if(read.method == 'fread')
  #   colnames(input.data) <- study$renamed.File.Columns.sorted
  # else
  #   colnames(input.data) <- study$renamed.File.Columns





  if(nrow(input.data) > 0)
  {
    dims <- dim(input.data)
    print_and_log(sprintf('Input file loaded (%s x %s)!',thousand_sep(dims[1]),dims[2]),'info')
  }
  else{
    print_and_log('Empty input data! file is ignored.','warning')
    return(NULL)
  }


  # get some information before anything changes
  input.data <- input_data_preAnalysis_report(input.data)


  ### do some simple processings, like toupper(),...
  input.data = input_data_preprocessing(input.data)
  # input.data <- input.data
  print_and_log(paste(read.method, 'was used for reading this file.'),'info' , display=.QC$config$debug$verbose)

  return(input.data)
}

input_data_preprocessing <- function(input.data) {


  input.data[, EFFECT_ALL := toupper(EFFECT_ALL)]
  input.data[, OTHER_ALL := toupper(OTHER_ALL)]

  if('CHR' %in% names(input.data)) # if chr column exists
    if(!is.numeric(input.data$CHR)) # if chr is not numeric
      input.data[, CHR := toupper(CHR)]

  return(input.data)
}


input_data_preAnalysis_report <- function(input.data)
{

  if('CHR' %in% names(input.data))
  {

    print_and_log('Variants distribution on each chromosome in input file...','info',display=.QC$config$debug$verbose)
    print_and_log(kable(input.data[, .N ,by=CHR][order(as.numeric(CHR))],format = "rst"),
                  'info',
                  cat= FALSE,
                  display= .QC$config$debug$verbose)
  }


  ##-----

  if('IMPUTED' %in% names(input.data))
  {

    tbl <- input.data[, .N ,keyby=IMPUTED]

    print_and_log('Imputation distribution in input file...','info',display=.QC$config$debug$verbose)
    print_and_log(kable(tbl,format = "rst"),
                  'info',
                  cat= FALSE,
                  display= .QC$config$debug$verbose)
  }

  return(input.data)
}
