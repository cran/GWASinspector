#' Unifies column separators in all data rows
#'
#' It is common that column separator is different between the header and the rest of the file. This can be problematic for reading the file correctly and should be treated before loading the data.
#' This function is based on Rtools and *awk* function to read and save a new file with equal column separators between all rows of data.
#'
#' @param input.file path to the input file
#' @param output.file path of the outfile path; if not defined , a tag is added to input file
#' @param sep what character to be used add separator; tab character is default
#' @return A new text file is generated from the old one with a consistent column separator character.
#'
reformat.columns <- function(input.file, output.file = NULL, sep='\t')
{
  if(!check.awk())
  {
    stop('This function requires gawk command!\nRTools not found or PATH variable does not contain it! see: https://cran.rstudio.com/bin/windows/Rtools/',call. = FALSE)
  }



  if(is.null(input.file) || input.file =='' || !file.exists(input.file))
    stop('input file not found!', call. = FALSE)

  input.file.name <- basename(input.file)
  input.file.path <- dirname(input.file)

  if(is.null(output.file))
    output.file <- sprintf('%s/%s_%s',input.file.path,'formatted',input.file.name)


  # only the header row is required
  data <- fread(file = input.file,
                     sep = 'auto',
                     nrows = 0,
                     header = TRUE)

  # how many columns are found?
  col.count <- dim(data)[2]
  awk_command <- sprintf('gawk \'BEGIN{OFS="\t"}{print $%s}\' %s > %s',
                         paste(1:col.count,collapse = ',$'),
                         input.file,
                         output.file)

  shell(awk_command)
  print(sprintf('Reformatted file saved to : \'%s\' !',output.file))
  print(sprintf('Column Count : \'%s\' !',col.count))
}



## check if folders are writable
checkFolderPermission <- function(config) {

  ## check output folder permission
  if(file.access(config$paths$dir_output, 2 ) != 0)
    runStopCommand(sprintf("Algorithm can not save output files and plots at \'%s\'! check folder permission!",
                           config$paths$dir_output))



  ## check reference folder permission- required for saving alt ref file
  if(file.access(config$paths$dir_references, 2 ) != 0)
    runStopCommand(sprintf("Algorithm can not save reference files at \'%s\'! check folder permission!",
                           config$paths$dir_references))

}



# Make sure awk is installed. This is part of Rtools.
check.awk <- function() {

  # installed <- invisible(system('gawk --v') == 0)

  if(Sys.which('gawk') != '')
    return(TRUE)
  else
    return(FALSE)


}

# Make sure wc command is installed. This is part of Rtools.
check.wc <- function() {

  if(Sys.which('wc') != '')
    return(TRUE)
  else
    return(FALSE)

}


check.java <- function() {

  if(Sys.which('java') != '')
    return(TRUE)
  else
    return(FALSE)

}

# Make sure gzip command is installed. This is part of Rtools.
check.gzip <- function() {
  if(Sys.which('gzip') != '')
    return(TRUE)
  else
    return(FALSE)
}

# Make sure gzip command is installed. This is part of Rtools.
check.unzip <- function() {
  if(Sys.which('unzip') != '')
    return(TRUE)
  else
    return(FALSE)
}

check.xlsx.package <- function(existing.packages)
{
  if (is.element('xlsx', existing.packages))
    return(TRUE)
  else
    return(FALSE)
}

check.rsqlite.package <- function(existing.packages)
{
  if (is.element('RSQLite', existing.packages))
    return(TRUE)
  else
    return(FALSE)
}

get.R.version <- function(){

  name <- trimws(gsub('\\([^)]*\\)',
               R.Version()$version.string,
               replacement = '',
               ignore.case = TRUE) )

  arch <- ifelse(grepl(pattern = '64',x =  R.Version()$arch) , '64bit', '32bit')

  return(paste(name,arch,sep = ' - '))
}



get.parallel.core.count <- function()
{

  # return 1 , if core value is not numeric or is set to 1 or lower
    if(is.null(.QC$cpu.core) || !is.numeric(.QC$cpu.core) || .QC$cpu.core <= 1  )
    return(1)

  #===================================
  # selected core value is more than 1

  #Sys.info()['sysname']
  if(grepl(pattern = 'win' , tolower(.Platform$OS.type)))
  {
    print.and.log("Parallel processing is not avaialble in Windows OS! single core will be used.",'warning')
    return(1)

  }

  if(!.QC$parallel.package.exists)
  {
    print.and.log("Parallel package is not installes in R! single core will be used.",'warning')
    return(1)

  }


  pc.cores <- parallel::detectCores()


  if(.QC$cpu.core <= pc.cores)
    return(.QC$cpu.core)
  else
    return(pc.cores)

}


check.parallel.package <- function(existing.packages)
{
  if (is.element('parallel', existing.packages))
    return(TRUE)
  else
    return(FALSE)
}


check.kableExtra.package <- function(existing.packages)
{
  if (is.element('kableExtra', existing.packages))
  {
    # requireNamespace("kableExtra")
    # require(kableExtra)
    # added to global environmentas as new variable if the package exists
   # `%>%` <- NULL
   # `%>%` <<- kableExtra::`%>%`
    return(TRUE)
  }
  else
    return(FALSE)
}


check.rJava.package <- function(existing.packages)
{
  if (is.element('rJava', existing.packages))
    return(TRUE)
  else
    return(FALSE)
}


check.ggplot2.version <- function(existing.packages)
{
  if (is.element('ggplot2', existing.packages))
    return(paste(packageVersion('ggplot2'),collapse = '.'))
  else
    return('not installed')
}


check.pandoc <- function()
{
  if (rmarkdown::pandoc_available())
    return(TRUE)
  else
    return(FALSE)
}

get.OS <- function()
{
  sys <- Sys.info()
  return(paste(sys['sysname'], sys['release']))
}


removeDuplicatedLines <- function(input.data) {
  print.and.log('looking for duplicated lines ...','info')

  dup_lines <- which(duplicated(input.data))
  .QC$thisStudy$dup_lines_count <- length(dup_lines)

  if(.QC$thisStudy$dup_lines_count > 0)
  {

    tbl <- input.data[dup_lines, .N ,keyby=CHR]

    print.and.log('duplicated lines distribution in input file...','info',display=.QC$config$debug$verbose)
    print.and.log(kable(tbl,format = "rst"),
                  'info',
                  cat= FALSE,
                  display= .QC$config$debug$verbose)


    print.and.log(sprintf('Duplicated lines found in input file: %s lines were removed.',
                          format(.QC$thisStudy$dup_lines_count,big.mark = ',',scientific = FALSE)),
                  'warning',display=.QC$config$debug$verbose)
    input.data <- input.data[!dup_lines,]
  }

  return(input.data)
}
