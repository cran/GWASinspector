#' Copies the template header translation table to the local machine
#'
#' This template file is used to translate a dataset`s column names (the header) into the standard names used by GWASinspector.
#' The file contains a two-column table, with the left column containing the standard column-names and the right the alternatives.
#' Both the standard and alternative columns must be fully capitalized. This is a text file which includes most common variable/header names and can be edited according to user specifications.
#' The default filename is \strong{alt_headers.txt}.
#'
#' @param dir.path Path to the folder for saving a header-translation table file.
#' @return Copies a sample header-translation table in the specified folder.
#'
get_headerTranslation <- function(dir.path)
{
  file.name = 'alt_headers.txt'

  if(missing(dir.path))
    stop('No directory specified.', call. = FALSE)

  if(!is(dir.path,"character"))
    stop('Argument is not in a correct format. A string parameter is required.', call. = FALSE)
  else
    dir.path <- gsub('/+$', '' , dir.path)



  # check if folder exists
  if(!dir.exists(dir.path))
    stop(sprintf('Directory not found \'%s\'',dir.path), call. = FALSE)

  # check if folder is writable
  if(file.access(dir.path, 2) != 0)
    stop(sprintf('Permission denied to write at \'%s\'',dir.path), call. = FALSE)

  # check if header file exists in package
  headerFile <- system.file("extdata", "alt_headers.txt", package = "GWASinspector")

  if(!file.exists(headerFile))
    stop('Header translation file not found in the package!', call. = FALSE)



  file.name = 'alt_headers.txt'
  file.path = file.path(dir.path , file.name)


  header.file.exists  <- TRUE

  while (header.file.exists) {

    if (file.exists(file.path))
    {
      file.name = sprintf('alt_headers_%s.txt', sample(1:100000, 1))
      file.path = file.path(dir.path , file.name)
    } else
      header.file.exists <- FALSE

  }


    saveResult <- file.copy(from = headerFile,
                            to = file.path,
                            overwrite = FALSE)


    # double check if file is saved
    if(saveResult & file.exists(file.path))
      message(sprintf('Header translation file saved: \'%s\'!',file.path))
    else
      warning(sprintf('Could not save sample header translation file at: \'%s\'!',file.path), call. = FALSE)


}
