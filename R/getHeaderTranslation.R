#' Save a sample header translation table file
#'
#' This template file is used to translate a dataset`s column names (the header) into the standard names used by GWASinspector.
#' The file contains a two-column table, with the left column containing the standard column-names and the right the alternatives.
#' Both the standard and alternative columns must be fully capitalized.  This is a text file which includes most common variable/header names and can be edited according to user specifications.
#' The default name of this file is *alt_headers.txt*. configuration file should be edited if this name is changed by user (**header_translations** property).
#'
#' @param dir.path Path to the folder for saving a header-translation table file.
#' @return Saves a sample header-translation table file in the specified folder.
#' @examples
#' get.headerTranslation(tempdir())
#'
get.headerTranslation <- function(dir.path = NULL)
{
  file.name = 'alt_headers.txt'

  if(is.null(dir.path))
{
    stop('The directory for saving the sample alternate header file can not be empty.', call. = FALSE)
}
  else
  {
    dir.path <- gsub('/+$', '' , dir.path)
    file.path = file.path(dir.path , file.name)
  }



  # check if folder exists
  if(!dir.exists(dir.path))
    stop(sprintf('Directory not found \'%s\'',dir.path), call. = FALSE)

  # check if folder is writable
  if(file.access(dir.path, 2) != 0)
    stop(sprintf('Permission denied to write at \'%s\'',dir.path), call. = FALSE)

  # check if config file exists in package
  configFile <- system.file("extdata", "alt_headers.txt", package = "GWASinspector")

  if(!file.exists(configFile))
    stop('Header translation file not found in the package!', call. = FALSE)


  # save file
  # exit if file already exists
  if(!file.exists(file.path))
  {
    saveResult <- file.copy(from = configFile,
                            to = file.path,
                            overwrite = FALSE)

    # double check if file is saved
    if(saveResult & file.exists(file.path))
      message(sprintf('Header translation file saved: \'%s\'!',file.path))
    else
      warning(sprintf('Could not save sample header translation file at: \'%s\'!',file.path), call. = FALSE)
  }
  else
  {
    warning(sprintf('Header translation table file already exists: \'%s\'!',file.path), call. = FALSE)
    warning('Exited without overwriting the file!', call. = FALSE)
  }


}
