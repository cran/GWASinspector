#' Save a sample configuration file
#'
#' Save a sample configuration file for running GWASinspector package. This templates should be edited and then used for running the QC.
#' User can save multiple copies to be used for different sets of files.
#' Default name is *config.ini*, which can be changed while saving the file or afterwards.
#'
#' @param dir.path Path to the folder for saving a sample configuration file.
#' @return Saves a sample configuration file (config.ini) in the specified folder.
#' @examples
#' get.config(tempdir())
#'
get.config <- function(dir.path = NULL)
{

  file.name = 'config.ini'

  if(is.null(dir.path))
  {
    stop('The directory for saving the sample configuration file can not be empty.', call. = FALSE)
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
  configFile <- system.file("extdata", "config.ini", package = "GWASinspector")

  if(!file.exists(configFile))
    stop('Config file not found in the package!', call. = FALSE)


  # save file
  # exit if file already exists
  if(!file.exists(file.path))
  {
    saveResult <- file.copy(from = configFile,
                            to = file.path,
                            overwrite = FALSE)

    # double check if file is saved
    if(saveResult & file.exists(file.path))
      message(sprintf('Sample config file saved: \'%s\'!',file.path))
    else
      warning(sprintf('Could not save sample config file at: \'%s\'!',file.path), call. = FALSE)
  }
  else
  {
    warning(sprintf('Config file already exists: \'%s\'!',file.path), call. = FALSE)
    warning('Exited without overwriting the file!' , call. = FALSE)
  }


}
