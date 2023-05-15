#' Copies the template configuration file to the local machine
#'
#' This templates should be edited and then used for setting up and running the QC pipeline.
#' Default filename is \strong{config.ini}.
#'
#' @param dir.path Path to the folder for saving a sample configuration file.
#' @return Copies a sample configuration file (config.ini) in the specified folder.
#'
get_config <- function(dir.path)
{
  if(missing(dir.path))
  {
    stop('No directory specified.',
         call. = FALSE)
  }

  if(!is(dir.path,"character"))
    stop('Argument is not in a correct format. A string parameter is required.', call. = FALSE)
  else
    dir.path <- gsub('/+$', '' , dir.path)

  # check if folder exists
  if (!dir.exists(dir.path))
    stop(sprintf('Directory not found \'%s\'', dir.path), call. = FALSE)

  # check if folder is writable
  if (file.access(dir.path, 2) != 0)
    stop(sprintf('Permission denied to write at \'%s\'', dir.path),
         call. = FALSE)

  # check if config file exists in package
  configFile <-
    system.file("extdata", "config.ini", package = "GWASinspector")

  if (!file.exists(configFile))
    stop('Config file not found in the package!', call. = FALSE)


  # save file
  # exit if file already exists

  file.name = 'config.ini'
  file.path = file.path(dir.path , file.name)

  config.file.exists  <- TRUE

  while (config.file.exists) {
    if (file.exists(file.path))
    {
      file.name = sprintf('config_%s.ini', sample(1:100000, 1))
      file.path = file.path(dir.path , file.name)
    } else
      config.file.exists <- FALSE

  }


  saveResult <- file.copy(from = configFile,
                          to = file.path,
                          overwrite = FALSE)

  # double check if file is saved
  if (saveResult &
      file.exists(file.path))
  {
    message(sprintf('Sample config file saved: \'%s\'!', file.path))
    return(invisible(file.path))
  }
  else
    warning(sprintf('Could not save sample config file at: \'%s\'!', file.path),
            call. = FALSE)


}
