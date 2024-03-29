setupLogOptions <- function(config,
                            log.file.path = NULL) {
  ## set log file path
  if(is.null(log.file.path))
  log.file.path<-sprintf('%s/%s_log.txt',
                         config$paths$dir_output,
                         config$paths$filename_output_tag)

  if(file.exists(log.file.path))
    file.remove(log.file.path)

  flog.appender(appender.file(log.file.path), name='GWASinspector_logger')
  flog.threshold(INFO)
  return(log.file.path)
}


print_and_log <- function(message,
                          level = "info",
                          cat = TRUE,
                          display = TRUE) {


  ## conver to FALSE if called from outside the paclage pipeline
  .QC$verbose <- ifelse(is.null(.QC$verbose),FALSE,.QC$verbose)

  if (.QC$verbose & cat & display) {
    cat(paste("-", toupper(level), " : ", message, sep = " "), fill = TRUE)
  }
  else if (.QC$verbose & display) {
    print(message, quote = FALSE)
  }

#  tryCatch({
    if (level == "info") {
      flog.info(msg = message, name = "GWASinspector_logger")
    } else if (level == "warning") {
      flog.warn(msg = message, name = "GWASinspector_logger")
    } else if (level == "fatal") {
      flog.fatal(msg = message, name = "GWASinspector_logger")
      runStopCommand()
    }
  # },
  # error = function(err) {}
  # )
}


runStopCommand<-function(message = NULL, call=FALSE)
{
  # this is hanldled by on.exit() command on main function
  # removeFunctionVariablesFromRAM() #terminationFunctions.r
  if(is.null(message))
    stop("Function is interrupted due to an error! check the log file for details.", call. = call)
  else
    stop(message, call. = call)
}







writeTXTreport <- function(message){

  write(message,
        file= .QC$thisStudy$txt.report.path,
        append=TRUE)

}

writeFileComparisonTXTreport <- function(message){

  write(message,
        file= .QC$config$paths$txt.report,
        append=TRUE)

}
