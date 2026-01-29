setupLogOptions <- function(config,
                            log.file.path = NULL) {
  ## set log file path
  if(is.null(log.file.path))
  log.file.path<-sprintf('%s/%s_log.txt',
                         config$paths$dir_output,
                         config$paths$filename_output_tag)

  if(file.exists(log.file.path))
    file.remove(log.file.path)

  log_appender(appender_file(log.file.path), namespace='GWASinspector_logger')
  log_threshold(INFO, namespace='GWASinspector_logger')
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
      log_info(message, namespace = "GWASinspector_logger")
    } else if (level == "warning") {
      log_warn(message, namespace = "GWASinspector_logger")
    } else if (level == "fatal") {
      log_fatal(message, namespace = "GWASinspector_logger")
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
