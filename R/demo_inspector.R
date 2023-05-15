#' Runs an example QC
#'
#' This function runs the QC algorithm on a fabricated GWAS result file.
#'
#' @param result.dir character. Path to the output folder for saving QC result files
#' @return QC reports from running the algorithm on a sample GWAS file are generated and saved in the specified folder.
#'
demo_inspector <- function(result.dir)
{

  if(missing(result.dir))
    stop("Function arguments are not set.",call. = FALSE)

  if(!dir.exists(result.dir))
    stop("Directory does not exist.",call. = FALSE)

  # get embedded file
  ex.config.file <- system.file("extdata", "config.ini", package = "GWASinspector")

  if(is.null(ex.config.file) || !file.exists(ex.config.file))
    stop("Sample file not found. Re-install the package")


  user.options<-getDefaultSystemOptions()
  changeROptions() ## change R Options for better performance. THIS WILL BE REVERSED AT THE END OF QC RUN!
  # reset to default setting and empty the created environmet on exit (including errors)
  on.exit({
    removeFunctionVariablesFromRAM() #terminationFunctions.R
    resetDefaultSystemOptions(user.options) # sytem restore
  })



  inspector <- new("Inspector")
  .QC$StudyList <- new("StudyList")
  .QC$file.counter <- 1
  .QC$config<-checkConfigFile(ex.config.file)
  .QC$config <- example_config(result.dir, .QC$config)

  start.time = proc.time()
  inspector@start_time <- Sys.time()

  .QC$config$new_items$starttime <- Sys.time()
  .QC$config$test.run <- FALSE
  .QC$alt.reference.data <- data.table(numeric(0))
  .QC$reference.data.effect <- data.table()


  log.file.path <- setupLogOptions(.QC$config)

  ##==============
  print_and_log("Inspection started!")

  check_tools()

  printConfigVariables(.QC$config)

  ##==============
  .QC$headerKV<-tryCatch(getFileHeaderKV(),
                         error= function(err)
                         {
                           print_and_log(paste('Error in creating header table.',err$message),'fatal')
                         }
  )


  cat('\n---------- [analyzing input files] ----------\n',fill = TRUE)
  .QC$qc.study.list <- lapply(.QC$config$paths$filename,
                             create_file_specific_config)


  .QC$reference.data <- uploadReferenceFile()
  checkReferenceFileIntegrity()

  .QC$qc.study.list <- lapply(.QC$qc.study.list,
                              process_each_file)

  invisible(gc())

  .QC$config$new_items$endtime <- Sys.time()

  create_report_files()

  inspector@end_time <- .QC$config$new_items$endtime

  print_and_log(sprintf('\nFinished analyzing %s files!',length(.QC$qc.study.list)))
  print_and_log(sprintf("Run time: %s",timetaken(start.time)))# END LOG

}
