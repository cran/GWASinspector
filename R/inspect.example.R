#' Runs an example QC
#'
#' This function runs the QC algorithm on a fabricated GWAS result file.
#'
#' @param result.dir character. Path to the output folder for saving QC result files
#' @return QC reports from running the algorithm on a sample GWAS file are generated and saved in the specified folder.
#' @examples
#' \donttest{
#' inspect.example(tempdir())
#' }
inspect.example<-function(result.dir = NULL)
{
  if(is.null(result.dir) || !dir.exists(result.dir))
    stop("Please provide an existing folder for saving the output files.",call. = FALSE)


  # reset to default setting and empty the created environmet on exit (including errors)
  on.exit({
    removeFunctionVariablesFromRAM() #terminationFunctions.R
    resetDefaultSystemOptions(user.options) # sytem restore
  })

  user.options<-getDefaultSystemOptions()
  changeROptions() ## change R Options for better performance. THIS WILL BE REVERSED AT THE END OF QC RUN!


  .QC$config <- example.config(result.dir)
  .QC$config<-checkConfigFile(.QC$config)
  start.time = proc.time()
  .QC$config$new_items$starttime <- Sys.time()
  .QC$config$test.run <- FALSE
  .QC$alt.reference.data <- data.table(numeric(0))
  .QC$reference.data.effect <- data.table()


  log.file.path <- setupLogOptions(.QC$config)

  ##==============
  print.and.log("Inspection started!",'info')

  check.tools()

  printConfigVariables(.QC$config)

  ##==============
  .QC$headerKV<-tryCatch(getFileHeaderKV(),
                         error= function(err)
                         {
                           print.and.log(paste('Error in cerating header table.',err$message),'fatal')
                         }
  )


  message('\n---------- [analyzing input files] ----------\n')
  .QC$qc.study.list <- lapply(.QC$config$paths$filename,
                             create.file.specific.config)


  .QC$reference.data <- uploadReferenceFile()
  checkReferenceFileIntegrity()

  .QC$qc.study.list <- lapply(.QC$qc.study.list,
                              process.each.file)

  invisible(gc())

  .QC$config$new_items$endtime <- Sys.time()

  create.report.files()


  print.and.log(sprintf('\nFinished analyzing %s files!',length(.QC$qc.study.list)))
  print.and.log(sprintf("Run time: %s",timetaken(start.time)))# END LOG

}
