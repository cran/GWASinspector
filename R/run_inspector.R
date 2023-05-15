#' Runs the QC pipeline on a set of GWAS result files
#'
#' This is the main function of the package for running the QC algorithm on a set of GWAS result files.
#' It requires an object of class \linkS4class{Inspector} which should be created by \code{\link{setup_inspector}}.
#' Check the package vignette and tutorial for more details on this topic.
#'
#' @param inspector An instance of \linkS4class{Inspector} class. Check \code{\link{setup_inspector}} for more details.
#' @param verbose logical. If FALSE, no messages will show up in the terminal and are only saved in the log file.
#' @param test.run logical. If TRUE, only the first 1000 lines of each data file are loaded and analyzed;
#' plots and saving the cleaned output dataset are skipped. Default value is FALSE.
#' @return Reports from running the algorithm on a single or a series of GWAS result files are generated and saved.
#'
run_inspector <- function(inspector, verbose = TRUE, test.run=FALSE)
{

  if(missing(inspector))
    stop('Function arguments are not set.',call. = FALSE)

  if (!is(inspector, "Inspector"))
    stop("Object must be of class Inspector.", call. = FALSE)

  if(!validate_Inspector(inspector, printWarnings = FALSE))
    stop("Function interrupted.", call. = FALSE)

.QC$verbose <- verbose



  # reset to default setting and empty the createdenvironmet on exit (including errors)
  on.exit({
    removeFunctionVariablesFromRAM() #terminationFunctions.R
    resetDefaultSystemOptions(user.options) # sytem restore
  })
  #############################################################


  #### creating a restore point ####
  user.options<-getDefaultSystemOptions()
  ## change R Options for better performance. THIS WILL BE REVERSED AT THE END OF QC RUN!
  changeROptions()

  #### get config
  ## TODO create it, not borrow it
  .QC$config<-make_config(inspector)

  start.time <-  proc.time()
  .QC$config$new_items$starttime <- Sys.time()

  .QC$StudyList <- new("StudyList")


  inspector@start_time <- .QC$config$new_items$starttime

  #### 1
  if(.QC$verbose)
  {
    cat('\n===============================================', fill = TRUE)
    cat(sprintf('=========== %s v.%s ===========',
                    .QC$package.name,
                    .QC$script.version), fill = TRUE)
    cat('===============================================', fill = TRUE)
  }




  log.file.path <- setupLogOptions(.QC$config)
  ##==============
  print_and_log(sprintf('%s (v.%s) - Inspection started!',.QC$package.name, .QC$script.version),
                'info',
                display=.QC$config$debug$verbose)

  print_and_log(sprintf("Log file saved at \'%s\'",log.file.path))

  check_tools() #suppFunctions.R
  .QC$headerKV<-getFileHeaderKV()

  ##print and log some config variable
  printConfigVariables(.QC$config)


    ## upload reference file
  ## =====================================
  if(.QC$verbose)
    cat('\n---------- [uploading allele frequency reference file] ----------', fill = TRUE)

  # store data from a search and use it for future searches
  .QC$stored.reference.data <- NULL
  # store hIDs that were looked in the database
  .QC$searched.hIDs <- NULL

  # the following variable can be a data table or a database object
  # this is decided from the file extension
  .QC$reference.data <- uploadReferenceFile()
  checkReferenceFileIntegrity() ## check reference columns and values , exit if not satisfactory - add SOURCE , DATE_ADDED columns





  ## upload alternative reference file
  ## =====================================
  # if user has defined a file name which doesnot exist => this file will be created and filled and updated
  # if user has defined a file name that exist => it will be loaded and checked for correct columns; it will be updated during algorithm
  .QC$alt.reference.data <- data.table(numeric(0)) # create an empty data table for keeping unfound variables in each study and save as alt_ref

  if(!is.null(.QC$config$supplementaryFiles$allele_ref_alt) &&
     !is.na(.QC$config$supplementaryFiles$allele_ref_alt) &&
     .QC$config$alt_ref_file_exists)
  { # load the file if exists
    if(.QC$verbose)
      cat('\n---------- [uploading allele frequency alternative reference file] ----------', fill = TRUE)


    .QC$alt.reference.data <- tryCatch(uploadAltReferenceFile(),
                                       error= function(err)
                                       {
                                         print_and_log(paste('Error loading alternative reference dataset: ',err$message),'warning')
                                         return(data.table())
                                       }
    )


    ## check reference columns and values , return an empty data table if not satisfactory
    tryCatch(checkAltReferenceFileIntegrity(),
             error= function(err)
             {
               print_and_log(paste('Error verifying alternative reference dataset: ',err$message),'warning')
               return(data.table(numeric(0)))
             })
  }

  alt.reference.data.rowcount <- nrow(.QC$alt.reference.data)  # check row count if alt reference file. this is used to update this file if row count has increased in algorithm




  ## upload Beta (effect) reference file
  ## =====================================
  .QC$reference.data.effect <- data.table()

  if(!is.na(.QC$config$supplementaryFiles$beta_ref_std))
  {
    if(.QC$verbose)
      cat('\n---------- [uploading Beta (Effect) reference file] ----------', fill = TRUE)

    .QC$reference.data.effect <- tryCatch(uploadBetaReferenceFile(),
                                          error= function(err)
                                          {
                                            print_and_log(paste('Error loading effect size reference dataset: ',err$message),'warning')
                                            return(data.table())
                                          }
    )

    ## check reference columns , does NOT break algorithm if file is invalid
    # nrow == 0 if not valid
    ## FIXME consult about column names
    tryCatch(checkBetaReferenceFileIntegrity(),
             error= function(err)
             {
               print_and_log(paste('Error verifying effect size reference dataset: ',err$message),'warning')
               return(data.table())
             }
    )
  }



  ## TODO move it to validation process
  if(.QC$verbose)
    cat('\n---------- [analyzing input files] ----------\n', fill = TRUE)

  set_test_run_variables(test.run)

  .QC$qc.study.list <- lapply(inspector@input_files,
                              create_file_specific_config)


  # remove problematic file from list
  .QC$qc.study.list[sapply(.QC$qc.study.list, is.null)] <- NULL

  if(length(.QC$qc.study.list) == 0 )
    print_and_log('All Files Removed From Further Analysis During Header Checking!','fatal')




  .QC$file.counter <- 1

  # display study files, missing columns and line number to user
  # ask if algorithm should be done
  # verify_files_with_user(.QC$qc.study.list)

  print_and_log('\n','info',cat=TRUE)

  ## upload files one by one
  # process and check for missing ...
  # match
  # plot
  ## =====================================
  #print_and_log(mem_used(),'info',cat= FALSE)

  .QC$qc.study.list <- lapply(.QC$qc.study.list,
                              process_each_file)


  # remove null files from list
  .QC$qc.study.list[sapply(.QC$qc.study.list, is.null)] <- NULL

  if(length(.QC$qc.study.list) == 0 )
    print_and_log('All Files Removed From Further Analysis During Processing!','fatal')


  .QC$config$new_items$endtime <- Sys.time()
  inspector@end_time <- .QC$config$new_items$endtime

  # compare different studies together if there are more than one
  # draw precision , skew-kurt and effect plots
  ## =====================================
  multi_file_comparison()



  ## save alternate reference data
  # update if row number of alt ref file has increased from beginning
  ## =====================================
  if(nrow(.QC$alt.reference.data) > 0 & nrow(.QC$alt.reference.data) > alt.reference.data.rowcount)
  {
    tryCatch(save_alternate_reference(),
             error= function(err)
             {
               print_and_log(paste('Error saving alternate reference file: ',err$message),'warning')
             }
    )

  }


  ## data cleaning ...
  # remove refrence and alt reference data from RAM
  # remove effect plot of each study from RAM
  # ===========================
  if(exists("reference.data",envir = .QC))
    if(is(.QC$reference.data, "SQLiteConnection"))
      RSQLite::dbDisconnect(.QC$reference.data)

  rm(reference.data , envir = .QC)
  rm(alt.reference.data , envir = .QC)
  rm(reference.data.effect , envir = .QC)

  for(i in 1:length(.QC$qc.study.list))
  {
    .QC$qc.study.list[[i]]$effect.plot <- NULL
  }

  invisible(gc())




  #print_and_log(mem_used(),'info',cat= FALSE)
  # REPORT html and excel
  ## ==============================================
  create_report_files()


  ### post QC processes
  # 1 find significant variants



  #### -------- finishing algorithm ---------- ###
  ## =====================================



  #print_and_log(mem_used(),'info',cat= FALSE)
  ## 9- termination functions
  print_and_log(sprintf('\nFinished analyzing %s files!',length(.QC$qc.study.list)))
  print_and_log(sprintf("Run time: %s",timetaken(start.time)))# END LOG


  inspector@StudyList <- .QC$StudyList
  return(invisible(inspector))
}
