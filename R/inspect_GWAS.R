#' Runs the QC algorithm
#'
#' This is the main function for running the algorithm on GWAS result files.
#'
#' @param config.file character. Path to a configuration (.ini) file to configurate the QC. For a sample configuration file, see get.config(). For details on how to edit the config file, see the tutorial.
#' @param user.verification logical. If TRUE, the algorithm will pause and ask the user to verify that it has selected the correct input files for QC. Default vaule is FALSE.
#' @param test.run logical. if TRUE, only the first 1000 lines of each data file are loaded and analysed, and no plots or final dataset is produced.
#' @return QC reports from running the algorithm on a single or a series of GWAS result files are generated and saved.
#'
inspect<-function(config.file = NULL, user.verification = FALSE, test.run = FALSE){


  #### 0 check if config file exists ####
  ## =====================================
  if(is.null(config.file))
  {
    runStopCommand('Configuration file not provided.')

  } else if(!file.exists(config.file))
  {
    runStopCommand(sprintf('Configuration file not found! check if the path is correct: %s', config.file))

  }

  ## 1
  # =============================================
  #### creating a restore point ####
  # 'terminationFunctions.R'
  user.options<-getDefaultSystemOptions()
  changeROptions() ## change R Options for better performance. THIS WILL BE REVERSED AT THE END OF QC RUN!

  # check if use.verification parameter os a logical parameter or not
  if(!is.logical(user.verification))
    runStopCommand('user verification parameter should be a logical value! (TRUE/FALSE)')


  # reset to default setting and empty the createdenvironmet on exit (including errors)
  on.exit({
    removeFunctionVariablesFromRAM() #terminationFunctions.R
    resetDefaultSystemOptions(user.options) # sytem restore
  })


  ## start Algorithm
  # =============================================

  message('\n=============================================')
  message(sprintf('=========== %s v.%s ===========',
              .QC$package.name,
              .QC$script.version))
  message('=============================================')



  .QC$config<-read.ini(filepath = config.file)



  message('\n---------- [validating Config file] ----------\n')

  ## check config file

  # check if file is in correct format
  if(!checkConfigFileSections(.QC$config))
    runStopCommand('Config file is not in a correct format! run get.config() to obtain a template.')

  ## stops the process if paths or directories not found
  .QC$config<-checkConfigFile(.QC$config) # checkConfigFile.R


  ## =====================================
  start.time = proc.time()
  .QC$config$new_items$starttime <- Sys.time()
  # assign(x = "config" , value = config , envir = .QC)


  log.file.path <- setupLogOptions(.QC$config)
  ##==============
  print.and.log("Inspection started!",
                'info')


  #print.and.log(mem_used(),'info',cat= FALSE)
  # check if awk and wc and gzip commands are installed and accessible.
  # awk is used for synchronizing each row separator character
  # wc is used for counting file lines (check if all lines of file is read successfully)
  # gzip is used for reading zip files
  # TODO this is part of RTools and user is urged to install it.
  cores <- 1
  .QC$cpu.core <- cores # checked and changed in check.tools() function
  check.tools() #suppFunctions.R


  print.and.log(sprintf("Log file saved at \'%s\'",log.file.path),
                'info')

  ##print and log some config variable
  printConfigVariables(.QC$config)


  ## reading header table from alt_header txt file
  # ==============================================
  .QC$headerKV<-tryCatch(getFileHeaderKV(),
                         error= function(err)
                         {
                           print.and.log(paste('Error in cerating header table.',err$message),'fatal')
                         }
  )

  ## initiate study file specific variables
  # including plot and file paths and saved data file paths and report variables
  ## =====================================
  # print.and.log('===============================================','info')
  #
  # print.and.log('---------- [analyzing input files] ----------\n','info')

  message('\n---------- [analyzing input files] ----------\n')

  #print.and.log(mem_used(),'info',cat= FALSE)

  # change config variables if this is a test run
  # new varibles are created
  # .QC$config$test.run <- TRUE
  # .QC$config$test.row.count <- 1000
  set.test.run.variables(test.run)




  if(.QC$cpu.core == 1)
    .QC$qc.study.list <- lapply(.QC$config$paths$filename,
                                create.file.specific.config)
  else
    .QC$qc.study.list <- parallel::mclapply(X= .QC$config$paths$filename,
                                            FUN = create.file.specific.config,
                                            mc.cores = .QC$cpu.core)




  # remove problematic file from list
  .QC$qc.study.list[sapply(.QC$qc.study.list, is.null)] <- NULL

  if(length(.QC$qc.study.list) == 0 )
    print.and.log('All Files Removed From Further Analysis During Header Checking!','fatal')




  .QC$file.counter <- 1

  # display study files, missing columns and line number to user
  # ask if algorithm should be done
  verify.files.with.user(.QC$qc.study.list, user.verification)





  ## upload reference file
  ## =====================================
  message('\n\n---------- [uploading allele frequency reference file] ----------')

  # the following variable can be a data table or a database object
  # this is decided from the file extension
  .QC$reference.data <- uploadReferenceFile()
  checkReferenceFileIntegrity() ## check reference columns and values , exit if not satisfactory - add SOURCE , DATE_ADDED columns





  ## upload alternative reference file
  ## =====================================
  # if user has defined a file name which doesnot exist => this file will be created and filled and updated
  # if user has defined a file name that exist => it will be loaded and checked for correct columns; it will be updated during algorithm
  .QC$alt.reference.data <- data.table(numeric(0)) # create an empty data table for keeping unfound variables in each study and save as alt_ref

  if(!is.na(.QC$config$supplementaryFiles$allele_ref_alt) & .QC$config$alt_ref_file_exists)
  { # load the file if exists
    message('\n\n---------- [uploading allele frequency alternative reference file] ----------')


    .QC$alt.reference.data <- tryCatch(uploadAltReferenceFile(),
                                       error= function(err)
                                       {
                                         print.and.log(paste('Error loading alternative reference dataset: ',err$message),'warning')
                                         return(data.table())
                                       }
    )


    ## check reference columns and values , return an empty data table if not satisfactory
    tryCatch(checkAltReferenceFileIntegrity(),
             error= function(err)
             {
               print.and.log(paste('Error verifying alternative reference dataset: ',err$message),'warning')
               return(data.table(numeric(0)))
             })
  }

  alt.reference.data.rowcount <- nrow(.QC$alt.reference.data)  # check row count if alt reference file. this is used to update this file if row count has increased in algorithm




  ## upload Beta (effect) reference file
  ## =====================================
  .QC$reference.data.effect <- data.table()

  if(!is.na(.QC$config$supplementaryFiles$beta_ref_std))
  {
    message('\n\n---------- [uploading Beta (Effect) reference file] ----------')
    .QC$reference.data.effect <- tryCatch(uploadBetaReferenceFile(),
                                          error= function(err)
                                          {
                                            print.and.log(paste('Error loading effect size reference dataset: ',err$message),'warning')
                                            return(data.table())
                                          }
    )

    ## check reference columns , does NOT break algorithm if file is invalid
    # nrow == 0 if not valid
    ## FIXME consult about column names
    tryCatch(checkBetaReferenceFileIntegrity(),
             error= function(err)
             {
               print.and.log(paste('Error verifying effect size reference dataset: ',err$message),'warning')
               return(data.table())
             }
    )
  }


  print.and.log('\n','info',cat=TRUE)

  ## upload files one by one
  # process and check for missing ...
  # match
  # plot
  ## =====================================
  #print.and.log(mem_used(),'info',cat= FALSE)

  if(.QC$cpu.core == 1)
    .QC$qc.study.list <- lapply(.QC$qc.study.list,
                                process.each.file)
  else
    .QC$qc.study.list <- parallel::mclapply(X = .QC$qc.study.list,
                                            FUN = process.each.file,
                                            mc.cores = .QC$cpu.core)


  # remove null files from list
  .QC$qc.study.list[sapply(.QC$qc.study.list, is.null)] <- NULL

  if(length(.QC$qc.study.list) == 0 )
    print.and.log('All Files Removed From Further Analysis During Processing!','fatal')


  # compare different studies together if there are more than one
  # draw precision , skew-kurt and effect plots
  ## =====================================
  multi.file.comparison()



  ## save alternate reference data
  # update if row number of alt ref file has increased from beginning
  ## =====================================
  if(nrow(.QC$alt.reference.data) > 0 & nrow(.QC$alt.reference.data) > alt.reference.data.rowcount)
  {
    tryCatch(save.alternate.reference(),
             error= function(err)
             {
               print.and.log(paste('Error saving alternate reference file: ',err$message),'warning')
             }
    )

  }


  ## data cleaning ...
  # remove refrence and alt reference data from RAM
  # remove effect plot of each study from RAM
  # ===========================

  rm(reference.data , envir = .QC)
  rm(alt.reference.data , envir = .QC)
  rm(reference.data.effect , envir = .QC)

  for(i in 1:length(.QC$qc.study.list))
  {
    .QC$qc.study.list[[i]]$effect.plot <- NULL
  }

  invisible(gc())

  .QC$config$new_items$endtime <- Sys.time()

  ## FOR DEBUGGING
  ## FIXME DELETE
  # saveRDS(.QC$qc.study.list , paste0(.QC$config$paths$dir_output,'/studyVariables.rds'))
  # saveRDS(.QC$config , paste0(.QC$config$paths$dir_output,'/configVariables.rds'))


  #print.and.log(mem_used(),'info',cat= FALSE)
  # REPORT html and excel
  ## ==============================================
  create.report.files()


  ### post QC processes
  # 1 find significant variants



  #### -------- finishing algorithm ---------- ###
  ## =====================================



  #print.and.log(mem_used(),'info',cat= FALSE)
  ## 9- termination functions
  print.and.log(sprintf('\nFinished analyzing %s files!',length(.QC$qc.study.list)))
  print.and.log(sprintf("Run time: %s",timetaken(start.time)))# END LOG



  ## end Algorithm using on.exit() at top of page


}
