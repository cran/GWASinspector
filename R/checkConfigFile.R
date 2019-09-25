# this file includes functions for checking and validating config file
##   included functions are:
# 1- primaryCheckConfigFile => checking only input file and directory paths
# 2- checkConfigFile => checking and validating all options in config file
# 3- evaluateListsAsStrings => changes strings in config file to R vectors


checkConfigFile <- function(config) {
  ### initiate config variables
  # newly added variable to config (NOT IN CONFIG FILE)

  # check if alt-referece file exists or not
  # it will be loaded if exists and updated and saved at the end of algorithm
  # it will be created if doesnot exist and saved at the end of algorithm
  config$alt_ref_file_exists <- FALSE

  # check if there is WRITE persmission on output folder (files and report and plots are saved)
  config$output.path.writable <- FALSE

  # check if there is WRITE persmission on reference folder ( alt reference file is saved)
  config$reference.path.writable <- FALSE
  #=====================#


  ## 0
  #### check INPUT DIRECTORY

  dir <- config$paths$dir_data

  if (is.null(dir)) {
    runStopCommand("Input directory is not set in config file!")

  } else if (!dir.exists(dir)) {
    runStopCommand(sprintf("Input directory \'%s\'wrong or not found!", dir))

  }


  ## ======================================================================
  ## 1
  ## if output DIR is not found in config file, set it equal to input DIR
  if (is.empty(config$paths$dir_output)) {
    config$paths$dir_output <- config$paths$dir_data
  }
  else if (!dir.exists(config$paths$dir_output)) {
    ## try create output DIR if set but not found
    outputDirCreated <-
      TRUE ## to double check if folder is created successfully
    outputDirCreated <- tryCatch(
      dir.create(config$paths$dir_output),

      error = function(err) {
        runStopCommand(
          sprintf(
            "Could not create outpur directory at \'%s\'! see below error:\n%s",
            config$paths$dir_output,
            err
          )
        )
        return(FALSE)
      }
    )

    ## double check if output_dir was created correctly
    if (outputDirCreated != TRUE)
      runStopCommand('Could not create output Directory!')
  }



  ## =======================================================

  ## if reference DIR is not found in config file, set it equal to input DIR
  ## TODO - OPTION: remove dir_reference from config file (reference file should be given in full path)
  if (is.empty(config$paths$dir_references)) {
    config$paths$dir_references <- config$paths$dir_data
  }
  else if (!dir.exists(config$paths$dir_references)) {
    runStopCommand("Reference directory not found!")
  }


  ## ========================================================
  ## 2
  ##change \ to / in paths
  config$paths$dir_data <-
    gsub(
      pattern = '\\\\',
      replacement = '/' ,
      x = config$paths$dir_data
    )
  config$paths$dir_output <-
    gsub(
      pattern = '\\\\',
      replacement = '/' ,
      x = config$paths$dir_output
    )
  config$paths$dir_references <-
    gsub(
      pattern = '\\\\',
      replacement = '/' ,
      x = config$paths$dir_references
    )

  ## double check to remove / at the end of path e.g '/var/user/user1/' to 'var/user/user1'
  config$paths$dir_data <-
    sub(
      pattern = '/$',
      replacement = '' ,
      x = config$paths$dir_data
    )
  config$paths$dir_output <-
    sub(
      pattern = '/$',
      replacement = '' ,
      x = config$paths$dir_output
    )
  config$paths$dir_references <-
    sub(
      pattern = '/$',
      replacement = '' ,
      x = config$paths$dir_references
    )



  ##==========================================================
  ## 3
  #### INPUT FILE
  ##check if it is set in config file
  if (is.empty(config$paths$filename))
    runStopCommand('Input file not specified in the config file!')

  ## create full path of inputfiles
  input.file.names <- list.files(
    path = config$paths$dir_data,
    pattern = config$paths$filename,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE,
    ignore.case = TRUE
  )

  ## file should have a valid extension
  input.file.names <-
    input.file.names[file_ext(input.file.names) %in% c('gz', 'zip', 'txt', 'dat', 'csv')]

  ## check if file exists
  if (length(input.file.names) == 0)
  {
    runStopCommand(
      sprintf(
        "Input file not found at \'%s\'! (['gz|zip|txt|dat|csv'] extensions are accepted)!",
        config$paths$dir_data
      )
    )
  } else{
    # save a copy of pattern that is in config file
    config$paths$original.filename <- config$paths$filename

    # put the found file list in this variable
    config$paths$filename <- input.file.names
  }


  # check if output folder is empty or not. warn user that existing files will be overwritten
  ##==========================================================
  ## create template of input files
  # like study1|study2|study3
  found.file.template <-
    paste(sapply(config$paths$filename , function(x)
      tools::file_path_sans_ext(basename(x))) ,
      collapse = '|')

  existing.files <- list.files(
    path = config$paths$dir_output,
    pattern = found.file.template,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE,
    ignore.case = TRUE
  )

  config$new_items$non.empty.output.folder <-
    ifelse(length(existing.files) > 0 , TRUE , FALSE)


  ##==========================================================

  # 4
  #### HEADER file

  ## if header file is not mentioned, upload it from package source
  # if(is.empty(config$supplementaryFiles$header_translations)){
  #
  #   # check if package default header file is present and accessible. exit if not found
  #   default.alt.header.file <- system.file("extdata", "alt_headers.txt", package = "GWASinspector")
  #
  #   if(file.exists(default.alt.header.file))
  #     config$supplementaryFiles$header_translations <- default.alt.header.file
  #   else
  #     runStopCommand('Alternate Header file is not defined and the default file is also not found in package!')
  # }
  # else
  # {
  #   ## create full path of header file
  #   config$supplementaryFiles$header_translations <- sprintf('%s/%s',
  #                                                            config$paths$dir_references,
  #                                                            config$supplementaryFiles$header_translations)
  #
  #   ## check if header file exists
  #   if(!file.exists(config$supplementaryFiles$header_translations)){
  #     runStopCommand(sprintf("Alternative Header file not found at \'%s\'!",config$supplementaryFiles$header_translations))
  #   }
  # }


  ## STOP the QC if header translation file is not mentioned or not found
  if (is.empty(config$supplementaryFiles$header_translations))
  {
    runStopCommand("Alternative Header file not specified in the configuration file!")
  }
  else {
    ## create full path of header file
    config$supplementaryFiles$header_translations <-
      sprintf(
        '%s/%s',
        config$paths$dir_references,
        config$supplementaryFiles$header_translations
      )

    ## check if header file exists
    if (!file.exists(config$supplementaryFiles$header_translations)) {
      runStopCommand(
        sprintf(
          "Alternative Header file not found at \'%s\'!",
          config$supplementaryFiles$header_translations
        )
      )
    }
  }


  # 5
  ###reference file
  if (is.empty(config$supplementaryFiles$allele_ref_std)) {
    runStopCommand("Reference file not set!")
  }
  else {
    ## create full path of reference file
    config$supplementaryFiles$allele_ref_std <- sprintf(
      '%s/%s',
      config$paths$dir_references,
      config$supplementaryFiles$allele_ref_std
    )

    ## check if reference file exists
    if (!file.exists(config$supplementaryFiles$allele_ref_std)) {
      runStopCommand(
        sprintf(
          "Reference file not found at %s!",
          config$supplementaryFiles$allele_ref_std
        )
      )
    }

    ## stop if a database file with sqlite extension is specified but required package is not installed
    if (tools::file_ext(config$supplementaryFiles$allele_ref_std) == 'sqlite' &&
        !check.rsqlite.package(installed.packages()[, 1]))
      runStopCommand(
        'RSQLite package is not installed on this computer!\nyou can use the table version of reference file (RData/RDS/txt).'
      )
  }


  config$supplementaryFiles$allele_ref_std_population <-
    checkConfigParameters(
      toupper(config$supplementaryFiles$allele_ref_std_population),
      'list',
      'EUR' ,
      c('COMMON', 'EAS', 'AMR', 'AFR', 'EUR', 'SAS')
    )



  ### alternative reference file.

  # if it is not selected by user, the process of matching with alt-ref is skipped.
  if (is.empty(config$supplementaryFiles$allele_ref_alt)) {
    config$supplementaryFiles$allele_ref_alt <- NA
  } else{
    # if this file is defined by uesr, it will be uploaded and updated after each study file is analyzed.
    # if it does not exist , it will be created after first study file, updated with next studies in the row and saved at exit

    ## create full path of alt reference file
    config$supplementaryFiles$allele_ref_alt <- sprintf(
      '%s/%s',
      config$paths$dir_references,
      config$supplementaryFiles$allele_ref_alt
    )

    if (file.exists(config$supplementaryFiles$allele_ref_alt))
      config$alt_ref_file_exists <- TRUE
  }


  ### starandard reference file for EFFECT (BETA)
  if (is.empty(config$supplementaryFiles$beta_ref_std)) {
    config$supplementaryFiles$beta_ref_std <- NA
  } else{
    ## create full path of alt reference file
    config$supplementaryFiles$beta_ref_std <- sprintf('%s/%s',
                                                      config$paths$dir_references,
                                                      config$supplementaryFiles$beta_ref_std)

    ## check if beta reference file exists
    if (!file.exists(config$supplementaryFiles$beta_ref_std)) {
      runStopCommand(
        sprintf(
          "Beta (effect) reference file not found at \'%s\'!",
          config$supplementaryFiles$beta_ref_std
        )
      )
    }
  }


  # 6
  ###read parameters
  ##set NA string set from config
  if (is.empty(config$input_parameters$na.string))
    config$input_parameters$na.string <- c("NA", "nan", "NaN", ".")
  else
    config$input_parameters$na.string <-
    checkConfigParameters(config$input_parameters$na.string,
                          'text',
                          c("NA", "nan", "NaN", "."))

  ## 7
  ##set column separator from config file
  if (is.empty(config$input_parameters$column_separator))
    config$input_parameters$column_separator <- 'auto'
  else
    config$input_parameters$column_separator <-
    checkConfigParameters(config$input_parameters$column_separator,
                          'text',
                          'auto')

  ##set effect type string  from config
  if (is.empty(config$input_parameters$effect_type))
    config$input_parameters$effect_type <- 'BETA'
  else
    config$input_parameters$effect_type <-
    checkConfigParameters(
      toupper(config$input_parameters$effect_type),
      'list',
      default = 'BETA',
      range = c("BETA", "OR")
    )

  # this item is used for display in report and plots
  # if effect = OR , it will be converted to BEta by ln(OR)
  # in this case, ln(OR) is displayed in plots and reports instead of BETA
  config$input_parameters$effect_type_string <-
    ifelse(config$input_parameters$effect_type == 'BETA' ,
           'BETA',
           'Ln(OR)')



  ## output tag infront of inputfilename

  config$paths$filename_output_tag <-
    checkConfigParameters(config$paths$filename_output_tag,
                          'text',
                          'QC')


  config$plot_specs$plot_title <-
    checkConfigParameters(config$plot_specs$plot_title,
                          'text',
                          'none')

  ## 8
  ## Imputed alleles - T or F
  if (is.empty(config$input_parameters$imputed_T))
    config$input_parameters$imputed_T <- "TRUE|T|YES|Y"
  else
    config$input_parameters$imputed_T <-
    paste(evaluateListsAsStrings('imputed_T', config$input_parameters$imputed_T),
          collapse = '|')

  if (is.empty(config$input_parameters$imputed_F))
    config$input_parameters$imputed_F <- "FALSE|F|NO|N"
  else
    config$input_parameters$imputed_F <-
    paste(evaluateListsAsStrings('imputed_F', config$input_parameters$imputed_F),
          collapse = '|')

  ##IMP_QUALITY
  config$filters$minimal_impQ_value  <-
    checkConfigParameters(config$filters$minimal_impQ_value, 'numeric',-0.5)
  config$filters$maximal_impQ_value  <-
    checkConfigParameters(config$filters$maximal_impQ_value, 'numeric', 1.5)

  ##HQ Filters ####
  #####

  config$filters$HQfilter_FRQ  <-
    checkConfigParameters(config$filters$HQfilter_FRQ, 'numeric', 0.01)
  config$filters$HQfilter_HWE  <-
    checkConfigParameters(config$filters$HQfilter_HWE, 'numeric', 1e-6)
  config$filters$HQfilter_cal  <-
    checkConfigParameters(config$filters$HQfilter_cal, 'numeric', 0.95)
  config$filters$HQfilter_imp  <-
    checkConfigParameters(config$filters$HQfilter_imp, 'numeric', 0.3)


  # config$input_parameters$calculate_missing_p

  config$input_parameters$calculate_missing_p  <-
    checkConfigParameters(config$input_parameters$calculate_missing_p,
                          'logical', FALSE)

  #remove chromosomes from QC
  config$remove_chromosomes$remove_X  <-
    checkConfigParameters(config$remove_chromosomes$remove_X,
                          'logical', FALSE)
  config$remove_chromosomes$remove_Y  <-
    checkConfigParameters(config$remove_chromosomes$remove_Y,
                          'logical', FALSE)
  config$remove_chromosomes$remove_M  <-
    checkConfigParameters(config$remove_chromosomes$remove_M,
                          'logical', FALSE)
  config$remove_chromosomes$remove_XY  <-
    checkConfigParameters(config$remove_chromosomes$remove_XY,
                          'logical', FALSE)

  config$filters$threshold_diffEAF  <-
    checkConfigParameters(config$filters$threshold_diffEAF, 'numeric', 0.15)

  ##Threshoold for creating histograms and plots,
  #  < 1 is converted to percent and then to real value in histogram plot function



  # cut off P-value for manhattan and p-p plot
  config$plot_specs$plot_cutoff_p  <-
    checkConfigParameters(config$plot_specs$plot_cutoff_p, 'numeric', 0.01)


  config$plot_specs$make_plots  <-
    checkConfigParameters(config$plot_specs$make_plots, 'logical', FALSE)



  ##output items
  config$output_parameters$save_final_dataset <-
    checkConfigParameters(config$output_parameters$save_final_dataset,
                          'logical',
                          FALSE)

  config$output_parameters$gzip_final_dataset <-
    checkConfigParameters(config$output_parameters$gzip_final_dataset,
                          'logical',
                          FALSE)

  config$output_parameters$html_report <-
    checkConfigParameters(config$output_parameters$html_report,
                          'logical',
                          TRUE)


  config$output_parameters$object_file <-
    checkConfigParameters(config$output_parameters$object_file,
                          'logical',
                          TRUE)

  # converted to upppercase for consistency
  config$output_parameters$out_header <-
    checkConfigParameters(
      toupper(config$output_parameters$out_header),
      'list',
      'STANDARD' ,
      c('GENABEL', 'GWAMA', 'PLINK', 'META', 'STANDARD', 'GCTA')
    )

  #
  config$output_parameters$out_sep <-
    checkConfigParameters(config$output_parameters$out_sep,
                          'list',
                          '\t' ,
                          c('\t', ';', '|', ','))
  #
  config$output_parameters$out_na <-
    checkConfigParameters(config$output_parameters$out_na,
                          'text',
                          'NA')

  #
  config$output_parameters$out_dec <-
    checkConfigParameters(config$output_parameters$out_dec,
                          'list',
                          '.' ,
                          c('.'))


  #
  config$plot_specs$graphic_device <-
    checkConfigParameters(
      tolower(config$plot_specs$graphic_device),
      'list',
      'png' ,
      c('png', 'jpeg', 'tiff')
    )



  # debug options
  # this items are removed from user file and kept for debugger
  config$debug$verbose  <-
    checkConfigParameters(config$debug$verbose,
                          'logical', FALSE)

  config$debug$save_pre_modification_file  <-
    checkConfigParameters(config$debug$save_pre_modification_file,
                          'logical', FALSE)

  config$debug$reduced.AF.plot  <-
    checkConfigParameters(config$debug$reduced.AF.plot,
                          'logical', TRUE)

  config$debug$test_row_count  <-
    checkConfigParameters(config$debug$test_row_count,
                          'numeric', 1000)

  ## END oF VALIDATION


  # ======
  # save plots as png or jpeg or jpeg
  # if neither exists , set make_plot as false
  # config$graphic.device is set in the following function
  # the reslut value for graphic.device is T/F

  # following parameters are set:
  # config$graphic.device
  # config$plot_specs$make_plots is set to FALSE if device is not available
  # .QC$graphic.device 'png' 'jpeg' 'tiff'
  # .QC$img.extension  '.png' '.jpeg' '.tiff'
  if (config$plot_specs$make_plots)
    config <- set.graphic.device(config)
  else
  {
    # set as default values if user does not want the plots
    config$graphic.device <- TRUE
    .QC$graphic.device  <- 'png'
    .QC$img.extension <- '.png'
  }
  # ============


  ## multi file comparison plot paths
  config$paths$effsizePlotPath <-
    sprintf('%s/Checkgraph_effect-size%s',
            config$paths$dir_output,
            .QC$img.extension)

  config$paths$precisionPlotPath <-
    sprintf('%s/Checkgraph_precision%s',
            config$paths$dir_output,
            .QC$img.extension)

  config$paths$skew_kurt <- sprintf('%s/Checkgraph_skew_kurt%s',
                                    config$paths$dir_output,
                                    .QC$img.extension)

  config$paths$html.report <- sprintf(
    '%s/%s_%s',
    config$paths$dir_output,
    config$paths$filename_output_tag,
    'Report.html'
  )

  config$paths$txt.report <- sprintf(
    '%s/%s_%s',
    config$paths$dir_output,
    config$paths$filename_output_tag,
    'Report.txt'
  )

  config$paths$xlsx.report <- sprintf(
    '%s/%s_%s',
    config$paths$dir_output,
    config$paths$filename_output_tag,
    'Report.xlsx'
  )

  config$paths$save_pre_modification_file <- sprintf(
    '%s/%s_%s',
    config$paths$dir_output,
    config$paths$filename_output_tag,
    'early_file.txt'
  )

  ## check folder permission
  checkFolderPermission(config)

  return(config)
}


checkConfigFileSections <- function(config)
{
  if (all(names(config) %in% c(
      'paths',
      'supplementaryFiles',
      'input_parameters',
      'output_parameters',
      'remove_chromosomes',
      'plot_specs',
      'filters',
      'debug')))
  return(TRUE)
  else
    return(FALSE)

}



evaluateListsAsStrings <- function(input.name, input.string)
{
  output.list <- tryCatch(
    eval(parse(text = input.string)),
    warning = function(err) {
      print.and.log(
        sprintf(
          'error evaluating \'%s\' parameter in config file! should be written as a vector.',
          input.name
        ),
        'fatal'
      )
      return(NULL)
    },
    error = function(err) {
      print.and.log(
        sprintf(
          'error evaluating \'%s\' parameter in config file! should be written as a vector.',
          input.name
        ),
        'fatal'
      )
      return(NULL)
    }
  )

  ## double check for correct conversion
  if (is.vector(output.list))
    return(output.list)
  else
    print.and.log(
      sprintf(
        '\'%s\' parameter should be written as a vector in config file!',
        input.name
      ),
      'fatal'
    )

}



## TODO not complete yet
printConfigVariables <- function(config) {
  print.and.log(sprintf('Input directory: \'%s\'', config$paths$dir_data),
                'info')

  print.and.log(sprintf('Output directory: \'%s\'', config$paths$dir_output),
                'info')

  print.and.log(sprintf('References directory: \'%s\'', config$paths$dir_references),
                'info')

  print.and.log(
    sprintf(
      'Alternate header file: \'%s\'',
      config$supplementaryFiles$header_translations
    ),
    'info'
  )

  print.and.log(
    sprintf(
      'Allele Frequency Reference file: \'%s\'',
      config$supplementaryFiles$allele_ref_std
    ),
    'info'
  )

  if (!is.na(config$supplementaryFiles$allele_ref_alt))
    print.and.log(
      sprintf(
        'Allele Frequency Alternative Reference file: \'%s\'',
        config$supplementaryFiles$allele_ref_alt
      ),
      'info'
    )

  if (!is.na(config$supplementaryFiles$beta_ref_std))
    print.and.log(
      sprintf(
        'Effect Size Reference file: \'%s\'',
        config$supplementaryFiles$beta_ref_std
      ),
      'info'
    )


  print.and.log(paste('Operating CPU cores:', .QC$cpu.core), 'info')


  if (!.QC$config$graphic.device)
    print.and.log('No graphic devices are available. Plotting will be skipped!',
                  'warning')


}



## if the input paramete can be converted to correct type, return the values ,else set to default value
checkConfigParameters <-
  function(parameter, type , default , range = NULL) {
    ## convert to uppercase for consistency
    # if(is.character(parameter))
    #    parameter <- toupper(parameter)



    if (type == 'numeric')
    {
      if (is.null(parameter) ||
          parameter == '' ||
          is.na(as.numeric(parameter)))
        output <- default
      else
        output  <- as.numeric(parameter)
    }

    if (type == 'logical')
    {
      if (is.null(parameter) ||
          parameter == '' ||
          is.na(as.logical(parameter)))
        output <- default
      else
        output  <- as.logical(parameter)
    }


    if (type == 'text')
    {
      if (is.null(parameter) ||
          parameter == '')
        output <- default
      else{
        parameter <- gsub('"', '', parameter) # remove " from variables
        output  <- parameter
      }
    }


    if (type == 'list')
    {
      if (is.null(parameter) ||
          length(parameter) == 0 || parameter %notin% range)
        output <- default
      else{
        parameter <- gsub('"', '', parameter)# remove " from variables
        output <- parameter
      }
    }


    return(output)
  }



set.test.run.variables <- function(test.run)
{
  if (test.run)
  {
    .QC$config$output_parameters$save_final_dataset <- FALSE
    .QC$config$plot_specs$make_plots <- FALSE
    .QC$config$test.run <- TRUE

  }
  else
  {
    .QC$config$test.run <- FALSE
  }

}
