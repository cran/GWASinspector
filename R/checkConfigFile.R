# this file includes functions for checking and validating config file
##   included functions are:
# 1- primaryCheckConfigFile => checking only input file and directory paths
# 2- checkConfigFile => checking and validating all options in config file
# 3- evaluateListsAsStrings => changes strings in config file to R vectors


checkConfigFile <- function(config.file) {
  if (!file.exists(config.file))
    runStopCommand(sprintf('Configuration file not found at: %s', config.file))
  else
    config <- read.ini(filepath = config.file)


  # check if file is in correct format
  if (!checkConfigFileSections(config))
    runStopCommand('Config file is not in a correct format! run get_config() to obtain a template.')

  ### initiate config variables
  # newly added variable to config (NOT IN CONFIG FILE)


  #=====================#


  if (is_empty(config$paths$filename))
    config$paths$filename <- ".+"

  if (is_empty(config$paths$dir_data))
    runStopCommand("Data directory not set!")

  if (is_empty(config$paths$dir_output))
    runStopCommand("Output directory not set!")


  if (is_empty(config$paths$dir_references))
    runStopCommand("Reference directory not set!")

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



  ##### file names

  input.file.names <- list.files(
    path = config$paths$dir_data,
    pattern = config$paths$filename,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE,
    ignore.case = TRUE)


  # check file order string
  if (is_empty(config$input_parameters$file_order_string))
    config$input_parameters$file_order_string <- ""
  else
    config$input_parameters$file_order_string <-
    paste(
      evaluateListsAsStrings('file_order_string', config$input_parameters$file_order_string),
      collapse = '|'
    )

  ## order the files
  orderIndex <- grep(config$input_parameters$file_order_string,
                     input.file.names,
                     ignore.case = T)

  if(length(orderIndex) > 0)
    input.file.names <- c(input.file.names[orderIndex],
                          input.file.names[-orderIndex])

  ## keep valid extensions
  config$paths$input_files <-input.file.names[file_ext(input.file.names) %in% c('gz', 'zip', 'txt', 'dat', 'csv', 'bz2')]



  # check if output folder is empty or not. warn user that existing files will be overwritten
  ##==========================================================
  ## create template of input files
  # like study1|study2|study3


  # TODO
  # found.file.template <-
  #   paste(sapply(config$paths$filename , function(x)
  #     tools::file_path_sans_ext(basename(x))) ,
  #     collapse = '|')
  #
  # existing.files <- list.files(
  #   path = config$paths$dir_output,
  #   pattern = found.file.template,
  #   full.names = TRUE,
  #   recursive = FALSE,
  #   all.files = FALSE,
  #   ignore.case = TRUE
  # )

  # config$new_items$non.empty.output.folder <-
  #   ifelse(length(existing.files) > 0 , TRUE , FALSE)


  ##==========================================================

  # 4
  #### HEADER file

  ## STOP the QC if header translation file is not mentioned or not found
  if (is_empty(config$supplementaryFiles$header_translations))
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

  }


  # 5
  ###reference file
  if (is_empty(config$supplementaryFiles$allele_ref_std)) {
    runStopCommand("Reference file not set!")
  }
  else {
    ## create full path of reference file
    config$supplementaryFiles$allele_ref_std <- sprintf(
      '%s/%s',
      config$paths$dir_references,
      config$supplementaryFiles$allele_ref_std
    )
  }


  # config$supplementaryFiles$allele_ref_std_population <-
  #   checkConfigParameters(
  #     toupper(config$supplementaryFiles$allele_ref_std_population),
  #     'list',
  #     config$supplementaryFiles$allele_ref_std_population , # return the input string if not one of the below
  #     c('COMMON', 'EAS', 'AMR', 'AFR', 'EUR', 'SAS')
  #   )

  ## check if population string is correct
  if(grepl(pattern = "sqlite",x = config$supplementaryFiles$allele_ref_std )
     && is_empty(config$supplementaryFiles$allele_ref_std_population))
    runStopCommand("allele_ref_std_population parameter can not be empty.")

  if(!is.element(toupper(config$supplementaryFiles$allele_ref_std_population),
                 c('COMMON', 'EAS', 'AMR', 'AFR', 'EUR', 'SAS')))
    runStopCommand(sprintf("Population \"%s\" not available.",config$supplementaryFiles$allele_ref_std_population))
  else
    config$supplementaryFiles$allele_ref_std_population <- toupper(config$supplementaryFiles$allele_ref_std_population)



  ### alternative reference file.

  # if it is not selected by user, the process of matching with alt-ref is skipped.
  if (is_empty(config$supplementaryFiles$allele_ref_alt)) {
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


  }


  ### starandard reference file for EFFECT (BETA)
  if (is_empty(config$supplementaryFiles$beta_ref_std)) {
    config$supplementaryFiles$beta_ref_std <- NA
  } else{
    ## create full path of alt reference file
    config$supplementaryFiles$beta_ref_std <- sprintf('%s/%s',
                                                      config$paths$dir_references,
                                                      config$supplementaryFiles$beta_ref_std)

  }


  # 6
  ###read parameters
  ##set NA string set from config
  if (is_empty(config$input_parameters$na.string))
    config$input_parameters$na.string <- c("NA", "nan", "NaN", ".")
  else
    config$input_parameters$na.string <-
    checkConfigParameters(config$input_parameters$na.string,
                          'text',
                          c("NA", "nan", "NaN", "."))

  ## 7
  ##set column separator from config file
  if (is_empty(config$input_parameters$column_separator))
    config$input_parameters$column_separator <- 'auto'
  else
    config$input_parameters$column_separator <-
    checkConfigParameters(config$input_parameters$column_separator,
                          'text',
                          'auto')

  ##set effect type string  from config
  if (is_empty(config$input_parameters$effect_type))
    config$input_parameters$effect_type <- 'BETA'
  else
    config$input_parameters$effect_type <-
    checkConfigParameters(
      toupper(config$input_parameters$effect_type),
      'list',
      default = 'BETA',
      range = c("BETA", "OR")
    )





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
  if (is_empty(config$input_parameters$imputed_T))
    config$input_parameters$imputed_T <- "TRUE|T|YES|Y"
  else
    config$input_parameters$imputed_T <-
    paste(
      evaluateListsAsStrings('imputed_T', config$input_parameters$imputed_T),
      collapse = '|'
    )

  if (is_empty(config$input_parameters$imputed_F))
    config$input_parameters$imputed_F <- "FALSE|F|NO|N"
  else
    config$input_parameters$imputed_F <-
    paste(
      evaluateListsAsStrings('imputed_F', config$input_parameters$imputed_F),
      collapse = '|'
    )


  ##IMP_QUALITY
  config$filters$minimal_impQ_value  <-
    checkConfigParameters(config$filters$minimal_impQ_value, 'numeric', -0.5)
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
                          'logical',
                          FALSE)

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

  config$output_parameters$save_as_effectSize_reference <-
    checkConfigParameters(config$output_parameters$save_as_effectSize_reference,
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

  config$output_parameters$add_column_multiallelic <-
	checkConfigParameters(config$output_parameters$add_column_multiallelic,
                          'logical',
                          FALSE)

  config$output_parameters$add_column_HQ <-
	checkConfigParameters(config$output_parameters$add_column_HQ,
                          'logical',
                          FALSE)

  config$output_parameters$add_column_rsid <-
    checkConfigParameters(config$output_parameters$add_column_rsid,
                          'logical',
                          FALSE)

  config$output_parameters$add_column_hid <-
    checkConfigParameters(config$output_parameters$add_column_hid,
                          'logical',
                          FALSE)

  config$output_parameters$add_column_AF <-
    checkConfigParameters(config$output_parameters$add_column_AF,
                          'logical',
                          FALSE)

  config$output_parameters$add_column_AFmismatch <-
    checkConfigParameters(config$output_parameters$add_column_AFmismatch,
                          'logical',
                          FALSE)

   config$output_parameters$ordered <-
    checkConfigParameters(config$output_parameters$ordered,
                          'logical',
                          FALSE)

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






  return(config)
}




checkConfigFileSections <- function(config)
{
  required.sections <-  c(
    'paths',
    'supplementaryFiles',
    'input_parameters',
    'output_parameters',
    'remove_chromosomes',
    'plot_specs',
    'filters'
  )

  missing.sections <- which( required.sections %notin% names(config) )

  if (length(missing.sections) == 0)
    return(TRUE)
  else
  {
    cat("\nError in configuration file.")
    cat("\nMissing sections:", paste(required.sections[missing.sections],collapse = '  |  '))
    return(FALSE)
  }

}



evaluateListsAsStrings <- function(input.name, input.string)
{
  output.list <- tryCatch(
    eval(parse(text = input.string)),
    warning = function(err) {
      print_and_log(
        sprintf(
          'error evaluating \'%s\' parameter in config file! should be written as a vector.',
          input.name
        ),
        'fatal'
      )
      return(NULL)
    },
    error = function(err) {
      print_and_log(
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
    print_and_log(
      sprintf(
        '\'%s\' parameter should be written as a vector in config file!',
        input.name
      ),
      'fatal'
    )

}



## TODO not complete yet
printConfigVariables <- function(config) {
  print_and_log(sprintf('Input directory: \'%s\'', config$paths$dir_data),
                'info')

  print_and_log(sprintf('Output directory: \'%s\'', config$paths$dir_output),
                'info')

  print_and_log(sprintf('References directory: \'%s\'', config$paths$dir_references),
                'info')

  print_and_log(
    sprintf(
      'Alternate header file: \'%s\'',
      config$supplementaryFiles$header_translations
    ),
    'info'
  )

  print_and_log(
    sprintf(
      'Allele Frequency Reference file: \'%s\'',
      config$supplementaryFiles$allele_ref_std
    ),
    'info'
  )

  if (!is.na(config$supplementaryFiles$allele_ref_alt))
    print_and_log(
      sprintf(
        'Allele Frequency Alternative Reference file: \'%s\'',
        config$supplementaryFiles$allele_ref_alt
      ),
      'info'
    )

  if (!is.na(config$supplementaryFiles$beta_ref_std))
    print_and_log(
      sprintf(
        'Effect Size Reference file: \'%s\'',
        config$supplementaryFiles$beta_ref_std
      ),
      'info'
    )


  if (!.QC$config$graphic.device)
    print_and_log('No graphic devices are available. Plotting will be skipped!',
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
          is.na(as.logical(toupper(parameter))))
        output <- default
      else
        output  <- as.logical(toupper(parameter))
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




set_test_run_variables <- function(test.run)
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



make_config <- function(object)
{
  config <- list(
    paths = list(
      filename = object@paths$filename,
      filename_output_tag = object@paths$filename_output_tag,
      dir_data = object@paths$dir_data,
      dir_output = object@paths$dir_output,
      dir_references = object@paths$dir_references,
      input_files = object@input_files
    ),
    supplementaryFiles = list(
      header_translations = object@supplementaryFiles$header_translations,
      allele_ref_std =  object@supplementaryFiles$allele_ref_std,
      allele_ref_std_population =  object@supplementaryFiles$allele_ref_std_population,
      allele_ref_alt =  object@supplementaryFiles$allele_ref_alt,
      beta_ref_std = object@supplementaryFiles$beta_ref_std
    ),
    input_parameters = list(
      effect_type = object@input_parameters$effect_type,
      column_separator = object@input_parameters$column_separator,
      na.string = object@input_parameters$na.string,
      imputed_T = object@input_parameters$imputed_T,
      imputed_F = object@input_parameters$imputed_F,
      calculate_missing_p = object@input_parameters$calculate_missing_p
    ),
    output_parameters = list(
      save_final_dataset = object@output_parameters$save_final_dataset,
      save_as_effectSize_reference = object@output_parameters$save_as_effectSize_reference,
      gzip_final_dataset = object@output_parameters$gzip_final_dataset,
      out_header = object@output_parameters$out_header,
      out_sep = object@output_parameters$out_sep,
      out_na = object@output_parameters$out_na,
      out_dec = object@output_parameters$out_dec,
      html_report = object@output_parameters$html_report,
      object_file = object@output_parameters$object_file,
      add_column_multiallelic = object@output_parameters$add_column_multiallelic,
      add_column_HQ= object@output_parameters$add_column_HQ,
      add_column_AFmismatch = object@output_parameters$add_column_AFmismatch,
      add_column_rsid = object@output_parameters$add_column_rsid,
      add_column_hid = object@output_parameters$add_column_hid,
      add_column_AF = object@output_parameters$add_column_AF,
      ordered = object@output_parameters$ordered
    ),
    remove_chromosomes = list(
      remove_X = object@remove_chromosomes$remove_X,
      remove_Y = object@remove_chromosomes$remove_Y,
      remove_XY = object@remove_chromosomes$remove_XY,
      remove_M = object@remove_chromosomes$remove_M
    ),
    plot_specs = list(
      make_plots = object@plot_specs$make_plots,
      plot_cutoff_p = object@plot_specs$plot_cutoff_p,
      graphic_device = object@plot_specs$graphic_device,
      plot_title = object@plot_specs$plot_title
    ),
    filters = list(
      HQfilter_FRQ = object@filters$HQfilter_FRQ,
      HQfilter_HWE = object@filters$HQfilter_HWE,
      HQfilter_cal = object@filters$HQfilter_cal,
      HQfilter_imp = object@filters$HQfilter_imp,
      threshold_diffEAF = object@filters$threshold_diffEAF,
      minimal_impQ_value = object@filters$minimal_impQ_value,
      maximal_impQ_value = object@filters$maximal_impQ_value
    ),
    debug = list(
      verbose = object@debug$verbose,
      save_pre_modification_file = object@debug$save_pre_modification_file,
      reduced.AF.plot = object@debug$reduced.AF.plot,
      test_row_count = object@debug$test_row_count
    )
  )


  # this item is used for display in report and plots
  # if effect = OR , it will be converted to BEta by ln(OR)
  # in this case, ln(OR) is displayed in plots and reports instead of BETA
  config$input_parameters$effect_type_string <-
    ifelse(config$input_parameters$effect_type == 'BETA' ,
           'BETA',
           'Ln(OR)')


  if (!is.na(config$supplementaryFiles$allele_ref_alt) &&
      !is_empty(config$supplementaryFiles$allele_ref_alt) &&
      file.exists(config$supplementaryFiles$allele_ref_alt))
    config$alt_ref_file_exists <- TRUE
  else
    config$alt_ref_file_exists <- FALSE



  if (config$plot_specs$make_plots)
    config <- set_graphic_device(config)
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

  return(config)
}
