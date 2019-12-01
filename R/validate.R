validate.Inspector <- function(object, printWarnings = TRUE)
{
  if (!is(object, "Inspector"))
    stop("Object must be of class Inspector", call. = FALSE)

  ##
  errors <- character()
  warnings <- character()

  ## data directory
  if (!dir.exists(object@paths$dir_data))
    errors <- c(errors, "Data directory was not found.")




  ## output directory


  if (!dir.exists(object@paths$dir_output)) {
    tryCatch(
      dir.create(object@paths$dir_output),
      error = function(err) {
        msg <- sprintf("Could not create outpur directory at \'%s\'! %s",
                       object@paths$dir_output,
                       err)
        errors <- c(errors, msg)
      }
    )
  }

  if (file.access(object@paths$dir_output, 2) != 0)
    errors <-
    c(errors,
      sprintf("Permission denied: \'%s\'", object@paths$dir_output))


  ## reference directory
  if (!dir.exists(object@paths$dir_references))
    errors <- c(errors, "Reference directory not found.")
  else if (file.access(object@paths$dir_references, 2) != 0)
    errors <-
    c(errors,
      sprintf("Permission denied: \'%s\'", object@paths$dir_references))


  ## filenames

  if (nrow(data.frame(object@input_files)) == 0)
    errors <-
    c(errors, "Input file list is empty.")




  ## check if header file exists
  if (!file.exists(object@supplementaryFiles$header_translations))
    errors <-
    c(
      errors,
      sprintf(
        "Alternative Header file not found at \'%s\'!",
        object@supplementaryFiles$header_translations
      )
    )
  else
  {
    test_out <- tryCatch({
      headerTable <- read.table(
        file = object@supplementaryFiles$header_translations,
        sep = '',
        header = FALSE,
        stringsAsFactors = FALSE
      )

      if (ncol(headerTable) != 2L)
        errors <-
          c(errors,
            sprintf(
              '\'headerTable\' should have two columns but has %s!',
              ncol(headerTable)
            ))

      if (any(duplicated(headerTable[, 2])))
        errors <-
          c(errors, 'Duplicated items found in header table!')

    },
    error = function(err)
    {
      errors <- c(errors, sprintf('Error in reading header file: %s', err$message))
    })
  }



  ## check if reference file exists
  if (!file.exists(object@supplementaryFiles$allele_ref_std))
    errors <-
    c(
      errors,
      sprintf(
        "Reference file not found at \'%s\'!",
        object@supplementaryFiles$allele_ref_std
      )
    )


  ## check if beta reference file exists
  if (!is.na(object@supplementaryFiles$beta_ref_std) &&
      !file.exists(object@supplementaryFiles$beta_ref_std))
    errors <- c(
      errors,
      sprintf(
        "Beta (effect) reference file not found at \'%s\'!",
        object@supplementaryFiles$beta_ref_std
      )
    )


  ##       WARNINGS
  if(object@output_parameters$save_final_dataset == FALSE)
    warnings <- c(warnings, "Cleaned output file will not be saved.")

  if(object@plot_specs$make_plots == FALSE)
    warnings <- c(warnings, "Plots are skipped.")

  if(object@output_parameters$html_report == FALSE)
    warnings <- c(warnings, "HTML report will not be saved.")

  if(object@remove_chromosomes$remove_X == TRUE)
    warnings <- c(warnings, "Variants on chromosome X will be deleted.")

  if(object@remove_chromosomes$remove_XY == TRUE)
    warnings <- c(warnings, "Variants on chromosome XY will be deleted.")

  if(object@remove_chromosomes$remove_Y == TRUE)
    warnings <- c(warnings, "Variants on chromosome Y will be deleted.")

  if(object@remove_chromosomes$remove_M == TRUE)
    warnings <- c(warnings, "Variants on chromosome MT will be deleted.")



  ######################################################################

  if(object@StudyList@studyCount > 0)
    object@StudyList <- new("StudyList")  # empty the existing StudyList


  if (printWarnings && length(warnings) > 0)
  {
   cat("   Warnings")
    for (warning in warnings)
      cat("\n\t -", warning)
  }


  ## display the results
  if (length(errors) == 0)
  {
    return(TRUE)
  }
  else
  {
    cat("\n   Errors")
    for (error in errors)
      cat("\n\t -", error)
    cat("\n   Object was not validated.\n")
    return(FALSE)
  }



}
