#' Importing a QC-configuration file into R
#'
#' To run a QC in GWASinspector, first generate a config file using \code{\link{get.config}}, and edit it to suit your requirements.
#' Next, use the function \code{\link{setup.inspector}} to check the configuration file and import it into R.
#' This will create an object of the inspector class, which can then be processed using \code{\link{run.inspector}}.
#'
#' @param config.file character. Path to a configuration (.ini) file. For a sample configuration file, see \code{\link{get.config}}.
#' @param validate logical. Whether to validate the object.
#' @return returns a new instance of \linkS4class{Inspector} class.
#' @examples
#' config.file <- get.config(tempdir())
#' job <- setup.inspector(config.file , validate = FALSE)
#' job
#'
setup.inspector <- function(config.file , validate = TRUE)
{


  if (missing(config.file))
    runStopCommand('Configuration file not set.')
  else
    configuration <- checkConfigFile(config.file)


  object <- new(
    "Inspector",
    paths = list(
      filename = configuration$paths$filename,
      filename_output_tag = configuration$paths$filename_output_tag,
      dir_data = configuration$paths$dir_data,
      dir_output = configuration$paths$dir_output,
      dir_references = configuration$paths$dir_references
    ),
    supplementaryFiles = list(
      header_translations = configuration$supplementaryFiles$header_translations,
      allele_ref_std =  configuration$supplementaryFiles$allele_ref_std,
      allele_ref_std_population =  configuration$supplementaryFiles$allele_ref_std_population,
      allele_ref_alt =  configuration$supplementaryFiles$allele_ref_alt,
      beta_ref_std = configuration$supplementaryFiles$beta_ref_std
    ),
    input_parameters = list(
      effect_type = configuration$input_parameters$effect_type,
      column_separator = configuration$input_parameters$column_separator,
      na.string = configuration$input_parameters$na.string,
      imputed_T = paste(configuration$input_parameters$imputed_T,collapse = '|'),
      imputed_F = paste(configuration$input_parameters$imputed_F,collapse = '|'),
      calculate_missing_p = configuration$input_parameters$calculate_missing_p
    ),
    output_parameters = list(
      save_final_dataset = configuration$output_parameters$save_final_dataset,
      gzip_final_dataset = configuration$output_parameters$gzip_final_dataset,
      out_header = configuration$output_parameters$out_header,
      out_sep = configuration$output_parameters$out_sep,
      out_na = configuration$output_parameters$out_na,
      out_dec = configuration$output_parameters$out_dec,
      html_report = configuration$output_parameters$html_report,
      object_file =  configuration$output_parameters$object_file
    ),
    remove_chromosomes = list(
      remove_X = configuration$remove_chromosomes$remove_X,
      remove_Y = configuration$remove_chromosomes$remove_Y,
      remove_XY = configuration$remove_chromosomes$remove_XY,
      remove_M = configuration$remove_chromosomes$remove_M
    ),
    plot_specs = list(
      make_plots = configuration$plot_specs$make_plots,
      plot_cutoff_p = configuration$plot_specs$plot_cutoff_p,
      graphic_device = configuration$plot_specs$graphic_device,
      plot_title = configuration$plot_specs$plot_title
    ),
    filters = list(
      HQfilter_FRQ = configuration$filters$HQfilter_FRQ,
      HQfilter_HWE = configuration$filters$HQfilter_HWE,
      HQfilter_cal = configuration$filters$HQfilter_cal,
      HQfilter_imp = configuration$filters$HQfilter_imp,
      threshold_diffEAF = configuration$filters$threshold_diffEAF,
      minimal_impQ_value = configuration$filters$minimal_impQ_value,
      maximal_impQ_value = configuration$filters$maximal_impQ_value
    ),
    debug = list(
      verbose = configuration$debug$verbose,
      save_pre_modification_file = configuration$debug$save_pre_modification_file,
      reduced.AF.plot = configuration$debug$reduced.AF.plot,
      test_row_count = configuration$debug$test_row_count
    ),
    input_files = configuration$paths$input_files,
    created_at = Sys.time(),
    start_time = Sys.time(),
    end_time = Sys.time()
  )

  if(!validate) ## return the object without validating
    return(object)
  else if(validate.Inspector(object))
    return(object)

}
