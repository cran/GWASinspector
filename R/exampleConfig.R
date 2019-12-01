example.config <- function(result.dir, config)
{
  message('This is a example run for GWASinspector algorithm which is based on a sample fabricated GWAS result file.')
  message('New folder will be created for saving sample output files in your active directory.')
  message('You can check the log file for further information.')

  # Sys.sleep(3)


  # remove the trailing backslash from path
  result.dir <- gsub('/+$', '' , result.dir)

  new.dir <- paste(result.dir,'example_output',sep = '/')
  create.new.dir <- NULL
  create.new.dir <- tryCatch({
    if (!dir.exists(new.dir))
      dir.create(new.dir)
  }, error = function(err)
  {
    print.and.log(paste('Could not create output directory', err$message),
                  'fatal')
  })


  package.path <- find.package('GWASinspector',quiet = FALSE)
  extdata.path <- paste(package.path,'extdata',sep='/')


  config$paths$dir_data <- extdata.path
  config$paths$dir_references <- extdata.path
  config$paths$dir_output <- new.dir

  config$paths$filename <- file.path( config$paths$dir_data,'sample.txt.gz')


  config$supplementaryFiles$header_translations <- file.path(config$paths$dir_references ,'alt_headers.txt')
  config$supplementaryFiles$allele_ref_std <- file.path(config$paths$dir_references ,'HapMap_CEU_r28_b36_EDIT_v10c_sample.rds')
  config$supplementaryFiles$allele_ref_std_population <- 'common'

  config$output_parameters$object_file <- FALSE

  config$plot_specs$plot_cutoff_p <- 0.2
  config$plot_specs$plot_title <- 'Example Run'


  config$input_parameters$effect_type_string <- "BETA"
  config$graphic.device <- TRUE
  .QC$graphic.device  <- 'png'
  .QC$img.extension <- '.png'

  config$paths$xlsx.report <- sprintf(
    '%s/%s_%s',
    config$paths$dir_output,
    config$paths$filename_output_tag,
    'Report.xlsx'
  )

  config$paths$html.report <- sprintf(
    '%s/%s_%s',
    config$paths$dir_output,
    config$paths$filename_output_tag,
    'Report.html'
  )

  if(file.exists(config$paths$filename) &&
     file.exists(config$supplementaryFiles$header_translations) &&
     file.exists( config$supplementaryFiles$allele_ref_std))
    return(config)
  else
    stop("Sample files were not found. Re-install the package.")
}
