example.config <- function(result.dir)
{
  message('This is a example run for GWASinspector algorithm which is based on a sample fabricated GWAS result file.')
  message('New folder will be created for saving sample output files in your active directory.')
  message('You can check the log file for further information.')

  Sys.sleep(4)


  # remove the trailing backslash from path
  result.dir <- gsub('/+$', '' , result.dir)

  new.dir <- paste(result.dir,'example_output',sep = '/')


  package.path <- find.package('GWASinspector',quiet = FALSE)
  extdata.path <- paste(package.path,'extdata',sep='/')

  configFile <- system.file("extdata", "config.ini", package = "GWASinspector")
  config <- read.ini(filepath = configFile)

  config$paths$dir_data <- extdata.path
  config$paths$dir_references <- extdata.path
  config$paths$dir_output <- new.dir

  config$paths$filename <- 'sample.txt.gz'
  config$supplementaryFiles$allele_ref_std <- 'HapMap_CEU_r28_b36_EDIT_v10c_sample.rds'
  config$supplementaryFiles$allele_ref_std_population <- 'common'

  config$output_parameters$object_file <- FALSE

  config$plot_specs$plot_cutoff_p <- 0.2
  config$plot_specs$plot_title <- 'Example Run'

  return(config)
}
