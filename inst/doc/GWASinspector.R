## ----echo = FALSE-------------------------------------------------------------
  knitr::opts_chunk$set(
  eval=FALSE,
  results = "hide",
  collapse = TRUE, 
  comment = "#>",
  results = "asis"
)

## -----------------------------------------------------------------------------
#  # this will automatically download and install the dependencies.
#  install.packages("GWASinspector")

## -----------------------------------------------------------------------------
#  # get the installation function from our website:
#  source('http://GWASinspector.com/references/install_GWASinspector.R')
#  
#  # this function will check R packages and install the dependencies from CRAN.
#  install.GWASinspector(package.path = 'path/to/packageFile.gz')
#  

## -----------------------------------------------------------------------------
#  library(GWASinspector)
#  

## -----------------------------------------------------------------------------
#  system_check()
#  

## -----------------------------------------------------------------------------
#  get_headerTranslation('/path/to/referenceFolder') # copy the file to selected folder
#  

## -----------------------------------------------------------------------------
#  get_config('/home/user') # copy the file to selected folder
#  

## -----------------------------------------------------------------------------
#  
#  ## load the package
#  library(GWASinspector)
#  
#  ## import the QC-configuration file
#  job <- setup_inspector("/home/user/config.ini")
#  
#  ## check the created instance
#  ## input result files that will be inspected are also displayed
#  job
#  
#  ## run the algorithm
#  job <- run_inspector(job)
#  
#  ## check the results
#  ## comprehensive report and result file are already saved in the output folder
#  result_inspector(job)
#  

## -----------------------------------------------------------------------------
#  library(GWASinspector)
#  demo_inspector('/sample_dir')

