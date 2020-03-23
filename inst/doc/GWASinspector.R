## ---- echo = FALSE-------------------------------------------------------
  knitr::opts_chunk$set(
  eval=FALSE,
  results = "hide",
  collapse = TRUE, 
  comment = "#>",
  results = "asis"
)

## ------------------------------------------------------------------------
#  # this will automatically download and install the dependencies.
#  install.packages("GWASinspector")

## ------------------------------------------------------------------------
#  # get the installation function from our website:
#  source('http://GWASinspector.com/references/install_GWASinspector.R')
#  
#  # this function will check R packages and install the dependencies from CRAN.
#  install.GWASinspector(package.path = 'path/to/packageFile.gz')
#  

## ------------------------------------------------------------------------
#  require(GWASinspector)
#  

## ------------------------------------------------------------------------
#  system.check()
#  

## ------------------------------------------------------------------------
#  get.headerTranslation('/path/to/referenceFolder') # copy the file to selected folder
#  

## ------------------------------------------------------------------------
#  get.config(getwd()) # copy the file to selected folder
#  

## ------------------------------------------------------------------------
#  
#  ## load the package
#  require(GWASinspector)
#  
#  ## import the QC-configuration file
#  job <- setup.inspector("/home/user/config.ini")
#  
#  ## check the created instance
#  ## input result files that will be inspected are also displayed
#  job
#  
#  ## run the algorithm
#  job <- run.inspector(job)
#  
#  ## check the results
#  ## comprehensive report and result file are already saved in the output folder
#  result.inspector(job)
#  

## ------------------------------------------------------------------------
#  require(GWASinspector)
#  demo.inspector('/sample_dir')

