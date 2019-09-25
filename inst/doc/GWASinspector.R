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

## ------------------------------------------------------------------------
#  system.check()

## ------------------------------------------------------------------------
#  get.headerTranslation.file('c:/path/to/referenceFolder') # copies the file to selected folder
#  

## ------------------------------------------------------------------------
#  setwd('c:/path/to/workingDirectory') # copies the file to selected folder
#  

## ------------------------------------------------------------------------
#  get.config(getwd()) # copies the file to the working directory
#  

## ------------------------------------------------------------------------
#  inspect('config.ini')

## ------------------------------------------------------------------------
#  library(GWASinspector)
#  inspect.example('/sample_dir')

