# use package specific environment instead of global environment
# all variables should be put here
# this environment is removed on exit
.QC <- new.env()


.QC$package.name <- 'GWASinspector'

.QC$package.description <- 'Comprehensive, efficient and easy to use quality control of genome-wide association study results'

.QC$script.version <- '1.4.6'

.QC$url <- "GWASinspector.com"

.QC$help <- 'Check out our website for more help and support'

.QC$vignetteLink <- 'Manual is available from vignette("GWASinspector")'
