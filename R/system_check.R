#' Checks which required and optional packages are available
#'
#' This functions checks the availability of all required and optional packages for GWASinspector.
#' The optional packages are not mandatory for running the algorithm; but will add useful functionalities.
#' @return System information and required functionalities for the QC algorythm are checked and reported as a data frame.
#' @examples system.check()
#'
system.check <- function()
{
  # message('\nThis table displays system information and checks if the required packages are available for QC process.\nCheck the package manual for detailed information.')

  check.tools()


  tbl <- data.frame(t(data.frame(
    'Platform' = .QC$OS,
    'System' = .QC$r.version,
    'QC script version' = .QC$script.version,
    'Package description' = '  ',
    "pandoc" = ifelse(.QC$pandoc.exists , paste(rmarkdown::pandoc_version(),collapse = '.') , 'not available'),
    "kableExtra" = ifelse(.QC$kableExtra.package.exists,   paste(packageVersion('kableExtra'),collapse = '.') , 'not available'),
    "parallel" = ifelse(.QC$parallel.package.exists,  paste(packageVersion('parallel'),collapse = '.') , 'not available'),
    "CPU cores" = ifelse(.QC$parallel.package.exists,  parallel::detectCores() , 'not available'),
    "ggplot2" = .QC$ggplot2.version,
    "xlsx" = ifelse(.QC$xlsx.package.exists,  paste(packageVersion('xlsx'),collapse = '.') , 'not available'),
    "rJava" = ifelse(.QC$rJava.package.exists,  paste(packageVersion('rJava'),collapse = '.') , 'not available'),
    "RSQLite" = ifelse(.QC$rsqlite.package.exists,  paste(packageVersion('RSQLite'),collapse = '.') , 'not available'),
    "Capabilities" = '  '
  )))

  colnames(tbl) <- 'version/availability'
  row.names(tbl) <- sapply(row.names(tbl), function(x) return(gsub('\\.',' ',x)))
  row.names(tbl)[5:12] <- sapply(row.names(tbl)[5:12] , function(x) return(paste('~',x)))

  # ==
  cap <- data.frame(capabilities(c('png','jpeg','tiff','cairo')))
  colnames(cap) <- 'version/availability'
  cap$`version/availability` <- as.character(cap$`version/availability`)
  row.names(cap) <- sapply(row.names(cap) , function(x) return(paste('~',x)))


  # function modified according to CRAN comments
  # returns the table instead of printing it
  return(rbind(tbl,cap))

  # # ==
  # print(knitr::kable(rbind(tbl,cap),format='rst'))
  #
  # # ==
  # if(.QC$java.exists){
  #   jav <- data.frame(system('java -version',intern = TRUE))
  #   print(knitr::kable(jav,format='rst',col.names = 'java'))
  # }
  # else
  #   message('\nWarning: java not found on system!\n')


}
