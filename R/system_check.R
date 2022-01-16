#' Checks which required and optional packages are available
#'
#' Checks if required and optional packages are installed on the system.
#' Although the optional packages do not contribute to the QC itself, having them available
#' will allow for Excel and HTML formatted log files, which are easier to read and interpret.
#'
#' @return System information and required functionalities for the QC algorithm are checked and reported as a data frame.
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
    "ggplot2" = .QC$ggplot2.version,
    "xlsx" = ifelse(.QC$xlsx.package.exists,  paste(packageVersion('xlsx'),collapse = '.') , 'not available'),
    "rJava" = ifelse(.QC$rJava.package.exists,  paste(packageVersion('rJava'),collapse = '.') , 'not available'),
    "RSQLite" = ifelse(.QC$rsqlite.package.exists,  paste(packageVersion('RSQLite'),collapse = '.') , 'not available'),
    "Capabilities" = '  '
  )))

  colnames(tbl) <- 'version/availability'
  row.names(tbl) <- sapply(row.names(tbl), function(x) return(gsub('\\.',' ',x)))
  row.names(tbl)[5:10] <- sapply(row.names(tbl)[5:10] , function(x) return(paste('~',x)))

  # ==
  cap <- data.frame(capabilities(c('png','jpeg','tiff','cairo')))
  colnames(cap) <- 'version/availability'
  cap$`version/availability` <- as.character(cap$`version/availability`)
  row.names(cap) <- sapply(row.names(cap) , function(x) return(paste('~',x)))


  # function modified according to CRAN comments
  # returns the table instead of printing it
  return(rbind(tbl,cap))


}
