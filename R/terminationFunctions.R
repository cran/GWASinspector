## this file cotains functions that are called at the very first or the e nd of the algorithm
# getDefaultSystemOptions => R options before QC algorithm starts
# changeROptions => set R options for better performance
# resetDefaultSystemOptions =? restore R options to user default values


getDefaultSystemOptions<-function()
{
  return(options())
}




removeFunctionVariablesFromRAM<-function(){

  print_and_log('cleaning workspace...','info')

  if(exists("reference.data",envir = .QC))
    if(is(.QC$reference.data, "SQLiteConnection"))
      RSQLite::dbDisconnect(.QC$reference.data)


  rm(.QC)

 #  dev.close.result <- tryCatch({
 #      invisible(graphics.off())
 #      while (!is.null(grDevices::dev.list()))  grDevices::dev.off()
 #      return(TRUE)
 #    },
 #    error = function(err){
 #      return(FALSE)
 #    }
 # )
 #
 #
 #  if(dev.close.result)
 #    print_and_log('Graphic devices are closed.','info')


  invisible(gc()) ## Free RAM

  if(.QC$verbose)
  {
    cat('\n=============================================', fill = TRUE)
    cat('============= FINISHED QC ===================', fill = TRUE)
    cat('=============================================', fill = TRUE)
  }

}

##change R options back to what it was
resetDefaultSystemOptions<-function(user.options)
{
  options(user.options)
}


# This is an issue with ggsave() from ggplot2 package
# https://github.com/tidyverse/ggplot2/issues/2787
removeRedundantPlotFile <- function()
{
  redundantFilePath = paste(.QC$config$paths$dir_output,"Rplots.pdf",sep='/')
  try(
    {
      if(file.exists(redundantFilePath))
        file.remove(redundantFilePath)
    },silent = TRUE)
}
