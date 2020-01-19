## this file cotains functions that are called at the very first or the e nd of the algorithm
# getDefaultSystemOptions => R options before QC algorithm starts
# changeROptions => set R options for better performance
# resetDefaultSystemOptions =? restore R options to user default values


getDefaultSystemOptions<-function()
{
  return(options())
}




removeFunctionVariablesFromRAM<-function(){

  print.and.log('cleaning workspace...','info')

  if(exists("reference.data",envir = .QC) &&  class(.QC$reference.data) == "SQLiteConnection")
    RSQLite::dbDisconnect(.QC$reference.data)


  rm(.QC)
  # rm(list = ls(envir = .QC))

  # if(!is.null(dev.list()))
  #   invisible(dev.off(dev.list()))

  dev.close.result <- tryCatch({
      invisible(graphics.off())
      invisible(dev.off())
      return(TRUE)
    },
    error = function(err){
      return(FALSE)
    }
 )


  if(dev.close.result)
    print.and.log('Graphic devices are closed.','info')


  invisible(gc()) ## Free RAM

  message('\n=============================================')
  message('============= FINISHED QC ===================')
  message('=============================================')

}

##change R options back to what it was
resetDefaultSystemOptions<-function(user.options)
{
  options(user.options)
}
