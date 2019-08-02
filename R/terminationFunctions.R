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

  # FIXME delete after debugging
  # rm(config, envir = .QC) ## removing config variable from memory
  #
  # rm(log.name.main, envir = .QC) ## removing logger variable from memory
  #
  # rm(qc.study.list, envir = .QC) ## removing items variable from memory
  #
  # rm(thisStudy, envir = .QC) ## removing items variable from memory
  #
  # rm(reference.data, envir = .QC) ## removing items variable from memory
  #
  # rm(alt.reference.data, envir = .QC) ## removing items variable from memory

  # FIXME delete after debugging
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
