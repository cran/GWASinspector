
uploadBetaReferenceFile<-function()
{

  betaReferenceFile <- .QC$config$supplementaryFiles$beta_ref_std
  file.extension<-tolower(file_ext(betaReferenceFile))

  #load the file based on its extension
  if(file.extension %in% c('csv','txt','dat')){
    beta_ref_std<-fread(betaReferenceFile,
                        data.table = TRUE)
  }else if(file.extension == "rdata"){
    load(betaReferenceFile)

    # make sure RDATA file contains the correct variable name
    if(!exists('beta_ref_std'))
      print.and.log('Beta (Effect) reference file should contain \'beta_ref_std\' variable! use RDS or txt file instead of RDATA!','warning')

  }else if(file.extension == 'rds')
  {
    beta_ref_std<-readRDS(betaReferenceFile)
  }else if (file.extension == "gz"){
    beta_ref_std <- read.table(gzfile(betaReferenceFile),
                               sep = "",
                               header = TRUE,
                               stringsAsFactors = FALSE)

    close(gzfile(betaReferenceFile))
  }

  # convert to data.table
  if(!is.data.table(beta_ref_std))
    beta_ref_std<-as.data.table(beta_ref_std)


  print.and.log(sprintf("Effect-size Reference file \'%s\' loaded (%s x %s)!",
                        betaReferenceFile,nrow(beta_ref_std),ncol(beta_ref_std)),
                'info')

  return(beta_ref_std)
}



##check if Beta reference file has all the crucial columns => stop if column is missing
# only ID and EFFECT columns are required.
checkBetaReferenceFileIntegrity <- function() {

  if(is.null(.QC$reference.data.effect) | nrow(.QC$reference.data.effect) == 0 )
  {
    print.and.log('Effect-Size referene file is empty!','warning')
    return(data.table())
  }


  #=============================

  ref.col.names <- colnames(.QC$reference.data.effect)

  ##FIXME which columns are required
  required.ref.col.names <- c("hID", "EFF_ALL" ,"NON_EFF_ALL", "EFFECT" )


  missing.ref.col.index <- which(required.ref.col.names %notin% ref.col.names)


  if(length(missing.ref.col.index) > 0){
    print.and.log(sprintf('Missing crucial column in Effect-Size reference file : \'%s\' !',
                          paste(required.ref.col.names[missing.ref.col.index],collapse = '|')),
                  'warning')
    print.and.log('Effect-size comparison is skipped!', 'warning')
    .QC$reference.data.effect <- data.table()
  }
  else
  {

    # changing column names for consistency with AF dataset for matching function
    names(.QC$reference.data.effect)[names(.QC$reference.data.effect) == 'EFF_ALL'] <- 'ALT'
    names(.QC$reference.data.effect)[names(.QC$reference.data.effect) == 'NON_EFF_ALL'] <- 'REF'

    print.and.log('Effect-Size Reference file validated!', 'info')
  }
}


