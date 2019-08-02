## functions related to Reference file
## uploadReferenceFile => uploads a csv or RDS or RDATA file (RDS is faster and more efficient)


uploadReferenceFile<-function()
{
  config <- .QC$config
  allele_ref_std<-NULL
  referenceFile<-config$supplementaryFiles$allele_ref_std
  file.extension<-tolower(file_ext(referenceFile))


  # if it is a database
  if(file.extension == "sqlite")
  {
    # connect to database
    allele_ref_stdDB <- getRSQLiteDatabase(referenceFile)

    # check if database has tables - break the process if table count == 0
    allele_ref_stdDB_table <- getRSQLiteDatabase.tableCount(allele_ref_stdDB)

    # check if selected subpopulation exists in DB
    getRSQLiteDatabase.SubPopulationExists(allele_ref_stdDB)

    # check if table has INDEL variants or not
   # .QC$reference.data.has.INDEL <- getRSQLiteDatabase.hasINDEL(allele_ref_stdDB)

    print.and.log(sprintf("Reference database found at \'%s\' : including %s tables!",
                          referenceFile,allele_ref_stdDB_table),
                  'info')

    return(allele_ref_stdDB)
  }


  # if it is not a database

  if(file.extension %in% c('csv','txt','dat')){
    allele_ref_std<-fread(referenceFile,
                          data.table = TRUE)
  }else if(file.extension == "rdata"){
    load(referenceFile)

    # make sure RDATA file contains the correct variable name
    if(!exists('allele_ref_std'))
      runStopCommand('Reference file should contain \'allele_ref_std\' variable! use RDS or txt file instead of RDATA!')


  }else if(file.extension == 'rds')
  {
    allele_ref_std<-readRDS(referenceFile)
  }

  # convert to data.table
  if(!is.data.table(allele_ref_std))
    allele_ref_std<-as.data.table(allele_ref_std)

  print.and.log(sprintf("Reference file \'%s\' loaded (%s x %s)!",
                        referenceFile,thousand.sep(nrow(allele_ref_std)),ncol(allele_ref_std)),
                'info')

  return(allele_ref_std)
}



##check if reference file has all the crucial columns => stop if column is missing
checkReferenceFileIntegrity <- function() {

  # 1
  # return if it is a database
  #TODO check table columns
  if(!is.data.table(.QC$reference.data))
    return(NULL)



  # 2
  # check columns if reference is a data table
  ref.col.names <- colnames(.QC$reference.data)

  ##TODO set this is a setting file
  #required.ref.col.names <- c("SNP","CHR", "POS" ,"ALT", "REF", "AF")
  required.ref.col.names <- c("hID", "REF" ,"ALT", "AF" )
  missing.ref.col.index <- which(required.ref.col.names %notin% ref.col.names)

  if(length(missing.ref.col.index) > 0)
    print.and.log(sprintf('Missing crucial column in reference file : \'%s\' !',
                          paste(required.ref.col.names[missing.ref.col.index],collapse = '|')),
                  'fatal')
  else
  {

    ## add two columns for consistency with alternate reference file
    # source column will be used for reporting how many items were found in std-ref or alt-ref
    # data_added is not required for this file. if it is missing it will be set to NA
    if('DATE_ADDED' %notin% ref.col.names)
      .QC$reference.data[,DATE_ADDED := NA ]


    ## if reference data contains this columns, convert to this value for syjchronicity
    ## if reference data does not contains this columns, create and fill with standard value
    ## IMPORTANT: this column is used for creating ALLELE FREQUENCY plot AF PLOT
    #if('SOURCE' %notin% ref.col.names)
    if('SOURCE' %notin% ref.col.names)
      .QC$reference.data[,SOURCE := 'Std_ref']


    # check if reference data has INDEL values
    #if(any(endsWith(.QC$reference.data$hID,'2')))
    # if(all(is.na(.QC$reference.data$TSA)))
    #   .QC$reference.data.has.INDEL <- FALSE
    # else
    #   .QC$reference.data.has.INDEL <- TRUE


    print.and.log('Reference file validated!',
                  'info')
  }
}
