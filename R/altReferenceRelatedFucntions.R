
uploadAltReferenceFile<-function()
{

  altReferenceFile <- .QC$config$supplementaryFiles$allele_ref_alt
  file.extension<-tolower(file_ext(altReferenceFile))

  #load the file based on its extension
  if(file.extension %in% c('csv','txt','dat')){
    allele_ref_alt_std<-fread(altReferenceFile,
                              data.table = TRUE)
  }else if(file.extension == "rdata"){
    load(altReferenceFile)

    # make sure RDATA file contains the correct variable name
    if(!exists('allele_ref_alt_std'))
      runStopCommand('Alternative Reference file should contain \'allele_ref_alt_std\' variable! use RDS or txt file instead of RDATA!')


  }else if(file.extension == 'rds')
  {
    allele_ref_alt_std <- readRDS(altReferenceFile)


  }else if(file.extension == "zip") {

    # name if txt file inside zip file
    # FIXME how to know the exact name!?!?
    embeded.file <- paste0(tools::file_path_sans_ext(basename(altReferenceFile)),'.txt')

    allele_ref_alt_std <- read.table(unz(description = altReferenceFile,
                                         filename = embeded.file),
                                     sep = "",
                                     header = TRUE,
                                     stringsAsFactors = FALSE)

    close(unz(
      description = altReferenceFile,
      filename = embeded.file)
    )


  } else if (file.extension == "gz"){
    allele_ref_alt_std <- read.table(gzfile(altReferenceFile),
                                     sep = "",
                                     header = TRUE,
                                     stringsAsFactors = FALSE)

    close(gzfile(altReferenceFile))
  }


  # convert to data.table and set KEY
  if(!is.data.table(allele_ref_alt_std))
    allele_ref_alt_std<-data.table::setDT(allele_ref_alt_std, key = "hID")


  # FIXME only hID must be used.
  # set key for fast access
  # if('hID' %in% names(allele_ref_alt_std))
  #   setkey(allele_ref_alt_std,"hID")
  # else
  #   setkey(allele_ref_alt_std,"ID")



  print.and.log(sprintf("Alternative Reference file \'%s\' loaded (%s x %s)!",
                        altReferenceFile,thousand.sep(nrow(allele_ref_alt_std)),ncol(allele_ref_alt_std)),
                'info')


  # if('DATE_ADDED' %notin% ref.col.names)
  #   .QC$reference.data[,DATE_ADDED := NA ]
  #
  # if('SOURCE' %notin% ref.col.names)
  #   .QC$reference.data[,SOURCE := 'Alt_ref']

  return(allele_ref_alt_std)
}



##check if reference file has all the crucial columns => stop if column is missing
checkAltReferenceFileIntegrity <- function() {
  ref.col.names <- colnames(.QC$alt.reference.data)

  ##TODO set this is a setting file
  required.ref.col.names <- c("hID","ID", "REF" ,"ALT", "AF" ,'DATE_ADDED','SOURCE')

  missing.ref.col.index <- which(required.ref.col.names %notin% ref.col.names)

  if(length(missing.ref.col.index) > 0)
  {
    print.and.log(sprintf('Missing crucial column in alternative reference file : \'%s\' !',
                          paste(required.ref.col.names[missing.ref.col.index],collapse = '|')),
                  'warning',display=.QC$config$debug$verbose)

    .QC$alt.reference.data <- data.table()
  }
  else
  {
    #check for duplicated hIDs
    checkDuplicatedHID_in_Alt_Ref()

    print.and.log('Alternative Reference file validated!',
                  'info')
  }
}


checkDuplicatedHID_in_Alt_Ref <- function()
{
  dups <- which(duplicated(.QC$alt.reference.data$hID)) # the firs item is kept, other duplicates are removed
  if(length(dups) > 0)
  {
    .QC$alt.reference.data <- .QC$alt.reference.data[!dups,]
    print.and.log(sprintf('%s duplicated items found in alternate reference and are removed!',
                          length(dups)),
                  'warning',display=.QC$config$debug$verbose)

  }
}



update.alternate.reference <- function(input.data) {

  if(any(c('MARKER','CHR','POSITION') %notin% colnames(input.data)))
  {
    print.and.log('Alternate allele frequency data set was not updated due to a missing column in dataset.','warning',display=.QC$config$debug$verbose)
    return(NULL)
  }



  # find variants that were not dounf in either references and have a valid allele frequency
  unknown.variants <- subset(input.data[is.na(REF) & !is.na(EFF_ALL_FREQ)],
                               select=c('MARKER','hID','EFFECT_ALL','OTHER_ALL','EFF_ALL_FREQ'))


  if(nrow(unknown.variants) > 0 )
  {
    # rename columns according to refrence file standard
    names(unknown.variants)[names(unknown.variants) == 'MARKER'] <- 'ID'
    names(unknown.variants)[names(unknown.variants) == 'EFFECT_ALL'] <- 'ALT'
    names(unknown.variants)[names(unknown.variants) == 'OTHER_ALL'] <- 'REF'
    names(unknown.variants)[names(unknown.variants) == 'EFF_ALL_FREQ'] <- 'AF'


    # TODO check if required for all columns
    # if(!is.numeric(unknown.variants$SNP))
    #   unknown.variants$SNP <- as.numeric(unknown.variants$SNP)


    # add data and source column
    unknown.variants[,DATE_ADDED := as.character(Sys.Date())]
    unknown.variants[,SOURCE := .QC$thisStudy$file.name]


    # bind new variants with previous alt reference file and save the data
    if(nrow(.QC$alt.reference.data ) > 0)
      .QC$alt.reference.data <- rbind(.QC$alt.reference.data , unknown.variants)
    else
      .QC$alt.reference.data <- unknown.variants


    # ==

    print.and.log(sprintf('Alternative Reference file is updated with %s rows', thousand.sep(nrow(unknown.variants))),
                  'info')
  }
}



save.alternate.reference <- function()
{
  altReferenceFile <- .QC$config$supplementaryFiles$allele_ref_alt
  file.extension <- tolower(file_ext(altReferenceFile))

  message('\n---------- [saving alternate reference file] ----------')

  #load the file based on its extension
  if(file.extension %in% c('csv','txt','dat')){
    saveDataSet(dataset = .QC$alt.reference.data,
                file.path = .QC$config$supplementaryFiles$allele_ref_alt,
                zipped = FALSE,
                order = FALSE)
  }else if(file.extension == "rdata"){
    # save function can not save part of environment
    # and should be presented as a new object
    allele_ref_alt_std <- .QC$alt.reference.data
    save(allele_ref_alt_std, file = .QC$config$supplementaryFiles$allele_ref_alt)
    rm(allele_ref_alt_std)
  }else if(file.extension == 'rds')
  {
    saveRDS(.QC$alt.reference.data, file = .QC$config$supplementaryFiles$allele_ref_alt, version = '2')


  }else if(file.extension %in% c('gz','zip')) {
    saveDataSet(dataset = .QC$alt.reference.data,
                file.path = .QC$config$supplementaryFiles$allele_ref_alt,
                zipped = TRUE,
                order = FALSE)
  }


  print.and.log('Alternate reference file is saved!','info')
}
