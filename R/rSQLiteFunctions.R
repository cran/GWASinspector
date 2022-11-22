getRSQLiteDatabase <- function(database.path){
  if(file.exists(database.path)){
    DB<-RSQLite::dbConnect(RSQLite::SQLite(),database.path)
    return(DB)
  }
  else{
    print.and.log(paste('database file not found at:',database.path),'fatal')
    return(NULL)
  }

}

getRSQLiteDatabase.tableCount <- function(database){

  tblCount <- length(RSQLite::dbListTables(database))

  if(tblCount > 0)
    return(tblCount)
  else
    print.and.log('Database is empty!','fatal')

}

getRSQLiteDatabase.hasINDEL <- function(database){
  # get the first table of database
  tbl <- RSQLite::dbListTables(database)[1]
  # check if it has INDELS
  INDEL.count <- RSQLite::dbGetQuery(database,sprintf('select TSA from %s where TSA in ("insertion","deletion") LIMIT 5',tbl))

  if(nrow(INDEL.count) == 5)
    return(TRUE)
  else
    return(FALSE)

}

getRSQLiteDatabase.SubPopulationExists <- function(database){

  population.Column = switch(.QC$config$supplementaryFiles$allele_ref_std_population,
                             'AMR'= 'AMR_AF',
                             'EUR'= 'EUR_AF',
                             'SAS'= 'SAS_AF',
                             'EAS' = 'EAS_AF',
                             'AFR'='AFR_AF',
                             'COMMON' = 'AF')


  # get table fields from the first table of database
  tblFields <- RSQLite::dbListFields(database, RSQLite::dbListTables(database)[1])


  if(!is.element(population.Column , tblFields))
  {

    # display available population data
    available.population <- c(character(0))

    if(is.element('AF' , tblFields))
      available.population <- cbind(available.population,'COMMON')

    if(is.element('EUR_AF' , tblFields))
      available.population <- cbind(available.population,'EUR')

    if(is.element('AMR_AF' , tblFields))
      available.population <- cbind(available.population,'AMR')

    if(is.element('SAS_AF' , tblFields))
      available.population <- cbind(available.population,'SAS')

    if(is.element('EAS_AF' , tblFields))
      available.population <- cbind(available.population,'EAS')

    if(is.element('AMR_AF' , tblFields))
      available.population <- cbind(available.population,'AMR')


    available.population <- paste(available.population,collapse = ' | ')
    available.population <- paste('Available population data:' , available.population)

    print.and.log(sprintf('Allele frequency data for \'%s\' population is not found in database! %s',
                          .QC$config$supplementaryFiles$allele_ref_std_population,
                          available.population),
                  'fatal')


    # print.and.log('COMMON super population code is selected by default.','warning')
    #
    # .QC$config$supplementaryFiles$allele_ref_std_population <- 'AF'
  }
  else
  {
    print.and.log(sprintf('Allele frequency data for \'%s\' population will be used.',
                          .QC$config$supplementaryFiles$allele_ref_std_population),
                  'info')

    .QC$config$supplementaryFiles$allele_ref_std_population <- population.Column

  }

}

getRSQLiteDatabase.tableNames <- function(database)
{

  tableNames <- RSQLite::dbListTables(database)
  return(as.data.table(tableNames))
}

compareInputfileWithReferenceDataBase <- function(input.data)
{

  if(!is.element('CHR',names(input.data)))
  {
    print.and.log('CHR column is required for matching with reference database!','warning')
    print.and.log('input file is ignored!','warning')
    return(NULL)
  }




  rn <- unlist(input.data$hID)

  if(.QC$config$output_parameters$add_column_rsid &&  ## get REF_RSID if user has selected in config file
     is.element("REF_RSID",RSQLite::dbListFields(.QC$reference.data,'variants')))
    {
      rs <- RSQLite::dbGetQuery(.QC$reference.data,
                              sprintf('SELECT hID,REF_RSID,REF,ALT , %s as AF FROM variants WHERE "hID" = :x' ,
                                      .QC$config$supplementaryFiles$allele_ref_std_population) ,
                              param = list(x = rn))
      } else {
        rs <- RSQLite::dbGetQuery(.QC$reference.data,
                              sprintf('SELECT hID,REF,ALT , %s as AF FROM variants WHERE "hID" = :x' ,
                                      .QC$config$supplementaryFiles$allele_ref_std_population) ,
                              param = list(x = rn))
    }


  ## merging data
  if(is.null(data.table::key(input.data)) || data.table::key(input.data) != 'hID')
    data.table::setkey(input.data,hID)

  rs <- data.table::setDT(rs , key = 'hID')

  ##input.data <- merge(x = input.data, y = rs, by.x = 'hID', by.y = 'hID',all.x=TRUE)
  # removed MERGE with data table join
  input.data <- rs[input.data]



  # multi allelic variants will be matched more than once and create duplicated lines.
  # only the first one is requiored and the rest should be removed.
  dup.allele <- which(duplicated(input.data,by=c('hID','EFFECT_ALL','OTHER_ALL','REF','ALT')))

  if(length(dup.allele) > 0)
  {
    #print.and.log(paste('found duplicates in reference matching :',length(dup.allele)),'warning')
    input.data <- input.data[!dup.allele]
  }


  .QC$stored.reference.data <- rs
  .QC$searched.hIDs <- rn

  rm(rs)
  rm(rn)




  # a row may match to more than one result due to multiple values in refrence
  # first check if there are any duplicated rows that AF of that population is zero. keep the one with  AF>0
  # FIXME second, check if there are still duplicated hIDs. delete all
  # 1
  # e.g  A|GT  0
  #      A|GTGT 0.1
  # duplicate.hID.zeroAF <- input.data[which(duplicated(input.data$hID) |
  #                                           duplicated(input.data$hID , fromLast = TRUE))][AF == 0 ,]$hID

  # ignore variable of the duplicated variants on the same position is set to NA becuase one of them is removed in the next step
  #input.data[duplicate.hID.zeroAF, ignore := as.numeric(NA)]

  # input.data <- input.data[!which(is.element(input.data$hID, duplicate.hID.zeroAF) & AF == 0), ]


  # remove varaints on the same chr:pos if it can not be found out
  # these are
  # 1- INS or DEL on the same position with different alleles
  # 2- different INDEL types on the same position (INS/DEL)
  input.data <- save.remove.ambiguous.variants(input.data)



  # find multi-allelic variants
  input.data[, MULTI_ALLELIC := ifelse(is.na(AF),
                                       NA,
                                       ifelse(grepl(',',AF),
                                              1,
                                              0))]



  ## try allele matching on multi-allelic variants
  # if(is.element('Yes',.QC$thisStudy$tables$multi_allele_count_preProcess$`Multi-allelic`))
  if(any(input.data$MULTI_ALLELIC == 1, na.rm = TRUE))
  {
    input.data[MULTI_ALLELIC == 1,
               c('ALT','AF') := clean.multi_alleles(EFFECT_ALL , OTHER_ALL, REF, ALT, AF) ,
               by = list(EFFECT_ALL , OTHER_ALL,REF, ALT,AF)]

    # some multi-allele INDEL AFs are all 0 and will be returned the same way due to missing alleles
    # e.g. AAC,AA,TT  0,0,0   ==> this AF can be converted to 0
    #input.data[VT == 2 & MULTI_ALLELIC == 1 &  grepl(',', AF) & all(strsplit(AF,',')[[1]] == "0") , AF := "0" ]
    input.data[VT == 2 & MULTI_ALLELIC == 1 &  grepl(',', AF) & !grepl('[1-9]',AF) , AF := "0" ]

    ## get frequency table for multi-allelic variants
    #.QC$thisStudy$tables$multi_allele_count_postProcess <- getMultiAlleleCountTbl(input.data,'AF')

  }




  # AF column may be character type due to remaining ',' => convert to numeric
  # AF of multi-allelics that could notbe matched are set to NA
  if(!is.numeric(input.data$AF))
    input.data[, AF := as.numeric(AF)]


  # FIXME do not convert ALT to NA because it is used to count unmatched multiallelic variants
  #input.data[is.na(AF) , `:=` (REF = NA , ALT = NA)]
  # input.data[is.na(AF) , REF := NA ]




  ## add column for consistency with table version
  input.data[,DATE_ADDED := NA ]

  ## add std_ref to found variants
  input.data[!is.na(REF), SOURCE := 'Std_ref' ]

  # ' REF , ALT , AF , DATE_ADDED , SOURCE' columns are added to input data

  return(input.data)
}

compareInputfileWithStoredReferenceFile <- function(input.data)
{
  # search variants that are new in this file and are not previously searched
  new.hIDs <- unlist(input.data[!is.element(hID , .QC$searched.hIDs),hID])

  # search new hIDs
  rs <- search.database(new.hIDs)

  # some information
  # show how many variant information is in memory
#  print.and.log(paste0('Variant info from previous search: ',nrow(.QC$stored.reference.data)),
#                'info',
#                display=.QC$config$debug$verbose)
#
#  print.and.log(paste0('Variants searched before: ',length(.QC$searched.hIDs)),
#                'info',
#                display=.QC$config$debug$verbose)
#
#  print.and.log(paste0('Variant not seen before: ',length(new.hIDs)),
#                'info',
#                display=.QC$config$debug$verbose)
#
#  print.and.log(paste0('New found variants: ',nrow(rs)),
#                'info',
#                display=.QC$config$debug$verbose)
  ## ==========

  # add new hIDs to the previously searched hID list
  .QC$searched.hIDs <- c(.QC$searched.hIDs , new.hIDs)

  # bind old and new found reference data
  if(nrow(rs) > 0)
  {
    .QC$stored.reference.data <- list(.QC$stored.reference.data , rs)
    .QC$stored.reference.data <- rbindlist(.QC$stored.reference.data)
  }

  # add key if not existing
  if(is.null(data.table::key(.QC$stored.reference.data)) || data.table::key(.QC$stored.reference.data) != 'hID')
    data.table::setkey(.QC$stored.reference.data , hID)

  if(is.null(data.table::key(input.data)) || data.table::key(input.data) != 'hID')
    data.table::setkey(input.data , hID)


  input.data <- .QC$stored.reference.data[input.data]

  # multi allelic variants will be matched more than once and create duplicated lines.
  # only the first one is required and the rest should be removed.
  dup.allele <- which(duplicated(input.data,by=c('hID','EFFECT_ALL','OTHER_ALL','REF','ALT')))

  if(length(dup.allele) > 0)
  {
    #print.and.log(paste('found duplicates in reference matching :',length(dup.allele)),'warning')
    input.data <- input.data[!dup.allele]
  }


  # a row may match to more than one result due to multiple values in reference
  # first check if there are any duplicated rows that AF of that population is zero. keep the one with  AF>0
  # FIXME second, check if there are still duplicated hIDs. delete all
  # 1
  # e.g  A|GT  0
  #      A|GTGT 0.1
  # duplicate.hID.zeroAF <- input.data[which(duplicated(input.data$hID) |
  #                                           duplicated(input.data$hID , fromLast = TRUE))][AF == 0 ,]$hID

  # ignore variable of the duplicated variants on the same position is set to NA becuase one of them is removed in the next step
  #input.data[duplicate.hID.zeroAF, ignore := as.numeric(NA)]

  # input.data <- input.data[!which(is.element(input.data$hID, duplicate.hID.zeroAF) & AF == 0), ]


  # remove variants on the same chr:pos if it can not be found out
  # these are
  # 1- INS or DEL on the same position with different alleles
  # 2- different INDEL types on the same position (INS/DEL)
  input.data <- save.remove.ambiguous.variants(input.data)



  # find multi-allelic variants
  input.data[, MULTI_ALLELIC := ifelse(is.na(AF),
                                       NA,
                                       ifelse(grepl(',',AF),
                                              1,
                                              0))]



  ## try allele matching on multi-allelic variants
  # if(is.element('Yes',.QC$thisStudy$tables$multi_allele_count_preProcess$`Multi-allelic`))
  if(any(input.data$MULTI_ALLELIC == 1, na.rm = TRUE))
  {
    input.data[MULTI_ALLELIC == 1,
               c('ALT','AF') := clean.multi_alleles(EFFECT_ALL , OTHER_ALL, REF, ALT, AF) ,
               by = list(EFFECT_ALL , OTHER_ALL,REF, ALT,AF)]

    # some multi-allele INDEL AFs are all 0 and will be returned the same way due to missing alleles
    # e.g. AAC,AA,TT  0,0,0   ==> this AF can be converted to 0
    #input.data[VT == 2 & MULTI_ALLELIC == 1 &  grepl(',', AF) & all(strsplit(AF,',')[[1]] == "0") , AF := "0" ]
    input.data[VT == 2 & MULTI_ALLELIC == 1 &  grepl(',', AF) & !grepl('[1-9]',AF) , AF := "0" ]

    ## get frequency table for multi-allelic variants
    #.QC$thisStudy$tables$multi_allele_count_postProcess <- getMultiAlleleCountTbl(input.data,'AF')

  }




  # AF column may be character type due to remaining ',' => convert to numeric
  # AF of multi-allelics that could not be matched are set to NA
  if(!is.numeric(input.data$AF))
    input.data[, AF := as.numeric(AF)]


  # FIXME do not convert ALT to NA because it is used to count unmatched multiallelic variants
  #input.data[is.na(AF) , `:=` (REF = NA , ALT = NA)]
  # input.data[is.na(AF) , REF := NA ]




  ## add column for consistency with table version
  input.data[,DATE_ADDED := NA ]

  ## add std_ref to found variants
  input.data[!is.na(REF), SOURCE := 'Std_ref' ]

  # ' REF , ALT , AF , DATE_ADDED , SOURCE' columns are added to input data

  return(input.data)
}


search.database <- function(variantIDs)
{

  if(.QC$config$output_parameters$add_column_rsid &&  ## get REF_RSID if user has selected in config file
     is.element("REF_RSID",RSQLite::dbListFields(.QC$reference.data,'variants')))
  {
    rs <- RSQLite::dbGetQuery(.QC$reference.data,
                              sprintf('SELECT hID,REF_RSID,REF,ALT , %s as AF FROM variants WHERE "hID" = :x' ,
                                      .QC$config$supplementaryFiles$allele_ref_std_population) ,
                              param = list(x = variantIDs))
  } else {
    rs <- RSQLite::dbGetQuery(.QC$reference.data,
                              sprintf('SELECT hID,REF,ALT , %s as AF FROM variants WHERE "hID" = :x' ,
                                      .QC$config$supplementaryFiles$allele_ref_std_population) ,
                              param = list(x = variantIDs))
  }


  rs <- data.table::setDT(rs,key = 'hID')

  return(rs)

}




merge.rs <- function(id,rs)
{
  # var <- rs[fmatch(id,rs$hID),]
  # return(list(var$REF,var$ALT,var$AF))
}
