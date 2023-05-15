compareInputfileWithReferenceData <- function(input.data)
{
  #  .QC$reference.data ==> this can be a data table or a database object
  #  .QC$alt.reference.daOta
  ## 1
  if(!is.null(.QC$stored.reference.data))
    input.data <- compareInputfileWithStoredReferenceFile(input.data)
  else if(is.data.table(.QC$reference.data))
    input.data <- compareInputfileWithReferenceFile(input.data)
  else
    input.data <- compareInputfileWithReferenceDataBase(input.data)


  ## 2
  ## check the unknown variants with alternative reference file
  ## this refrence file will be empty if is not set by user and this step is automatically skipped

  input.data <- tryCatch(compareInputfileWithAlternateReferenceFile(input.data),
                         error = function(x)
                         {
                           print_and_log(paste('Error in searching alternate reference database:',x$message),
                                         'warning')
                           return(input.data)
                         })

  return(input.data)
}

compareInputfileWithAlternateReferenceFile <- function(input.data)
{
  #check if alt ref data has any rows & there are unfound variants
  if(nrow(.QC$alt.reference.data) > 0 & input.data[is.na(REF),.N] > 0 )
  {
    print_and_log('Comparing input file with alternate reference file ...','info')
    # we only want to check unfound variants in the alt ref dataset
    # we only need original file columns, so , added columns from checking with reference file from previous step should be deleted
    if(.QC$thisStudy$hID.added)
    {

      # double check KEY property for merging data
      if(is.null(data.table::key(.QC$alt.reference.data)) || data.table::key(.QC$alt.reference.data) != 'hID')
        data.table::setkey(.QC$alt.reference.data , hID)

      if(is.null(data.table::key(input.data)) || data.table::key(input.data) != 'hID')
        data.table::setkey(input.data , hID)


      # FIXME change merge to data table join
      tmp.data <- merge(x = input.data[is.na(REF) ,!c("REF","ALT","SOURCE","DATE_ADDED","AF")],
                        y = .QC$alt.reference.data,
                        by.x = "hID",
                        by.y = "hID",
                        all.x = TRUE)
    }
    else
    {

      setkey(input.data,"MARKER")
      setkey(.QC$alt.reference.data,"ID")

      tmp.data <- merge(x = input.data[is.na(REF) ,!c("REF","ALT","SOURCE","DATE_ADDED","AF")],
                        y = subset(.QC$alt.reference.data,
                                   select = c('ID','REF','ALT','AF','DATE_ADDED','SOURCE')),
                        by.x = "MARKER",
                        by.y = "ID",
                        all.x = TRUE)
    }



    # reomve ID column from alternate reference file
    if(is.element('ID', names(tmp.data)))
      tmp.data[,ID := NULL] ## removed for consistency with reference data

    if(!is.element('MULTI_ALLELIC', names(tmp.data)))
      tmp.data[,MULTI_ALLELIC := as.numeric(NA)] ## added for consistency with reference data

    # if(!is.element('ignore', names(tmp.data)))
    #   tmp.data[,ignore := NA] ## added for consistency with reference data

    tmp.data[,DATE_ADDED := as.character(DATE_ADDED)]
    input.data[,DATE_ADDED := as.character(DATE_ADDED)]

    # bind matched data with reference set with  matched data from alt reference set
    input.data<-rbind(input.data[!is.na(REF),],tmp.data)
    rm(tmp.data)
  }

  return(input.data)
}


compareInputfileWithReferenceFile<-function(input.data)
{

   # merge using hID if it Exists
  if(.QC$thisStudy$hID.added)
  {

    #set key for fast access
    setkey(.QC$reference.data,"hID")
    setkey(input.data,"hID")

    input.data<-merge(x=input.data,
                      y=.QC$reference.data,
                      by.x="hID",
                      by.y="hID",
                      all.x=TRUE)


  }
  else # merge using rsID if hID does not Exist
  {

    setkey(input.data,"MARKER")
    setkey(.QC$reference.data,"ID")

    input.data<-merge(x=input.data,
                      y=subset(.QC$reference.data,
                               select = c('ID','REF','ALT','AF','DATE_ADDED','SOURCE')),
                      by.x="MARKER",
                      by.y="ID",
                      all.x=TRUE)

  }

  # find multi-allelic variants
  input.data[, MULTI_ALLELIC := ifelse(grepl(',',AF),1,0)]

  ## get frequency table for multi-allelic variants
  ## moved to after step 3
  ##.QC$thisStudy$tables$multi_allele_count_preProcess <- getMultiAlleleCountTbl(input.data,'AF')


  ## TODO the following sectoin is similar to lines 194-234 rSQLiteFunctions.R
  ## try allele matching on multi-allelic variants
  # if(is.element('Yes',.QC$thisStudy$tables$multi_allele_count_preProcess$`Multi-allelic`))
  if(any(input.data$MULTI_ALLELIC == 1, na.rm = TRUE))
  {
    input.data[MULTI_ALLELIC == 1,
               c('ALT','AF') := clean_multi_alleles(EFFECT_ALL , OTHER_ALL, REF, ALT, AF) ,
               by = list(EFFECT_ALL , OTHER_ALL,REF, ALT,AF)]

    # some multi-allele INDEL AFs are all 0 and will be returned the same way due to missing alleles
    # e.g. AAC,AA,TT  0,0,0   ==> this AF can be converted to 0
    #input.data[VT == 2 & MULTI_ALLELIC == 1 &  grepl(',', AF) & all(strsplit(AF,',')[[1]] == "0") , AF := "0" ]
    input.data[VT == 2 & MULTI_ALLELIC == 1 &  grepl(',', AF) & !grepl('[1-9]',AF) , AF := "0" ]

    ## get frequency table for multi-allelic variants
    #.QC$thisStudy$tables$multi_allele_count_postProcess <- getMultiAlleleCountTbl(input.data,'AF')

  }




  # AF column may be character type due to remaining ',' => convert to numeric
  # AF of multi-allelics than could notbe matched are set as NA
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



compareInputfileWithBetaReferenceFile<-function(input.data)
{


  m_b <- mean(input.data$EFFECT)
  sd_b <- sd(input.data$EFFECT)
  ci_b <- m_b + (5 * sd_b)


  if(is.element('hID', names(input.data)))
  {
    input.data <- input.data[,c('EFFECT_ALL', 'OTHER_ALL',  'EFFECT','HQ','hID','PVALUE')]
    setkey(.QC$reference.data.effect, "hID")
    #setkey(input.data, "hID")

    ## 1
    # merge data with reference file
    matched.data<-merge(x=.QC$reference.data.effect,
                        y=input.data,
                        by="hID",
                        all.x=TRUE)

  }
  else
  {
    print_and_log('hID is missing. Can not compare input file with reference dataset.','warning',display=.QC$config$debug$verbose)
    return(data.table())
  }


  matched.data <- matched.data[HQ == TRUE &
                                 !is.na(EFFECT.y) &
                                 EFFECT_ALL == ALT &
                                 OTHER_ALL == REF &
                                 EFFECT.y > (-1 * ci_b) &
                                 EFFECT.y < ci_b , ]



  # p-value of the variants in reference dataset is < 0.001
  # X: reference Beta
  # Y: input file Beta
  .QC$thisStudy$effect.rho_3 <- signif(cor(matched.data$EFFECT.x ,matched.data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_3.n <- nrow(matched.data)

  temp_data <- matched.data[PVALUE.x < 0.0001,]
  .QC$thisStudy$effect.rho_4 <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_4.n <- nrow(temp_data)

  temp_data <- matched.data[PVALUE.x < 0.00001,]
  .QC$thisStudy$effect.rho_5 <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_5.n <- nrow(temp_data)


  temp_data <- matched.data[PVALUE.x < 0.000001,]
  .QC$thisStudy$effect.rho_6 <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_6.n <- nrow(temp_data)


    # filtering the variants in input rsult files on P_value
  temp_data <- matched.data[PVALUE.y < 0.001,]
  .QC$thisStudy$effect.rho_3.y <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_3.n.y <- nrow(temp_data)


  temp_data <- matched.data[PVALUE.y < 0.0001,]
  .QC$thisStudy$effect.rho_4.y <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_4.n.y <- nrow(temp_data)

  temp_data <- matched.data[PVALUE.y < 0.00001,]
  .QC$thisStudy$effect.rho_5.y <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_5.n.y <- nrow(temp_data)


  temp_data <- matched.data[PVALUE.y < 0.000001,]
  .QC$thisStudy$effect.rho_6.y <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)
  .QC$thisStudy$effect.rho_6.n.y <- nrow(temp_data)


  .QC$thisStudy$tables$betaCor.tbl <- t(data.table( 'P-value < 0.001' = c(sprintf("%s (%s)",.QC$thisStudy$effect.rho_3 , .QC$thisStudy$effect.rho_3.n),
                                                                          sprintf("%s (%s)",.QC$thisStudy$effect.rho_3.y , .QC$thisStudy$effect.rho_3.n.y)),
                                                    'P-value < 0.0001' = c(sprintf("%s (%s)",.QC$thisStudy$effect.rho_4 , .QC$thisStudy$effect.rho_4.n),
                                                                           sprintf("%s (%s)",.QC$thisStudy$effect.rho_4.y , .QC$thisStudy$effect.rho_4.n.y)),
                                                    'P-value < 0.00001' = c(sprintf("%s (%s)",.QC$thisStudy$effect.rho_5 , .QC$thisStudy$effect.rho_5.n),
                                                                            sprintf("%s (%s)",.QC$thisStudy$effect.rho_5.y , .QC$thisStudy$effect.rho_5.n.y)),
                                                    'P-value < 0.000001' = c(sprintf("%s (%s)",.QC$thisStudy$effect.rho_6 , .QC$thisStudy$effect.rho_6.n),
                                                                             sprintf("%s (%s)",.QC$thisStudy$effect.rho_6.y , .QC$thisStudy$effect.rho_6.n.y))))

  colnames(.QC$thisStudy$tables$betaCor.tbl) <- c('r*','r**')

  return(matched.data[PVALUE.x < 0.0001,])
}
