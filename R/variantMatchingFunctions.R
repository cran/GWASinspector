compareInputfileWithReferenceData <- function(input.data)
{
  #  .QC$reference.data ==> this can be a data table or a database object
  #  .QC$alt.reference.data
  input.file.colNames <- names(input.data)

  ## 1
  if(is.data.table(.QC$reference.data))
    input.data <- compareInputfileWithReferenceFile(input.data)
  else
    input.data <- compareInputfileWithReferenceDataBase(input.data)


  ## 2
  ## check the unknown variants with alternative refernce file
  ## this refrence file will be empty if is not set by user and this step is automatically skipped

  input.data <- tryCatch(compareInputfileWithAlternateReferenceFile(input.data, input.file.colNames),
                         error = function(x)
                         {
                           print.and.log(paste('Error in searching alternate reference database:',x$message),
                                         'warning')
                           return(input.data)
                         })

  return(input.data)
}

compareInputfileWithAlternateReferenceFile <- function(input.data,input.file.colNames)
{
  #check if alt ref data has any rows & there are unfound variants
  if(nrow(.QC$alt.reference.data) > 0 &
     length(which(is.na(input.data$REF))) > 0 )
  {
    print.and.log('Comparing input file with alternate reference file ...','info')
    # we only want to check unfound variants in the alt ref dataset
    # we only need original file columns, so , added columns from checking with reference file from previous step should be deleted
    if(.QC$thisStudy$hID.added)
    {
      setkey(.QC$alt.reference.data,"hID")
      setkey(input.data,"hID")

      tmp.data <- merge(x = subset(input.data[is.na(REF) & MULTI_ALLELIC == 0 ,] ,
                                   select = input.file.colNames),
                        y = .QC$alt.reference.data,
                        by.x = "hID",
                        by.y = "hID",
                        all.x = TRUE)
    }
    else
    {

      setkey(input.data,"MARKER")
      setkey(.QC$alt.reference.data,"ID")

      tmp.data <- merge(x = subset(input.data[is.na(REF),] ,
                                   select = input.file.colNames),
                        y = subset(.QC$alt.reference.data,
                                   select = c('ID','REF','ALT','AF','DATE_ADDED','SOURCE')),
                        by.x = "MARKER",
                        by.y = "ID",
                        all.x = TRUE)
    }


    # if(!is.element('TSA', names(tmp.data)))
    #   tmp.data[,TSA := NA] ## added for consistency with reference data

    # reomve ID column from alternate reference file
    if(is.element('ID', names(tmp.data)))
      tmp.data[,ID := NULL] ## removed for consistency with reference data

    if(!is.element('MULTI_ALLELIC', names(tmp.data)))
      tmp.data[,MULTI_ALLELIC := as.numeric(NA)] ## added for consistency with reference data

    # if(!is.element('ignore', names(tmp.data)))
    #   tmp.data[,ignore := NA] ## added for consistency with reference data

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
  .QC$thisStudy$tables$multi_allele_count_preProcess <- getMultiAlleleCountTbl(input.data,'AF')

  return(input.data)
}



compareInputfileWithBetaReferenceFile<-function(input.data)
{


  m_b <- mean(input.data$EFFECT)
  sd_b <- sd(input.data$EFFECT)
  ci_b <- m_b + (5 * sd_b)


  if(is.element('hID', names(input.data)))
  {
    input.data <- input.data[,c('EFFECT_ALL', 'OTHER_ALL',  'EFFECT','HQ','hID')]
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
    print.and.log('hID is missing. Can not compare input file with reference dataset.','warning',display=.QC$config$debug$verbose)
    return(data.table())
  }


  matched.data <- matched.data[HQ == TRUE &
                                 !is.na(EFFECT.y) &
                                 EFFECT_ALL == ALT &
                                 OTHER_ALL == REF &
                                 EFFECT.y > (-1 * ci_b) &
                                 EFFECT.y < ci_b , ]



  # p-value of the variants in refernce dataset is < 0.001
  .QC$thisStudy$effect.rho_3 <- signif(cor(matched.data$EFFECT.x ,matched.data$EFFECT.y),4)


  temp_data <- matched.data[PVALUE < 0.0001,]
  .QC$thisStudy$effect.rho_4 <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)

  temp_data <- matched.data[PVALUE < 0.00001,]
  .QC$thisStudy$effect.rho_5 <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)


  temp_data <- matched.data[PVALUE < 0.000001,]
  .QC$thisStudy$effect.rho_6 <- signif(cor(temp_data$EFFECT.x ,temp_data$EFFECT.y),4)

  return(matched.data[PVALUE < 0.0001,])
}
