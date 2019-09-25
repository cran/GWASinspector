#### COLUMN NAMES
# CALLRATE
# CHR
# EFF_ALL_FREQ
# EFFECT
# EFFECT_ALL
# HWE_PVAL
# IMP_QUALITY
# IMPUTED
# MARKER
# N_TOTAL
# OTHER_ALL
# POSITION
# PVALUE
# STDERR
# STRAND




##CHR column should be as character type
process.column.CHR <- function(input.data){

  .QC$thisStudy$character.chromosome <- FALSE

  if(!is.numeric(input.data$CHR))
  {

    # check if column has character values
    # they are converted to numbers here, but should be deconverted to character before saving final dataset
    .QC$thisStudy$character.chromosome <- any(input.data$CHR %in% c('0X' ,'0Y','X','Y','XY','M','MT'))

    input.data[CHR == '0X', CHR := '23']
    input.data[CHR == 'X'  , CHR := '23']
    input.data[CHR == 'Y'  , CHR := '24']
    input.data[CHR == '0Y'  , CHR := '24']
    input.data[CHR == 'XY' , CHR := '25']
    input.data[CHR == 'M'  , CHR := '26']
    input.data[CHR == 'MT' , CHR := '26']


    #  input.data$CHR<-as.numeric(input.data$CHR)
    input.data[,CHR := as.numeric(CHR)]
  }

  # convert out of range chromosome values wiht NA
  invalid.items <- which(input.data$CHR < 1 | input.data$CHR > 26)

  .QC$thisStudy$column.INVALID.list$CHR <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items , CHR := NA]
  }

  .QC$thisStudy$column.NA.list$CHR <- which(is.na(input.data$CHR))


  return(input.data)
}

## convert chromosomes from numeric to character
deconvert.column.CHR <- function(input.data){

  if(!is.character(input.data$CHR))
    input.data[,CHR := as.character(CHR)]
  #input.data$CHR<-as.character(input.data$CHR)

  input.data[CHR == '23'  , CHR := 'X']
  input.data[CHR == '24'  , CHR := 'Y']
  input.data[CHR == '25' , CHR := 'XY']
  input.data[CHR == '26'  , CHR := 'M']

  return(input.data)
}

##########
process.column.EFFECT_ALL<- function(input.data){

  ## convert empty strings to NA
  input.data[trimws(EFFECT_ALL) == '' , EFFECT_ALL := NA]

  invalid.items <- which(!is.na(input.data$EFFECT_ALL) & !grepl(pattern = '^[0ATCGIDR9-]+', x = input.data$EFFECT_ALL ,perl = TRUE))

  .QC$thisStudy$column.INVALID.list$EFFECT_ALL <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, EFFECT_ALL := NA]
  }

  .QC$thisStudy$column.NA.list$EFFECT_ALL <- which(is.na(input.data$EFFECT_ALL))

  # chekc if allele has none Base characters
  if(any(is.element(input.data$EFFECT_ALL,c('D','I','-','0','R'))))
  {
    .QC$thisStudy$hanNoneBaseAlleles <- TRUE
    print.and.log('Input file has none base character for INDEL variants!','warning',display=.QC$config$debug$verbose)
  }

  return(input.data)
}

##########
process.column.OTHER_ALL<- function(input.data){

  ## convert empty strings to NA
  input.data[trimws(OTHER_ALL) == '' , OTHER_ALL := NA]

  invalid.items <- which(!is.na(input.data$OTHER_ALL) & !grepl(pattern = '^[0ATCGIDR9-]+', x = input.data$OTHER_ALL ,perl = TRUE))

  .QC$thisStudy$column.INVALID.list$OTHER_ALL <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, OTHER_ALL := NA]
  }

  .QC$thisStudy$column.NA.list$OTHER_ALL <- which(is.na(input.data$OTHER_ALL))

  return(input.data)
}


##########
process.column.IMPUTED<- function(input.data){
  config <- .QC$config
  # FIXME check if convert to uppercase is required for as.logical function
  # FIXME TRUE FALSE or 0 1


  # input.data$IMPUTED<-as.logical(input.data$IMPUTED) ## this is new conversion fron char to logical
  #converting from logical to numeric converts T to 1

  if(!is.numeric(input.data$IMPUTED)){

    input.data[,IMPUTED := toupper(IMPUTED)]

    input.data[!is.na(IMPUTED),  IMPUTED := gsub(pattern = config$input_parameters$imputed_T,
                                                 x = IMPUTED,
                                                 replacement = 1)]

    input.data[!is.na(IMPUTED),  IMPUTED := gsub(pattern = config$input_parameters$imputed_F,
                                                 x = IMPUTED,
                                                 replacement = 0)]


    # input.data$IMPUTED= gsub(pattern = config$input_parameters$imputed_T,
    #                          x = input.data$IMPUTED,
    #                          replacement = 1)
    # input.data$IMPUTED= gsub(pattern = config$input_parameters$imputed_F,
    #                          x = input.data$IMPUTED,
    #                          replacement = 0)

    #input.data$IMPUTED <- as.numeric(input.data$IMPUTED)
    input.data[,IMPUTED := as.numeric(IMPUTED)]
  }

  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  invalid.items <- which(input.data$IMPUTED != 1 & input.data$IMPUTED != 0)

  .QC$thisStudy$column.INVALID.list$IMPUTED <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, IMPUTED := NA]
  }


  .QC$thisStudy$column.NA.list$IMPUTED <- which(is.na(input.data$IMPUTED))



  return(input.data)
}


##########
process.column.STRAND<- function(input.data){

  input.data<-switchNegativeStrandsToPositive(input.data) ## variantModifierFUnction.R

  # convert empty strings to NA
  input.data[trimws(STRAND) == '' ,STRAND := NA]

  # all negative strands are converted to positive, so strand  is invalid if not NA or +
  invalid.items <- which(!is.na(input.data$STRAND) & input.data$STRAND != '+')
  .QC$thisStudy$column.INVALID.list$STRAND <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, STRAND := NA]
  }


  .QC$thisStudy$column.NA.list$STRAND <- which(is.na(input.data$STRAND))

  return(input.data)
}


##########
process.column.POSITION<- function(input.data){
  if(!is.numeric(input.data$POSITION))
    input.data[,POSITION := as.numeric(POSITION)]
  # input.data$POSITION<-as.numeric(input.data$POSITION)


  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  invalid.items <- which(input.data$POSITION <= 0)

  .QC$thisStudy$column.INVALID.list$POSITION <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, POSITION := NA]
  }

  ##
  .QC$thisStudy$column.NA.list$POSITION <- which(is.na(input.data$POSITION))

  return(input.data)
}

##########
process.column.EFFECT<- function(input.data){

  # ==== convert Beta to OR
  # odds_ratio = exp(Beta)
  # Beta = natural log(odds_taio)

  if(.QC$config$input_parameters$effect_type == 'BETA')
  {
    # change the name of beta column to EFFECT in input data
    names(input.data)[names(input.data) == 'BETA'] <- 'EFFECT'

    # change the name of beta column to EFFECT in renamed columns
    .QC$thisStudy$renamed.File.Columns.sorted[grepl(x = .QC$thisStudy$renamed.File.Columns.sorted, pattern = 'BETA')] =
      'EFFECT'
  }
  else
  {
    if(!is.numeric(input.data$OR))
      input.data[,OR := as.numeric(OR)]

    input.data[,EFFECT := log(OR)]

    # check if any of the EFFECT values are Inf and remove them
    if(any(is.infinite(input.data$EFFECT)))
    {
      infinite.ORs <- which(is.infinite(input.data$EFFECT))

      if(length(infinite.ORs) > 0)
      {
        print.and.log(paste0('variants with Infinite OR value are removed from input file: ',length(infinite.ORs)),'warning',display=.QC$config$debug$verbose)
        input.data[is.infinite(EFFECT), EFFECT := NA]
        # input.data <- input.data[!infinite.ORs,]
      }
    }
  }

  # ===

  if(!is.numeric(input.data$EFFECT))
    input.data[,EFFECT := as.numeric(EFFECT)]




  invalid.items <- which(input.data$EFFECT == -1 &
                           (input.data$PVALUE == -1 | input.data$STDERR == -1))

  .QC$thisStudy$column.INVALID.list$EFFECT <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, EFFECT := NA]
  }

  .QC$thisStudy$column.NA.list$EFFECT <- which(is.na(input.data$EFFECT))


  return(input.data)
}

##########
process.column.STDERR<- function(input.data){
  if(!is.numeric(input.data$STDERR))
    input.data[,STDERR := as.numeric(STDERR)]
  #input.data$STDERR<-as.numeric(input.data$STDERR)


  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  uncertain.items <- which(input.data$STDERR == 0) ## it may be due too poor rounding

  .QC$thisStudy$column.INVALID.list$zero.STDERR <- uncertain.items
  if(length(uncertain.items) > 0){
    input.data[uncertain.items, STDERR := NA]
  }


  invalid.items <- which(input.data$STDERR < 0) ## -1 is covered as < 0
  .QC$thisStudy$column.INVALID.list$STDERR <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, STDERR := NA]
  }

  .QC$thisStudy$column.NA.list$STDERR <- which(is.na(input.data$STDERR))


  return(input.data)
}


##########
process.column.PVALUE<- function(input.data){
  if(!is.numeric(input.data$PVALUE))
    input.data[,PVALUE := as.numeric(PVALUE)]
  # input.data$PVALUE<-as.numeric(input.data$PVALUE)


  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  uncertain.items <- which(input.data$PVALUE == -1 )
  .QC$thisStudy$column.INVALID.list$minusone.PVALUE <- uncertain.items

  if(length(uncertain.items) > 0){
    input.data[uncertain.items, PVALUE := NA]
  }


  invalid.items <- which(input.data$PVALUE > 1 | input.data$PVALUE <= 0)
  .QC$thisStudy$column.INVALID.list$PVALUE <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, PVALUE := NA]
  }

  .QC$thisStudy$column.NA.list$PVALUE <- which(is.na(input.data$PVALUE))

  return(input.data)
}

##########
process.column.EFF_ALL_FREQ<- function(input.data){
  if(!is.numeric(input.data$EFF_ALL_FREQ))
    input.data[,EFF_ALL_FREQ := as.numeric(EFF_ALL_FREQ)]
  #input.data$EFF_ALL_FREQ<-as.numeric(input.data$EFF_ALL_FREQ)


  .QC$thisStudy$column.INVALID.list$one.EFF_ALL_FREQ <- which(input.data$EFF_ALL_FREQ == 1)
  .QC$thisStudy$column.INVALID.list$zero.EFF_ALL_FREQ <- which(input.data$EFF_ALL_FREQ == 0)

  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  uncertain.items <- which(input.data$EFF_ALL_FREQ == -1 )
  .QC$thisStudy$column.INVALID.list$minusone.EFF_ALL_FREQ <- uncertain.items

  if(length(uncertain.items) > 0){
    input.data[uncertain.items, EFF_ALL_FREQ := NA]
  }



  invalid.items <- which(input.data$EFF_ALL_FREQ > 1 |
                           input.data$EFF_ALL_FREQ < 0) ## equal to 1 or 0 is monomorphic

  .QC$thisStudy$column.INVALID.list$EFF_ALL_FREQ <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, EFF_ALL_FREQ := NA]
  }


  .QC$thisStudy$column.NA.list$EFF_ALL_FREQ <- which(is.na(input.data$EFF_ALL_FREQ))

  return(input.data)
}

##########
process.column.HWE_PVAL<- function(input.data){
  if(!is.numeric(input.data$HWE_PVAL))
    input.data[,HWE_PVAL := as.numeric(HWE_PVAL)]
  #input.data$HWE_PVAL<-as.numeric(input.data$HWE_PVAL)



  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  uncertain.items <- which(input.data$HWE_PVAL == -1 )

  .QC$thisStudy$column.INVALID.list$minusone.HWE_PVAL <- uncertain.items

  if(length(uncertain.items) > 0){
    input.data[uncertain.items, HWE_PVAL := NA]

  }


  invalid.items <- which(input.data$HWE_PVAL > 1 |
                           input.data$HWE_PVAL <= 0)

  .QC$thisStudy$column.INVALID.list$HWE_PVAL <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, HWE_PVAL := NA]
  }


  .QC$thisStudy$column.NA.list$HWE_PVAL <- which(is.na(input.data$HWE_PVAL))

  # Fixed HWE P-value
  if(length(unique(input.data$HWE_PVAL)) == 1)
    .QC$thisStudy$fixed.hwep <- sprintf('YES (%s)' , input.data[1]$HWE_PVAL)
  else
    .QC$thisStudy$fixed.hwep <- 'No'


  return(input.data)
}

##########
process.column.IMP_QUALITY<- function(input.data){
  config <- .QC$config

  if(!is.numeric(input.data$IMP_QUALITY))
    input.data[,IMP_QUALITY := as.numeric(IMP_QUALITY)]
  #input.data$IMP_QUALITY<-as.numeric(input.data$IMP_QUALITY)


  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  invalid.items <- which(input.data$IMP_QUALITY <= config$filters$minimal_impQ_value |
                           input.data$IMP_QUALITY >= config$filters$maximal_impQ_value)

  .QC$thisStudy$column.INVALID.list$IMP_QUALITY <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, IMP_QUALITY := NA]
  }

  .QC$thisStudy$column.NA.list$IMP_QUALITY <- which(is.na(input.data$IMP_QUALITY))



  # Fixed imputation quality
  if(length(unique(input.data$IMP_QUALITY)) == 1)
    .QC$thisStudy$fixed.impq <-  sprintf('YES (%s)' , input.data[1]$IMP_QUALITY)
  else
    .QC$thisStudy$fixed.impq <- 'No'



  return(input.data)
}

##########
process.column.CALLRATE<- function(input.data){
  if(!is.numeric(input.data$CALLRATE))
    input.data[,CALLRATE := as.numeric(CALLRATE)]
  #input.data$CALLRATE<-as.numeric(input.data$CALLRATE)


  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  uncertain.items <- which(input.data$CALLRATE == -1 )
  .QC$thisStudy$column.INVALID.list$minusone.CALLRATE <- uncertain.items

  if(length(uncertain.items) > 0){
    input.data[uncertain.items, CALLRATE := NA]
  }


  invalid.items <- which(input.data$CALLRATE > 1 | input.data$CALLRATE < 0)

  .QC$thisStudy$column.INVALID.list$CALLRATE <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, CALLRATE := NA]
  }


  .QC$thisStudy$column.NA.list$CALLRATE <- which(is.na(input.data$CALLRATE))


  #  Fixed call rate
  if(length(unique(input.data$CALLRATE)) == 1)
    .QC$thisStudy$fixed.callrate <- sprintf('YES (%s)' , input.data[1]$CALLRATE)
  else
    .QC$thisStudy$fixed.callrate <- 'No'



  return(input.data)
}

##########
process.column.N_TOTAL<- function(input.data){
  if(!is.numeric(input.data$N_TOTAL))
    input.data[,N_TOTAL := as.numeric(N_TOTAL)]
  #input.data$N_TOTAL<-as.numeric(input.data$N_TOTAL)


  ## Check for inavlid or wrong or missing items => set them to NA
  ## used for report
  invalid.items <- which(input.data$N_TOTAL <= 0)

  .QC$thisStudy$column.INVALID.list$N_TOTAL <- invalid.items

  if(length(invalid.items) > 0){
    input.data[invalid.items, N_TOTAL := NA]
  }

  .QC$thisStudy$column.NA.list$N_TOTAL <- which(is.na(input.data$N_TOTAL))

  # maximum number of N
  .QC$thisStudy$MAX_N_TOTAL <- max(input.data$N_TOTAL,na.rm = TRUE)


  # Fixed sample size
  if(length(unique(input.data$N_TOTAL)) == 1)
    .QC$thisStudy$fixed.n_total <- sprintf('YES (%s)' , input.data[1]$N_TOTAL)
  else
    .QC$thisStudy$fixed.n_total <- 'No'


  return(input.data)
}


process.column.MARKER<- function(input.data){

  ## FIXME nothing set for invalid marker
  .QC$thisStudy$column.INVALID.list$MARKER <- numeric(0L)

  # TODO we are using hID so no need for rsID checking ?!?!
  .QC$thisStudy$column.NA.list$MARKER <- which(is.na(input.data$MARKER))

  return(input.data)
}
