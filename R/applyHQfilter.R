applyHQfilter <- function(input.data) {


  # HQ parameters
  HQfilter_FRQ <- as.numeric(.QC$config$filters$HQfilter_FRQ)
  HQfilter_HWE <- .QC$config$filters$HQfilter_HWE
  HQfilter_cal <- as.numeric(.QC$config$filters$HQfilter_cal)
  HQfilter_imp <- as.numeric(.QC$config$filters$HQfilter_imp)


  col.names <- colnames(input.data)
  row.count <- nrow(input.data)

  input.data[, HQ := TRUE]# SET ALL VARIANTS AS HIGH-QUALITY AT FIRST

  ## 1 - if all 4 columns are missing, all variants will stay HQ
  ## 2- if any of the 4 columns are missing they will not affect HQ/LQ selection
  ## 3- consider variants as LQ, if EFF_ALL_FREQ or IMP_QULAITY are NA
  ## 4- consider variants as HQ, if HWE_PVAL or CALLRATE are NA



  # check if allele frequencies are above threshold
  # pvalue threshold is bidirectional.
  if('EFF_ALL_FREQ' %in% col.names)
    input.data[ is.na(EFF_ALL_FREQ) |
                  EFF_ALL_FREQ < HQfilter_FRQ |
                  EFF_ALL_FREQ > 1 - HQfilter_FRQ ,
                HQ := FALSE]

  if('IMP_QUALITY' %in% col.names)
    input.data[ HQ == TRUE   & (is.na(IMP_QUALITY) | IMP_QUALITY < HQfilter_imp) ,
                HQ := FALSE]

  if('HWE_PVAL' %in% col.names)
    input.data[HQ == TRUE & !is.na(HWE_PVAL) & HWE_PVAL < HQfilter_HWE, # set as LQ if below threshold
               HQ := FALSE]

  if('CALLRATE' %in% col.names)
    input.data[HQ == TRUE & !is.na(CALLRATE) & CALLRATE < HQfilter_cal,# set as LQ if below threshold
               HQ := FALSE]




  # report the number

  HQtable <- table(input.data$HQ)
  HQ.count <- ifelse(is.na(HQtable['TRUE']) , 0 , HQtable['TRUE'])
  LQ.count <- ifelse(is.na(HQtable['FALSE']) , 0 , HQtable['FALSE'])


  .QC$thisStudy$HQ.count <- HQ.count
  .QC$thisStudy$LQ.count <- LQ.count




  return(input.data)
}
