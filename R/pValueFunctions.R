calculate.pvalue.correlation <- function(input.data) {

  # if the pvalue coumn is missing , correlation will not be calculated
  if(.QC$thisStudy$missing.PVALUE.column)
    return(input.data)


    # get non-NA pvalues
    data <- subset(input.data[,c('palindromic', 'PVALUE','PVALUE.calculated')],
                   !is.na(PVALUE) & !is.na(PVALUE.calculated) & PVALUE != 0 & PVALUE.calculated != 0)

    PVcor <- signif(cor(-log10(data$PVALUE),
                        -log10(data$PVALUE.calculated)),4)

    # PVcor.palindromic <- signif(cor(-log10(data[palindromic == TRUE]$PVALUE),
    #                                        -log10(data[palindromic == TRUE]$PVALUE.calculated)),3)

    .QC$thisStudy$rownum.PVcor <- nrow(data)
    .QC$thisStudy$PVcor <- PVcor
   # .QC$thisStudy$PVcor.palindromic <- PVcor.palindromic


}



calculate.PVALUE <- function(input.data)
{
  ## calculate pvalue based on stderr and beta

  ## this is used for
  ## 1- correlation checking
  ## 2- filling missing pvalues in data

  # STEPS
  # 1- calculate PVALUE
  # 2- check if column is totally missing ,create the columns and set as the calculated value


  input.data[,PVALUE.calculated := pchisq((EFFECT/STDERR)^2, 1, lower.tail=FALSE)]

  if(.QC$thisStudy$missing.PVALUE.column)
  {
    input.data$PVALUE <- input.data$PVALUE.calculated
    print.and.log('PVALUE column is created and filled from calculated values!','info')
  }




  return(input.data)
}


## this function is run at the end of algorithm before saving the final file.
## because missing pvalues should not be analyzed during QC
fill.missing.pvalues.from.calculated.pvalues <- function(input.data){

  if(.QC$config$parameters$calculate_missing_p){

    # check for missing PVALUES and repkkcae the missings
    input.data[is.na(PVALUE),
               PVALUE := PVALUE.calculated]
    print.and.log('Missing PVALUES replaced according to EFFECT and STDERR values!','info')
  }
  return(input.data)
}


## correct extreme values
#this is only done for plottings and does not affect the input data
correct.extreme.pvalues <- function(input.data){


  input.data[PVALUE < 10^-300, PVALUE := 10^-300]

  return(input.data)
}

correct.extreme.calculated.pvalues <- function(input.data){

  input.data[PVALUE.calculated < 10^-300, PVALUE.calculated := 10^-300]

  return(input.data)
}
