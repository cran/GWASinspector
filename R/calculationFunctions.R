calculate_af_correlation_std_ref <- function(input.data) {

  if('EFF_ALL_FREQ' %notin% colnames(input.data))
    return(NULL)

  # remove NA values from ALLELE FREQ
  # remove Multu allelic variants from dataset
  input.data <- input.data[!is.na(EFF_ALL_FREQ) & !is.na(AF)]

  if(nrow(input.data)  == 0)
    return(NULL)


  # only variants that are matched with standard reference are required
  cor.test.rho.all<- signif(cor(input.data[VT == 1 & SOURCE == 'Std_ref']$EFF_ALL_FREQ ,
                                input.data[VT == 1 & SOURCE == 'Std_ref']$AF),3)

  cor.test.rho.non.palindromic <- signif(cor(input.data[VT == 1 & palindromic == FALSE & SOURCE == 'Std_ref']$EFF_ALL_FREQ ,
                                             input.data[VT == 1 & palindromic == FALSE & SOURCE == 'Std_ref']$AF),3)

  cor.test.rho.palindromic <- signif(cor(input.data[VT == 1 & palindromic == TRUE  & SOURCE == 'Std_ref']$EFF_ALL_FREQ ,
                                         input.data[VT == 1 & palindromic == TRUE  & SOURCE == 'Std_ref']$AF),3)

  if(.QC$thisStudy$hasINDEL)
    .QC$thisStudy$AFcor.std_ref.indel <- signif(cor(input.data[VT == 2 & SOURCE == 'Std_ref']$EFF_ALL_FREQ ,
                                                    input.data[VT == 2 & SOURCE == 'Std_ref']$AF),3)


  .QC$thisStudy$AFcor.std_ref <- cor.test.rho.all
  .QC$thisStudy$AFcor.palindromic.std_ref <- cor.test.rho.palindromic
  .QC$thisStudy$AFcor.non.palindromic.std_ref <- cor.test.rho.non.palindromic


  # calculate AF correlation for each chromosome
  .QC$thisStudy$AFcor.std_ref.CHR <- sapply(split(input.data,input.data$CHR),
                                            function(x)
                                              c(as.character(x[1]$CHR),round(cor(x$AF,x$EFF_ALL_FREQ),3)))

	.QC$thisStudy$AFcor.std_ref.CHR <- as.data.table(.QC$thisStudy$AFcor.std_ref.CHR)
	.QC$thisStudy$AFcor.std_ref.CHR <- t(.QC$thisStudy$AFcor.std_ref.CHR)
	colnames(.QC$thisStudy$AFcor.std_ref.CHR) <- c("Chromosome","AF correlation")

}

calculate_af_correlation_alt_ref <- function(input.data) {

  if('EFF_ALL_FREQ' %notin% colnames(input.data))
    return(NULL)

  # remove NA values from ALLELE FREQ
  input.data <- input.data[!is.na(EFF_ALL_FREQ)]

  # only variants that are matched with alternate reference are required

  cor.test.rho.all<- signif(cor(input.data[VT == 1 & !is.na(SOURCE) &  SOURCE != 'Std_ref']$EFF_ALL_FREQ ,
                                input.data[VT == 1 & !is.na(SOURCE) &  SOURCE != 'Std_ref']$AF),3)

  cor.test.rho.non.palindromic <- signif(cor(input.data[VT == 1 & palindromic == FALSE & !is.na(SOURCE) & SOURCE != 'Std_ref']$EFF_ALL_FREQ ,
                                             input.data[VT == 1 & palindromic == FALSE & !is.na(SOURCE) & SOURCE != 'Std_ref']$AF),3)

  cor.test.rho.palindromic <- signif(cor(input.data[VT == 1 & palindromic == TRUE  & !is.na(SOURCE) & SOURCE != 'Std_ref']$EFF_ALL_FREQ ,
                                         input.data[VT == 1 & palindromic == TRUE  & !is.na(SOURCE) & SOURCE != 'Std_ref']$AF),3)

  if(.QC$thisStudy$hasINDEL)
    .QC$thisStudy$AFcor.alt_ref.indel<- signif(cor(input.data[VT == 2 & SOURCE != 'Std_ref']$EFF_ALL_FREQ ,
                                                   input.data[VT == 2 & SOURCE != 'Std_ref']$AF),3)


  .QC$thisStudy$AFcor.alt_ref <- cor.test.rho.all
  .QC$thisStudy$AFcor.palindromic.alt_ref <- cor.test.rho.palindromic
  .QC$thisStudy$AFcor.non.palindromic.alt_ref <- cor.test.rho.non.palindromic


}



calculateSkewness <- function(input.data){

  skewness <- signif(input.data[!is.na(EFFECT),
                                sum( (EFFECT  - mean(EFFECT))^3) / ((length(EFFECT)-1) * sd(EFFECT) ^ 3 )],3)

  .QC$thisStudy$skewness <- skewness



}

calculateSkewness_HQ <- function(input.data){

  skewness <- signif(input.data[HQ ==TRUE,
                                sum( (EFFECT  - mean(EFFECT))^3) / ((length(EFFECT)-1) * sd(EFFECT) ^ 3 )],3)


  .QC$thisStudy$skewness.HQ <- skewness



}


calculateKurtosis <- function(input.data){

  kurtosis <- signif(input.data[!is.na(EFFECT),
                                sum( (EFFECT  - mean(EFFECT))^4) / ((length(EFFECT)-1) * sd(EFFECT) ^ 4 )],3)

  .QC$thisStudy$kurtosis <- kurtosis




}

calculateKurtosis_HQ <- function(input.data){

  kurtosis <- signif(input.data[HQ ==TRUE,
                                sum( (EFFECT  - mean(EFFECT))^4) / ((length(EFFECT)-1) * sd(EFFECT) ^ 4 )],3)

  .QC$thisStudy$kurtosis.HQ <- kurtosis


}



calculateVischerStats <- function(input.data){

  if('EFF_ALL_FREQ' %notin% colnames(input.data))
    return(NULL)


  input.data <- input.data[!is.na(EFF_ALL_FREQ)]
  nMax = max(input.data$N_TOTAL ,na.rm = TRUE)
  Visschers.stat <- signif(input.data[,median(2 *
                                                input.data$EFF_ALL_FREQ *
                                                (1-input.data$EFF_ALL_FREQ) *
                                                nMax *
                                                (input.data$STDERR^2))], 3)


  .QC$thisStudy$Visschers.stat <- Visschers.stat

}

calculateVischerStats_HQ <- function(input.data){

  if('EFF_ALL_FREQ' %notin% colnames(input.data))
    return(NULL)

  input.data <- input.data[!is.na(EFF_ALL_FREQ)]

  nMax = max(input.data$N_TOTAL ,na.rm = TRUE)

  Visschers.stat <- signif(input.data[HQ == TRUE, median(2 *
                                                input.data$EFF_ALL_FREQ *
                                                (1-input.data$EFF_ALL_FREQ) *
                                                nMax *
                                                (input.data$STDERR^2))], 3)

  .QC$thisStudy$Visschers.stat.HQ <- Visschers.stat

}


calculateLambda <- function(input.data){

  lambda <- 'NA'
  lambda.imp <- 'NA'
  lambda.gen <- 'NA'

  if('PVALUE' %in% colnames(input.data))
    lambda <- signif(input.data[!is.na(PVALUE), median(qchisq(PVALUE, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)],3)

  if(all(c('IMPUTED','PVALUE') %in% colnames(input.data)))
  {
    lambda.imp <-  signif(input.data[!is.na(PVALUE) & IMPUTED == 1, median(qchisq(PVALUE, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)],3)
    lambda.gen <- signif(input.data[!is.na(PVALUE) & IMPUTED == 0, median(qchisq(PVALUE, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)],3)
  }

  .QC$thisStudy$lambda <- lambda
  .QC$thisStudy$lambda.imp <- lambda.imp
  .QC$thisStudy$lambda.gen <- lambda.gen


}
