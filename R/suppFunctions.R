## negates in function
## for returning items that are NOT in a list
'%notin%' <- Negate('%in%')

## for more comprehensible code
## used for dataset filtering
# input.data[not(palindromic)]
not <- is.na

## swithc alleles
## TODO multicharacter alleles e.g.  'TCC' => 'AGG'
switch_allele<-function(allele)
{
  switched.allele<-switch(as.character(allele),
                          'A'='T',
                          'T'='A',
                          'C'='G',
                          'G'='C',
                          allele)

  return(switched.allele)
}


switch_allele_vectorized<-function(allele)
{

  switched.allele <- ifelse(allele == "A", "T",
                            ifelse(allele == "T", "A",
                                   ifelse(allele == "C", "G",
                                          ifelse(allele == "G", "C",allele))))
  return(switched.allele)
}

switch_alleles_vectorized<-function(allele)
{

  ifelse(nchar(allele) == 1,
         switch_allele_vectorized(allele), ## SNPs
         paste(sapply(strsplit(allele,''), switch_allele_vectorized),collapse = "")) # IN/DELs

}



############

getDefaultSystemOptions<-function()
{
  return(options())
}



# set traceback to NULL , nothing after Stop()
# turn back with options(error=traceback) and options(warn=0)
changeROptions<-function(){
  options(error= NULL)
  options(warn= -1)

  options(java.parameters = "- Xmx1024m") # for xlsx package

  if(capabilities('cairo'))
    options(bitmapType='cairo')
  else
    print_and_log('cairo module not available! run system_check()','warning')
}


##change R options back to what it was
resetDefaultSystemOptions<-function(user.options)
{
  options(user.options)
}




## FIXME not used - plots are displayed even if 1 variant is present
calculatePlotThreshold <- function(input.data) {
  config <- .QC$config
  rowCount <- nrow(input.data)
  plot.itemcount.threshold <- abs(config$thresholds$use_threshold) ## get absolute value in case of negative values


  if(plot.itemcount.threshold <= 1){

    plot.itemcount.threshold <- as.integer(rowCount * plot.itemcount.threshold)

    return(plot.itemcount.threshold)

  }else if(plot.itemcount.threshold > 1 &&  plot.itemcount.threshold > rowCount){  #check if threshold is not more than file rows

    print_and_log('Plot threshold more than file size! value set to 10 percent.','warning',display=.QC$config$debug$verbose)

    plot.itemcount.threshold <- as.integer(rowCount / 10)

    return(plot.itemcount.threshold)

  }else{
    return(as.integer(plot.itemcount.threshold))
  }
}


calculatePercent <- function(sample.count,
                             total.count,
                             decimal.place=2,
                             number=FALSE,
                             pretty= FALSE)
{

  if(sample.count > total.count |
     sample.count < 0 |
     total.count <= 0 |
     !is.numeric(sample.count) |
     !is.numeric(total.count) |
     is.na(sample.count))
  {
    #print_and_log('Error calculating percent value!','fatal')
    return(NA)
  }

  p <- sample.count / total.count

  if(number)
    return(signif(p , decimal.place)) # exact value - used for comparisons
  else
  {
    if(pretty)
      return(sprintf('%s (%s)' , format(sample.count, big.mark="," , scientific = FALSE) ,paste(signif((p * 100),decimal.place), '%',sep = ''))) # converted to => n (p%) with thousand seprator
    else
      return(paste(signif((p * 100),decimal.place), '%',sep = '')) # converted to => p%
  }

}





check_tools <- function()
{
  existing.packages <- installed.packages()[,1]


  # fileFunctions.R

  # used for synchronizing column sep character in all rows
  .QC$awk.exists <- check_awk()

  # used for counting file lines. important to make sure if all rows has been read by fread.
  # number of lines and number of variants should match
  .QC$wc.exists <- check_wc()

  .QC$java.exists <- check_java()

  # used for reading zipped files.
  # fread does not have a built in reader for zip files
  .QC$gzip.exists <- check_gzip()
  .QC$unzip.exists <- check_unzip()


  .QC$xlsx.package.exists <- check_xlsx_package(existing.packages)

  .QC$rsqlite.package.exists <- check_rsqlite_package(existing.packages)

  .QC$kableExtra.package.exists <- check_kableExtra_package(existing.packages)

  .QC$rJava.package.exists <- check_rJava_package(existing.packages)

  .QC$r.version<- get_R_version()

  .QC$pandoc.exists <-  check_pandoc()

  .QC$ggplot2.version <- check_ggplot2_version(existing.packages)

  .QC$OS <- get_OS()
}



## used for checking config file input variables
# return false is parameter is not set by user or not found
is_empty <- function(parameter) {
  #if(!exists(deparse(substitute(parameter))) || is.null(parameter) || parameter == '')
  if(is.null(parameter) || parameter == '')
    return(TRUE)
  else
    return(FALSE)
}


thousand_sep <- function(input)
{
  tryCatch(
    return(format(as.numeric(input), big.mark="," , scientific = FALSE)),
    error = function(war) return(NA)
  )

}


increment <- function(value,factor)
{
  if(is.numeric(value) & is.numeric(factor))
    return(value + factor)
  else
    return(NA)
}


# convert triallelic alleles to bi-allelic according to allele frequency in a population
# A,C  0.1,0  => A 0.1
# A,C  0,0  => A,C 0
# A,C  0.1,0.5  => A,C 0.1,0.5
tri_to_bi_alleles <- function(alleles,freqs) {
  afs <- strsplit(freqs,',')
  alts <- strsplit(alleles,',')

  if(all(afs[[1]] == '0')  | all(afs[[1]] != '0'))
  {
    alleles
  }
  else
  {
    zero <- which(afs[[1]] == '0')
    paste(alts[[1]][-zero] , collapse = ',')
  }
}

tri_to_bi_freq <- function(alleles,freqs) {
  afs <- strsplit(freqs,',')
  alts <- strsplit(alleles,',')

  if(all(afs[[1]] == '0'))
    '0'
  else if (all(afs[[1]] != '0'))
    freqs
  else
  {

    zero <- which(afs[[1]] == '0')
    paste(afs[[1]][-zero] , collapse = ',')
  }
}


clean_multi_alleles <- function(A1,A2,REF,ALT,AF)
{
  tbl <- cbind(data.table('allele' = strsplit(ALT,',')[[1]]),
               data.table('freqs' = strsplit(AF,',')[[1]]))

  # default value - input arguments are returened if variant can not be handled
  output <- list(ALT,AF)

  tryCatch(
    {

      if(A1 == REF && is.element(A2,tbl$allele)) # perfect match
      {
        variant <- tbl[allele == A2,]
        output <- list(variant$allele, variant$freqs)
      }
      else if(A2 == REF && is.element(A1,tbl$allele)) # flipped variants
      {
        variant <- tbl[allele == A1,]
        output <- list(variant$allele, variant$freqs)
      }
      else if (A1 == switch_alleles_vectorized(REF) && is.element(switch_alleles_vectorized(A2),tbl$allele)) # flipped and switched variants
      {
        variant <- tbl[allele == switch_alleles_vectorized(A2),]
        output <- list(variant$allele, variant$freqs)
      }
      else if(A2 == switch_alleles_vectorized(REF) && is.element(switch_alleles_vectorized(A1),tbl$allele)) # switched variants
      {
        variant <- tbl[allele == switch_alleles_vectorized(A1),]
        output <- list(variant$allele, variant$freqs)
      }

      return(output)
    },
    error = function(err) {
      print_and_log(sprintf('Error in variant %s , %s , %s , %s',A1,A2,REF,ALT),'warning',display=.QC$config$debug$verbose)
      return(output)
    }
  )

}


getMultiAlleleCountTbl <- function(input.data) {

  # tbl <- input.data[,.N, keyby=.(VT,is.na(REF), MULTI_ALLELIC)]
  tbl <- input.data[,.N, keyby=.(VT,(SOURCE != "Std_ref" | is.na(SOURCE)), MULTI_ALLELIC)]

  tbl <- t(data.table(
    'Bi-allelic SNP' =  ifelse(length(tbl[VT == 1 & !(SOURCE) & !MULTI_ALLELIC,N]) == 0 ,
                               '0',
                               format(tbl[VT == 1 & !(SOURCE) & !MULTI_ALLELIC,N],big.mark = ',',scientific = FALSE)),
    'Multi-allelic SNP' =  ifelse(length(tbl[VT == 1 & !(SOURCE) & MULTI_ALLELIC,N]) == 0 ,
                                  '0' ,
                                  format(tbl[VT == 1 & !(SOURCE) & MULTI_ALLELIC,N],big.mark = ',',scientific = FALSE)),
    'SNPs not found in standard reference dataset' =  ifelse(length(tbl[VT == 1 & (SOURCE) ,N]) == 0 ,
                                                             '0' ,
                                                             format(tbl[VT == 1 & (SOURCE),N],big.mark = ',',scientific = FALSE)),
    'Bi-allelic INDEL' = ifelse(length(tbl[VT == 2 & !(SOURCE) & !MULTI_ALLELIC,N]) == 0 ,
                                '0' ,
                                format( tbl[VT == 2 & !(SOURCE) & !MULTI_ALLELIC,N],big.mark = ',',scientific = FALSE)),
    'Multi-allelic INDEL' =  ifelse(length(tbl[VT == 2 & !(SOURCE) & MULTI_ALLELIC,N]) == 0 ,
                                    '0' ,
                                    format(tbl[VT == 2 & !(SOURCE) & MULTI_ALLELIC,N],big.mark = ',',scientific = FALSE)),
    'INDELs not found in standard reference dataset' =  ifelse(length(tbl[VT == 2 & (SOURCE) ,N]) == 0 ,
                                                               '0' ,
                                                               format(tbl[VT == 2 & (SOURCE) ,N],big.mark = ',',scientific = FALSE))
  ))

  colnames(tbl) <- 'count'


  return(tbl)
}

##not used anymore
getVariantTypeCountTbl <- function(input.data)
{
  rowCount <- nrow(input.data)
  tbl <- input.data[,.(.N),keyby=VT]
  tbl$VT <- as.character(tbl$VT )
  tbl[VT == '1' , VT := 'SNP']
  tbl[VT == '2' , VT := 'INDEL']
  names(tbl) <- c('Variant type','Count')

  tbl[,Count := calculatePercent(Count, rowCount, pretty = TRUE)]
  return(tbl)
}
