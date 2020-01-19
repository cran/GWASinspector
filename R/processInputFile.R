processInputFile <- function(input.data) {



  # saved files
  #       a.	[filename_output]_SNPs_invalid_allele.txt
  #           Invalid other-allele values (max. 30)
  #       b.	[filename_output]_SNPs_duplicates.txt
  #           Duplicate SNP names
  #       c.	[filename_output]_SNPs_removed.txt
  #           Unusable SNPs
  #       d.	[filename_output]_SNPs_improbable_values.txt
  #           SNPs with invalid values (max. 1000 entries are saved)
  #

  #### step 1: checking crucial variables

  #### step 2: saving unusable variants [SNPs are removed from dataset]
  ## list of SNPs with missing or invalid crucial variables ---- “[filename]_SNPs_removed.txt


  #### step 3: checking non-crucial variables


  #### step 4: saving invalid variants  [SNPs are removed from dataset]
  # •	Monomorphic SNPs
  #   o	Identical alleles
  #   o	Missing or invalid effect/other allele  => [these are removed in previous step]
  #   o	Allele frequency = 1 or 0
  #      [filename]_SNPs_invalid_OTHER_ALL.txt in previous version and [filename]_SNPs_invalid_ALLELE.txt

  #   o	(duplicates)
  #     [filename]_SNPs_duplicates.txt

  #### step 5: saving imrprobable datas [SNPs are NOT removed from dataset]
  ## invalid values in non-crucial variables ---- [filename]_SNPs_improbable_values.txt





  #### step 1: checking crucial variables

  column.names <- colnames(input.data)
  .QC$thisStudy$input.data.rowcount <- nrow(input.data)



  input.data.backup  <-  input.data ##a copy of input data is kept for final report (without NA)
  input.data.backup <- as.data.table(input.data.backup)

  input.data <- tryCatch(processCrucialColumns(input.data), #is in this file
                         error = function(err) {
                           print.and.log(paste('Error in processing crucial columns:' , err$message), 'warning')
                           return(NULL)
                         }
  )

   ## ------------------------------------------------
  #### step 2: saving unusable variants - SNPs with missing or invalid crucial variables
  ## find rows with NA in crucial columns
  ## first 100 is saved
  ## data is removed from dataset
  ## [filename_output]_SNPs_removed.txt

  input.data <-tryCatch(save.and.remove.unusable.variants(input.data,input.data.backup), #saveFileFunctions.R
                        error = function(err) {
                          print.and.log(paste('Error in removing crucial columns:' , err$message), 'warning')
                          return(NULL)
                        }
  )

  if(is.null(input.data))
    return(NULL)

  # remove duplicate Marker names and save as separate file
  input.data <- tryCatch(removeDuplicateVariants(input.data),
                         error = function(err) {
                           print.and.log(paste('Error in processing duplicated columns:' , err$message), 'warning')
                           return(NULL)
                         }
  )

  if(is.null(input.data))
    return(NULL)

  invisible(gc())


  ##===============================
  remaining.rows<-nrow(input.data)
  .QC$thisStudy$rowcount.step1 <- remaining.rows

  if(remaining.rows == 0)
  {
    print.and.log('ALL ROWS WERE DELETED WHILE CHECKING CRUCIAL VAIRABLES! CHECK INPUT FILE FOR DATA INTEGRITY!',
                  'warning')
    return(NULL)
  }



  #########   END OF STEP 1 variant removal #########


  ## ------------------------------------------------
  #### step 3: checking non-crucial variables

  input.data.backup  <-  input.data ##a copy of input data is kept for final report (without NA)
  input.data.backup <- as.data.table(input.data.backup)

  invisible(gc())

  input.data <- tryCatch(processNonCrucialColumns(input.data),
                         error = function(err) {
                           print.and.log(paste('Error in processing non-crucial columns:' , err$message), 'warning')
                           return(NULL)
                         }
  )

  if(is.null(input.data))
    return(NULL)

  ## ---------------------------------------------
  #### step 5: saving imrprobable datas [SNPs are NOT removed from dataset]
  ## invalid values in non-crucial variables ---- [filename]_SNPs_improbable_values.txt

  save.NA.Dataset(input.data , input.data.backup)

  rm(input.data.backup) # file is saved and no need for this variable anymore
  invisible(gc())
  ## ------------------------------------------------
  #### step 4: saving invalid variants
  #### clean the file from monomorphic, duplicates ,chromosomes


  input.data<-tryCatch(removeMonomorphicVariants(input.data), ## variantModifierFUnction.R
                       error = function(err) {
                         print.and.log(paste('Error in processing monomorphic columns:' , err$message), 'warning')
                         return(NULL)
                       }
  )

  if(is.null(input.data))
    return(NULL)


  invisible(gc())

  # find INDELs and SNPs
  input.data <- tryCatch(variantDiscrimination(input.data),
                         error = function(err) {
                           # all variants are set as SNP if there was an error in  this phase
                           print.and.log(paste('Error in discriminating variant types:' , err$message), 'warning')
                           print.and.log('All variants are set as SNP.', 'warning')
                           input.data[, VT := 1]
                           return(input.data)
                         }
  )


  ## remove chormosomal snps based on user input
  input.data<-tryCatch(removeChromosomeVariants(input.data), # variantModifierFUnctions.R
                       error = function(err) {
                         print.and.log(paste('Error in processing chromosal variants:' , err$message), 'warning')
                         return(NULL)
                       }
  )

  if(is.null(input.data))
    return(NULL)

  #### check how many rows are left####

  remaining.rows<-nrow(input.data)
  .QC$thisStudy$rowcount.step2 <- remaining.rows

  if(remaining.rows == 0)
  {
    print.and.log('ALL ROWS WERE DELETED IN STEP2! CHECK INPUT FILE FOR DATA INTEGRITY!',
                  'warning')
    return(NULL)
  }


  #########   END OF STEP 2 variant removal #########



  return(input.data)
}



processNonCrucialColumns <- function(input.data) {

  column.names <- colnames(input.data)


  # THIS HAS BECOME A CRUCIAL COLUMN FOR hID
  # if('CHR' %in% column.names)
  #   input.data <- tryCatch(process.column.CHR(input.data),
  #                          error = function(err) {
  #                            print.and.log(paste('Error in processing CHR:' , err$message), 'warning')
  #                            return(NULL)
  #                          }
  #   )


  if('IMPUTED' %in% column.names)
    input.data <- tryCatch(process.column.IMPUTED(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing IMPUTED:' , err$message), 'warning')
                             return(NULL)
                           }
    )

  if('MARKER' %in% column.names)
    input.data <- tryCatch(process.column.MARKER(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing MARKER:' , err$message), 'warning')
                             return(NULL)
                           }
    )


  if('STRAND' %in% column.names)
    input.data <- tryCatch(process.column.STRAND(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing STRAND:' , err$message), 'warning')
                             return(NULL)
                           }
    )

  if('PVALUE' %in% column.names)
    input.data <- tryCatch(process.column.PVALUE(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing PVALUE:' , err$message), 'warning')
                             return(NULL)
                           }
    )

  if('EFF_ALL_FREQ' %in% column.names)
    input.data <- tryCatch(process.column.EFF_ALL_FREQ(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing EFF_ALL_FREQ:' , err$message), 'warning')
                             return(NULL)
                           }
    )

  if('HWE_PVAL' %in% column.names)
    input.data <- tryCatch(process.column.HWE_PVAL(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing HWE_PVAL:' , err$message), 'warning')
                             return(NULL)
                           }
    )

  if('IMP_QUALITY' %in% column.names)
    input.data <- tryCatch(process.column.IMP_QUALITY(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing IMP_QUALITY:' , err$message), 'warning')
                             return(NULL)
                           }
    )

  if('CALLRATE' %in% column.names)
    input.data <- tryCatch(process.column.CALLRATE(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing CALLRATE:' , err$message), 'warning')
                             return(NULL)
                           }
    )

  if('N_TOTAL' %in% column.names)
    input.data <- tryCatch(process.column.N_TOTAL(input.data),
                           error = function(err) {
                             print.and.log(paste('Error in processing N_TOTAL:' , err$message), 'warning')
                             return(NULL)
                           }
    )
  return(input.data)
}


processCrucialColumns <- function(input.data) {

  input.data <- process.column.EFFECT_ALL(input.data)
  input.data <- process.column.OTHER_ALL(input.data)
  input.data <- process.column.EFFECT(input.data)
  input.data <- process.column.STDERR(input.data)
 # input.data <- process.column.MARKER(input.data)
  input.data <- process.column.CHR(input.data)
  input.data <- process.column.POSITION(input.data)

  return (input.data)
}


variantDiscrimination <- function(input.data) {

  # these are SNP conditions
  input.data[is.element(EFFECT_ALL , c('A','G','C','T')) & is.element(OTHER_ALL , c('A','G','C','T'))
             , VT := 1]

  # the rest is INDEL
  input.data[is.na(VT), VT := 2]


  # check if input file has INDEL values
  if(input.data[VT == 2, .N] > 0)
    .QC$thisStudy$hasINDEL <- TRUE

  # check if reference data has INDEL values

  # if(.QC$thisStudy$hasINDEL & !.QC$reference.data.has.INDEL)
  #   print.and.log('Input file has INDEL variants but Reference dataset does not include any!','warning')


  if('CHR' %in% names(input.data)) # if chr column exists
  {
    # hID is added if CHR exists
    input.data <- add_hIDcolumn(input.data)
  }

  return(input.data)
}




## if it is called on input file => data.file is true
## if it is called on oither files like effect size reference file => data.file is FALSE
add_hIDcolumn <- function(input.data, data.file = TRUE){

  if('CHR' %notin% names(input.data)){ # can not make hID if there is not a chr column

    print.and.log('CHR column is missing! variant matching will be done by rsID.','warning')

  }
  else if(is.element('hID', names(input.data))){ # do not make hid if already exists in file

    print.and.log('hID column already existed in input file and was not generated!','warning')

    if(data.file)
      .QC$thisStudy$hID.added <- TRUE
  }
  else{
    # VT column is already added to file
    input.data[,hID := paste(CHR,POSITION,VT,sep = ':')]
     #input.data[,hID := sprintf('%s:%s_%s_%s',CHR,POSITION,EFFECT_ALL,OTHER_ALL)]


    if(data.file)
      .QC$thisStudy$hID.added <- TRUE
  }

  return(input.data)


}
