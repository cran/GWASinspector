## added columns to dataset
# 'found' in allele.match() => variant is found in ref set
# 'palindromic' in allele.match() =>
# 'match' in allele.match() => eff_all == minor and oth_all == major
# 'flip' in allele.match() => variant needs flipping
# 'switch' in allele.match() => variant needs switching
# 'wring' in allele.match() => alleles do not match with ref set
# 'PVALUE.calculated' in replacePVALUE() => pvalue calculated from STDERR and EFFECT
# 'HQ' in applyHQfilter() => tag HQ variants based on AF, hwe_p , IMP_Q , callrate
# 'unusable' processcolumns() => missing or invalid crucial column
# 'highDiffEAF' alleleMatch() => items with highDiff PVALUE - AF
# for.AF checks if function is run for allele frequency or EFFECT PLOT

allele.match <- function(matched.data) {


  setkey(matched.data,EFFECT_ALL,OTHER_ALL,ALT,REF)

  ##FIXME might differ according to refrence data
  matched.data[REF == '', REF := NA] ## convert empty value to NA
  matched.data[!is.na(REF), found := TRUE] ## found in reference data


  # Palindromic variants:
  # A-T
  # vs A-T : match (1)
  # vs T-A : flip (2)
  # anything else : wrong (3)

  # non-Palindromic variants:
  # A-G
  # vs A-G : match (1)
  # vs G-A: flip (2)
  # vs C-T: flip (5)
  # vs T-C : switch (4)
  # anything else : wrong (6)


  ## A-T or G-C are palindrimic
  matched.data[VT == 1 & EFFECT_ALL == switch.allele.vectorized(OTHER_ALL), palindromic := TRUE]

  ## (1) matched exactly like the reference file
  matched.data[(found) & (EFFECT_ALL == ALT & OTHER_ALL == REF) , match := TRUE]


  ## (2) A-G in file vs G-A in reference file should be flipped
  # no difference in palindrommic or non-palindromic
  matched.data[(found) & (EFFECT_ALL == REF & OTHER_ALL == ALT), flip := TRUE]



  ## matching INDELs
  # FIXME the following condition is put because if INDELs have non base chars the can not be matched and will be
  # removed as a mismatch
  matched.data[.QC$thisStudy$hanNoneBaseAlleles & (found) & VT == 2 , match := TRUE]


  # (4) non-palindromic variants should be swithched if
  # A-C VS T-G
  matched.data[(found) &
                 is.na(palindromic) &
                 is.na(match) &
                 is.na(flip) &
                 (EFFECT_ALL == switch.allele.vectorized(ALT) &
                    OTHER_ALL == switch.allele.vectorized(REF)),
               switch := TRUE]

  # (5) non-palindromic variants should be flipped if
  # A-C VS G-T
  matched.data[(found) &
                 is.na(palindromic) &
                 is.na(match) &
                 is.na(switch) &
                 is.na(flip) &
                 (EFFECT_ALL == switch.allele.vectorized(REF) &
                    OTHER_ALL == switch.allele.vectorized(ALT)),
               flip := TRUE]




  # (6) non ppalindromic variants that are not fliped or switched or matched should be wrong
  # A-C VS A-G

  matched.data[(found) &
                 VT == 1 &
                 is.na(match) &
                 is.na(switch) &
                 is.na(flip) ,
               # is.na(ignore),
               wrong := TRUE]


  matched.data[(found) &
                 VT == 2 &
                 is.na(switch) &
                 is.na(flip) &
                 EFFECT_ALL == 'R',
               flip := TRUE]

  matched.data[(found) &
                 VT == 2 &
                 is.na(match) &
                 is.na(switch) &
                 is.na(flip) &
                 #  is.na(ignore) &
                 !is.element(EFFECT_ALL,c('I','D','0','-','R')),
               wrong := TRUE]


  return(matched.data)

}



allele.match.effectPlot <- function(matched.data) {

  ## TODO no need for allele matching because this is already done in previous steps
  ## and reference dataset is also matched with reference dataset

  return(matched.data)


  setkey(matched.data,EFFECT_ALL,OTHER_ALL,ALT,REF)


  # Palindromic variants:
  # A-T
  # vs A-T : match (1)
  # vs T-A : flip (2)
  # anything else : wrong (3)

  # non-Palindromic variants:
  # A-G
  # vs A-G : match (1)
  # vs G-A: flip (2)
  # vs C-T: flip (5)
  # vs T-C : switch (4)
  # anything else : wrong (6)



  matched.data[EFFECT_ALL == switch.allele.vectorized(OTHER_ALL), palindromic := TRUE]


  ## (1) matched exactly like the reference file
  matched.data[EFFECT_ALL == ALT & OTHER_ALL == REF , match := TRUE]


  ## (2) A-G in file vs G-A in reference file should be flipped
  # no difference in palindrommic or non-palindromic
  matched.data[EFFECT_ALL == REF & OTHER_ALL == ALT, flip := TRUE]



  # (4) non-palindromic variants should be swithched if
  # A-C VS T-G
  matched.data[is.na(palindromic) &
                 is.na(match) &
                 is.na(flip) &
                 (EFFECT_ALL == switch.allele.vectorized(ALT) &
                    OTHER_ALL == switch.allele.vectorized(REF)),
               switch := TRUE]

  # (5) non-palindromic variants should be flipped if
  # A-C VS G-T
  matched.data[ is.na(palindromic) &
                 is.na(match) &
                 is.na(switch) &
                 is.na(flip) &
                 (EFFECT_ALL == switch.allele.vectorized(REF) &
                    OTHER_ALL == switch.allele.vectorized(ALT)),
               flip := TRUE]

 # matched.data[ flip == TRUE, EFFECT := -1 * EFFECT]

  return(matched.data)

}


switch.flip.variant <- function(matched.data){
  # switched variants should be :
  # only snp alleles are switched
  matched.data[VT == 1 & switch == TRUE, `:=` (EFFECT_ALL = ALT ,
                                               OTHER_ALL = REF)]

  # flipped variants are switched and EFFECT and frequency are changed
  # FIXME only SNPs alleles are flipped
  # FIXME indel alleles are not changed , but indel AF is flipped
  if('EFF_ALL_FREQ' %in% colnames(matched.data))
  {
    matched.data[ VT == 1 & flip == TRUE, `:=` (EFFECT_ALL = ALT
                                                ,OTHER_ALL = REF
                                                ,EFFECT = -1 * EFFECT
                                                ,EFF_ALL_FREQ = 1 - EFF_ALL_FREQ)]

    matched.data[ VT == 2 & flip == TRUE, `:=` (EFFECT = -1 * EFFECT
                                                ,EFF_ALL_FREQ = 1 - EFF_ALL_FREQ)]

  }
  else
  {
    matched.data[ VT == 1 &  flip == TRUE, `:=` (EFFECT_ALL = ALT
                                                 ,OTHER_ALL = REF
                                                 ,EFFECT = -1 * EFFECT)]

    matched.data[ VT == 2 &  flip == TRUE, EFFECT := -1 * EFFECT]

  }

  # variants that are not found in reference dataset are changed if frequenct above 0.5
  #FIXME this step was discussed and removed
  # becuase in new ref sets like 1000g , minor allele is not used
  # matched.data[is.na(matched.data$POS) & matched.data$EFF_ALL_FREQ > 0.5 , EFF_ALL_FREQ := 1 - EFF_ALL_FREQ]
  # matched.data[is.na(matched.data$POS) & matched.data$EFF_ALL_FREQ > 0.5 , EFFECT := -1 * EFFECT]}

  return(matched.data)
}

process.matched.data <- function(matched.data) {

  #report about unknown and mismatches
  reportAlleleMatchStat(matched.data)

  # remove mismatches from dataset
  matched.data <- save.remove.mismatch.Variants(matched.data)


  if(nrow(matched.data) == 0)
  {
    print.and.log('ALL ROWS WERE DELETED AFTER ALLELE MATCHING! CHECK INPUT FILE FOR DATA INTEGRITY (step 3)!',
                  'warning')

    return(NULL)
  }

  variable.statistics.pre.matching(matched.data)

  # only for debugging purposes
  # checks for .QC$config$debug$save_pre_modification_file == TRUE
  save.pre_modification.file(matched.data)

  matched.data <- switch.flip.variant(matched.data)

  variable.statistics.post.matching(matched.data)

  #report variables with highDIff
  check.diffEAF(matched.data)

  return(matched.data)
}


reportAlleleMatchStat <- function(matched.data) {
  study <- .QC$thisStudy
  nrow.total <- study$rowcount.step2

  # variables that are found in standard reference file
  study$found.rows.std <- matched.data[found == TRUE & SOURCE == 'Std_ref' , .N]
  study$switched.rows.std <- matched.data[switch == TRUE & SOURCE == 'Std_ref' , .N]
  study$flipped.rows.std <- matched.data[flip == TRUE & SOURCE == 'Std_ref' , .N]
  study$mismatched.rows.std <- matched.data[wrong == TRUE & SOURCE == 'Std_ref' , .N]
  study$multiAlleleVariants.rowcount <- matched.data[MULTI_ALLELIC == 1 & grepl(',',ALT) , .N]

  # variables that are found in alternate reference file
  study$found.rows.alt <- matched.data[found == TRUE & SOURCE != 'Std_ref' , .N]
  study$switched.rows.alt <- matched.data[switch == TRUE & SOURCE != 'Std_ref' , .N]
  study$flipped.rows.alt <- matched.data[flip == TRUE & SOURCE != 'Std_ref' , .N]
  study$mismatched.rows.alt <- matched.data[wrong == TRUE & SOURCE != 'Std_ref' , .N]

  # # pallindromics
  # study$palindromic.rows <- matched.data[palindromic == TRUE , .N]

  # variables that are not found in standard reference file
  study$not.found.rows.std <- nrow.total - study$found.rows.std

  # variables that are not found in either standard or alternate reference file
  study$not.found.rows.alt <- nrow.total - study$found.rows.std - study$found.rows.alt



  ## table of how many variants where foundi each refrence file
  study$tables$match.ref.table <- t(table(matched.data[!is.na(SOURCE)]$SOURCE))
  if(nrow(study$tables$match.ref.table) == 0 )
    study$tables$match.ref.table <- data.table(x = 0)



  .QC$thisStudy <- study
}


check.diffEAF <- function(input.data) {
  study <- .QC$thisStudy


  DIFFthreshold <- .QC$config$filters$threshold_diffEAF

  if('EFF_ALL_FREQ' %in% colnames(input.data))
  {
    input.data[(found) & abs(EFF_ALL_FREQ - AF) > DIFFthreshold , highDiffEAF := TRUE]

    palindormicHighDiffEAF <- length(which( input.data$palindromic &
                                              input.data$highDiffEAF))

    palindormicExtremeDiffEAF <- length(which(input.data$palindromic &
                                                ((input.data$EFF_ALL_FREQ > 0.65 & input.data$AF < 0.35) |
                                                   (input.data$EFF_ALL_FREQ < 0.35 & input.data$AF > 0.65))
    )
    )



    nonpalindormicHighDiffEAF <- length(which(is.na(input.data$palindromic) &
                                                input.data$highDiffEAF))
    ##--
    study$palindormicHighDiffEAF <- ifelse(palindormicHighDiffEAF > 0 , palindormicHighDiffEAF, 0 )
    study$palindormicExtremeDiffEAF <- ifelse(palindormicExtremeDiffEAF > 0 , palindormicExtremeDiffEAF , 0)
    study$nonpalindormicHighDiffEAF <- ifelse(nonpalindormicHighDiffEAF > 0 , nonpalindormicHighDiffEAF , 0)
  }
  else
  {
    study$palindormicHighDiffEAF <- 'NA'
    study$palindormicExtremeDiffEAF <- 'NA'
    study$nonpalindormicHighDiffEAF <- 'NA'
  }
  ##--


  .QC$thisStudy <- study
}


save.remove.mismatch.Variants <- function(input.data) {

  # remove multi allelic variants
  # variants that are found in database but cannot be matched due to many variants on the same positions
  multiAllelic.rows<-which(input.data[,MULTI_ALLELIC == 1 & grepl(',', ALT)])

  if(length(multiAllelic.rows) > 0)
  {
    multiAllelic.data <- subset(input.data[head(multiAllelic.rows,100),] ,
                                select = .QC$thisStudy$renamed.File.Columns.sorted)
    print.and.log(sprintf("%s multi allelic variants removed from dataset (step 3)!",length(multiAllelic.rows)),
                  'warning',display=.QC$config$debug$verbose)
    saveDataSet(multiAllelic.data,
                .QC$thisStudy$SNPs_multi_allelic.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec)

    input.data <- input.data[!multiAllelic.rows,]
  }



  #----------------mismatch variants are reported
  #----- ----------and removed from dataset



  mismatched.rows<-which(input.data[, (wrong)])

  ## remove reference file columns from final dataset
  if(length(mismatched.rows) >0){
    mismatched.data <- subset(input.data[head(mismatched.rows,100),] ,
                              select = .QC$thisStudy$renamed.File.Columns.sorted)
    print.and.log(sprintf("%s mismatched variants removed from dataset (step 3)!",length(mismatched.rows)),
                  'warning',display=.QC$config$debug$verbose)
    saveDataSet(mismatched.data,
                .QC$thisStudy$SNPs_mismatches.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec)

    input.data <- input.data[!mismatched.rows,] ## remove mismatched from dataset

  }

  .QC$thisStudy$rowcount.step3 <- nrow(input.data)


  return(input.data)
}


save.remove.ambiguous.variants <- function(input.data)
{

  if(!.QC$thisStudy$hanNoneBaseAlleles)
    return(input.data)

  # remove varaints on the same chr:pos if it can not be found out
  # these are
  # 1- INS or DEL on the same position with different alleles
  # 2- different INDEL types on the same position (INS/DEL)

  # duplicate.positions <- which(duplicated(input.data[VT ==2,]$hID) |
  #                                duplicated(input.data[VT ==2,]$hID, fromLast = TRUE))

  duplicate.positions <- which(duplicated(input.data$hID) |
                                 duplicated(input.data$hID, fromLast = TRUE))

  if(length(duplicate.positions) > 0)
  {
    duplicate.positions.variants <- input.data[duplicate.positions,]

    #only one line if a duplicated variant in enough. remove others
    duplicate.positions.variants <- duplicate.positions.variants[!(duplicated(hID)),]

    .QC$thisStudy$ambiguos.rows <- nrow(duplicate.positions.variants)
    print.and.log(sprintf("%s ambiguous variants removed from dataset (step 3)!",.QC$thisStudy$ambiguos.rows),
                  'warning',display=.QC$config$debug$verbose)

    saveDataSet(subset(head(duplicate.positions.variants,100) ,
                                select = .QC$thisStudy$renamed.File.Columns.sorted),
                .QC$thisStudy$SNPs_ambiguous.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec)

    input.data <- input.data[!duplicate.positions, ]

    rm(duplicate.positions.variants)
    rm(duplicate.positions)
  }

  return(input.data)
}
