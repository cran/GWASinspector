##this file includes functions for editing or removing variants
## removeVariants => calls the second functions based on config
## removeChromoseVariants => remove all variants of a chromosome based on config
## switchNegativeStrandsToPositive => change strand value and allele values for negative stranded rows
## removeInvalidVariants => remove raw from dataset because of an invalid data

removeChromosomeVariants<-function(input.data){

	config <- .QC$config

  ## delete chromose data
	if(any(config$remove_chromosomes$remove_X ,
	       config$remove_chromosomes$remove_Y,
	       config$remove_chromosomes$remove_XY,
	       config$remove_chromosomes$remove_M))
    input.data <- find_and_remove_ChromosomeVariant(input.data)


  return(input.data)
}


find_and_remove_ChromosomeVariant<-function(input.data){
	config <- .QC$config

	print_and_log("Removing specified chromosomes ...")

  ##REMOVE a chromosome from input file
  if(config$remove_chromosomes$remove_X == TRUE){
    r.index <- which(input.data$CHR == 23)
    .QC$thisStudy$x.chr.count.removed <- length(r.index)

    print_and_log(sprintf('data from chromosome X will be deleted! (%s rows)',
                          thousand_sep(.QC$thisStudy$x.chr.count.removed)),'warning',display=.QC$config$debug$verbose)
    input.data<-input.data[!r.index,]
  }

  if(config$remove_chromosomes$remove_Y == TRUE){
    r.index <- which(input.data$CHR == 24)
    .QC$thisStudy$y.chr.count.removed <- length(r.index)

    print_and_log(sprintf('data from chromosome Y will be deleted! (%s rows)',
                          thousand_sep(.QC$thisStudy$y.chr.count.removed)),'warning',display=.QC$config$debug$verbose)

    input.data<-input.data[!r.index,]
  }

  if(config$remove_chromosomes$remove_XY == TRUE){
    r.index <- which(input.data$CHR == 25)
    .QC$thisStudy$xy.chr.count.removed <- length(r.index)

    print_and_log(sprintf('data from chromosome XY will be deleted! (%s rows)',
                          thousand_sep(.QC$thisStudy$xy.chr.count.removed)),'warning',display=.QC$config$debug$verbose)
    input.data<-input.data[!r.index,]
  }

  if(config$remove_chromosomes$remove_M == TRUE){
    r.index <- which(input.data$CHR == 26)
    .QC$thisStudy$m.chr.count.removed <- length(r.index)

    print_and_log(sprintf('data from chromosome M will be deleted! (%s rows)',
                          thousand_sep(.QC$thisStudy$m.chr.count.removed)),'warning',display=.QC$config$debug$verbose)
    input.data<-input.data[!r.index,]
  }
  ####

  return(input.data)
}

##TODO could be converted to search by reference to be faster
switchNegativeStrandsToPositive<-function(input.data)
{

  negative.strand.index<-which(input.data$STRAND == '-')


  if(length(negative.strand.index) > 0)
  {

    input.data[negative.strand.index,'EFFECT_ALL'] <- apply(input.data[negative.strand.index,'EFFECT_ALL'],
                                                          1,
                                                          function(x) switch_allele(x))

    input.data[negative.strand.index,'OTHER_ALL'] <- apply(input.data[negative.strand.index,'OTHER_ALL'],
                                                         1,
                                                         function(x) switch_allele(x))

    input.data[negative.strand.index,'STRAND'] <- '+'

    print_and_log(sprintf('\'%s\' Negative strands found and switched!',thousand_sep(length(negative.strand.index))),
                  'info')
    .QC$thisStudy$neg.strand.count <- length(negative.strand.index)
  }

  return(input.data)
}



removeMonomorphicVariants <- function(input.data){

  ## find monomorphic allele indexes

  EFFECT_ALL.equal.OTHER_ALL <- which((input.data$EFFECT_ALL == input.data$OTHER_ALL))

  ## get union of lists to be removed
  monomorphic.alleles <- Reduce(union, list(.QC$thisStudy$column.INVALID.list$one.EFF_ALL_FREQ,
                                            .QC$thisStudy$column.INVALID.list$zero.EFF_ALL_FREQ,
                                            EFFECT_ALL.equal.OTHER_ALL)
  )


  if(length(monomorphic.alleles) > 0)
  {


    .QC$thisStudy$monomorphic.count <- length(monomorphic.alleles)
    # _SNPs_monomorphic.txt
    saveDataSet(input.data[head(monomorphic.alleles,100),],
                .QC$thisStudy$SNPs_monomorphic.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec,
				ordered = .QC$config$output_parameters$ordered)

    print_and_log(sprintf('\'%s\' monomorphic variants removed from file (step 2)!',
                          thousand_sep(length(monomorphic.alleles))),
                  'warning',display=.QC$config$debug$verbose)

    input.data<-input.data[!monomorphic.alleles,]
  }

  return(input.data)
}


removeDuplicateVariants <- function(input.data)
{

	setkey(input.data,CHR,POSITION,EFFECT_ALL,OTHER_ALL)

	dup.allele <- which(
                        # duplicated(input.data$hID) | duplicated(input.data$hID, fromLast = TRUE)
                        duplicated(input.data, by = key(input.data)) |
                        duplicated(input.data, by = key(input.data), fromLast = TRUE)
                    )

  if(length(dup.allele) > 0)
  {

    # save duplicate rows to output folder
    saveDataSet(input.data[head(dup.allele,100), .QC$thisStudy$renamed.File.Columns.sorted , with = FALSE],
                .QC$thisStudy$SNPs_duplicates.path,
                columnSeparator = .QC$config$output_parameters$out_sep,
                naValue = .QC$config$output_parameters$out_na,
                decValue = .QC$config$output_parameters$out_dec,
				ordered = .QC$config$output_parameters$ordered)



    tbl <- input.data[dup.allele, .N ,keyby=CHR]

    print_and_log('duplicated variants distribution in input file...','info',display=.QC$config$debug$verbose)
    print_and_log(kable(tbl,format = "rst"),
                  'info',
                  cat= FALSE,
                  display= .QC$config$debug$verbose)




    # remove duplicates from dataset
    input.data<-input.data[!dup.allele,]

    .QC$thisStudy$duplicate.count <- length(dup.allele)

    print_and_log(sprintf('\'%s\' duplicated variants removed from file.',
                          thousand_sep(length(dup.allele))),
                  'warning',display=.QC$config$debug$verbose)
  }

  return(input.data)
}

removeDuplicateVariants_postMatching <- function(input.data)
{

	setkey(input.data,CHR,POSITION,EFFECT_ALL,OTHER_ALL)

	dup.allele <- which(
                        # duplicated(input.data$hID) | duplicated(input.data$hID, fromLast = TRUE)
                        duplicated(input.data, by = key(input.data)) |
                        duplicated(input.data, by = key(input.data), fromLast = TRUE)
                    )

  if(length(dup.allele) > 0)
  {

    # save duplicate rows to output folder
#     saveDataSet(input.data[head(dup.allele,100), .QC$thisStudy$renamed.File.Columns.sorted , with = FALSE],
#                 .QC$thisStudy$SNPs_duplicates_postMatch,
#                 columnSeparator = .QC$config$output_parameters$out_sep,
#                 naValue = .QC$config$output_parameters$out_na,
#                 decValue = .QC$config$output_parameters$out_dec,
# 				ordered = .QC$config$output_parameters$ordered)


    # remove duplicates from dataset
    input.data<-input.data[!dup.allele,]

    print_and_log(sprintf('\'%s\' duplicated variants removed from file after matching with reference dataset.',
                          thousand_sep(length(dup.allele))),
                  'warning',display=.QC$config$debug$verbose)
  }

  return(input.data)
}

