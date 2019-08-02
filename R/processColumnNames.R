
changeColumnNames <- function(input.colnames,selected.header.format){


  if(selected.header.format == 'STANDARD')
    return(input.colnames)

  ###


  ## 1
  if(selected.header.format == 'GENABEL')
    new.col.names <- sapply(input.colnames,changeColumnNameGenABEL)


  ## 2
  else if(selected.header.format == 'GWAMA')
    new.col.names <- sapply(input.colnames,changeColumnNameGWAMA)


  ## 3
  else if(selected.header.format == 'PLINK')
    new.col.names <- sapply(input.colnames,changeColumnNamePLINK)


  ## 4
  else if(selected.header.format == 'META')
    new.col.names <- sapply(input.colnames,changeColumnNameMETA)


  ## 5
  else if(selected.header.format == 'GCTA')
    new.col.names <- sapply(input.colnames,changeColumnNameGCTA)

  return(new.col.names)
}



changeColumnNameGenABEL <- function(input.colName)
{
  value <- switch(toupper(input.colName),
                  "MARKER"="name",
                  "CHR"="chromosome",
                  "POSITION"="position",
                  "STRAND"="strand",
                  "EFFECT_ALL"="allele1",
                  "OTHER_ALL"="allele2",
                  "EFF_ALL_FREQ"="effallelefreq",
                  "N_TOTAL"="n",
                  "EFFECT"="beta",
                  "STDERR"="sebeta",
                  "PVALUE"="p",
                  "HWE_PVAL"="pexhwe",
                  "CALLRATE"= "call",
                  #"EFFECT_ALL"="effallele",
                  "PGC"="pgc", ##### not in our data
                  "LAMBDA.ESTIMATE"="lambda.estimate",  ##### not in our data
                  "LAMBDA.SE"="lambda.se",  ##### not in our data
                  "BUILD"="build",  ##### not in our data
                  input.colName
  )

  return(value)
}



changeColumnNameGWAMA <- function(input.colName)
{
  value <- switch(toupper(input.colName),
                  "MARKER"="MARKER",
                  "CHR"="CHR",
                  "POSITION"="POSITION",
                  "STRAND"="STRAND",
                  "EFFECT_ALL"="EA",
                  "OTHER_ALL"="NEA",
                  "EFF_ALL_FREQ"="EAF",
                  "N_TOTAL"="N",
                  "EFFECT"="BETA",
                  "STDERR"="SE",
                  "PVALUE"="P",
                  "CALLRATE"= "N",
                  "IMPUTED"="IMPUTED",
                  "IMP_QUALITY"="IMP_QUALITY",
                  input.colName
  )
}


changeColumnNamePLINK <- function(input.colName)
{
  value <- switch(toupper(input.colName),
                  "MARKER"="SNP",
                  "CHR"="CHR",
                  "POSITION"="BP",
                  "STRAND"="STRAND",
                  "EFFECT_ALL"="A1",
                  "OTHER_ALL"="A2",
                  "EFF_ALL_FREQ"="EFF_ALL_FREQ",
                  "N_TOTAL"="N",
                  "EFFECT"="BETA",
                  "STDERR"="SE",
                  "PVALUE"="P",
                  "CALLRATE"= "N",
                  "IMPUTED"="IMPUTED",
                  "IMP_QUALITY"="IMP_QUALITY",
                  input.colName
  )
}


changeColumnNameMETA <- function(input.colName)
{
  value <- switch(toupper(input.colName),
                  "MARKER"="rsid",
                  "CHR"="chr",
                  "POSITION"="pos",
                  "STRAND"="strand",
                  "EFFECT_ALL"="allele_B", ## FIXME A or B !!!
                  "OTHER_ALL"="allele_A",
                  "EFF_ALL_FREQ"="EFF_ALL_FREQ",
                  "N_TOTAL"="N",
                  "EFFECT"="beta",
                  "STDERR"="se",
                  "PVALUE"="P_value",
                  "CALLRATE"= "N",
                  "IMPUTED"="imputed",
                  "INFO"="info",
                  input.colName
  )
}


changeColumnNameGCTA <- function(input.colName)
{
  value <- switch(toupper(input.colName),
                  "MARKER"="SNP",
                  "EFFECT_ALL"="A1",
                  "OTHER_ALL"="A2",
                  "EFF_ALL_FREQ"="freq",
                  "N_TOTAL"="N",
                  "EFFECT"="b",
                  "STDERR"="se",
                  "PVALUE"="p",
                  input.colName
  )

  return(value)
}


## this functtion is not used anymore because order is set since loading the file
changeColumnOrder <- function(input.data) {
  col.order <- c( 'MARKER','CHR','POSITION','EFFECT_ALL','OTHER_ALL',
                  'STRAND','EFFECT','STDERR','PVALUE','EFF_ALL_FREQ',
                  'HWE_PVAL','CALLRATE','N_TOTAL','IMPUTED','IMP_QUALITY')

  col.count <- ncol(input.data)
  col.index.all <- c(1:col.count)
  # put the known headers in first columns
  col.index.wanted <- na.omit(match( col.order, names(input.data)))
  col.index.unknown <- col.index.all[-col.index.wanted]
  # put unknown headers after known columns
  col.index.sorted <- c(col.index.wanted ,  col.index.unknown)
  # sort columns
  input.data <- input.data[,col.index.sorted,with=FALSE]

  return(input.data)
}
