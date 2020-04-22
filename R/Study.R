create.Study <- function(study)
{

  object <- new(
    "Study",
    File = list(
      file.path = study$file.path,
      file.name = study$file.name,
      file.extension = study$file.extension,
      file.line.count = study$file.line.count,
      file.line.count = study$file.endsWithNewLine,
      dup_lines_count = study$dup_lines_count,
      original.File.Columns = study$original.File.Columns,
      renamed.File.Columns.sorted = study$renamed.File.Columns.sorted
    ),
    Counts = list(
      input.data.rowcount = study$input.data.rowcount,
      duplicate.count = study$duplicate.count,
      rowcount.step1 = study$rowcount.step1,
      rowcount.step2 = study$rowcount.step2,
      rowcount.step3 = study$rowcount.step3,
      found.rows = study$found.rows.std + study$found.rows.alt,
      mismatched.rows = study$mismatched.rows.std + study$mismatched.rows.alt,
      ambiguos.rows = study$ambiguos.rows,
      switched.rows = study$switched.rows.std + study$switched.rows.alt,
      flipped.rows = study$flipped.rows.std + study$flipped.rows.alt,
      monomorphic.count = study$monomorphic.count,
      palindromic.rows = study$palindromic.rows,
      non.palindromic.rows = study$non.palindromic.rows,
      neg.strand.count = study$neg.strand.count,
      not.found.rows = study$not.found.rows.std + study$not.found.rows.alt,
      multiAlleleVariants.rowcount = study$multiAlleleVariants.rowcount,
      HQ.count = study$HQ.count,
      LQ.count = study$LQ.count
    ),
    Correlations = list(
      AFcor.alt_ref = study$AFcor.alt_ref,
      AFcor.alt_ref.indel = study$AFcor.alt_ref.indel,
      AFcor.non.palindromic.alt_ref = study$AFcor.non.palindromic.alt_ref,
      AFcor.non.palindromic.std_ref = study$AFcor.non.palindromic.std_ref,
      AFcor.palindromic.alt_ref = study$AFcor.palindromic.alt_ref,
      AFcor.palindromic.std_ref = study$AFcor.palindromic.std_ref,
      AFcor.std_ref = study$AFcor.std_ref,
      AFcor.std_ref.CHR = study$AFcor.std_ref.CHR,
      AFcor.std_ref.indel = study$AFcor.std_ref.indel,
      PVcor = study$PVcor,
      PVcor.palindromic = study$PVcor.palindromic
    ),
    Statistics = list(
      lambda = study$lambda,
      lambda.gen = study$lambda.gen,
      lambda.imp = study$lambda.imp,
      kurtosis = study$kurtosis,
      kurtosis.HQ = study$kurtosis.HQ,
      Visschers.stat = study$Visschers.stat,
      Visschers.stat.HQ = study$Visschers.stat.HQ,
      skewness = study$skewness,
      skewness.HQ = study$skewness.HQ,
      fixed.callrate = study$fixed.callrate,
      fixed.hwep = study$fixed.hwep,
      fixed.impq = study$fixed.impq,
      fixed.n_total = study$fixed.n_total,
      MAX_N_TOTAL = study$MAX_N_TOTAL,
      MAX_N_CASES = study$MAX_N_CASES,
      hasINDEL = study$hasINDEL
    ),
    starttime = study$starttime,
    endtime = study$endtime,
    Successful_run = TRUE
    )

  return(object)
}
