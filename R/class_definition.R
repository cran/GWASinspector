
methods::setOldClass("proc_time")


#' An S4 class to represent a list of inspected GWAS study result files.
#'
#' This class is embedded in the \linkS4class{Inspector} class and should not be initiated separately.
#'
#' @slot studyList A list of GWAS study result files. Each member of this list is of class \code{\linkS4class{Study}}.
#' @slot studyCount A numeric value indicating how many items of class \code{\linkS4class{Study}} are included.
#' @docType class
StudyList <- setClass(
  "StudyList",
  slots = list(
    studyList = "list",
    studyCount = "numeric"
  ),
  prototype = list(
    studyList = list(),
    studyCount = 0
  )
)


#' An S4 class to represent the methods and parameters for inspecting a list of GWAS study result files.
#'
#' An object of this class is created by \code{\link{setup.inspector}} function. Each section of the
#' configuration file is represented as a list of attributes in this object.
#'
#' @slot paths A list of parameters which indicate \strong{Paths} section from configuration file.
#' @slot supplementaryFiles A list of parameters which indicate \strong{supplementaryFiles} section from configuration file.
#' @slot input_parameters A list of parameters which indicate \strong{input_parameters} section from configuration file.
#' @slot output_parameters A list of parameters which indicate \strong{output_parameters} section from configuration file.
#' @slot remove_chromosomes A list of parameters which indicate \strong{remove_chromosomes} section from configuration file.
#' @slot plot_specs A list of parameters which indicate \strong{plot_specs} section from configuration file.
#' @slot filters A list of parameters which indicate \strong{filters} section from configuration file.
#' @slot debug A list of parameters which indicate \strong{debug} section from configuration file.
#' @slot input_files A list of files that will be inspected during the run.
#' @slot created_at The time that object was created.
#' @slot start_time The time that object was run.
#' @slot end_time The time that run was finished.
#' @slot StudyList An object of \linkS4class{StudyList} class.
#' @docType class
Inspector <- setClass(
  "Inspector",
  slots = list(
    ## config file
    paths = "list",
    supplementaryFiles = "list",
    input_parameters = "list",
    output_parameters = "list",
    remove_chromosomes = "list",
    plot_specs = "list",
    filters = "list",
    debug="list",
    ## Additional items,
    input_files = "vector",
    created_at = "POSIXt",
    start_time = "POSIXt",
    end_time = "POSIXt",
    StudyList= "StudyList"
  ),
  prototype = list(
    paths = list(
      filename = ".+",
      filename_output_tag = "QC",
      dir_data = "D:/",
      dir_output = "D:/output",
      dir_references = "D:/",
      input_files = c()
    ),
    supplementaryFiles = list(
      header_translations = "alt_headers.txt",
      allele_ref_std = "",
      allele_ref_std_population = "",
      allele_ref_alt = "",
      beta_ref_std = ""
    ),
    input_parameters = list(
      effect_type = "character",
      column_separator = "\t",
      na.strings = c(),
      imputed_T = "",
      imputed_F = "",
      calculate_missing_p = FALSE
    ),
    output_parameters = list(
      save_final_dataset = FALSE,
      gzip_final_dataset = TRUE,
      out_header = "standard",
      out_sep = "\t",
      out_na = NA,
      out_dec = ".",
      html_report = TRUE,
      object_file = TRUE
    ),
    remove_chromosomes = list(
      remove_X = FALSE,
      remove_Y = FALSE,
      remove_XY = FALSE,
      remove_M = FALSE
    ),
    plot_specs = list(
      make_plots = TRUE,
      plot_cutoff_p = 0.01,
      graphic_device = "",
      plot_title = ""
    ),
    filters = list(
      HQfilter_FRQ = 0.01,
      HQfilter_HWE = 1e-6,
      HQfilter_cal = 0.95,
      HQfilter_imp = 0.3,
      threshold_diffEAF = 0.15,
      minimal_impQ_value = -0.5,
      maximal_impQ_value = 1.5
    ),
    debug=list(
      verbose = FALSE,
      save_pre_modification_file = FALSE,
      reduced.AF.plot = TRUE,
      test_row_count = 1000
    ),
    created_at = Sys.time(),
    start_time = Sys.time(),
    end_time = Sys.time(),
    StudyList = new("StudyList")
  )
)


#' An S4 class to represent an inspected GWAS study result file.
#'
#' This class is embedded in the \linkS4class{StudyList} class and should not be initiated separately.
#'
#' @slot File A list representing GWAS result file specifications
#' @slot Counts A list representing different variant counts from the GWAS result file.
#' @slot Correlations A list representing different allele frequency and P-value correlations in the GWAS result file.
#' @slot Statistics A list representing summary statistics from the GWAS result file.
#' @slot Successful_run A logical value indicating whether the run was successful or not.
#' @slot starttime The time that file inspection started.
#' @slot endtime The time that file inspection ended.
#' @docType class
Study <- setClass(
  "Study",
  slots = list(
    starttime = "POSIXt",
    endtime = "POSIXt",
    File = "list",
    Counts = "list",
    Correlations = "list",
    Statistics = "list",
    Successful_run = "logical"
  ),
  prototype = list(
    File = list(
      file.path = character(0),
      file.name = character(0),
      file.extension = character(0),
      file.line.count = numeric(0),
      dup_lines_count = numeric(0),
      original.File.Columns = c(character(0)),
      renamed.File.Columns.sorted = c(character(0))
    ),
    Counts = list(
      input.data.rowcount = numeric(0),
      duplicate.count = numeric(0),
      rowcount.step1 = numeric(0),
      rowcount.step2 = numeric(0),
      rowcount.step3 = numeric(0),
      found.rows = numeric(0),
      mismatched.rows = numeric(0),
      ambiguos.rows = numeric(0),
      switched.rows = numeric(0),
      flipped.rows = numeric(0),
      monomorphic.count = numeric(0),
      palindromic.rows = numeric(0),
      non.palindromic.rows = numeric(0),
      neg.strand.count = numeric(0),
      not.found.rows = numeric(0),
      multiAlleleVariants.rowcount = numeric(0),
      HQ.count = numeric(0),
      LQ.count = numeric(0)
    ),
    Correlations = list(
      AFcor.alt_ref = numeric(0),
      AFcor.alt_ref.indel = numeric(0),
      AFcor.non.palindromic.alt_ref = numeric(0),
      AFcor.non.palindromic.std_ref = numeric(0),
      AFcor.palindromic.alt_ref = numeric(0),
      AFcor.palindromic.std_ref = numeric(0),
      AFcor.std_ref = numeric(0),
      AFcor.std_ref.CHR = numeric(0),
      AFcor.std_ref.indel = numeric(0),
      PVcor = numeric(0),
      PVcor.palindromic = numeric(0)
    ),
    Statistics = list(
      lambda = numeric(0),
      lambda.gen = numeric(0),
      lambda.imp = numeric(0),
      kurtosis = numeric(0),
      kurtosis.HQ = numeric(0),
      Visschers.stat = numeric(0),
      Visschers.stat.HQ = numeric(0),
      skewness = numeric(0),
      skewness.HQ = numeric(0),
      fixed.callrate = numeric(0),
      fixed.hwep = numeric(0),
      fixed.impq = numeric(0),
      fixed.n_total = numeric(0),
      MAX_N_TOTAL = numeric(0),
      N.max = numeric(0),
      hasINDEL = FALSE
    ),
    starttime = Sys.time(),
    endtime =  Sys.time(),
    Successful_run = FALSE
  )
)




#####################################################

setMethod(
  "show",
  "Inspector",
  function(object) {
    if (!is(object, "Inspector")) {
      stop("Object must be of class Inspector", call. = FALSE)
    }

    cat("\nInspector object specifications:\n\n")
    cat("   Paths")
    cat("\n\tFile pattern:\t", object@paths$filename)
    cat("\n\tOutput tag:\t", object@paths$filename_output_tag)
    cat("\n\tData directory:\t", object@paths$dir_data)
    cat("\n\tOutput directory:\t", object@paths$dir_output)
    cat("\n\tReference directory:\t", object@paths$dir_references)



    cat("\n\n   Supplementary Files")
    cat("\n\tHeader translation:\t", basename(object@supplementaryFiles$header_translations))
    cat("\n\tStandard allele reference:\t", basename(object@supplementaryFiles$allele_ref_std))
    cat("\n\tSelected population:\t", object@supplementaryFiles$allele_ref_std_population)

    if (!is.na(object@supplementaryFiles$allele_ref_alt)) {
      cat("\n\talternate allele reference:\t", basename(object@supplementaryFiles$allele_ref_alt))
    }
    if (!is.na(object@supplementaryFiles$beta_ref_std)) {
      cat("\n\tBeta reference:\t", basename(object@supplementaryFiles$beta_ref_std))
    }


    cat("\n\n   Input parameters")
    cat("\n\tcolumn separator:\t", gsub(x = object@input_parameters$column_separator, pattern = "	", replacement = "\t", fixed = T))
    cat("\n\tNA string:\t", object@input_parameters$na.string)
    cat("\n\tImputed string:\t", object@input_parameters$imputed_T)
    cat("\n\tNon-imputed string:\t", object@input_parameters$imputed_F)
    cat("\n\tEffect-type:\t", object@input_parameters$effect_type)


    cat("\n\n   Output parameters")
    cat("\n\tSave an output clean file:\t", object@output_parameters$save_final_dataset)
    cat("\n\tCompressing the result file :\t", object@output_parameters$gzip_final_dataset)
    cat("\n\tResult file header:\t", object@output_parameters$out_header)
    cat("\n\tColumn separator:\t", gsub(x = "\t", pattern = "	", replacement = "\t", fixed = T))
    cat("\n\tDecimal value:\t", object@output_parameters$out_dec)
    cat("\n\tNA string:\t", object@output_parameters$out_na)
    cat("\n\tSaving Html report file:\t", object@output_parameters$html_report)
    cat("\n\tSaving study object file:\t", object@output_parameters$object_file)



    cat("\n\n   Remove chromosome variants")
    cat("\n\tX:\t", object@remove_chromosomes$remove_X)
    cat("\n\tY:\t", object@remove_chromosomes$remove_Y)
    cat("\n\tXY:\t", object@remove_chromosomes$remove_XY)
    cat("\n\tM:\t", object@remove_chromosomes$remove_M)



    cat("\n\n   Plot specifications")
    cat("\n\tSave plots:\t", object@plot_specs$make_plots)
    cat("\n\tGraphic device:\t", object@plot_specs$graphic_device)
    cat("\n\tPlot title:\t", object@plot_specs$plot_title)
    cat("\n\tP-value cut off:\t", object@plot_specs$plot_cutoff_p)


    cat("\n\n   HQ variant filters")
    cat("\n\tSave plots:\t", object@filters$HQfilter_FRQ)
    cat("\n\tGraphic device:\t", object@filters$HQfilter_HWE)
    cat("\n\tPlot title:\t", object@filters$HQfilter_cal)
    cat("\n\tP-value cut off:\t", object@filters$HQfilter_imp)



    cat("\n\n   Filters")
    cat("\n\tAF threshold for HQ variants:\t", object@filters$HQfilter_FRQ)
    cat("\n\tHWE P_value threshold for HQ variants:\t", object@filters$HQfilter_HWE)
    cat("\n\tCall rate threshold for HQ variants:\t", object@filters$HQfilter_cal)
    cat("\n\tImp. qual. threshold for HQ variants:\t", object@filters$HQfilter_imp)
    cat("\n\tthreshold for the difference between reported and reference AF:\t", object@filters$threshold_diffEAF)
    cat("\n\tMinimal possible Imp. qual. :\t", object@filters$minimal_impQ_value)
    cat("\n\tMaximal possible Imp. qual. :\t", object@filters$maximal_impQ_value)


    cat("\n\n   Input files")
    if(length(object@input_files) > 0)
    {
      for(i in 1:length(object@input_files))
         cat("\n\t" , i , "-" ,object@input_files[[i]])
    }else
      cat("\n\tempty")


    cat("\n\n   Results")
    if(object@StudyList@studyCount > 0)
    {
      cat(sprintf("\n\t%s report objects are included. use \"result.inspector()\" function for detail.", object@StudyList@studyCount))
      cat(sprintf("\n\tLast run: %s", object@start_time))
    }else
      cat("\n\tQC has not been performed yet.")


    cat("\n\n")
  }
)


setMethod(
  "show",
  "Study",
  function(object) {
    if (!is(object, "Study")) {
      stop("Object must be of class Study", call. = FALSE)
    }

    cat("This is a \"Study\" class object from GWASinspector package. Check the manual for more information.")
    cat("\nUse @ sign for a specific attribute.")
  }
)


setMethod(
  "show",
  "StudyList",
  function(object) {
    if (!is(object, "StudyList")) {
      stop("Object must be of class StudyList", call. = FALSE)
    }
  }
)
