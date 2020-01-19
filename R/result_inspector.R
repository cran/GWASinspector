#' Display the reports of running the Inspector algorithm
#'
#' This function displays a brief report about the results of running the Inspector algorithm on a list of GWAS result files.
#' The full report including plots, cleaned files and summary statistics are generated and saved in the output folder during the algorithm run.
#'
#' @param inspector An instance of \linkS4class{Inspector} class. Check \code{\link{setup.inspector}} for more details.
#' @return A data.table containing a brief report about the results.
#'
result.inspector <- function(inspector)
{

  if(missing(inspector))
    stop('Function arguments are not set.',call. = FALSE)

  if (!is(inspector, "Inspector"))
    stop("Object must be of class Inspector", call. = FALSE)

  if (inspector@StudyList@studyCount == 0)
    stop("Inspector is not run yet.", call. = FALSE)



  dt <- data.table(`#` = numeric(0),
                   Name = character(0),
                   Result = logical(0),
                   Rows = numeric(0),
                   "Output rows" = numeric(0),
                   "Found in Ref" = numeric(0),
                   "AF correlation (std)" = numeric(0),
                   "P-value correlation" = numeric(0),
                   "Lambda" = numeric(0))

  for (i in 1:inspector@StudyList@studyCount)
  {
    study <- inspector@StudyList@studyList[[i]]

    success <- FALSE
    if(inspector@StudyList@studyList[[i]]@Successful_run)
      success <- TRUE

    dt <- rbind(dt, data.table(`#` = i,
                               Name = basename(study@File$file.path),
                               Result = success,
                Rows = study@Counts$input.data.rowcount,
                "Output rows" = study@Counts$rowcount.step3,
                "Found in Ref" = study@Counts$found.rows,
                "AF correlation (std)" =study@Correlations$AFcor.std_ref,
                "P-value correlation" = study@Correlations$PVcor,
                "Lambda" = study@Statistics$lambda
                )
    )
  }


  print(knitr::kable(dt,format = 'rst',row.names = FALSE,align = 'l',caption = 'f'))

  return(invisible(dt))
}
