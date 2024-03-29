% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_definition.R
\docType{class}
\name{Inspector-class}
\alias{Inspector-class}
\alias{Inspector}
\title{An S4 class to represent the methods and parameters for inspecting a list of GWAS study result files.}
\description{
An object of this class is created by \code{\link{setup_inspector}} function. Each section of the
configuration file is represented as a list of attributes in this object.
}
\section{Slots}{

\describe{
\item{\code{paths}}{A list of parameters which indicate \strong{Paths} section from configuration file.}

\item{\code{supplementaryFiles}}{A list of parameters which indicate \strong{supplementaryFiles} section from configuration file.}

\item{\code{input_parameters}}{A list of parameters which indicate \strong{input_parameters} section from configuration file.}

\item{\code{output_parameters}}{A list of parameters which indicate \strong{output_parameters} section from configuration file.}

\item{\code{remove_chromosomes}}{A list of parameters which indicate \strong{remove_chromosomes} section from configuration file.}

\item{\code{plot_specs}}{A list of parameters which indicate \strong{plot_specs} section from configuration file.}

\item{\code{filters}}{A list of parameters which indicate \strong{filters} section from configuration file.}

\item{\code{debug}}{A list of parameters which indicate \strong{debug} section from configuration file.}

\item{\code{input_files}}{A list of files that will be inspected during the run.}

\item{\code{created_at}}{The time that object was created.}

\item{\code{start_time}}{The time that object was run.}

\item{\code{end_time}}{The time that run was finished.}

\item{\code{StudyList}}{An object of \linkS4class{StudyList} class.}
}}

