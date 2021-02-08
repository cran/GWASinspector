
drawPlots <- function(processed.data) {

  study <- .QC$thisStudy
  config <- .QC$config


  ##--------------
  # draw plots if it set in config file
  if(config$plot_specs$make_plots){

    ## pvalue correlation plot will not be produced if pvalue column was missing from input file
    ## because it is replaced with calculated pvalues and correlation coefficient will be 1
    if(!study$missing.PVALUE.column)
      # pValueCorPlot.R
      tryCatch({

# 		    plotScatterSmooth.observedP.vs.expectedP(input.data = processed.data,
#                                           plot_cutoff_p = config$plot_specs$plot_cutoff_p,
#                                           PVcor = study$PVcor,
#                                           pvalCorSmPlotPath = study$pvalCorSmPlotPath,
#                                           plot.subtitle = study$plot.title )

        if(study$PVcor > 0.95)
          plot.observedP.vs.expectedP(input.data = processed.data,
                                           plot_cutoff_p = config$plot_specs$plot_cutoff_p,
                                           PVcor = study$PVcor,
                                           pvalCorPlotPath = study$pvalCorPlotPath,
                                           plot.subtitle = study$plot.title )
        else
          plot.observedP.vs.expectedP_dual(input.data = processed.data,
                                      plot_cutoff_p = config$plot_specs$plot_cutoff_p,
                                      PVcor = study$PVcor,
                                      pvalCorPlotPath = study$pvalCorPlotPath,
                                      plot.subtitle = study$plot.title )

          },
               error = function(err)
                 print.and.log(paste('error in pvalue plot.',err$message),'warning',display=.QC$config$debug$verbose)
      )



    # mafPlotFUnction.R
    # only variants that are matched with standard reference are required
    # tryCatch(plot.DataMAF.vs.RefMAF(processed.data[SOURCE == 'Std_ref' & !is.na(AF)],
    #                                 study$stdMafPlotPath,
    #                                 study$AFcor.std_ref,
    #                                 study$AFcor.palindromic.std_ref,
    #                                 study$AFcor.std_ref.indel,
    #                                 paste0(study$plot.title , ' (Standard Reference)' )),
    #          error = function(err)
    #            print.and.log(paste('error in AF plot.',err$message),'warning',display=.QC$config$debug$verbose)
    # )

	tryCatch(plotScatterSmooth.DataMAF.vs.RefMAF(processed.data[SOURCE == 'Std_ref' & !is.na(AF)],
                                    study$stdMafSmPlotPath,
                                    study$AFcor.std_ref,
                                    study$AFcor.palindromic.std_ref,
                                    study$AFcor.std_ref.indel,
                                    study$plot.title ),
             error = function(err)
               print.and.log(paste('error in AF plot.',err$message),'warning',display=.QC$config$debug$verbose)
    )

    # only variants that are matched with alternative reference are required
    if(!is.na(config$supplementaryFiles$allele_ref_alt) & processed.data[!is.na(SOURCE) & SOURCE != 'Std_ref' , .N] > 0)
      tryCatch( plot.DataMAF.vs.RefMAF(processed.data[!is.na(SOURCE) & SOURCE != 'Std_ref'],
                                       study$altMafPlotPath,
                                       study$AFcor.alt_ref,
                                       study$AFcor.palindromic.alt_ref,
                                       study$AFcor.alt_ref.indel,
                                       paste0(study$plot.title, ' (Alternative Reference)')),
                error = function(err)
                  print.and.log(paste('error in AF-Alt plot.',err$message),'warning',display=.QC$config$debug$verbose)
      )

    # histPlotFUnction.R

    tryCatch(plot.variable.frequency.histograms(processed.data,
                                                study,
                                                study$plot.title ),
             error = function(err)
               print.and.log(paste('error in histograms plot.',err$message),'warning',display=.QC$config$debug$verbose)
    )



    ## data contains calculated pvalue; so there is no need to define stderr and beta columns
    ## manh_plot.R
    tryCatch(man.plot(processed.data,
                      chr = 'CHR',
                      pvalue = 'PVALUE',
                      position ='POSITION',
                      fileName = study$manPlotPath,
                      plot.subtitle = study$plot.title,
                      plot.title = 'Manhattan Plot',
                      p.threshold = config$plot_specs$plot_cutoff_p,
                      sig.threshold.log = -log10(5*10^-8),
                      check.columns = FALSE),
             error = function(err)
               print.and.log(paste('error in manhattan plot.',err$message),'warning',display=.QC$config$debug$verbose)
    )


    ## QQplotFunction.R
    tryCatch(QQ_plots(processed.data ,
                      config$plot_specs$plot_cutoff_p,
                      study$plot.title ),
             error = function(err)
               print.and.log(paste('error in QQ plot.',err$message),'warning',display=.QC$config$debug$verbose)
    )

  }
  else
  {
    print.and.log('Plots are skipped!','warning',display=.QC$config$debug$verbose)
  }
}
