effect_size_ref_make <- function(dataset, db_path)
{
  if ('PVALUE' %in% names(dataset))
  {
    savedFile <- tryCatch({
      dataset <-  dataset[PVALUE <= 0.001, list(hID, EFFECT_ALL, OTHER_ALL, EFFECT, PVALUE)]
      names(dataset) <-  c("hID", "EFF_ALL", "NON_EFF_ALL", "EFFECT", "PVALUE")

      saveRDS(object = dataset,
              file = db_path,
              version = '2')

      print.and.log(
        sprintf('Effect-size reference file saved with %s rows.', nrow(dataset)),
        'info',
        display = .QC$config$debug$verbose
      )
    }, error = function(err)
    {
      print.and.log(paste( 'Error occured when saving the effect-size reference!', err$message),
        'warning')
    })

  } else{
    print.and.log('PVALUE not found in dataset to save as effect-size reference file!','warning',
                  display = .QC$config$debug$verbose)
  }
}
