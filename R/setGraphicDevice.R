set.graphic.device <- function(config)
{
  # return if no capability exists
  if( !capabilities('jpeg') & !capabilities('png') & !capabilities('tiff'))
  {
    config$plot_specs$make_plots <- FALSE
    config$graphic.device <- FALSE
    return(config)
  }



  # check if selected device is available
  if(capabilities('png') & config$plot_specs$graphic_device == 'png')
  {
    set.graphic.device_png()
  }
  else if(capabilities('jpeg') & config$plot_specs$graphic_device == 'jpeg')
  {
    set.graphic.device_jpeg()
  }
  else if(capabilities('tiff') & config$plot_specs$graphic_device == 'tiff')
  {
    set.graphic.device_tiff()
  }
  else if(capabilities('png'))
  {
    set.graphic.device_png()
  }
  else if(capabilities('jpeg'))
  {
    set.graphic.device_jpeg()
  }
  else if(capabilities('tiff'))
  {
    set.graphic.device_tiff()
  }

  config$graphic.device <- TRUE
  return(config)
}


set.graphic.device_png <- function()
{
  .QC$graphic.device ='png'
  .QC$img.extension = '.png'
}

set.graphic.device_jpeg <- function()
{
  .QC$graphic.device ='jpeg'
  .QC$img.extension = '.jpeg'
}

set.graphic.device_tiff <- function()
{
  .QC$graphic.device ='tiff'
  .QC$img.extension = '.tiff'
}
