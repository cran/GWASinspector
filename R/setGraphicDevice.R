set_graphic_device <- function(config)
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
    set_graphic_device_png()
  }
  else if(capabilities('jpeg') & config$plot_specs$graphic_device == 'jpeg')
  {
    set_graphic_device_jpeg()
  }
  else if(capabilities('tiff') & config$plot_specs$graphic_device == 'tiff')
  {
    set_graphic_device_tiff()
  }
  else if(capabilities('png'))
  {
    set_graphic_device_png()
  }
  else if(capabilities('jpeg'))
  {
    set_graphic_device_jpeg()
  }
  else if(capabilities('tiff'))
  {
    set_graphic_device_tiff()
  }

  config$graphic.device <- TRUE
  return(config)
}


set_graphic_device_png <- function()
{
  .QC$graphic.device ='png'
  .QC$img.extension = '.png'
}

set_graphic_device_jpeg <- function()
{
  .QC$graphic.device ='jpeg'
  .QC$img.extension = '.jpeg'
}

set_graphic_device_tiff <- function()
{
  .QC$graphic.device ='tiff'
  .QC$img.extension = '.tiff'
}
