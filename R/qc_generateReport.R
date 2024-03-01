#'
#'
#'
#' Used https://github.com/csoneson/alevinQC/tree/devel and https://stackoverflow.com/a/37122513 as inspiration for producing an html QC report
#'
#' @import rmarkdown
#'
#'

qc_generateReport = function(list_of_qc_things, overwrite = FALSE, output_format = NULL, output_file = "QC_Report.Rmd", output_dir="./"){
  
  #Determine the template
  templateFile = system.file("rmd/gimapQCTemplate.Rmd", package="gimap")
  
  if (file.exists(templateFile)) {
    if (file.exists(output_file) & !overwrite){
      stop("there is already an output .Rmd file", output_file,
           ". Please remove or rename this file, or choose another output_file name.", 
           call. = FALSE)
    } else {
      file.copy(from = templateFile, to = output_file, overwrite = overwrite)
    }
    
  } else {
    stop("The Rmd template file ", templateFile, " does not exist", call. = FALSE)
  }
  
  #Process the Arguments
  args                <- list()
  args$input          <- templateFile
  args$output_dir     <- output_dir
  args$output_format  <- output_format
  args$output_file    <- output_file
  
  #Run the render 
  outputFileName = do.call('render', args=args)
  invisible(outputFileName)
}