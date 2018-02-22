#' Save Converted mzML
#'
#' @param input a list of `peaks` and `header` which has been created using the `convert_mzml` function
#' @param filename a character string of the filename to save converted files to.
#'
#' @export
#' @author Tom Wilson \email{tpw2@@aber.ac.uk}

save_mzml <- function(input,filename)
{
    filename <- paste0(filename, '.mzML')

    mzR::writeMSData(
    object = input$peaks,
    file = filename,
    header = input$header,
    outformat = 'mzml'
  )

return(invisible(NULL))

}
