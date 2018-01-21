#' Save Converted mzML
#'
#' @param peaks a list of \code{peaks} and \code{header} which has been created using the \code{convert_mzml} function
#' @param filename a character string of the filename to save converted files to.
#'
#' @export
#' @author Tom Wilson \email{tpw2@@aber.ac.uk}

save_mzml <- function(peaks,filename)
{
    filename <- paste0(filename, '.mzML')

    mzR::writeMSData(
    object = peaks$peaks,
    file = filename,
    header = peaks$header,
    outformat = 'mzml'
  )

return(invisible(NULL))

}
