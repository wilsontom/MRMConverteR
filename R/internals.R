#' Clean Scan Index
#'
#' @param x a vector of raw scan indexes
#' @keywords internal
#'

clean_index <- function(x)
{
  xa <- gsub('X..SRM.SIC.', '', x)
  xb <- gsub('SRM.SIC.', '', xa)
  return(xb)
}

#' Split Scan Index
#'
#' @param x a vector of cleaned scan indexes
#' @keywords internal
#'

split_index <- function(x)
{
  if (length(x) < 1) {
  }
  xa <- strsplit(x, '\\.')

  xlen <- vapply(xa, length, length(xa))

  if (any(xlen != 4)) {
    stop('some text here')
  }
  xms1 <- lapply(xa, function(x)
    (paste0(x[[1]], '.', x[[2]])))
  xms2 <- lapply(xa, function(x)
    (paste0(x[[3]], '.', x[[4]])))

  xmz <- data.frame(ms1 = unlist(xms1), ms2 = unlist(xms2))

  return(xmz)

}

#' Extract Q1 m/z value
#'
#' @param x a vector of MS1-MS2 pairs
#' @keywords internal
#'


strip_ms1 <- function(x)
{
  xms <- strsplit(as.character(x), '-')

  xlen <- vapply(xms, length, length(xms))

  if (any(xlen != 2)) {
    stop('some text here')
  }

  xms1 <- unlist(lapply(xms, function(x)
    (x[[1]])))

  return(as.numeric(xms1))

}

#' mzR Header names
#'
#' @keywords internal
#'

header_names <- function()
{
  fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                    package = "msdata")
  ms_fl <- mzR::openMSfile(fl, backend = "pwiz")

  hdr <- mzR::header(ms_fl)
  return(names(hdr))
}




