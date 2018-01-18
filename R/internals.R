globalVariables(c('.', 'rt', 'int', 'product', 'polarity', 'mz'))


#' @keywords internal
#'

extract_polarity <- function(x)
{
  x <- as.character(x)

  if(length(grep('X..', x)) == 1){
    p <- -1
  }else{
    p <- 1
  }
  if(x == 'TIC'){
    p <- 0
  }
  return(p)
}

#' @keywords internal
extract_precursor <- function(x)
{
  if (x == 'TIC') {
    return(0)
  } else{
    x <- gsub('X..SRM.SIC.|SRM.SIC.', '', x)
    xa <- strsplit(x, '\\.')
    return(as.numeric(paste0(xa[[1]][1], '.', xa[[1]][2])))
  }
}

#' @keywords internal
extract_prodcut <- function(x)
{
  if (x == 'TIC') {
    return(0)
  } else{
    x <- gsub('X..SRM.SIC.|SRM.SIC.', '', x)
    xa <- strsplit(x, '\\.')
    return(as.numeric(paste0(xa[[1]][3], '.', xa[[1]][4])))
  }
}



#' @keywords internal
header_names <- function(x){

	c("seqNum",                 "acquisitionNum",
	  "msLevel",                "polarity",
	  "peaksCount",             "totIonCurrent",
	  "retentionTime",          "basePeakMZ",
	  "basePeakIntensity",      "collisionEnergy",
	  "ionisationEnergy",       "lowMZ",
	  "highMZ",                 "precursorScanNum",
	  "precursorMZ",            "precursorCharge",
	  "precursorIntensity",     "mergedScan",
	  "mergedResultScanNum",    "mergedResultStartScanNum",
	  "mergedResultEndScanNum", "injectionTime",
	  "spectrumId");
}


