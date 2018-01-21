globalVariables(c('.', 'rt', 'int', 'product', 'polarity', 'mz'))


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


