#' Convert MRM-MS files to LC-MS style files
#'
#' Convert MRM-MS \code{.mzML} files which contain a series of time-intensity pairs (chromatograms) to a series of LC-MS scan style \code{m/z}-intensity pairs.
#'
#' Once the data is represented as a matrix of \code{m/z} and intensity for each scan; the data can be further processed by open source tools such as \code{xcms}.
#'
#' @param input a character stirng of the absolute file path of a MRM-MS \code{.mzML} file
#' @param outpath a character string of the location to save converted files to
#'
#' @export
#' @importFrom stats aggregate
#' @author Tom Wilson \email{tpw2@@aber.ac.uk}


convert <- function(input, outpath)
{
  if (tools::file_ext(input) != 'mzML') {
    stop(deparse(substitute(input)), ' is not a .mzML file', call. = FALSE)
  }

  open_raw <- mzR::openMSfile(input, backend = 'pwiz')

  if (length(open_raw) != 0) {
    warning('length is > 0')
  }

  chrom_raw <- mzR::chromatograms(open_raw)

  ind_name <- lapply(chrom_raw, function(x)
    (names(x)))

  if (ind_name[[1]][2] == 'TIC') {
    chrom_raw[[1]] <- NULL
  }


  for (i in seq_along(chrom_raw)) {
    chrom_raw[[i]] <-
      data.frame(chrom_raw[[i]], mz = rep(names(chrom_raw[[i]])[2]))
  }


  mz_raw <- lapply(chrom_raw, function(x)
    (clean_index(x$mz)))

  mz_split <- lapply(mz_raw, split_index)

  chrom_clean <- NULL
  for (i in seq_along(chrom_raw)) {
    chrom_clean[[i]] <-
      data.frame(
        rt = chrom_raw[[i]]$time,
        mz = mz_split[[i]]$ms1,
        int = chrom_raw[[i]][, 2],
        mzagg = paste0(mz_split[[i]]$ms1, '-', mz_split[[i]]$ms2)
      )
  }

  # assume a cycle time of 1.0 Seconds

  for (i in seq_along(chrom_clean)) {
    chrom_clean[[i]]$rt <- round(chrom_clean[[i]]$rt, digits = 1)
  }

  # scan_time <- as.numeric(names(chrom_split))

  chrom_all <- do.call('rbind', chrom_clean)

  chrom_split <- split(chrom_all, chrom_all$rt)

  scan_time <- as.numeric(names(chrom_split))

  # sum int across agg peak
  chrom_agg <- NULL
  for (i in seq_along(chrom_split)) {
    chrom_agg[[i]] <-
      aggregate(chrom_split[[i]]$int, by = list(chrom_split[[i]]$mzagg), sum)
    names(chrom_agg[[i]]) <- c('mz', 'int')
  }


  for (i in seq_along(chrom_agg)) {
    chrom_agg[[i]]$mz <- strip_ms2(chrom_agg[[i]]$mz)
  }



  chrom_agg <- lapply(chrom_agg, function(x)
    (x[order(x$mz), ]))


  hdnm <- header_names()

  hd_tmp <-
    data.frame(matrix(nrow = length(chrom_agg), ncol = length(hdnm)))
  names(hd_tmp) <- hdnm

  hd_tmp[, 'seqNum'] <-
    hd_tmp[, 'acquisitionNum'] <-
    seq(from = 1,
        to = nrow(hd_tmp),
        by = 1)
  hd_tmp[, 'msLevel'] <- rep(1)
  hd_tmp[, 'polarity'] <- rep(1)

  peak_count <- vapply(chrom_agg, nrow, length(chrom_agg))

  hd_tmp[, 'peaksCount'] <- peak_count

  tic <- unlist(lapply(chrom_agg, function(x)
    (sum(x$int))))
  bpi <- unlist(lapply(chrom_agg, function(x)
    (max(x$int))))
  bpmz <-
    unlist(lapply(chrom_agg, function(x)
      (x$mz[which(x$int == max(x$int))])))

  hd_tmp[, 'totIonCurrent'] <- as.numeric(tic)
  hd_tmp[, 'basePeakIntensity'] <- as.numeric(bpi)
  hd_tmp[, 'basePeakMZ'] <- as.numeric(bpmz)

  hd_tmp[, 'retentionTime'] <- scan_time

  hd_tmp[is.na(hd_tmp)] <- 0

  if ('filterString' %in% names(hd_tmp)) {
    hd_tmp[, 'filterString'] <- as.character(hd_tmp[, 'filterString'])
  }
  destfile <- paste0(outpath, '/', 'convert_', basename(input))
  destfile <- gsub('.mzML', '.mzML', destfile)

  for (i in seq_along(chrom_agg)) {
    chrom_agg[[i]] <- as.matrix(chrom_agg[[i]])
  }

  mzR::writeMSData(
    data = chrom_agg,
    filename = destfile,
    header = hd_tmp,
    outformat = 'mzml'
  )

}
