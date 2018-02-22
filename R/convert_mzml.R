#' Convert MRM-MS files to LC-MS style files
#'
#' Data aquired using Liquid Chromatography Triple Quadrupole Mass Spectromety (LC-QqQ-MS) operating in Selective Reaction Monitoiring (SRM) or Multiple Reaction Monitoring (MRM) mode
#' is stored in `mzML` files as a series of chromatrograms (retention time vs intensity). `convert_mzml` converts `mzML` chromatogram data into spectra; which is how data acquired using
#' LC-MS is stored. Once SRM/MRM data is represented as a series of spectra (_m/z_ vs intensity) then tools such as `xcms` can be used to process and analyse the data. 
#'
#' @param mzML a character string of the absolute file path of a SRM/MRM-MS `.mzML` file
#' @return a list of two elements
#'   - __peaks__ a list of spectra 
#'   - __header__ a `data.frame` of spectra header data
#'
#' @export
#' @importFrom stats aggregate
#' @importFrom dplyr %>% select bind_rows group_by summarise ungroup arrange rename tibble mutate
#' @author Tom Wilson \email{tpw2@@aber.ac.uk}
#' @examples
#' library(msdata)
#' librar(mzR)
#' 
#'
#'
#'
#'
#'
#'
#'



convert_mzml <- function(mzML)
{
  if (tools::file_ext(mzML) != 'mzML') {
    stop(deparse(substitute(mzML)), ' is not a .mzML file', call. = FALSE)
  }

  open_raw <- mzR::openMSfile(mzML, backend = 'pwiz')

  if (length(open_raw) != 0) {
    warning('length is > 0')
  }

  chrom_raw <- mzR::chromatograms(open_raw)


  chrom_header <- mzR::chromatogramHeader(open_raw)

  filter_name <- chrom_header[,'chromatogramId']


  filter_index <-
    tibble(
      polarity = chrom_header[,'polarity'],
      precursor = chrom_header[,'precursorIsolationWindowTargetMZ'],
      product = chrom_header[,'productIsolationWindowTargetMZ'],
      chromidx = chrom_header[,'chromatogramIndex']
    )


  idn_out <-
    which(filter_index[, 'polarity'] == -1 &
            is.na(filter_index[, 'precursor']) &
            is.na(filter_index[, 'product']))

  if (length(idn_out) > 0) {
    filter_index <- filter_index[-idn_out,]
    chrom_raw[[idn_out]] <- NULL
  }

  chrom_tmp <- NULL
  for (i in seq_along(chrom_raw)) {
    chrom_tmp[[i]] <-
      data.frame(
        rt = chrom_raw[[i]]$time,
        int = chrom_raw[[i]][, 2],
        polarity = rep(filter_index[i, 'polarity']),
        precursor = rep(filter_index[i, 'precursor']),
        product = rep(filter_index[i, 'product'])
      )
  }

  chrom_all <- chrom_tmp %>% bind_rows() %>% rename(., mz = product)
  chrom_all[, 'rt'] <- round(chrom_all[, 'rt'], digits = 1)

  chrom_agg <-
    chrom_all %>% group_by(., polarity, rt, mz) %>% summarise(., int = sum(int)) %>% ungroup(.) %>%
    mutate(., spf = paste0(.$polarity, ':', .$rt)) %>% arrange(., rt)

  uuid <-
    data.frame(spf = unique(chrom_agg[, 'spf']),
               id = seq(from = 1, to = length(unique(chrom_agg[, 'spf']$spf))))


  midx <- match(chrom_agg$spf, uuid[, 'spf'])
  chrom_agg <- chrom_agg %>% mutate(., spid = uuid$id[midx])


  chrom_split <- split(chrom_agg, ordered(chrom_agg$spid))

  chrom_pol <-
    purrr::map(chrom_split, ~ {
      unique(.$polarity)
    }) %>% unlist(.)

  chrom_rt <-  purrr::map(chrom_split, ~ {
    unique(.$rt)
  }) %>% unlist(.)
  peaks_count <- purrr::map(chrom_split, ~ {
    nrow(.)
  }) %>% unlist(.)

  bpi <- purrr::map(chrom_split, ~ {
    max(.$int)
  }) %>% unlist(.)
  tic <- purrr::map(chrom_split, ~ {
    sum(.$int)
  }) %>% unlist(.)
  bpmz <-
    purrr::map(chrom_split, ~ {
      .$mz[which(.$int == max(.$int))]
    }) %>% unlist(.)


  # # this is too slow
  peaks <-
    purrr::map(chrom_split, ~ {
      select(., -c(polarity, rt, spf, spid)) %>% as.matrix(.)
    })

  names(peaks) <- NULL

  hdnm <- header_names()

  hd_tmp <-
    data.frame(matrix(nrow = length(peaks), ncol = length(hdnm)))
  names(hd_tmp) <- hdnm

  hd_tmp[, 'acquisitionNum'] <-
    seq(from = 1,
        to = nrow(hd_tmp),
        by = 1)

  hd_tmp[, 'seqNum'] <- hd_tmp[, 'acquisitionNum']

  hd_tmp[, 'msLevel'] <- rep(1, nrow(hd_tmp))
  hd_tmp[, 'polarity'] <- chrom_pol
  hd_tmp[, 'retentionTime'] <- chrom_rt
  hd_tmp[, 'peaksCount'] <- peaks_count
  hd_tmp[, 'totIonCurrent'] <- as.numeric(tic)
  hd_tmp[, 'basePeakIntensity'] <- as.numeric(bpi)
  hd_tmp[, 'basePeakMZ'] <- as.numeric(bpmz)

  hd_tmp[is.na(hd_tmp)] <- 0

  if ('filterString' %in% names(hd_tmp)) {
    hd_tmp[, 'filterString'] <- as.character(hd_tmp[, 'filterString'])
  }

  return(list(peaks = peaks, header = hd_tmp))

}
