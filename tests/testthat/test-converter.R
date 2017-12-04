context('MRMConverteR')

test_that('MRMConverteR', {
  example_file_local <-
    system.file('extdata/example_qqq.mzML', package = 'MRMConverteR')

  mrm_raw <- mzR::openMSfile(example_file_local)

  expect_true(isS4(mrm_raw))
  expect_that(length(mrm_raw), equals(0))
  expect_that(mzR::nChrom(mrm_raw), is_more_than(1))

  mrm_chroms <- mzR::chromatograms(mrm_raw)

  expect_true(is.list(mrm_chroms))
  expect_error(mzR::header(mrm_raw))


  chrom_names <-
    unlist(lapply(mrm_chroms, function(x)
      (names(x)[[2]])))[-1]

  expect_true(is.vector(MRMConverteR:::clean_index(chrom_names)))

  chrom_names_clean <- MRMConverteR:::clean_index(chrom_names)
  chrom_split <- MRMConverteR:::split_index(chrom_names_clean)

  expect_true(is.data.frame(chrom_split))
  expect_true(ncol(chrom_split) == 2)
  expect_true(nrow(chrom_split) == length(chrom_names_clean))

  check_size <- chrom_split[, 'ms1'] > chrom_split[, 'ms2']

  expect_true(all(check_size == TRUE))

  mzr_hd <- MRMConverteR:::header_names()

  expect_true(is.vector(mzr_hd))
  expect_true(length(mzr_hd) == length(unique(mzr_hd)))

  fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                    package = "msdata")

  expect_error(convert(f1, ''))

  })
