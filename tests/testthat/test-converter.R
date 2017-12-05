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

  xa <- names(mrm_chroms[[10]])[2]
  xb <- names(mrm_chroms[[1]])[2]
  expect_true(is.numeric(MRMConverteR:::extract_polarity(xa)))
  expect_true(is.numeric(MRMConverteR:::extract_precursor(xa)))
  expect_true(is.numeric(MRMConverteR:::extract_prodcut(xa)))

  expect_that(MRMConverteR:::extract_polarity(xb), equals(0))
  expect_that(MRMConverteR:::extract_precursor(xb), equals(0))
  expect_that(MRMConverteR:::extract_prodcut(xb), equals(0))

  mzr_hd <- MRMConverteR:::header_names()

  expect_true(is.vector(mzr_hd))
  expect_true(length(mzr_hd) == length(unique(mzr_hd)))

  fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                    package = "msdata")

  expect_error(convert(f1, ''))

  ps <- convert(example_file_local, NULL, return = TRUE)

  expect_true(is.list(ps))

  psdim <- vapply(ps, ncol, length(ps))
  expect_true(all(psdim == 2))

  })
