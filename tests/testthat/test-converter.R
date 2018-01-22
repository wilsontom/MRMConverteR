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

  mrm_conv <- convert_mzml(example_file_local)

  expect_true(is.list(mrm_conv))
  expect_true(length(mrm_conv) == 2)
  expect_true(length(mrm_conv$peaks) == nrow(mrm_conv$header))

  expect_true(is.list(mrm_conv$peaks))
  expect_true(is.data.frame(mrm_conv$header))



  })
