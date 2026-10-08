test_that("ageing_error_matrix_to_fims() formats probabilities for FIMS", {
  probabilities <- matrix(
    c(0.8, 0.2, 0, 0.1, 0.7, 0.2, 0, 0.3, 0.7),
    nrow = 3,
    byrow = TRUE
  )

  result <- ageing_error_matrix_to_fims(probabilities, ages = 0:2)

  expect_named(
    result,
    c("type", "fleet", "age", "length", "timing", "observed", "unit", "uncertainty")
  )
  expect_equal(result[["type"]], rep("ageing_error", 9))
  expect_true(all(is.na(result[["fleet"]])))
  expect_equal(result[["age"]], rep(0:2, times = 3))
  expect_true(all(is.na(result[["length"]])))
  expect_true(all(is.na(result[["timing"]])))
  expect_equal(result[["observed"]], as.vector(t(probabilities)))
  expect_equal(result[["unit"]], rep("proportion", 9))
  expect_equal(result[["uncertainty"]], as.character(rep(0:2, each = 3)))
})

test_that("plot_output() returns one FIMS matrix per reader", {
  probabilities <- array(0, dim = c(2, 3, 3))
  probabilities[1, , ] <- matrix(
    c(0.8, 0.2, 0, 0.1, 0.7, 0.2, 0, 0.3, 0.7),
    nrow = 3,
    byrow = TRUE
  )
  probabilities[2, , ] <- diag(3)
  report <- list(
    AgeErrOut = probabilities,
    Aprob = matrix(c(0.2, 0.5, 0.3), nrow = 1),
    TheSD = matrix(c(0, 0.2, 0.3, 0, 0.2, 0.3), nrow = 2, byrow = TRUE),
    TheBias = matrix(0, nrow = 2, ncol = 3)
  )

  result <- plot_output(
    Data = matrix(c(1, 1, 2), nrow = 1),
    IDataSet = 1,
    MaxAge = 2,
    Report = report,
    subplot = integer(0),
    SaveDir = tempdir()
  )

  expect_named(result[["ageing_error_fims"]], c("reader1", "reader2"))
  expect_equal(
    result[["ageing_error_fims"]][["reader1"]][["observed"]],
    as.vector(t(probabilities[1, , ]))
  )
  expect_equal(
    result[["ageing_error_fims"]][["reader2"]][["observed"]],
    as.vector(t(probabilities[2, , ]))
  )
})
