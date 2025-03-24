# create initial data set created with the commands

# set.seed(123)
# data_test <- data.frame(reader1 = sample(5:10, 20, replace = TRUE)) |>
#   dplyr::mutate(reader2 = round(reader1 * rnorm(20, 1, 0.05)),
#                 reader3 = round(reader1 * rnorm(20, 1, 0.05)))
# data_test$reader2[18:20] <- NA
# data_test$reader3[15:17] <- NA

data_test <- data.frame(
  reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
  reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
  reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
)
temp_dir <- file.path(tempdir(), "test_AgeingError")
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

# testing tally_repeats()
test_that("Can tally repeats using tally_repeats()", {
  data_tallied <- data_test |> tally_repeats()
  testthat::expect_true(all(data_tallied[["counts"]] == c(1, 2, 3, 3, 2, 1, 1, 1, 1, 1, 1, 3)))
})

# testing create_data_file()
test_that("Can create a data file using create_data_file()", {
  data_tallied <- data_test |> tally_repeats()
  data_file <- create_data_file(data_tallied, dir = temp_dir, file_name = "test.dat")

  testthat::expect_true(file.exists(data_file))
  data_file_read <- readLines(data_file)
  testthat::expect_true(data_file_read[1] == "Range_of_ages")
  testthat::expect_true(data_file_read[6] == "3 # readers")
  testthat::expect_true(data_file_read[20] == "3 10 10 10")
  testthat::expect_true(length(data_file_read) == 20)
})

test_that("Can create a data file create_data_file() without tally_repeats()", {
  temp_dir <- file.path(tempdir(), "test_AgeingError")
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)
  }
  data_file <- create_data_file(data_test, dir = temp_dir, file_name = "test2.dat")

  testthat::expect_true(file.exists(data_file))
  data_file_read <- readLines(data_file)
  testthat::expect_true(data_file_read[1] == "Range_of_ages")
  testthat::expect_true(data_file_read[6] == "3 # readers")
  testthat::expect_true(data_file_read[20] == "3 10 10 10")
  testthat::expect_true(length(data_file_read) == 20)
})

# testing load_data()
test_that("Can load a data file using load_data()", {
  data_loaded <- load_data(
    DataFile = file.path(temp_dir, "test.dat"),
    NDataSet = 1,
    verbose = TRUE,
    EchoFile = file.path(temp_dir, "test_echo.out")
  )

  testthat::expect_true(length(data_loaded) == 16)
  testthat::expect_true(data_loaded$NReaders == 3)
  data_tallied <- data_test |> tally_repeats()
  testthat::expect_true(all(as.matrix(data_tallied) - data_loaded$TheData[1, , ] == 0))
})

test_that("Can load example data sets using load_data()", {
  example_path <- system.file("extdata", package = "AgeingError")
  data_files <- list.files(example_path, pattern = "dat", full.names = TRUE)
  for (i in 1:length(data_files)) {
    data_loaded <- load_data(
      DataFile = data_files[i],
      NDataSet = 1,
      verbose = TRUE,
      EchoFile = file.path(temp_dir, "test_echo.out")
    )
    testthat::expect_true(length(data_loaded) == 16)
    if (basename(data_files[i]) == "Sable.dat") {
      testthat::expect_true(data_loaded$NReaders == 6)
      testthat::expect_true(all(tail(data_loaded$TheData[1, , ], 1) == c(6, -999, -999, -999, -999, 2, 2)))
    }
  }
})

# test DoApplyAgeError()
test_that("Can run DoApplyAgeError()", {
  example_path <- system.file("extdata", package = "AgeingError")
  dat <- load_data(DataFile = file.path(example_path, "WHS2.dat"))
  spc <- load_specs(SpecsFile = file.path(example_path, "WHS2.spc"), DataSpecs = dat)

  model <- AgeingError::DoApplyAgeError(
    DataSpecs = dat,
    ModelSpecsInp = spc,
    SaveDir = temp_dir,
    verbose = TRUE
  )

  testthat::expect_true(file.exists(file.path(temp_dir, "AgeingError.lda")))
  testthat::expect_equal(model$fitv, 3576.323, tolerance = 0.001)
})

# test ProcessResults()
test_that("Can run ProcessResults()", {
  Output <- AgeingError::ProcessResults(
    SaveDir = temp_dir,
    CalcEff = FALSE,
    verbose = FALSE
  )
  testthat::expect_equal(Output$ModelSelection$BIC, -7096.8, tolerance = 0.1)
  testthat::expect_equal(
    Output$ErrorAndBiasArray["CV", "Age 9", "Reader 1"],
    Output$ErrorAndBiasArray["CV", "Age 9", "Reader 2"]
  )
  testthat::expect_equal(Output$ErrorAndBiasArray["CV", "Age 9", "Reader 1"], 0.4575283, tolerance = 0.0001)
  testthat::expect_true(file.exists(file.path(temp_dir, "AgeingError.rpt")))
  testthat::expect_true(file.exists(file.path(temp_dir, "AgeingError-DataSet-1TrueVsReadsByReader.png")))
})
