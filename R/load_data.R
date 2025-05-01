#' Deprecated function renamed to load_data()
#'
#' @param ... Any arguments associated with the deprecated function
#' @description
#' `r lifecycle::badge("deprecated")`
#' CreateData() has been renamed as [load_data()] to better clarify its purpose
#' @author Ian G. Taylor
#' @export
#' @seealso [load_data()]
CreateData <-
  function(...) {
    lifecycle::deprecate_stop(
      when = "2.1.0",
      what = "CreateData()",
      with = "load_data()"
    )
  }

#' Load the formatted ageing error data
#'
#' @param DataFile Filename for input data
#' @param NDataSet Number of data sets within `DataFile`
#' @param verbose Return messages to the console (in addition to any output to
#' `EchoFile`)
#' @param EchoFile A file path to a file that will be created or appended to if
#'   it already exists to store information about your data inputs. The default
#'   is `''`, which leads to output being printed to the screen rather than
#'   saved in a file. An example of a user-defined input would be
#'   `'EchoTMB.out'`.
#' @export
#' @author Andre E. Punt
load_data <- function(DataFile = "data.dat",
                      NDataSet = 1,
                      verbose = FALSE,
                      EchoFile = "") {
  # Put a first line in the EchoFile if the file does not already exist
  if (!file.exists(EchoFile)) {
    write("", file = EchoFile)
  }
  MatchTable <- function(Table,
                         Char1 = NULL,
                         Char2 = NULL,
                         Char3 = NULL,
                         Char4 = NULL,
                         Char5 = NULL) {
    ii <- rep(T, length(Table[, 1]))
    if (!is.null(Char1)) {
      ii <- ii & (Table[, 1] == Char1)
    }
    if (!is.null(Char2)) {
      ii <- ii & (Table[, 2] == Char2)
    }
    if (!is.null(Char3)) {
      ii <- ii & (Table[, 3] == Char3)
    }
    if (!is.null(Char4)) {
      ii <- ii & (Table[, 4] == Char4)
    }
    if (!is.null(Char5)) {
      ii <- ii & (Table[, 5] == Char5)
    }
    ii <- seq(1:length(Table[, 1]))[ii]
    return(ii)
  }

  # Read in the data file
  Data <- read.table(
    DataFile,
    comment.char = "?",
    fill = TRUE,
    blank.lines.skip = TRUE,
    stringsAsFactors = FALSE,
    col.names = 1:100
  )

  # Extract the minimum and maximum ages
  Index <- MatchTable(Data, Char1 = "Range_of_ages")
  MinAge <- as.numeric(Data[Index + 1, 1])
  MaxAge <- as.numeric(Data[Index + 1, 2])

  # Details of the readers
  IndexVals <- rep(0, NDataSet)
  Npnt <- rep(0, NDataSet)
  NReaders <- rep(0, NDataSet)
  MinusA <- rep(0, NDataSet)
  PlusA <- rep(0, NDataSet)
  RefAge <- rep(0, NDataSet)
  MaxReader <- 1
  for (Idataset in 1:NDataSet) {
    SearchTerm <- paste("Data_set_", Idataset, sep = "")
    IndexVals[Idataset] <- MatchTable(Data, Char1 = SearchTerm)
    Npnt[Idataset] <- as.numeric(Data[IndexVals[Idataset] + 1, 1])
    NReaders[Idataset] <- as.numeric(Data[IndexVals[Idataset] + 2, 1])
    MinusA[Idataset] <- as.numeric(Data[IndexVals[Idataset] + 3, 1])
    PlusA[Idataset] <- as.numeric(Data[IndexVals[Idataset] + 3, 2])
    if (PlusA[Idataset] < MinusA[Idataset]) {
      cli::cli_abort("Plus age must be greater than or equal to the minus age")
    }
    RefAge[Idataset] <- as.numeric(Data[IndexVals[Idataset] + 3, 3])
    if (RefAge[Idataset] <= MinusA[Idataset]) {
      cli::cli_abort("Reference age must be greater than the minimum age")
    }
    if (RefAge[Idataset] >= PlusA[Idataset]) {
      cli::cli_abort("Reference age must be less than the maximum age")
    }
    Readers <- as.numeric(Data[IndexVals[Idataset] + 4, 1:NReaders[Idataset]])
    MaxReader <- max(MaxReader, Readers)
    if (verbose) {
      cli::cli_alert_info("readers {Readers}")
    }
  }
  ReadPnt <- matrix(0, nrow = NDataSet, ncol = MaxReader)
  for (Idataset in 1:NDataSet) {
    ReadPnt[Idataset, 1:NReaders[Idataset]] <- as.numeric(
      Data[IndexVals[Idataset] + 4, 1:NReaders[Idataset]]
    )
  }

  # Now extract the data
  TheData <- array(-999, dim = c(NDataSet, max(Npnt), MaxReader + 1))
  Ipnt <- 0
  for (Idataset in 1:NDataSet) {
    for (Iline in 1:Npnt[Idataset]) {
      TheData[Idataset, Iline, 1:(NReaders[Idataset] + 1)] <- as.numeric(
        Data[(IndexVals[Idataset] + 4 + Iline), 1:(NReaders[Idataset] + 1)]
      )
    }
    if (verbose) {
      cli::cli_alert_info("Last line of data set {Idataset} is {TheData[Idataset, Npnt[Idataset], 1:(NReaders[Idataset] + 1)]}")
    }
  }

  # Do checks on the data set
  MaAge <- -1
  MiAge <- 1000
  NegVals <- 0
  for (IDataS in 1:NDataSet) {
    for (Ipnt in 1:Npnt[IDataS]) {
      for (Ireader in 1:NReaders[IDataS]) {
        if (TheData[IDataS, Ipnt, Ireader + 1] >= 0) {
          if (TheData[IDataS, Ipnt, Ireader + 1] > MaAge) {
            MaAge <- TheData[IDataS, Ipnt, Ireader + 1]
          }
          if (TheData[IDataS, Ipnt, Ireader + 1] < MiAge) {
            MiAge <- TheData[IDataS, Ipnt, Ireader + 1]
          }
        } else {
          NegVals <- 1
        }
      }
    }
  }
  if (NegVals == 1) {
    cli::cli_alert_warning("There are some missing data; the effective sample size calculation may be dubious")
  }

  # Create a tabular summary of the data
  write("Structure of the data set", EchoFile, append = TRUE)
  write("Data set # Entries Reader boolean", EchoFile, append = TRUE)
  ReaderStruc <- matrix(0, nrow = 1000, ncol = MaxReader + 2)
  Presense <- rep(0, MaxReader)
  ReaderSumm <- matrix(0, nrow = NDataSet, ncol = 3)
  NrowStruc <- 0
  for (IDataS in 1:NDataSet) {
    for (II in 1:Npnt[IDataS]) {
      Presense <- rep(0, MaxReader)
      for (Ireader in 1:NReaders[IDataS]) {
        if (TheData[IDataS, II, Ireader + 1] >= 0) {
          Presense[Ireader] <- ReadPnt[IDataS, Ireader]
        }
      }
      Ifound <- 0
      if (NrowStruc > 0) {
        for (JJ in 1:NrowStruc) {
          if (ReaderStruc[JJ, 1] == IDataS) {
            Jfound <- 1
            for (Ireader in 1:NReaders[IDataS]) {
              if (Presense[Ireader] != ReaderStruc[JJ, Ireader + 2]) {
                Jfound <- 0
              }
            }
            if (Jfound == 1) {
              Ifound <- JJ
            }
          }
        }
      }
      if (Ifound == 0) {
        NrowStruc <- NrowStruc + 1
        ReaderStruc[NrowStruc, 1] <- IDataS
        ReaderStruc[NrowStruc, 2] <- TheData[IDataS, II, 1]
        ReaderSumm[IDataS, 3] <- ReaderSumm[IDataS, 3] + 1
        for (Ireader in 1:NReaders[IDataS]) {
          ReaderStruc[NrowStruc, Ireader + 2] <- Presense[Ireader]
        }
      } else {
        ReaderStruc[Ifound, 2] <- ReaderStruc[Ifound, 2] + TheData[
          IDataS,
          II, 1
        ]
      }
    }
  }
  cli::cli_alert_info("Number of rows in NrowStruc {NrowStruc} = {NrowStruc}")
  ReaderStruc <- matrix(
    ReaderStruc[1:NrowStruc, ],
    nrow = NrowStruc,
    ncol = length(ReaderStruc[1, ])
  )
  print("ReaderStruc")
  print(ReaderStruc)
  for (II in 1:NrowStruc) {
    write(ReaderStruc[II, ], EchoFile, append = TRUE, ncolumns = MaxReader + 2)
  }
  print("ReaderSumm")
  print(ReaderSumm)

  # Create a reader summary; ReaderSumm[IDataS,2] is the maximum number of
  # readers for a given combination of readers
  for (II in 1:NrowStruc) {
    IDataS <- ReaderStruc[II, 1]
    ReaderSumm[IDataS, 1] <- IDataS
    MaxReaderOpt <- 0
    for (Ireader in 1:NReaders[IDataS]) {
      if (ReaderStruc[NrowStruc, Ireader + 2] > 0) {
        MaxReaderOpt <- MaxReaderOpt + 1
      }
    }
    if (MaxReaderOpt > ReaderSumm[IDataS, 2]) {
      ReaderSumm[IDataS, 2] <- MaxReaderOpt
    }
  }
  write("ReaderSumm", EchoFile, append = TRUE)
  write(t(ReaderSumm), EchoFile, append = TRUE, ncolumns = 3)
  print("ReaderSumm")
  print(ReaderSumm)

  # Outputs to screen
  write(
    paste("Number of reads by data set:       ", NReaders),
    EchoFile,
    append = TRUE
  )
  write(
    paste("Minimum and Maximum Ages:          ", MiAge, " ", MaAge),
    EchoFile,
    append = TRUE
  )
  write(
    "",
    EchoFile,
    append = TRUE
  )

  # Record total sample size and specify effective Ns
  EffNOpt <- rep(0, NDataSet)
  TotalN <- rep(0, NDataSet)
  EffN <- rep(0, NDataSet)
  for (IDataSet in 1:NDataSet) {
    TotalN[IDataSet] <- 0
    for (Ipnt in 1:Npnt[IDataSet]) {
      TotalN[IDataSet] <- TotalN[IDataSet] + TheData[IDataSet, Ipnt, 1]
    }
  }
  for (IDataSet in 1:NDataSet) {
    if (EffNOpt[IDataSet] <= 0) {
      EffN[IDataSet] <- TotalN[IDataSet]
    } else {
      EffN[IDataSet] <- EffNOpt[IDataSet]
    }
  }

  # Check for duplicates and condense as needed
  OneProblem <- 0
  for (IDataSet in 1:NDataSet) {
    Problem <- 0
    for (II in 2:Npnt[IDataSet]) {
      for (JJ in 1:(II - 1)) {
        if (TheData[IDataSet, JJ, 1] > 0) {
          Ifound <- 0
          for (Iread in 1:NReaders[IDataSet]) {
            if (TheData[IDataSet, JJ, Iread +
              1] != TheData[IDataSet, II, Iread + 1]) {
              Ifound <- 1
            }
          }
          if (Ifound == 0) {
            cli::cli_alert_warning("Lines {II} and {JJ} have the same ages")
            TheData[IDataSet, JJ, 1] <- TheData[IDataSet, JJ, 1] + TheData[
              IDataSet,
              II, 1
            ]
            TheData[IDataSet, II, 1] <- -1
            Problem <- 1
          }
        }
      } # JJ
    }
    if (Problem == 1) {
      cli::cli_alert_warning("Duplicate entries found for data set {IDataSet}; corrected data set in Echo.File")
      NLineOut <- 0
      for (II in 1:Npnt[IDataSet]) {
        if (TheData[IDataSet, II, 1] > 0) {
          NLineOut <- NLineOut + 1
          write(TheData[IDataSet, II, ], EchoFile, append = TRUE)
        }
      }
      write(paste("New lines ", NLineOut), EchoFile, append = TRUE)
      OneProblem <- 1
    }
  }

  ## Counter for storage
  MaxCells <- 1
  for (IDataS in 1:NDataSet) {
    Ncells <- 1
    for (Ireader in 1:ReaderSumm[IDataS, 2]) {
      Ncells <- Ncells * (MaxAge + 1)
    }
    if (Ncells > MaxCells) {
      MaxCells <- Ncells
    }
  }
  write(paste("total cells ", MaxCells), EchoFile, append = TRUE)

  Outs <- NULL
  Outs$MinAge <- MinAge
  Outs$MaxAge <- MaxAge
  Outs$NDataSet <- NDataSet
  Outs$MaxReader <- MaxReader
  Outs$TheData <- TheData
  Outs$Npnt <- Npnt
  Outs$ReadPnt <- ReadPnt
  Outs$NReaders <- NReaders
  Outs$MinusA <- MinusA
  Outs$PlusA <- PlusA
  Outs$RefAge <- RefAge
  Outs$ReaderSumm <- ReaderSumm
  Outs$ReaderStruc <- ReaderStruc
  Outs$MaxCells <- MaxCells
  Outs$TotalN <- TotalN
  Outs$EffN <- EffN
  if (verbose) {
    print(str(Outs))
  }
  return(Outs)
}
