#' Process results of the ageing error estimation
#'
#' @inheritParams DoApplyAgeError
#' @param CalcEff Calculate effective sample sizes (TRUE/FALSE)
#' @export
#' @author Andre E. Punt
ProcessResults <- function(Species = "AgeingError",
                           SaveDir = getwd(),
                           CalcEff = FALSE,
                           verbose = FALSE) {
  SaveFile <- file.path(
    SaveDir,
    paste0("/", Species, ".lda")
  )
  if (verbose) {
    print(SaveFile)
  }
  ReportFile <- file.path(
    SaveDir,
    paste0("/", Species, ".rpt")
  )
  if (verbose) {
    print(ReportFile)
  }

  load(SaveFile)
  if (verbose) {
    print(str(SaveAll))
  }
  NDataSet <- SaveAll$data$NDataSet
  MaxReader <- SaveAll$data$MaxReader

  write(
    paste(SaveAll$report$f, " ", SaveAll$report$Obj_fun),
    ReportFile
  )

  write(
    paste("Total number of readers:", MaxReader),
    ReportFile,
    append = TRUE
  )
  write(
    paste("Number of data sets:", NDataSet),
    ReportFile,
    append = TRUE
  )
  write(
    paste("Bias options by reader:", SaveAll$data$BiasOpt),
    ReportFile,
    append = TRUE
  )
  write(
    paste("Sigma options by reader:", SaveAll$data$SigOpt),
    ReportFile,
    append = TRUE
  )

  write(
    paste("Total objective function:", SaveAll$report$f),
    ReportFile,
    append = TRUE
  )
  Index <- which(abs(SaveAll$gradient) == max(abs(SaveAll$gradient)))
  write(
    paste("maximum gradient:", SaveAll$gradient[Index]),
    ReportFile,
    append = TRUE
  )

  write(
    paste("Number of readers: ", MaxReader),
    ReportFile,
    append = TRUE
  )
  write(
    paste("Range of ages: ", SaveAll$data$MinAge, " - ", SaveAll$data$MaxAge),
    ReportFile,
    append = TRUE
  )
  for (Ireader in 1:NDataSet) {
    write(
      paste(
        "Reader #", Ireader, " Minus/Plus ages: ",
        SaveAll$data$MinusA[Ireader], " / ", SaveAll$data$PlusA[Ireader]
      ),
      ReportFile,
      append = TRUE
    )
  }
  write(paste("Number of data sets:", NDataSet), ReportFile, append = TRUE)
  write(
    paste("Number of lines of data per data set: ", SaveAll$data$Npnt),
    ReportFile,
    append = TRUE
  )
  write(
    paste("Number of data points per data set:", SaveAll$data$TotalN),
    ReportFile,
    append = TRUE
  )
  write(
    "\n",
    ReportFile,
    append = TRUE
  )
  write(
    "# Bias options",
    ReportFile,
    append = TRUE
  )
  write(
    "#   -X: Mirrored with pattrn X",
    ReportFile,
    append = TRUE
  )
  write(
    "#   0: Unbiased [0 parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   1: Linearly proportional to age [1 parameter]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   2: Michaelis-Menten [3 parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    "\n# SD options",
    ReportFile,
    append = TRUE
  )
  write(
    "#   1: Constant CV [1 parameter]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   2: SD Michaelis-Menten [3 parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   3: CV Michaelis-Menten [3 parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   4: Known perfectly [0 parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   5: SD spline function of age [variable parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   6: SD linear piecewise function of age [
      variable parameters
      ]", ReportFile,
    append = TRUE
  )
  write(
    "#   7: SD linear function of age [2 parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    "#   8: CV linear function of age [2 parameters]",
    ReportFile,
    append = TRUE
  )
  write(
    paste("Reader BiasType SigmaType"),
    ReportFile,
    append = TRUE
  )
  for (Ireader in 1:MaxReader) {
    write(
      paste(
        Ireader, " ",
        SaveAll$data$BiasOpt[Ireader], " ",
        SaveAll$data$SigOpt[Ireader]
      ),
      ReportFile,
      append = TRUE
    )
  }
  write("", ReportFile, append = TRUE)

  # Bias and variance
  write("Reader Age CV SD Expected age", ReportFile, append = TRUE)
  for (Ireader in 1:MaxReader) {
    for (Age in 0:SaveAll$data$MaxAge) {
      if (Age > 1) {
        CVV <- SaveAll$report$TheSD[Ireader, Age + 1] / (Age * 1)
      } else {
        CVV <- SaveAll$report$TheSD[Ireader, Age + 1]
      }
      write(
        paste(
          Ireader, " ", Age, " ", round(CVV, 5), " ",
          round(SaveAll$report$TheSD[Ireader, Age + 1], 5), " ",
          round(SaveAll$report$TheBias[Ireader, Age + 1], 5)
        ),
        ReportFile,
        append = TRUE
      )
    }
  }
  write("", ReportFile, append = TRUE)

  # Estimated age-structure
  write("Estimated age-structure by data set", ReportFile, append = TRUE)
  write("===================================", ReportFile, append = TRUE)
  HeadString <- "Age "
  for (IDataSet in 1:NDataSet) {
    HeadString <- paste(
      HeadString,
      "Data set#",
      IDataSet
    )
  }
  write(HeadString, ReportFile, append = TRUE)
  for (Age in 0:SaveAll$data$MaxAge) {
    HeadString <- "Age "
    for (IDataSet in 1:NDataSet) {
      HeadString <- paste(
        HeadString,
        SaveAll$report$Aprob[IDataSet, Age + 1]
      )
    }
    write(HeadString, ReportFile, append = TRUE)
  }
  write("", ReportFile, append = TRUE)

  # Age-reading error matrices
  write("Final age-reading error matrices", ReportFile, append = TRUE)
  for (Ireader in 1:MaxReader) {
    write(paste("Matrix for reader# ", Ireader), ReportFile, append = TRUE)
    write(
      t(SaveAll$report$AgeErrOut[Ireader, , ]),
      ncolumns = SaveAll$data$MaxAge + 1,
      ReportFile,
      append = TRUE
    )
  }

  if (!is.null(SaveAll$sdreport)) {
    write("\nVariable estimate SD", ReportFile, append = TRUE)
    names <- row.names(summary(SaveAll$sdreport))
    write(t(cbind(names, summary(SaveAll$sdreport))), ReportFile,
      append = TRUE,
      ncolumns = 3
    )
  }

  # find the effective sample sizes
  if (CalcEff == TRUE) {
    ReaderSumm <- SaveAll$data$ReaderSumm
    ReaderStruc <- SaveAll$data$ReaderStruc
    write("Compute the effective sample sizes", ReportFile, append = TRUE)
    write("==================================", ReportFile, append = TRUE)
    GroupPointer <- 0
    Ages <- rep(0, MaxReader)
    ProbStore <- rep(0, SaveAll$data$MaxCells)
    for (IDataSet in 1:NDataSet) {
      write(paste("Data set: ", IDataSet), ReportFile, append = TRUE)
      write(
        paste("Data_set Group Group Line Readers Obs Obs_Numbers Pred_Numbers"),
        ReportFile,
        append = TRUE
      )
      Top <- 0
      Bot <- 0
      for (Kgroup in 1:ReaderSumm[IDataSet, 3]) {
        GroupPointer <- GroupPointer + 1

        # Find the total number of combinations of ages
        # (MaxAge+1)**number_of_readers
        Ncells <- 1
        for (Iread in 1:ReaderSumm[IDataSet, 2]) {
          Ncells <- Ncells * (SaveAll$data$MaxAge +
            1)
        }

        # Move through each possible combination of ages between 0
        # and MaxAge
        TotalProb <- 0
        for (II in 1:Ncells) {
          # Find the ages for this 'cell'
          Ndiv <- II
          for (Iread in 1:SaveAll$data$Nread[IDataSet]) Ages[Iread] <- -1
          for (Iread in 1:SaveAll$data$Nread[IDataSet]) {
            if (ReaderStruc[GroupPointer, Iread + 2] > 0) {
              DivJ <- 1
              if (Iread < SaveAll$data$Nread[IDataSet]) {
                for (Jread in SaveAll$data$Nread[IDataSet]:(Iread + 1)) {
                  if (Jread > 0 & Jread <= SaveAll$data$Nread[IDataSet]) {
                    if (ReaderStruc[GroupPointer, Jread + 2] > 0) {
                      DivJ <- DivJ * (SaveAll$data$MaxAge + 1)
                    }
                  }
                }
              }
              DivI <- floor((Ndiv - 1) / DivJ)
              Ages[Iread] <- DivI
              Ndiv <- Ndiv - DivI * DivJ
            }
          }

          # Find the probability for this cell, i.e. the probability
          # of an ageing reading of Ages(1)&Ages(2)&...
          if (SaveAll$data$Iperfect[IDataSet] == 0) {
            Prob1 <- 0
            for (Age1 in SaveAll$data$MinAge:SaveAll$data$MaxAge) {
              # Prior probability * product over readers
              Prob2 <- SaveAll$report$Aprob[IDataSet, Age1 + 1]
              for (Ireader in 1:SaveAll$data$Nread[IDataSet]) {
                Jreader <- SaveAll$data$ReadPnt[IDataSet, Ireader]
                AgeA <- Ages[Ireader]
                if (AgeA >= 0) {
                  Prob2 <- Prob2 * SaveAll$report$AgeErrOut[Jreader, Age1 +
                    1, AgeA + 1]
                }
              }
              Prob1 <- Prob1 + Prob2
            }
            ProbStore[II] <- Prob1
            TotalProb <- TotalProb + Prob1
          } else {
            # Product over readers
            Age1 <- Ages[1]
            Prob1 <- 1
            for (Ireader in 2:SaveAll$data$Nread[IDataSet]) {
              Jreader <- SaveAll$data$ReadPnt[IDataSet, Ireader]
              AgeA <- Ages[Ireader]
              if (Age1 >= 0 & AgeA >= 0) {
                Prob1 <- Prob1 * SaveAll$report$AgeErrOut[Ireader, Age1 +
                  1, AgeA + 1]
              }
            }
            ProbStore[II] <- Prob1
            TotalProb <- TotalProb + Prob1
          }
        } # for (II in 1:Ncells)

        # Now compute
        writeout <- NULL
        for (II in 1:Ncells) {
          # Find the ages for this 'cell'
          Ndiv <- II
          for (Iread in 1:SaveAll$data$Nread[IDataSet]) Ages[Iread] <- -1
          for (Iread in 1:SaveAll$data$Nread[IDataSet]) {
            if (ReaderStruc[GroupPointer, Iread + 2] > 0) {
              DivJ <- 1
              if (Iread < SaveAll$data$Nread[IDataSet]) {
                for (Jread in SaveAll$data$Nread[IDataSet]:(Iread + 1)) {
                  if (Jread > 0 & Jread <= SaveAll$data$Nread[IDataSet]) {
                    if (ReaderStruc[GroupPointer, Jread + 2] > 0) {
                      DivJ <- DivJ * (SaveAll$data$MaxAge + 1)
                    }
                  }
                }
              }
              DivI <- floor((Ndiv - 1) / DivJ)
              Ages[Iread] <- DivI
              Ndiv <- Ndiv - DivI * DivJ
            }
          }

          # Check for a match
          Jfound <- 0
          Pobs <- 0
          for (JJ in 1:SaveAll$data$Npnt[IDataSet]) {
            # Set Ifound to 1 if the ages don't match properly
            Ifound <- 0
            for (Iread in 1:SaveAll$data$Nread[IDataSet]) {
              if (Ages[Iread] != SaveAll$data$TheData[IDataSet, JJ, Iread +
                1] & SaveAll$data$TheData[IDataSet, JJ, Iread + 1] >= 0) {
                Ifound <- 1
              }
            }
            # We have a match so store the results
            if (Ifound == 0) {
              Pobs <- SaveAll$data$TheData[IDataSet, JJ, 1] / ReaderStruc[
                GroupPointer,
                2
              ]
              Jfound <- 1
              KK <- JJ
            }
          }

          # Find the probability for this cell
          Pest <- ProbStore[II] / TotalProb + 1e-190
          if (Jfound == 1) {
            HeadString <- paste(
              "Data Point: ", IDataSet, " ", GroupPointer,
              " ", Kgroup, " ", KK, " "
            )
            for (Iread in 1:SaveAll$data$Nread[IDataSet]) {
              HeadString <- paste(
                HeadString,
                Ages[Iread], " "
              )
            }
            HeadString <- paste(
              HeadString, SaveAll$data$TheData[
                IDataSet,
                KK, 1
              ], " ", round(
                Pobs * SaveAll$data$TotalN[IDataSet],
                7
              ), " ", round(Pest * SaveAll$data$TotalN[IDataSet], 12),
              " "
            )
            writeout <- rbind(writeout, HeadString)
          }

          # Compute the effective sample size muliplier
          Top <- Top + Pest * (1 - Pest) / ReaderStruc[GroupPointer, 2]
          Bot <- Bot + (Pest - Pobs)^2
        }
        write(writeout, ReportFile, append = TRUE)
      } # Kgroup
      EffPred <- Top / Bot * SaveAll$data$TotalN[IDataSet]
      write("Data_Set Predicted_EFF Assumed_Eff Sample_size", ReportFile,
        append = TRUE
      )
      write(paste(
        IDataSet, " ", EffPred, " ", SaveAll$data$EffN[IDataSet],
        " ", SaveAll$data$TotalN[IDataSet]
      ), ReportFile, append = TRUE)
      write("", ReportFile, append = TRUE)
    }
  } # IDataSet
  Npars <- length(SaveAll$sdreport$par.fixed)

  # find the effective sample sizes
  for (IDataSet in 1:NDataSet) {
    Data <- SaveAll$data$TheData[IDataSet, 1:SaveAll$data$Npnt[IDataSet], ]
    Output <- plot_output(
      Data = Data,
      IDataSet = IDataSet,
      MaxAge = SaveAll$data$MaxAge,
      Report = SaveAll$report,
      Nparameters = Npars,
      LogLike = SaveAll$report$Obj_fun,
      subplot = 1:3,
      ReaderNames = NULL,
      Species = Species,
      SaveDir = SaveDir
    )
  }
  # print(str(Output))
  return(Output)
}
