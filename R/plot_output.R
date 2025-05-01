#' Plot output
#'
#' Plots age comparisons and results from the fitted Ageing Error model
#'
#' @param Data Input data matrix
#' @param IDataSet Index of the data set used in creating the filename
#' @param MaxAge Maximum estimated age
#' @param Report Results from fitting the model
#' @param subplot Vector of which plots to create.
#' @param Nparameters Number of parameters
#' @param LogLike Negative log likelihood from fitting the model
#' @param ReaderNames Vector with names of each reader, defaults to
#'   'Reader1', 'Reader2', etc. if left at the default argument of `NULL`.
#'   If you pass a vector of strings, the vector must be the same length as
#'   `NCOL(Data) - 1`.
#' @param Species String used at beginning of the output files
#' @param SaveDir Directory for fitted model
#' @param verbose Report messages as function runs.
#' @param ... Additional arguments passed to [ageing_comparison()].
#' @return Returns AIC, AICc, and BIC for fitted model.
#'
#' @references Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
#' Quantifying age-reading error for use in fisheries stock assessments,
#' with application to species in Australias southern and eastern scalefish
#' and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.
#'
#' @author James T. Thorson, Ian G. Taylor
#'
#' @export
#'
plot_output <- function(Data,
                        IDataSet,
                        MaxAge,
                        Report,
                        subplot = 1:3,
                        Nparameters = 0,
                        LogLike = 0,
                        ReaderNames = NULL,
                        Species = "AgeingError",
                        SaveDir = getwd(),
                        verbose = FALSE,
                        ...) {
  fs::dir_create(SaveDir)

  # Interpret inputs
  Nreaders <- ncol(Data) - 1
  Ages <- Nages <- MaxAge + 1

  # Reader names
  if (is.null(ReaderNames)) {
    ReaderNames <- paste0("Reader", 1:Nreaders)
  }

  # Age-reading error matrices: dimensions are Reader, TrueAge, EstAge
  MisclassArray <- array(
    NA,
    dim = c(Nreaders, Ages, Ages),
    dimnames = list(
      paste("Reader", 1:Nreaders),
      paste("TrueAge", 0:MaxAge),
      paste("EstAge", 0:MaxAge)
    )
  )
  for (i in 1:Nreaders) {
    MisclassArray[i, , ] <- Report$AgeErrOut[i, , ]
  }

  # Estimated age-structure
  AgeStruct <- cbind(0:MaxAge, t(Report$Aprob))

  # Reader CV, SD and Bias
  Temp <- matrix(0, nrow = 5 * Nages * Nreaders, ncol = 5)
  for (Ireader in 1:Nreaders) {
    yrange <- (Ireader - 1) * Nages
    Temp[yrange + 1:Nages, 1] <- Ireader
    Temp[yrange + 1:Nages, 2] <- 0:MaxAge
    Temp[yrange + 1:Nages, 3] <- Report$TheSD[Ireader, ] / c(1, 1:MaxAge)
    Temp[yrange + 1:Nages, 4] <- Report$TheSD[Ireader, ]
    Temp[yrange + 1:Nages, 5] <- Report$TheBias[Ireader, ]
  }
  Temp <- t(Temp)

  ErrorAndBiasArray <- array(
    as.numeric(Temp),
    dim = c(5, Nages, Nreaders),
    dimnames = list(
      c("Reader", "True_Age", "CV", "SD", "Expected_age"),
      paste("Age", 0:MaxAge),
      paste("Reader", 1:Nreaders)
    )
  )
  if (verbose) {
    print(str(ErrorAndBiasArray))
  }
  # Estimate unobserved age for each otolith This is done by assigning each
  # otolith to the age which has maximum posterior probability (i.e. the
  # conditional mode, as is typically done for random effects)
  AgeProbs <- array(
    NA,
    dim = c(nrow(Data), Ages),
    dimnames = list(
      paste("Otolith", 1:nrow(Data)),
      paste("TrueAge", 0:MaxAge)
    )
  )
  OtI <- AgeI <- ReadI <- 1
  for (OtI in 1:nrow(Data)) {
    for (AgeI in 1:Ages) {
      AgeProbs[OtI, AgeI] <- 1
      for (ReadI in 1:Nreaders) {
        if (Data[OtI, ReadI + 1] != -999) {
          AgeRead <- Data[OtI, ReadI + 1]
          AgeProbs[OtI, AgeI] <- AgeProbs[OtI, AgeI] * (MisclassArray[
            ReadI,
            AgeI, AgeRead + 1
          ])^Data[OtI, 1]
        } # end check for value other than -999
      } # end loop over readers
    } # end loop over ages
  } # end loop over rows of data

  # Remove MaxAge before calculating 'TrueAge' because the MaxAge is a
  # plus-group, and ends up with maximum probability for most ages in the
  # upper tail ANDRE - Found another issue - should be age-1 because the
  # first age is zero
  TrueAge <- apply(
    AgeProbs,
    MARGIN = 1,
    FUN = function(Vec) {
      order(Vec[-length(Vec)], decreasing = TRUE)[1]
    }
  ) - 1

  DataExpanded <- Data[rep(1:nrow(Data), Data[, 1]), -1]
  DataExpanded[DataExpanded == -999] <- NA

  # Plot comparison of data for each pair of readers
  if (1 %in% subplot) {
    # make plots of input data for each reader pair
    for (ireader in 1:(Nreaders - 1)) {
      for (jreader in (ireader + 1):Nreaders) {
        ageing_comparison(
          xvec = DataExpanded[, ireader],
          yvec = DataExpanded[, jreader],
          xlab = ReaderNames[ireader],
          ylab = ReaderNames[jreader],
          maxage = max(DataExpanded, na.rm = TRUE),
          hist = FALSE,
          png = TRUE,
          SaveFile = SaveDir,
          filename = paste0(
            Species,
            "-DataSet-",
            IDataSet,
            "-",
            ReaderNames[ireader],
            "Vs",
            ReaderNames[jreader],
            ".png"
          ),
          verbose = verbose,
          ...
        )
      }
    }
  } # end check for whether subplot 1 was requested

  #  Plot estimated age structure
  if (2 %in% subplot) {
    png(
      file.path(
        SaveDir,
        paste0(
          Species,
          "-DataSet-",
          IDataSet,
          "EstimatedVsObservedAgeStructure.png"
        )
      ),
      width = 6,
      height = 6,
      units = "in",
      res = 200
    )
    par(mar = c(3, 3, 2, 0), mgp = c(1.5, 0.25, 0), tck = -0.02, oma = c(
      0,
      0, 0, 0
    ) + 0.1)
    plot(
      x = AgeStruct[, 1],
      y = AgeStruct[, 2],
      type = "s",
      lwd = 2,
      xlab = "Age",
      ylab = "Prop",
      main = "Estimated=Black, Observed=Red"
    )
    hist(
      as.matrix(DataExpanded),
      add = TRUE,
      freq = FALSE,
      breaks = seq(0, MaxAge, by = 1),
      col = rgb(red = 1, green = 0, blue = 0, alpha = 0.3)
    )
    dev.off()
  } # end check for whether subplot was requested

  # Plot true age against different age reads
  if (3 %in% subplot) {
    Ncol <- ceiling(sqrt(Nreaders))
    Nrow <- ceiling(Nreaders / Ncol)
    png(
      file.path(
        SaveDir,
        paste0(
          Species,
          "-DataSet-",
          IDataSet,
          "TrueVsReadsByReader.png"
        )
      ),
      width = Ncol * 3,
      height = Nrow * 3,
      units = "in",
      res = 200
    )
    par(
      mfrow = c(Nrow, Ncol),
      mar = c(3, 3, 2, 0),
      mgp = c(1.5, 0.25, 0),
      tck = -0.02,
      oma = c(0, 0, 5, 0) + 0.1
    )
    for (ReadI in 1:Nreaders) {
      Main <- ReaderNames[ReadI]

      # Add 0.5 to match convention in Punt model that otoliths are
      # read half way through year
      Temp <- cbind(TrueAge, Data[, ReadI + 1] + 0.5)
      # Exclude rows with no read for this reader
      Temp <- Temp[which(Data[, ReadI + 1] != -999), ]
      plot(
        x = Temp[, 1],
        y = Temp[, 2],
        ylim = c(0, MaxAge),
        xlim = c(0, MaxAge),
        col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
        xlab = "Mode predicted age | parameters",
        ylab = "Read age",
        lwd = 2,
        main = Main,
        pch = 21,
        cex = 0.2
      )
      lines(x = c(0, MaxAge), y = c(0, MaxAge), lwd = 1, lty = "dashed")
      lines(
        x = ErrorAndBiasArray["True_Age", , ReadI],
        y = ErrorAndBiasArray["Expected_age", , ReadI],
        type = "l",
        col = "red",
        lwd = 1
      )
      lines(
        x = ErrorAndBiasArray["True_Age", , ReadI],
        y = ErrorAndBiasArray["SD", , ReadI],
        type = "l",
        col = "blue",
        lwd = 1
      )
      lines(
        x = ErrorAndBiasArray["True_Age", , ReadI],
        y = ErrorAndBiasArray["Expected_age", , ReadI] +
          2 * ErrorAndBiasArray["SD", , ReadI],
        type = "l",
        col = "red",
        lwd = 1,
        lty = "dashed"
      )
      lines(
        x = ErrorAndBiasArray["True_Age", , ReadI],
        y = ErrorAndBiasArray["Expected_age", , ReadI] -
          2 * ErrorAndBiasArray["SD", , ReadI],
        type = "l",
        col = "red",
        lwd = 1,
        lty = "dashed"
      )
    }
    mtext(
      side = 3,
      outer = TRUE,
      text = paste0(
        "Reads(dot), Sd(blue), expected_read(red solid line),\n",
        " and 95% CI for expected_read(red dotted line)"
      ),
      line = 1
    )
    dev.off()
  } # end check for whether subplot was requested

  ## AIC
  Nll <- LogLike
  Df <- Nparameters
  n <- sum(ifelse(Data[, -1] == -999, 0, 1))
  Aic <- 2 * Nll + 2 * Df
  Aicc <- Aic + 2 * Df * (Df + 1) / (n - Df - 1)
  Bic <- 2 * Nll + Df * log(n)

  # Write definitions to file
  for (ReadI in 1:Nreaders) {
    Main <- ReaderNames[ReadI]
    write.csv(
      ErrorAndBiasArray[, , ReadI],
      file = file.path(
        SaveDir,
        paste0(
          Species,
          "_SS3_format_",
          Main,
          ".csv",
          sep = ""
        )
      )
    )
  }


  # Return stuff
  ModelSelection <- list(AIC = Aic, AICc = Aicc, BIC = Bic)
  Output <- list(
    ModelSelection = ModelSelection,
    ErrorAndBiasArray = ErrorAndBiasArray
  )
  return(Output)
}
