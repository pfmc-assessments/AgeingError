#' Step-wise model selection
#'
#' Run step-wise model selection to facilitate the exploration of several
#' modeling configurations using Akaike information criterion (AIC).
#'
#' @details AIC seems like an appropriate method to select among possible
#' values for `PlusAge`, i.e., the last row of `SearchMat`, because `PlusAge`
#' determines the number of estimated fixed-effect hyperparameters that are
#' used to define the true proportion-at-age hyperdistribution. This
#' hyperdistribution is in turn used as a prior when integrating across a true
#' age associated with each otolith. This true age, which is a latent effect,
#' can be interpreted as a random effect with one for each observation. So, the
#' use of AIC to select among parameterizations of the fixed effects defining
#' this hyperdistribution is customary (Pinheiro and Bates, 2009). This was
#' tested for sablefish, where AIC lead to a true proportion at age that was
#' biologically plausible.
#' @inheritParams RunFn
#' @param SearchMat A data frame explaining stepwise model selection options.
#'   The final column must be named `label` and contain a unique label for each
#'   parameter row (e.g., `Error_Reader1`, ..., `Bias_Reader1`, ..., `MinusAge`,
#'   `PlusAge`). All preceding columns are candidate numeric options to search.
#'
#'   There should be one label for each readers error and one for each readers
#'   bias + 2 labels, one for `MinusAge`, i.e., the age where the proportion at
#'   age begins to decrease exponentially with decreasing age, and one for
#'   `PlusAge`, i.e., the age where the proportion-at-age begins to decrease
#'   exponentially with increasing age.
#'
#'   Each element of a given options row is a possible value to search across
#'   for that reader. So, the number of option columns of `SearchMat` will be
#'   the maximum number of options that you want to include. Think of it as
#'   several vectors stacked row-wise where shorter rows are filled in with
#'   `NA` values. If reader two only has two options that the analyst wants to
#'   search over the remainder of the columns should be filled with `NA` values
#'   for that row.
#' @param InformationCriterion A string specifying the type of information
#'   criterion that should be used to choose the best model. The default is to
#'   use AIC, though AIC corrected for small sample sizes and BIC are also
#'   available.
#' @param SelectAges A logical input specifying if the boundaries should be
#'   based on `MinusAge` and `PlusAge`. The default is `TRUE`.
#'
#' @references
#' Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
#' Quantifying age-reading error for use in fisheries stock assessments,
#' with application to species in Australia's southern and eastern scalefish
#' and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.
#'
#' Pinheiro, J.C., and Bates, D. 2009. Mixed-Effects Models in S and S-PLUS.
#' Springer, Germany.
#' @author James T. Thorson
#'
#' @export
#' @seealso
#' * `run()` will run a single model, where this function runs multiple models.
#' * `plot_output()` is called internally by `run()` when processing results.
#'
#' @examples
#'
#' \dontrun{
#' example(run)
#' ##### Run the model (MAY TAKE 5-10 MINUTES)
#' fileloc <- file.path(tempdir(), "age")
#' dir.create(fileloc, showWarnings = FALSE)
#' write_files(
#'   dat = AgeReads2,
#'   dir = fileloc,
#'   minage = MinAge,
#'   maxage = MaxAge,
#'   refage = 10,
#'   minusage = 1,
#'   plusage = 30,
#'   biasopt = BiasOpt,
#'   sigopt = SigOpt,
#'   knotages = KnotAges
#' )
#' out <- run(directory = fileloc)
#' out$output$ModelSelection
#'
#' ##### Stepwise selection
#'
#' # Parameters
#' MaxAge <- ceiling(max(AgeReads2) / 10) * 10
#' MinAge <- 1
#'
#' ##### Stepwise selection
#' StartMinusAge <- 1
#' StartPlusAge <- 30
#'
#' # Define data frame explaining stepwise model selection options
#' # One row for each reader + 2 rows for
#' # PlusAge (age where the proportion-at-age begins to
#' # decrease exponentially with increasing age) and
#' # MinusAge (the age where the proportion-at-age begins to
#' # decrease exponentially with decreasing age)
#' # Each element of a given row is a possible value to search
#' # across for that reader
#' SearchMat <- as.data.frame(array(NA, dim = c(Nreaders * 2 + 2, 7)))
#' names(SearchMat) <- paste("Option", 1:7)
#' SearchMat$label <- c(
#'   paste("Error_Reader", 1:Nreaders),
#'   paste("Bias_Reader", 1:Nreaders), "MinusAge", "PlusAge"
#' )
#' # Readers 1 and 3 search across options 1-3 for ERROR
#' SearchMat[c(1, 3), c("Option 1", "Option 2", "Option 3")] <- rep(1, 2) %o% c(1, 2, 3)
#' # Reader 2 mirrors reader 1
#' SearchMat[2, "Option 1"] <- -1
#' # Reader 4 mirrors reader 3
#' SearchMat[4, "Option 1"] <- -3
#' # Reader 1 has no BIAS
#' SearchMat[5, "Option 1"] <- 0
#' # Reader 2 mirrors reader 1
#' SearchMat[6, "Option 1"] <- -1
#' # Reader 3 search across options 0-2 for BIAS
#' SearchMat[7, c("Option 1", "Option 2", "Option 3")] <- c(1, 2, 0)
#' # Reader 4 mirrors reader 3
#' SearchMat[8, "Option 1"] <- -3
#' # MinusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
#' SearchMat[9, paste("Option", 1:7)] <- c(
#'   StartMinusAge,
#'   StartMinusAge - 10,
#'   StartMinusAge - 4,
#'   StartMinusAge - 1,
#'   StartMinusAge + 1,
#'   StartMinusAge + 4,
#'   StartMinusAge + 10
#' )
#' SearchMat[9, paste("Option", 1:7)] <- ifelse(SearchMat[9, paste("Option", 1:7)] < MinAge,
#'   NA, SearchMat[9, paste("Option", 1:7)]
#' )
#' # PlusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
#' SearchMat[10, paste("Option", 1:7)] <- c(
#'   StartPlusAge,
#'   StartPlusAge - 10,
#'   StartPlusAge - 4,
#'   StartPlusAge - 1,
#'   StartPlusAge + 1,
#'   StartPlusAge + 4,
#'   StartPlusAge + 10
#' )
#' SearchMat[10, paste("Option", 1:7)] <- ifelse(SearchMat[10, paste("Option", 1:7)] > MaxAge,
#'   NA, SearchMat[10, paste("Option", 1:7)]
#' )
#'
#' # Run model selection
#' # This outputs a series of files
#' # 1. "Stepwise - Model loop X.txt" --
#' #   Shows the AIC/BIC/AICc value for all different combinations
#' #   of parameters arising from changing one parameter at a time
#' #   according to SearchMat during loop X
#' # 2. "Stepwise - Record.txt" --
#' #   The Xth row of IcRecord shows the record of the
#' #   Information Criterion for all trials in loop X,
#' #   while the Xth row of StateRecord shows the current selected values
#' #   for all parameters at the end of loop X
#' # 3. Standard plots for each loop
#' # WARNING: One run of this stepwise model building example can take
#' # 8+ hours, and should be run overnight
#' stepwise(
#'   SearchMat = SearchMat, Data = AgeReads2,
#'   NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge,
#'   RefAge = 10, MaxSd = 40, MaxExpectedAge = MaxAge + 10,
#'   SaveFile = fileloc, InformationCriterion = c("AIC", "AICc", "BIC")[3]
#' )
#' }
#'
stepwise <- function(SearchMat,
                       Data,
                       NDataSets,
                       KnotAges,
                       MinAge,
                       MaxAge,
                       RefAge,
                       MaxSd,
                       MaxExpectedAge,
                       SaveFile,
                       EffSampleSize = 0,
                       Intern = TRUE,
                       InformationCriterion = c("AIC", "AICc", "BIC"),
                       SelectAges = TRUE) {

  # stepwise currently assumes one aggregated data set per run. The helper
  # file-writing functions can encode multiple sets, but this workflow has not
  # been generalized or tested for that case yet.
  if (NDataSets != 1) {
    cli::cli_abort("stepwise() currently supports only NDataSets = 1 with the TMB workflow.")
  }

  InformationCriterion <- match.arg(InformationCriterion)

  if (is.matrix(SearchMat)) {
    SearchMat <- as.data.frame(SearchMat)
  }

  if (!is.data.frame(SearchMat)) {
    cli::cli_abort("SearchMat must be a data frame.")
  }

  if (!"label" %in% names(SearchMat)) {
    cli::cli_abort("SearchMat must include a final column named 'label'.")
  }

  if (tail(names(SearchMat), 1) != "label") {
    cli::cli_abort("SearchMat must include 'label' as the final column.")
  }

  if (ncol(SearchMat) < 2) {
    cli::cli_abort("SearchMat must include at least one option column and a final label column.")
  }

  RowLabels <- as.character(SearchMat$label)
  if (anyNA(RowLabels) || any(trimws(RowLabels) == "")) {
    cli::cli_abort("SearchMat labels must be non-missing, non-empty strings.")
  }

  if (anyDuplicated(RowLabels)) {
    cli::cli_abort("SearchMat labels must be unique.")
  }

  OptionCols <- names(SearchMat)[names(SearchMat) != "label"]
  SearchMat <- as.matrix(SearchMat[, OptionCols, drop = FALSE])
  storage.mode(SearchMat) <- "numeric"
  rownames(SearchMat) <- RowLabels

  if (!is.data.frame(Data)) {
    Data <- as.data.frame(Data)
  }

  # Downstream code expects the historical missing-value sentinel rather than
  # NA. Do this once up front so every trial uses identical cleaned input data.
  Data[is.na(Data)] <- -999

  # Infer reader count from the standard data format:
  # first column = count, remaining columns = reader observations.
  Nreaders <- ncol(Data) - 1

  # SearchMat rows are expected in this order:
  # [reader sigmas][reader biases][MinusAge][PlusAge].
  if (nrow(SearchMat) != (2 * Nreaders + 2)) {
    cli::cli_abort(
      "SearchMat must have 2 * Nreaders + 2 rows, where Nreaders = ncol(Data) - 1."
    )
  }

  # Keep all stepwise artifacts under SaveFile.
  fs::dir_create(SaveFile)

  ErrorRows <- grep("^Error_Reader[0-9]+$", RowLabels)
  BiasRows <- grep("^Bias_Reader[0-9]+$", RowLabels)
  MinusRow <- match("MinusAge", RowLabels)
  PlusRow <- match("PlusAge", RowLabels)

  if (length(ErrorRows) != Nreaders || length(BiasRows) != Nreaders || is.na(MinusRow) || is.na(PlusRow)) {
    cli::cli_abort(
      "SearchMat labels must include Error_Reader1..N, Bias_Reader1..N, MinusAge, and PlusAge."
    )
  }

  ErrorRows <- ErrorRows[order(as.integer(sub("^Error_Reader", "", RowLabels[ErrorRows])))]
  BiasRows <- BiasRows[order(as.integer(sub("^Bias_Reader", "", RowLabels[BiasRows])))]

  # Current best parameter vector starts from the first option in each row.
  ParamVecOpt <- SearchMat[, 1]
  NsearchCols <- ncol(SearchMat)

  # Rebuild age-boundary search options after each loop.
  # The first element is always the current best age, then nearby values from
  # the legacy kernel (0, -10, -4, -1, +1, +4, +10) are used where possible.
  # Output is padded or truncated to match the number of SearchMat columns.
  make_age_options <- function(current_age, min_age = -Inf, max_age = Inf, select_ages = TRUE) {
    if (!select_ages) {
      return(c(current_age, rep(NA_real_, NsearchCols - 1)))
    }

    kernel <- c(0, -10, -4, -1, 1, 4, 10)
    options <- current_age + kernel
    options[options < min_age | options > max_age] <- NA_real_

    if (NsearchCols <= length(options)) {
      return(options[seq_len(NsearchCols)])
    }

    c(options, rep(NA_real_, NsearchCols - length(options)))
  }

  Stop <- FALSE
  IcRecord <- NULL
  StateRecord <- NULL
  OuterIndex <- 0

  # Continue searching until Stop==TRUE
  while (Stop == FALSE) {
    # Per-loop bookkeeping objects.
    OuterIndex <- OuterIndex + 1
    Index <- 0
    IcVec <- NULL
    ParamMat <- NULL
    Reports <- list()
    ParamVecOptPreviouslyEstimates <- FALSE

    # Evaluate all one-step neighbors of the current parameter vector.
    # For each row of SearchMat, try each non-NA option while holding all
    # other parameters fixed.
    for (VarI in seq_len(nrow(SearchMat))) {
      for (ValueI in seq_along(stats::na.omit(SearchMat[VarI, ]))) {
        # Update the current vector of parameters
        ParamVecCurrent <- ParamVecOpt
        ParamVecCurrent[VarI] <- stats::na.omit(SearchMat[VarI, ])[ValueI]

        # Run each unique candidate once per loop. This avoids duplicate model
        # fits when the selected value is repeated in SearchMat.
        if ((all(ParamVecCurrent == ParamVecOpt) && ParamVecOptPreviouslyEstimates == FALSE) || !all(ParamVecCurrent == ParamVecOpt)) {
          # If running the current optimum, change so that it won't run again this loop
          if (all(ParamVecCurrent == ParamVecOpt)) ParamVecOptPreviouslyEstimates <- TRUE

          # Use one working directory per trial; files are overwritten each run.
          RunFile <- file.path(SaveFile, "Run")
          fs::dir_create(RunFile)

          # Increment Index
          Index <- Index + 1
          print(paste("Loop=", OuterIndex, " Run=", Index, " StartTime=", date(), sep = ""))

          # Split full parameter vector into components expected by write_files().
          SigOpt <- as.numeric(ParamVecCurrent[ErrorRows])
          BiasOpt <- as.numeric(ParamVecCurrent[BiasRows])
          MinusAge <- ParamVecCurrent[MinusRow]
          PlusAge <- ParamVecCurrent[PlusRow]

          # Build fresh .dat/.spc files and run one TMB fit for this candidate.
          write_files(
            dat = Data,
            dir = RunFile,
            file_dat = "data.dat",
            file_specs = "data.spc",
            minage = MinAge,
            maxage = MaxAge,
            refage = RefAge,
            minusage = MinusAge,
            plusage = PlusAge,
            biasopt = BiasOpt,
            sigopt = SigOpt,
            knotages = KnotAges
          )

          Out <- run(
            directory = RunFile,
            file_data = "data.dat",
            file_specs = "data.spc"
          )

          # Pull model selection metrics from run() output.
          Aic <- Out$output$ModelSelection$AIC
          Aicc <- Out$output$ModelSelection$AICc
          Bic <- Out$output$ModelSelection$BIC
          if (InformationCriterion == "AIC") IcVec <- c(IcVec, Aic)
          if (InformationCriterion == "AICc") IcVec <- c(IcVec, Aicc)
          if (InformationCriterion == "BIC") IcVec <- c(IcVec, Bic)

          # Store tested parameter vectors in the same order as IcVec.
          ParamMat <- rbind(ParamMat, ParamVecCurrent)
          utils::write.table(cbind(IcVec, ParamMat),
            file.path(SaveFile, paste0("Stepwise - Model loop ", OuterIndex, ".txt")),
            sep = "\t", row.names = FALSE
          )

          # Keep report from the selected run for loop-level output
          ReportPath <- file.path(RunFile, "AgeingError.rpt")
          if (file.exists(ReportPath)) {
            Reports[[Index]] <- readLines(ReportPath)
          }
        } # End if-statement for only running ParamVecOpt once per loop
      }
    } # End loop accross VarI and ValueI

    # Append this loop's criterion values and current selected state.
    IcRecord <- rbind(IcRecord, IcVec)
    StateRecord <- rbind(StateRecord, ParamVecOpt)
    utils::capture.output(
      list(
        IcRecord = IcRecord,
        StateRecord = StateRecord
      ),
      file = file.path(SaveFile, "Stepwise - Record.txt")
    )

    # Select the best candidate under the chosen criterion (smaller is better).
    Min <- which.min(IcVec)
    # Stop once no one-step neighbor improves the objective.
    if (all(ParamMat[Min, ] == ParamVecOpt)) Stop <- TRUE
    ParamVecOpt <- ParamMat[Min, ]

    # Refresh MinusAge options around the selected value for the next loop.
    CurrentMinusAge <- ParamVecOpt[MinusRow]
    SearchMat[MinusRow, ] <- make_age_options(
      current_age = CurrentMinusAge,
      min_age = MinAge,
      max_age = Inf,
      select_ages = SelectAges
    )

    # Refresh PlusAge options around the selected value for the next loop.
    CurrentPlusAge <- ParamVecOpt[PlusRow]
    SearchMat[PlusRow, ] <- make_age_options(
      current_age = CurrentPlusAge,
      min_age = -Inf,
      max_age = MaxAge,
      select_ages = SelectAges
    )
    # Persist the report from the selected candidate for this loop.
    if (length(Reports) >= Min && !is.null(Reports[[Min]])) {
      writeLines(Reports[[Min]], con = file.path(SaveFile, "AgeingError.rpt"))
    }

    # Copy final files from the best run to the SaveFile directory.
    RunFile <- file.path(SaveFile, "Run")
    FilesToCopy <- list.files(
      RunFile,
      pattern = "^AgeingError(\\.|-|_SS3_format_)",
      full.names = TRUE
    )
    if (length(FilesToCopy) > 0) {
      file.copy(FilesToCopy, to = SaveFile, overwrite = TRUE)
    }
  } # End while statement

  # Return full trace of the search path for downstream inspection.
  invisible(
    list(
      IcRecord = IcRecord,
      StateRecord = StateRecord,
      BestParameters = ParamVecOpt,
      Label = RowLabels,
      InformationCriterion = InformationCriterion
    )
  )
} # End stepwise


#' Deprecated function replaced by stepwise()
#'
#' @param ... Any arguments associated with the deprecated function
#' @description
#' `r lifecycle::badge("deprecated")`
#' StepwiseFn() has been replaced by [stepwise()]
#' @author James T. Thorson
#' @export
#' @seealso [stepwise()]
StepwiseFn <- function(...) {
  lifecycle::deprecate_warn(
    when = "2.2.1",
    what = "StepwiseFn()",
    with = "stepwise()"
  )
  stepwise(...)
}
