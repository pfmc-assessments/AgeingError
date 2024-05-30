#' Plot comparison of double age readings
#'
#' Plot with circles proportional to how many double readings
#' fell in each pair of coordinates
#'
#' @param xvec vector of values from reader A
#' @param yvec vector of values from reader B
#' @param scale.pts Documentation needed.
#' @param col.pts color for points
#' @param col.hist color for histograms
#' @param counts include text within each bubble showing count of values?
#' @param maxage maximum age to include in the plot (doesn't yet work well)
#' @param hist include a histogram along each axis?
#' @param hist.frac maximum value of histograms as fraction of maxage
#' @param xlab label for xvec
#' @param ylab label for yvec
#' @param title Optional title to add at top of plot
#' @param png Save plot to PNG file?
#' @param filename File name for PNG file.
#' @param SaveFile directory where plot will be saved.
#' NULL value will make it go to working directory.
#' @param verbose Report messages as function runs.
#' @author Ian G. Taylor
#' @export

ageing_comparison <- function(xvec, yvec, scale.pts = 2,
                              col.pts = grDevices::grey(.1, alpha = .5),
                              col.hist = grDevices::rgb(0, 0, .5, alpha = .7),
                              counts = TRUE, maxage = NULL,
                              hist = TRUE, hist.frac = .1,
                              xlab = "Age reader A",
                              ylab = "Age reader B",
                              title = NULL,
                              png = FALSE,
                              filename = "ageing_comparison.png",
                              SaveFile = NULL,
                              verbose = TRUE) {
  # get counts of each pair
  df1 <- as.data.frame(table(xvec, yvec), stringsAsFactors = FALSE)
  df1$xvec <- as.numeric(df1$xvec)
  df1$yvec <- as.numeric(df1$yvec)
  # remove rows with count of zero
  df1 <- df1[df1[, 3] != 0, ]
  if (length(df1[, 1]) == 0) {
    return()
  }
  # print(df1)
  # get axis limits
  if (is.null(maxage)) {
    maxage <- ceiling(max(xvec, yvec, na.rm = TRUE))
  }
  # set regular vs. internal axes depending on whether min is 0
  axs <- "i"
  if (min(xvec, yvec, na.rm = TRUE) == 0) {
    axs <- "r"
  }
  # open PNG file if requested
  if (png) {
    if (is.null(SaveFile)) {
      SaveFile <- getwd()
    }
    if (verbose) {
      message("writing image to ", file.path(SaveFile, filename))
    }
    grDevices::png(file.path(SaveFile, filename),
      width = 6.5, height = 6.5, pointsize = 10,
      res = 300, units = "in"
    )
  }
  # make empty plot
  plot(0,
    type = "n", xlim = c(0, maxage + 1), ylim = c(0, maxage + 1),
    xlab = xlab, ylab = ylab, axes = F, xaxs = axs, yaxs = axs
  )
  # add 1 to 1 line
  graphics::abline(0, 1, col = 1)

  # add histograms along the sides if requested
  # note: this system won't work
  if (hist) {
    hist.x <- graphics::hist(xvec, breaks = 0:(maxage + 1) - 0.5, plot = FALSE)
    hist.y <- graphics::hist(yvec, breaks = 0:(maxage + 1) - 0.5, plot = FALSE)
    scale.hist <- hist.frac * maxage /
      max(hist.x$counts, hist.y$counts, na.rm = TRUE)
    for (a in 0:maxage) {
      graphics::rect(
        xleft = a - 0.5,
        xright = a + 0.5,
        ybottom = 0,
        ytop = scale.hist * hist.x$counts[which(hist.x$mids == a)],
        col = col.hist,
        border = FALSE
      )
      graphics::rect(
        xleft = 0,
        xright = scale.hist * hist.y$counts[which(hist.y$mids == a)],
        ybottom = a - 0.5,
        ytop = a + 0.5,
        col = col.hist,
        border = FALSE
      )
    }
  }
  # add axes and box around the figure
  graphics::axis(1)
  graphics::axis(2)
  graphics::grid()
  # add lines at 0 if axes have padding around 0
  if (axs == "r") {
    graphics::rect(
      xleft = 0, ybottom = 0,
      xright = graphics::par()$usr[2], ytop = graphics::par()$usr[4]
    )
  } else {
    graphics::box()
  }
  # add points
  graphics::points(df1[, 1:2],
    col = col.pts,
    pch = 16, cex = scale.pts * sqrt(df1[, 3])
  )
  # add counts as text
  if (counts) {
    graphics::text(df1[, 1:2],
      col = "white",
      lab = paste(df1[, 3]), cex = scale.pts / 3
    )
  }
  # add title
  graphics::title(title)
  if (png) {
    grDevices::dev.off()
  }
  invisible(df1)
}
