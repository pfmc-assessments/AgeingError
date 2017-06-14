#' Plot comparison of double age readings
#'
#' Plot with circles proportional to how many double readings
#' fell in each pair of coordinates
#'
#' @param xvec vector of values from reader A
#' @param yvec vector of values from reader B
#' @param col.pts color for points
#' @param col.hist color for histograms
#' @param counts include text within each bubble showing count of values?
#' @param maxage maximum age to include in the plot (doesn't yet work well)
#' @param hist include a histogram along each axis?
#' @param hist.frac maximum value of histograms as fraction of maxage
#' @param xlab label for xvec
#' @param ylab label for yvec
#' @param title Optional title to add at top of plot
#' @author Ian G. Taylor
#' @export

ageing_comparison <- function(xvec, yvec, scale.pts=2, 
                              col.pts=grey(.1,alpha=.5),
                              col.hist=rgb(0,0,.5,alpha=.7),
                              counts=TRUE, maxage=NULL,
                              hist=TRUE, hist.frac=.1,
                              xlab="Age reader A",
                              ylab="Age reader B",
                              title=NULL){
  # get counts of each pair
  df1 <- as.data.frame(table(xvec,yvec), stringsAsFactors=FALSE)
  df1$xvec <- as.numeric(df1$xvec)
  df1$yvec <- as.numeric(df1$yvec)
  # remove rows with count of zero
  df1 <- df1[df1[,3]!=0,]
  #print(df1)
  # get axis limits
  if(is.null(maxage)){
    maxage <- ceiling(max(xvec, yvec, na.rm=TRUE))
  }
  # make empty plot
  plot(0,type='n', xlim=c(0,maxage+1), ylim=c(0,maxage+1),
       xlab=xlab, ylab=ylab, axes=F, xaxs='i', yaxs='i')
  # add 1 to 1 line
  abline(0, 1, col=1)

  # add histograms along the sides if requested
  # note: this system won't work 
  if(hist){
    hist.x <- hist(xvec, breaks=0:maxage, plot=FALSE)
    hist.y <- hist(yvec, breaks=0:maxage, plot=FALSE)
    scale.hist <- hist.frac*maxage/max(hist.x$counts, hist.y$counts, na.rm=TRUE)
    for(i in 1:maxage){
      rect(xleft=hist.x$breaks[i],               ybottom=0,
           xright=hist.x$breaks[i+1],            ytop=scale.hist*hist.x$counts[i+1],
           col=col.hist, border=FALSE)
      rect(xleft=0,                              ybottom=hist.y$breaks[i],
           xright=scale.hist*hist.y$counts[i+1], ytop=hist.y$breaks[i+1],
           col=col.hist, border=FALSE)
    }
  }
  # add axes and box around the figure
  axis(1)
  axis(2)
  grid()
  box()
  # add points
  points(df1[,1:2], col=col.pts, pch=16, cex=scale.pts*sqrt(df1[,3]))
  # add counts as text
  if(counts){
    text(df1[,1:2], col='white', lab=paste(df1[,3]), cex=scale.pts/3)
  }
  # add title
  title(title)
}
