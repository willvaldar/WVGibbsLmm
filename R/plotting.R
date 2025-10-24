#' Make a caterpillar plot of highest posterior density intervals from coda-format MCMC output
#'
#' Plots narrow and wide HPD intervals (default 50% and 95%), and shows both the posterior mean
#' (vertical dash) and the posterior medians (narrow gap; Tufte-inspired).
#' 
#' @return A data frame of all data plotted along with y axis locations
plot.hpd <- function(coda.object, wanted=colnames(coda.object),
		prob.wide=0.95,
		prob.narrow=0.50,
                xlab="HPD interval",
		names=NULL, ...)
{
	if (is.integer(wanted)) which.wanted <- wanted
	else which.wanted <- match(wanted, colnames(coda.object))
	num.wanted  <- length(which.wanted)

        d <- list(
          name = wanted,
          mu = colMeans(coda.object)[which.wanted],
          median = apply(coda.object, 2, median)[which.wanted],
          hpd.wide = coda::HPDinterval(coda.object, prob=prob.wide)[which.wanted,],
          hpd.narrow = coda::HPDinterval(coda.object, prob=prob.narrow)[which.wanted,]
        )
#	browser()
	if (!is.null(names)) d$name <- rep(names, length.out=length(wanted))
	ypos <- plot.ci(midvals = d$mu,
                        narrow.intervals = d$hpd.narrow,
                        wide.intervals = d$hpd.wide,
                        names = d$name,
                        xlab=xlab,
                        ...)
        d$ypos <- ypos
        points(d$median, ypos, pch="|", col="white")
        invisible(d)
} 

#' Make a caterpillar plot of confidence intervals
#'
#' Copied from original Diploffect.INLA distribution
#' 
plot.ci <- function(midvals, narrow.intervals, wide.intervals,
                    names = 1:length(midvals),
                    add = FALSE,
                    xlab = "Estimate", xlab.line = 2.5, xlab.cex = 1,
                    x.cex = 1, x.padj = 0, x.labels = TRUE,
                    ylab = "", y.cex = 1, yaxis = TRUE,
                    name.margin = 6, name.line = 4,
                    pch.midvals = 19, point.cex = 1,
                    col = rep("black", length(midvals)),
                    col.midvals = col,
                    include.top.axis = TRUE,
                    shift = 0,
                    type = "p",
                    use.this.lim = NULL,
                    main = "", main.cex = 1, main.line = 1,
                    ...) {
  
  nvals <- length(midvals)
  col.midvals <- rep(col.midvals, length.out=nvals)
  y.pos <- (1:nvals)-0.5
  if(!add){
    if(is.null(use.this.lim)){
      lim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
      lim <- c(-1,1) * diff(lim)*0.1 + lim
    }
    if(!is.null(use.this.lim)){
      lim <- use.this.lim
    }

    mar <- c(5, name.margin, 4, 2)+0.1
    oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
    plot(lim, c(0,nvals+0.5), type="n", axes=FALSE, ylab=ylab, xlab="", main=NA, ...)
    title(xlab=xlab, line=xlab.line, cex.main=xlab.cex)
    title(main=main, line=main.line, cex.main=main.cex)
    axis(1)
    if(include.top.axis){ axis(3, line=-1) }
    if(yaxis){
      axis(2, at=y.pos, labels=rev(names), las=1, lty=0, hadj=0, line=name.line)
    }
  }
  if("p"==type){
    for(i in 1:nvals){
      pos <- nvals-i + 0.5 + shift
      lines(wide.intervals[i,], rep(pos,2), col=col[i])
      lines(narrow.intervals[i,], rep(pos,2), lwd=3, col=col[i])
      points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i])
    }
  }
  invisible(rev(y.pos))
}

