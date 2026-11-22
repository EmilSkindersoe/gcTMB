#' Pairwise comparisons of a parameter between growth curves
#'
#' Compares a parameter of the population curves (random effects at zero) between every pair of rows
#' of \code{newdata}, e.g. the treatment levels: the difference of the parameter on the scale of its
#' link, its standard error from the covariance of the fixed effects, a Wald test of no difference,
#' and a confidence interval. With \code{b1type="rate"} in \code{\link{gcTMB}}, b1 is the maximum
#' growth rate, and \code{pairwise(fit, newdata, "b1")} compares the maximum growth rates.
#'
#' The tests use t quantiles with small-sample inference (\code{\link{gcTMBcontrol}}), and normal
#' quantiles otherwise. By default (\code{adjust="tukey"}), the p-values and the confidence
#' intervals are adjusted by the single-step method, the generalisation of Tukey's honestly significant
#' differences to any covariance of the estimates (Hothorn, Bretz and Westfall, 2008, the default of
#' the multcomp package): they use the distribution of the largest absolute test statistic of all the
#' pairs, given the correlations of the differences, so that the intervals cover all the differences
#' simultaneously with probability \code{level}. This distribution is simulated, with \code{nsim}
#' draws, which gives the adjusted p-values to about 0.002 for the default; the simulation has its own
#' seed and leaves the random number generator as it was. The other methods of \code{\link{p.adjust}},
#' e.g. "holm", adjust the p-values only, and the confidence intervals are then not adjusted.
#' For a parameter with the log link, the differences are log ratios, and the ratios with their
#' confidence limits are given as well. Whether all the curves have the same value of the parameter
#' is tested by a likelihood ratio test of the models with and without the covariates of the
#' parameter, see \code{\link{anova.gcTMB}}.
#'
#' @param x A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param newdata Data frame. The values of the covariates in fform, one row for each curve to compare;
#'   the rows are named as in \code{\link{timeto}}.
#' @param parm Character string. The parameter to compare, "b1" to "b5"
#' @param adjust Character string. Adjustment for the number of pairs: "tukey" (default) for the
#'   single-step adjustment of the p-values and simultaneous confidence intervals, a method of
#'   \code{\link{p.adjust}} for adjusted p-values only, e.g. "holm", or "none"
#' @param level Numeric between 0 and 1. The level of the confidence intervals
#' @param nsim Integer. Number of simulations of the largest absolute test statistic for \code{adjust="tukey"}
#'
#' @return A data frame of class pairwise.gcTMB with a row for each pair: the estimated difference, its
#'   standard error, the test statistic, the p-value and the confidence limits, with the ratio and its
#'   limits for the log link.
#' @references Hothorn, T., Bretz, F. and Westfall, P. (2008). Simultaneous inference in general
#'   parametric models. Biometrical Journal 50, 346-363.
#' @examples
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' # b1 as the maximum growth rate, on the log scale: differences of log rates
#' fit <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25, b1type="rate", link=c(b1="log"))
#' pairwise(fit, newdata=data.frame(trt=levels(FabaHeights25$trt)), parm="b1")
#' @export
pairwise <- function(x, newdata, parm="b1", adjust="tukey", level=0.95, nsim=1e5){
  if (!inherits(x, "gcTMB")) stop("x must be a gcTMB object.")
  k <- match(parm, paste("b", 1:5, sep=""))
  if (length(parm) != 1 || is.na(k)) stop("parm must be one of b1, b2, b3, b4 and b5.")
  adjust <- match.arg(adjust, c("tukey", p.adjust.methods))
  # the linear function of the coefficients for each row, named as the curves of timeto
  L <- newdatacurves(x, newdata, groups=FALSE)
  nc <- length(L)
  if (nc < 2) stop("newdata must have at least two rows to compare.")
  Lk <- t(vapply(L, function(Li) Li[k, ], numeric(ncol(L[[1]]))))
  free <- which(!x$fix)
  if (all(Lk[, free, drop=FALSE] == 0)) stop(parm, " has no estimated coefficients.")
  pv <- x$start
  nb <- length(free)
  pv[free] <- x$sdrl$par.fixed[seq_len(nb)]
  vc <- fixedcov(x)[seq_len(nb), seq_len(nb), drop=FALSE]

  pairs <- utils::combn(nc, 2)
  D <- Lk[pairs[1, ], , drop=FALSE] - Lk[pairs[2, ], , drop=FALSE]
  if (all(D == 0)) stop(parm, " is the same for all rows of newdata: it has no covariates that differ between them.")
  est <- as.vector(D %*% pv)
  Df <- D[, free, drop=FALSE]
  se <- sqrt(pmax(rowSums((Df %*% vc) * Df), 0))
  stat <- est/se
  stat[se == 0] <- NA
  if (adjust == "tukey"){
    # single-step: the distribution of the largest absolute statistic, with the correlations of the
    # differences
    ok <- which(se > 0)
    V <- Df[ok, , drop=FALSE] %*% vc %*% t(Df[ok, , drop=FALSE])
    mx <- maxabsdraws(stats::cov2cor(V), x$rdf, nsim)
    p <- rep(NA_real_, length(stat))
    p[ok] <- vapply(abs(stat[ok]), function(s) mean(mx >= s), numeric(1))
    z <- unname(stats::quantile(mx, level))
  } else {
    p <- stats::p.adjust(if (is.null(x$rdf)) 2*pnorm(-abs(stat)) else 2*pt(-abs(stat), x$rdf), adjust)
    z <- zquant(x, level)
  }
  res <- data.frame(Estimate=est, Std.Error=se, statistic=stat, p.value=p, lower=est - z*se, upper=est + z*se)
  names(res)[3] <- if (is.null(x$rdf)) "z value" else "t value"
  lnk <- x$link[[parm]]
  if (identical(lnk, "log")) res[c("ratio", "ratio.lower", "ratio.upper")] <- exp(res[c("Estimate", "lower", "upper")])
  rownames(res) <- paste(names(L)[pairs[1, ]], "-", names(L)[pairs[2, ]])
  attr(res, "parm") <- parm
  attr(res, "meaning") <- parmmeaning(x, k)
  attr(res, "link") <- lnk
  attr(res, "adjust") <- adjust
  attr(res, "level") <- level
  attr(res, "critical") <- z
  attr(res, "rdf") <- x$rdf
  class(res) <- c("pairwise.gcTMB", "data.frame")
  res
}

# nsim draws of the largest absolute value of standard normal statistics with correlation matrix R
# (possibly singular, e.g. all the pairs of three or more groups), or of t statistics on df degrees of
# freedom (df not NULL), with a common scale; by its own seed, leaving the random number generator as it was
maxabsdraws <- function(R, df=NULL, nsim=1e5){
  if (exists(".Random.seed", envir=globalenv(), inherits=FALSE)){
    oldseed <- get(".Random.seed", envir=globalenv(), inherits=FALSE)
    on.exit(assign(".Random.seed", oldseed, envir=globalenv()))
  } else {
    on.exit(if (exists(".Random.seed", envir=globalenv(), inherits=FALSE)) rm(".Random.seed", envir=globalenv()))
  }
  set.seed(1)
  e <- eigen(R, symmetric=TRUE)
  keep <- e$values > 1e-10 * max(e$values)
  A <- e$vectors[, keep, drop=FALSE] %*% diag(sqrt(e$values[keep]), sum(keep))
  Z <- matrix(stats::rnorm(nsim * sum(keep)), nsim) %*% t(A)
  if (!is.null(df)) Z <- Z / sqrt(stats::rchisq(nsim, df) / df)
  apply(abs(Z), 1, max)
}

# meaning of parameter k of a fit, for printing
parmmeaning <- function(x, k){
  if (k == 1) return(if (fitrate(x)) "maximum growth rate" else "asymptote")
  # chapmanrichards and korf: the rate or coefficient b2, the shape b3 and the intercept b4
  if (hasintercept(modelindex(x$model))) return(c(NA, gcmodel(modelindex(x$model))$b2, gcmodel(modelindex(x$model))$b3, "intercept", NA)[k])
  c(NA, "location in time", "scale", "threshold", "asymmetry")[k]
}

#' @rdname pairwise
#' @param digits Integer. Number of significant digits
#' @param ... Not used.
#' @export
print.pairwise.gcTMB <- function(x, digits=3, ...){
  parm <- attr(x, "parm")
  lnk <- attr(x, "link")
  cat("Pairwise differences of ", parm, " (", attr(x, "meaning"), ")",
      if (!identical(lnk, "identity")) paste(" on the scale of its ", lnk, " link", sep=""), ":\n\n", sep="")
  d <- x
  class(d) <- "data.frame"
  d$p.value <- format.pval(d$p.value, digits=3)
  print(d, digits=digits)
  cat("\n", if (is.null(attr(x, "rdf"))) "Wald z tests" else paste("Wald t tests on", attr(x, "rdf"), "degrees of freedom"),
      "; ", switch(attr(x, "adjust"), none="unadjusted p-values and ",
                   tukey="p-values and confidence limits adjusted by the single-step (Tukey) method, ",
                   paste("p-values adjusted by the ", attr(x, "adjust"), " method, unadjusted ", sep="")),
      100*attr(x, "level"), "% confidence limits",
      if (identical(lnk, "log")) ", and ratios", ".\n", sep="")
  invisible(x)
}
