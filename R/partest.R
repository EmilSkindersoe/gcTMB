#' Tests of the differences of a parameter between growth curves
#'
#' Compares a parameter of the population curves (random effects at zero) between the rows of
#' \code{newdata}, e.g. the treatment levels: every row against a reference row
#' (\code{contrasts="dunnett"}, the default), or every pair of rows (\code{contrasts="tukey"}). This is
#' the function for the hypothesis tests of the differences between curves; the summary of the fit
#' (\code{\link{summary.gcTMB}}) gives only the estimates and their standard errors. For each
#' comparison it gives the
#' difference of the parameter on the scale of its link, its standard error from the covariance of the
#' fixed effects, a Wald test of no difference, and a confidence interval. With \code{b1type="rate"} in
#' \code{\link{gcTMB}}, b1 is the maximum growth rate, and \code{partest(fit, newdata, "b1")} compares
#' the maximum growth rates.
#'
#' The tests use t quantiles with small-sample inference (\code{\link{gcTMBcontrol}}), and normal
#' quantiles otherwise. The contrasts and the adjustment for the number of comparisons are separate
#' choices, as in the multcomp package (Hothorn, Bretz and Westfall, 2008): the contrasts say which
#' differences are tested, and the adjustment how their p-values and confidence intervals account for
#' the number of them.
#'
#' By default (\code{adjust="single-step"}), the p-values and the confidence intervals are adjusted by
#' the single-step method, the default adjustment of \code{summary.glht} in multcomp: they use the
#' distribution of the largest absolute test statistic of all the comparisons, given the correlations
#' of the differences, so that the intervals cover all the differences simultaneously with probability
#' \code{level}. For the comparisons with a reference, this is Dunnett's test, and for all the pairs,
#' Tukey's honestly significant differences, both generalised to any covariance of the estimates. Dunnett's
#' test controls the familywise error rate of the comparisons with the reference only, so its critical
#' value is smaller than Tukey's for the same rows. The distribution is simulated, with \code{nsim}
#' draws, which gives the adjusted p-values to about 0.002 for the default; the simulation has its own
#' seed and leaves the random number generator as it was. The other methods of \code{\link{p.adjust}},
#' e.g. "holm", adjust the p-values only, and the confidence intervals are then not adjusted; "none"
#' gives the unadjusted Wald tests.
#'
#' For a parameter with the log link, the differences are log ratios, and the ratios with their
#' confidence limits are given as well. Whether all the curves have the same value of the parameter
#' is tested by a likelihood ratio test of the models with and without the covariates of the
#' parameter, see \code{\link{anova.gcTMB}}. With treatment contrasts in fform and the levels in
#' newdata in the order of the factor, the differences to the first level are the coefficients of the
#' parameter in the summary of the fit.
#'
#' @param x A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param newdata Data frame. The values of the covariates in fform, one row for each curve to compare;
#'   the rows are named as in \code{\link{timeto}}.
#' @param parm Character string. The parameter to compare, "b1" to "b5"
#' @param contrasts Character string. The comparisons: "dunnett" (default) for each row of newdata
#'   against the reference row \code{ref}, or "tukey" for all the pairs of rows
#' @param ref Integer or character string. The reference row of newdata for \code{contrasts="dunnett"}:
#'   its position, or its name as in the row names of the result; the first row by default. Not used
#'   for \code{contrasts="tukey"}
#' @param adjust Character string. Adjustment for the number of comparisons: "single-step" (default)
#'   for the single-step adjustment of the p-values and simultaneous confidence intervals (Tukey's or
#'   Dunnett's method, by \code{contrasts}), a method of
#'   \code{\link{p.adjust}} for adjusted p-values only, e.g. "holm", or "none"
#' @param level Numeric between 0 and 1. The level of the confidence intervals
#' @param nsim Integer. Number of simulations of the largest absolute test statistic for
#'   \code{adjust="single-step"}
#'
#' @return A data frame of class partest.gcTMB with a row for each comparison, named "a - b": the
#'   estimated difference a - b, its standard error, the test statistic, the p-value and the
#'   confidence limits, with the ratio and its limits for the log link.
#' @references Dunnett, C. W. (1955). A multiple comparison procedure for comparing several treatments
#'   with a control. Journal of the American Statistical Association 50, 1096-1121.
#'
#'   Hothorn, T., Bretz, F. and Westfall, P. (2008). Simultaneous inference in general parametric
#'   models. Biometrical Journal 50, 346-363.
#' @examples
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' # b1 as the maximum growth rate, on the log scale: differences of log rates
#' fit <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25, b1type="rate", link=c(b1="log"))
#' nd <- data.frame(trt=levels(FabaHeights25$trt))
#' # each treatment against the first level, by Dunnett's method
#' partest(fit, newdata=nd, parm="b1")
#' # all the pairs, by Tukey's method
#' partest(fit, newdata=nd, parm="b1", contrasts="tukey")
#' @export
partest <- function(x, newdata, parm="b1", contrasts=c("dunnett", "tukey"), ref=1, adjust="single-step", level=0.95, nsim=1e5){
  if (!inherits(x, "gcTMB")) stop("x must be a gcTMB object.")
  k <- match(parm, paste("b", 1:5, sep=""))
  if (length(parm) != 1 || is.na(k)) stop("parm must be one of b1, b2, b3, b4 and b5.")
  contrasts <- match.arg(contrasts)
  adjust <- match.arg(adjust, c("single-step", "tukey", p.adjust.methods))
  if (adjust == "tukey") adjust <- "single-step"
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

  if (contrasts == "tukey"){
    pairs <- utils::combn(nc, 2)
  } else {
    r <- if (is.character(ref)) match(ref, names(L)) else ref
    if (length(ref) != 1 || is.na(r) || !(r %in% seq_len(nc)))
      stop("ref must be the position or the name of a row of newdata: one of ", paste(names(L), collapse=", "), ".")
    pairs <- rbind(setdiff(seq_len(nc), r), r)
  }
  D <-Lk[pairs[1, ], , drop=FALSE] - Lk[pairs[2, ], , drop=FALSE]
  if (all(D == 0)) stop(parm, " is the same for all rows of newdata: it has no covariates that differ between them.")
  est <- as.vector(D %*% pv)
  Df <- D[, free, drop=FALSE]
  se <- sqrt(pmax(rowSums((Df %*% vc) * Df), 0))
  stat <- est/se
  stat[se == 0] <- NA
  if (adjust == "single-step"){
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
  attr(res, "contrasts") <- contrasts
  attr(res, "adjust") <- adjust
  attr(res, "level") <- level
  attr(res, "critical") <- z
  attr(res, "rdf") <- x$rdf
  class(res) <- c("partest.gcTMB", "data.frame")
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

#' @rdname partest
#' @param digits Integer. Number of significant digits
#' @param ... Not used.
#' @export
print.partest.gcTMB <- function(x, digits=3, ...){
  parm <- attr(x, "parm")
  lnk <- attr(x, "link")
  dun <- identical(attr(x, "contrasts"), "dunnett")
  cat(if (dun) "Differences to the reference" else "Pairwise differences", " of ", parm, " (", attr(x, "meaning"), ")",
      if (!identical(lnk, "identity")) paste(" on the scale of its ", lnk, " link", sep=""), ":\n\n", sep="")
  d <- x
  class(d) <- "data.frame"
  d$p.value <- format.pval(d$p.value, digits=3)
  print(d, digits=digits)
  cat("\n", if (is.null(attr(x, "rdf"))) "Wald z tests" else paste("Wald t tests on", attr(x, "rdf"), "degrees of freedom"),
      "; ", switch(attr(x, "adjust"), none="unadjusted p-values and ",
                   "single-step"=paste("p-values and confidence limits adjusted by the single-step (", if (dun) "Dunnett" else "Tukey", ") method, ", sep=""),
                   paste("p-values adjusted by the ", attr(x, "adjust"), " method, unadjusted ", sep="")),
      100*attr(x, "level"), "% confidence limits",
      if (identical(lnk, "log")) ", and ratios", ".\n", sep="")
  invisible(x)
}
