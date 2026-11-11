#' Diagnostic plots for gcTMB models
#'
#' Residuals against the fitted values, residuals against time, a normal QQ plot of the
#' residuals, and normal QQ plots of the predicted random effects, drawn as ggplot2 graphs on one
#' page.
#'
#' The default randomised quantile residuals are standard normal under the model for all families,
#' so the QQ plot is compared with the identity line; for the discrete families they are drawn
#' uniformly within the jump of the distribution function at each observation, and \code{seed}
#' makes the plots reproducible. For \code{type="response"}, the QQ plot line goes through the
#' quartiles. A lowess smoother (red) shows trends of the residuals against the fitted values and
#' time.
#'
#' The QQ plots of the random effects (\code{which=4}) show the predicted random effects of each
#' parameter and grouping factor. They are shrunk towards zero, so their spread is smaller than the
#' estimated standard deviation, and with few groups they say little about the normality of the
#' random effects.
#'
#' @param x A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param which Integer vector, a subset of 1:4. The plots to draw: 1 residuals against fitted values, 2 residuals against time, 3
#'   normal QQ plot of the residuals, 4 normal QQ plots of the predicted random effects (models with
#'   random effects only)
#' @param type Character string. Residual type, see \code{\link{residuals.gcTMB}}: "quantile" (default), "pearson" or
#'   "response"
#' @param random Logical. If TRUE, residuals and fitted values include the predicted random effects
#' @param ncol Integer or NULL. Number of columns of the page
#' @param seed Integer or NULL. Seed for the randomised quantile residuals
#' @param draw Logical. If FALSE, the graphs are returned without drawing them
#' @param ... Not used.
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, rform=b1 ~ 1 | plant_ID, data=FabaHeights25)
#' plot(fit, which=1:4, seed=1)
#' # a single graph, changed
#' p <- plot(fit, which=2, draw=FALSE)
#' p$time + ggplot2::labs(x="Days after sowing")
#'
#' @return a list of ggplot objects, invisibly, named fitted, time, qq and ranef, for changing or
#'   drawing single graphs
#' @export
plot.gcTMB <- function(x, which=1:3, type=c("quantile", "pearson", "response"), random=TRUE, ncol=NULL, seed=NULL, draw=TRUE, ...){
  type <- match.arg(type)
  if (!all(which %in% 1:4)) stop("which must be a subset of 1:4.")
  if (4 %in% which && length(x$reterms) == 0){
    if (length(which) == 1) stop("the model has no random effects.")
    which <- setdiff(which, 4)
  }
  if (!is.null(seed)) set.seed(seed)

  res <- residuals(x, type=type, random=random)
  d <- data.frame(fitted=fitted(x, random=random), time=as.vector(x$obj$env$data$x), res=res)
  ylab <- switch(type, quantile="Quantile residuals", pearson="Pearson residuals", response="Residuals")

  plots <- list()
  if (1 %in% which){
    plots$fitted <- residplot(d, "fitted", ylab) + ggplot2::labs(x="Fitted values", title="Residuals vs fitted")
  }
  if (2 %in% which){
    plots$time <- residplot(d, "time", ylab) + ggplot2::labs(x="Time", title="Residuals vs time")
  }
  if (3 %in% which){
    plots$qq <- qqplot(data.frame(res=res), identity=type != "response") +
      ggplot2::labs(x="Theoretical quantiles", y=ylab, title="Normal QQ plot")
  }
  if (4 %in% which){
    plots$ranef <- qqplot(ranefdata(x), identity=FALSE) +
      ggplot2::facet_wrap(~ term, scales="free_y") +
      ggplot2::labs(x="Theoretical quantiles", y="Predicted random effects", title="Random effects")
  }
  if (draw) drawplots(plots, ncol)
  invisible(plots)
}

# residuals against a variable, with a lowess smoother
residplot <- function(d, var, ylab){
  d$xv <- d[[var]]
  p <- ggplot2::ggplot(d, ggplot2::aes(x=xv, y=res)) +
    ggplot2::geom_hline(yintercept=0, colour="grey50") +
    ggplot2::geom_point(alpha=0.7)
  if (length(unique(d$xv)) >= 3){
    lw <- as.data.frame(lowess(d$xv, d$res))
    p <- p + ggplot2::geom_line(data=lw, ggplot2::aes(x=x, y=y), colour="red")
  }
  p + ggplot2::labs(y=ylab) + ggplot2::theme_bw()
}

# normal QQ plot of res, with the identity line or a line through the quartiles; by term if d has one
qqplot <- function(d, identity=TRUE){
  grp <- if (is.null(d$term)) rep(1, nrow(d)) else d$term
  d$theo <- NA
  for (g in unique(grp)){
    i <- which(grp == g)
    d$theo[i[order(d$res[i])]] <- qnorm(ppoints(length(i)))
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x=theo, y=res))
  if (identity){
    p <- p + ggplot2::geom_abline(intercept=0, slope=1, colour="red")
  } else {
    ql <- do.call(rbind, lapply(unique(grp), function(g){
      r <- d$res[grp == g]
      slope <- diff(quantile(r, c(0.25, 0.75), names=FALSE)) / diff(qnorm(c(0.25, 0.75)))
      data.frame(term=g, slope=slope, intercept=quantile(r, 0.25, names=FALSE) - slope*qnorm(0.25))
    }))
    if (is.null(d$term)) ql$term <- NULL
    p <- p + ggplot2::geom_abline(data=ql, ggplot2::aes(intercept=intercept, slope=slope), colour="red")
  }
  p + ggplot2::geom_point(alpha=0.7) + ggplot2::theme_bw()
}

# predicted random effects, one panel for each parameter (and covariate of random slopes) and term
ranefdata <- function(object){
  u <- object$sdrl$par.random
  uoff <- 0
  do.call(rbind, lapply(object$reterms, function(tm){
    nl <- length(tm$levels)
    out <- do.call(rbind, lapply(seq_along(tm$pars), function(k){
      data.frame(term=paste(relabels(tm)[k], "|", tm$group), res=u[uoff + (k - 1)*nl + seq_len(nl)])
    }))
    uoff <<- uoff + nl*length(tm$pars)
    out
  }))
}

# draw ggplot objects on one page, in ncol columns
drawplots <- function(plots, ncol=NULL){
  n <- length(plots)
  if (is.null(ncol)) ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n/ncol)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow, ncol)))
  for (i in seq_len(n)){
    print(plots[[i]], vp=grid::viewport(layout.pos.row=(i - 1) %/% ncol + 1, layout.pos.col=(i - 1) %% ncol + 1))
  }
  grid::popViewport()
}

# columns of the plot data used in the aesthetics
utils::globalVariables(c("xv", "res", "x", "y", "theo", "intercept", "slope"))
