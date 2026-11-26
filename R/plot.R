#' Plots of gcTMB models: validation and fitted curves
#'
#' Validation plots of the residuals (\code{type="validation1"}, the default): residuals against the
#' fitted values, residuals against x (e.g. time) and a normal QQ plot of the residuals, with the data
#' and the fitted curves drawn last on the same page, without a legend. Validation plots of the random
#' effects (\code{type="validation2"}, models with random effects only): normal QQ plots of the
#' predicted random effects. Fitted curves (\code{type="curves"}): for models with random effects, two
#' panels with one legend, the data with the curves of the individuals, and the population curves of
#' the fixed effects alone; for models without random effects, a single graph of the data with the
#' fitted curves. Several types are drawn on one page, e.g. \code{type=c("validation1", "validation2")};
#' with "validation1", "curves" adds nothing, as the curves already end the page. The graphs are
#' ggplot2 graphs.
#'
#' In models with random effects, the \code{outliers} observations with the largest absolute
#' residuals, the ones closest to being outliers, are labelled in the graphs of
#' \code{type="validation1"} by their individual: their level of the finest grouping factor of the
#' random effects (the one with the most levels, e.g. the plant of plants nested in plots). Models
#' without random effects have no labels.
#'
#' By default the gaussian and lognormal families show the normalized residuals, which are independent
#' standard normal under the model also with a correlation structure, and the other families show the
#' randomised quantile residuals. The randomised quantile residuals are standard normal under the model for all families,
#' so the QQ plot is compared with the identity line; for the discrete families they are drawn
#' uniformly within the jump of the distribution function at each observation, and \code{seed}
#' makes the plots reproducible. For \code{residuals="response"}, the QQ plot line goes through the
#' quartiles. The Pearson and response residuals of the families other than gaussian and lognormal
#' are not normally distributed under the model, so their QQ plot is drawn with a warning. A lowess smoother (red) shows trends of the residuals against the fitted values and
#' time.
#'
#' The simulations of \code{\link{simulate.gcTMB}} can also be used for simulation-based residuals,
#' e.g. with \code{createDHARMa()} of the package DHARMa; with \code{random=FALSE}, the simulations
#' are conditional on the predicted random effects, as the residuals of these plots.
#'
#' The QQ plots of the random effects (\code{type="validation2"}) show the predicted random effects of each
#' parameter and grouping factor. They are shrunk towards zero, so their spread is smaller than the
#' estimated standard deviation, and with few groups they say little about the normality of the
#' random effects.
#'
#' The curves plot shows the data with the fitted curves at \code{times}, by default over the range
#' of the observed times: for models with random effects, the curve of each individual (each
#' combination of the grouping variables of the random effects) with its predicted random effects,
#' and the population curves, with the random effects at zero, as thick lines for each combination
#' of the covariates of the fixed effects (\code{fform}). With \code{type="curves"}, the population
#' curves are in a panel of their own, next to the data and the curves of the individuals, and both
#' panels share the colours and the legend. Points and curves are coloured by these
#' covariates. A single numeric covariate is shown on a continuous colour scale, with the population
#' curves at its 10\%, 50\% and 90\% quantiles. The covariates of an individual are taken from its
#' first observation.
#'
#' @param x A gcTMB object, the result of \code{\link{gcTMB}}.
#' @param type Character vector, one or more of "validation1" (default) for the validation plots of
#'   the residuals with the data and the fitted curves last, without a legend; "validation2" for the
#'   normal QQ plots of the predicted random effects; "curves" for the data with the curves of the
#'   individuals and the population curves in two panels with one legend
#' @param which Integer vector, a subset of 1:3. The validation plots of \code{type="validation1"}
#'   to draw before the curves: 1 residuals against fitted values, 2 residuals against x, 3 normal QQ
#'   plot of the residuals
#' @param outliers Non-negative integer. The number of observations with the largest absolute
#'   residuals labelled by their individual in the graphs of \code{type="validation1"} (models with
#'   random effects only); 0 for none. Default 1. The same argument of \code{\link{viewcurves}} counts
#'   the individuals with the steepest increase and decrease between successive observations.
#' @param residuals Character string or NULL. Residual type of the validation plots, see
#'   \code{\link{residuals.gcTMB}}: "quantile", "pearson", "response" or "normalized" (decorrelated by the
#'   correlation structure of the fit, see \code{\link{corStruct}}). NULL (default) for "normalized" with the
#'   gaussian and lognormal families and "quantile" with the other families. For the lognormal family the
#'   quantile and pearson residuals are both the standardised residuals of log(y), so the normalized
#'   residuals differ from them only by the decorrelation and suit the scatter plots and the QQ plot alike
#' @param random Logical. If TRUE, residuals and fitted values of the validation plots include the
#'   predicted random effects
#' @param times Numeric vector or NULL. The times at which the curves are drawn, e.g.
#'   \code{seq(0, 61, by=0.1)} to draw them from time 0; NULL (default) for 200 times over the range of the
#'   observed times
#' @param ncol Integer or NULL. Number of columns of the page; NULL (default) for about as many columns as rows, the square root of the number of graphs rounded up
#' @param seed Integer or NULL. Seed for the randomised quantile residuals, which for the discrete
#'   families are random; it has no effect for the continuous families. NULL (default) sets no seed: the draws continue from the current state of the random number generator
#' @param draw Logical. If FALSE, the graphs are returned without drawing them
#' @param ... Not used.
#'
#' @examples
#' fit <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25)
#' # validation plots of the residuals, with the curves last
#' plot(fit)
#' # QQ plots of the predicted random effects
#' plot(fit, type="validation2")
#' # the data with the curves of the plants, and the population curve, from sowing
#' plot(fit, type="curves", times=seq(0, 61, by=0.1))
#' # by treatment
#' FabaHeights25$trt <- factor(FabaHeights25$Treatment)
#' fitt <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
#' plot(fitt, type="curves")
#' # a single graph, changed
#' p <- plot(fit, which=2, draw=FALSE)
#' p$time + ggplot2::labs(x="Days after sowing")
#' # the three plants with the largest residuals labelled
#' plot(fit, outliers=3)
#' # counts: the randomised quantile residuals are reproducible with a seed
#' fitp <- gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE))
#' plot(fitp, seed=1)
#'
#' @return a list of ggplot objects, invisibly, for changing or drawing single graphs: fitted, time,
#'   and qq for "validation1", ranef for "validation2", and curves for the curves plot (with
#'   \code{type="curves"} and random effects, one graph with the two panels as facets)
#' @export
plot.gcTMB <- function(x, type=c("validation1", "validation2", "curves"), which=1:3, outliers=1, residuals=NULL, random=TRUE,
                       times=NULL, ncol=NULL, seed=NULL, draw=TRUE, ...){
  type <- if (missing(type)) "validation1" else unique(match.arg(type, several.ok=TRUE))
  if (!is.numeric(outliers) || length(outliers) != 1 || is.na(outliers) || outliers < 0 || outliers != round(outliers)){
    stop("outliers must be a non-negative integer.")
  }
  # default: normalized residuals for the gaussian and lognormal families, quantile residuals otherwise
  if (is.null(residuals)) residuals <- if (x$family %in% c("gaussian", "lognormal")) "normalized" else "quantile"
  rtype <- match.arg(residuals, c("quantile", "pearson", "response", "normalized"))
  hasre <- length(x$reterms) > 0
  if ("validation2" %in% type && !hasre){
    if (length(type) == 1) stop("the model has no random effects.")
    type <- setdiff(type, "validation2")
  }
  plots <- list()
  d <- NULL
  if ("validation1" %in% type){
    if (!all(which %in% 1:3)) stop("which must be a subset of 1:3; the QQ plots of the random effects are type=\"validation2\".")
    if (!is.null(seed)) set.seed(seed)

    res <- stats::residuals(x, type=rtype, random=random)
    d <- data.frame(fitted=fitted(x, random=random), time=as.vector(x$obj$env$data$x), res=res)
    # the outliers observations with the largest absolute residuals, labelled in every graph by their
    # individual: their level of the finest grouping factor of the random effects (models with random
    # effects only)
    d$label <- NA_character_
    if (hasre && outliers > 0){
      fine <- which.max(vapply(x$reterms, function(tm) length(tm$levels), numeric(1)))
      top <- utils::head(order(abs(res), decreasing=TRUE, na.last=NA), outliers)
      d$label[top] <- obsgroups(x, fine)[top]
    }
    ylab <- switch(rtype, quantile="Quantile residuals", pearson="Pearson residuals", normalized="Normalized residuals", response="Residuals")

    if (1 %in% which){
      plots$fitted <- residplot(d, "fitted", ylab) + ggplot2::labs(x="Fitted values", title="Residuals vs fitted")
    }
    if (2 %in% which){
      xname <- all.vars(x$formula[[3]])[1]
      plots$time <- residplot(d, "time", ylab) + ggplot2::labs(x=xname, title=paste("Residuals vs", xname))
    }
    if (3 %in% which){
      # only the quantile residuals (and the standardised log residuals of the lognormal) are normal
      # under the model for families other than the gaussian
      if (rtype != "quantile" && !(x$family %in% c("gaussian", "lognormal"))){
        warning("the ", rtype, " residuals of the ", x$family, " family are not normally distributed under the model, ",
                "so the normal QQ plot can depart from the line even for a correct model; use residuals=\"quantile\".", call.=FALSE)
      }
      plots$qq <- qqplot(d[, c("res", "label")], identity=rtype != "response") +
        ggplot2::labs(x="Theoretical quantiles", y=ylab, title="Normal QQ plot")
    }
  }
  if ("validation2" %in% type){
    plots$ranef <- qqplot(ranefdata(x), identity=FALSE) +
      ggplot2::facet_wrap(~ term, scales="free_y") +
      ggplot2::labs(x="Theoretical quantiles", y="Predicted random effects", title="Random effects")
  }
  # the curves last: with the validation plots a single graph without a legend; alone, the individuals
  # with the data and the population curves side by side, with one legend; none with validation2 alone
  if ("validation1" %in% type){
    plots$curves <- curveplot(x, times, label=d$label) + ggplot2::theme(legend.position="none")
  } else if ("curves" %in% type) plots$curves <- curveplot(x, times, split=TRUE)
  if (draw){
    if (length(plots) == 1) print(plots[[1]]) else drawplots(plots, ncol)
  }
  invisible(plots)
}

# the data with the fitted curves: for models with random effects the curve of each individual (each
# combination of the grouping variables) with its predicted random effects, and the population curves
# (random effects zero) of each combination of the covariates of the fixed effects, at the times, by
# default over the range of the observed times. Points and curves are coloured by the covariates of the fixed effects.
# split=TRUE (models with random effects): two panels with one legend and colour scale, the data with the
# curves of the individuals, and the population curves alone.
# label: NULL, or a character vector over the observations, NA except for the one labelled in the graph
curveplot <- function(x, times=NULL, ngrid=200, split=FALSE, label=NULL){
  tvar <- all.vars(x$formula[[3]])
  if (length(tvar) != 1) stop("the time of the model must be a single variable for fit=TRUE.")
  mf <- model.frame(x$formula, data=x$data, na.action=na.omit)
  y <- model.response(mf)
  # binomial and beta-binomial counts cbind(successes, failures): proportions, the scale of the curves
  if (is.matrix(y)) y <- y[, 1]/rowSums(y)
  obs <- x$data[rownames(mf), , drop=FALSE]
  d <- data.frame(time=obs[[tvar]], y=as.vector(y))
  d$label <- if (!is.null(label) && length(label) == nrow(d)) label else NA_character_
  if (is.null(times)) times <- seq(min(d$time), max(d$time), length.out=ngrid)
  if (!is.numeric(times) || length(times) < 2 || anyNA(times)) stop("times must be a numeric vector of at least two times without missing values.")
  grid <- data.frame(sort(unique(times)))
  names(grid) <- tvar

  covvars <- unique(unlist(lapply(x$fterms, function(ft) all.vars(ft$terms))))
  grpvars <- unique(unlist(lapply(x$reterms, function(tm) all.vars(str2lang(tm$name)[[3]]))))
  slopevars <- unique(unlist(lapply(x$reterms, function(tm) all.vars(tm$xterms$terms))))
  # covariate combination of each row, for the colours
  covlab <- function(nd){
    if (length(covvars) == 0) return(rep("all", nrow(nd)))
    do.call(paste, c(lapply(covvars, function(v) as.character(nd[[v]])), sep=", "))
  }
  d$cov <- covlab(obs)

  # a single numeric covariate: a continuous colour scale, and the population curves at its 10%, 50%
  # and 90% quantiles
  numcov <- length(covvars) == 1 && is.numeric(obs[[covvars]])
  if (numcov){
    covlab <- function(nd) nd[[covvars]]
    d$cov <- covlab(obs)
  }

  # population curves
  if (numcov){
    pc <- data.frame(signif(unique(stats::quantile(obs[[covvars]], c(0.1, 0.5, 0.9), na.rm=TRUE, names=FALSE)), 3))
    names(pc) <- covvars
  } else if (length(covvars) > 0){
    pc <- unique(obs[, covvars, drop=FALSE])
    pc <- pc[stats::complete.cases(pc), , drop=FALSE]
  } else pc <- data.frame(row.names=1)
  pc$curve <- seq_len(nrow(pc))
  pg <- merge(pc, grid)
  pg$y <- predict(x, newdata=pg, random=FALSE)
  pg$time <- pg[[tvar]]
  pg$cov <- covlab(pg)

  # curves of the individuals, with the covariates of their first observation
  ig <- NULL
  if (length(x$reterms) > 0){
    vars <- unique(c(grpvars, covvars, slopevars))
    ind <- obs[, vars, drop=FALSE]
    key <- do.call(paste, c(lapply(grpvars, function(v) as.character(ind[[v]])), sep="\r"))
    ind <- ind[!duplicated(key), , drop=FALSE]
    ind <- ind[stats::complete.cases(ind), , drop=FALSE]
    ind$curve <- seq_len(nrow(ind))
    ig <- merge(ind, grid)
    ig$y <- predict(x, newdata=ig, random=TRUE)
    ig$time <- ig[[tvar]]
    ig$cov <- covlab(ig)
  }

  yname <- deparse(x$formula[[2]])
  split <- split && !is.null(ig)
  if (split){
    panels <- c("Individuals and data", "Population")
    d$panel <- factor(panels[1], levels=panels)
    ig$panel <- factor(panels[1], levels=panels)
    pg$panel <- factor(panels[2], levels=panels)
  }
  if (length(covvars) > 0){
    p <- ggplot2::ggplot(d, ggplot2::aes(x=time, y=y, colour=cov)) +
      ggplot2::geom_point(alpha=0.6)
    if (!is.null(ig)) p <- p + ggplot2::geom_line(data=ig, ggplot2::aes(group=curve), alpha=0.35)
    p <- p + ggplot2::geom_line(data=pg, ggplot2::aes(group=curve), linewidth=1.2) +
      ggplot2::labs(colour=paste(covvars, collapse=", "))
    if (numcov) p <- p + ggplot2::scale_colour_viridis_c(end=0.9)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(x=time, y=y)) +
      ggplot2::geom_point(alpha=0.6, colour="grey40")
    if (!is.null(ig)) p <- p + ggplot2::geom_line(data=ig, ggplot2::aes(group=curve), colour="grey50", alpha=0.49)
    p <- p + ggplot2::geom_line(data=pg, ggplot2::aes(group=curve), colour="blue", linewidth=1.2)
  }
  # the label of the observation with the largest residual, on top of the curves
  p <- p + labeltext(d, "time", "y")
  if (split) p <- p + ggplot2::facet_wrap(~ panel)
  p + ggplot2::labs(x=tvar, y=yname,
                    title=if (!is.null(ig)) "Fitted curves: individuals and population" else if (nrow(pc) == 1) "Fitted curve" else "Fitted curves") +
    ggplot2::theme_bw()
}

# residuals against a variable, with a lowess smoother
residplot <- function(d, var, ylab){
  d$xv <- d[[var]]
  p <- ggplot2::ggplot(d, ggplot2::aes(x=xv, y=res)) +
    ggplot2::geom_hline(yintercept=0, colour="grey50") +
    ggplot2::geom_point(alpha=0.7) +
    labeltext(d, "xv", "res")
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
  p + ggplot2::geom_point(alpha=0.7) + labeltext(d, "theo", "res") + ggplot2::theme_bw()
}

# the text of d$label next to its point (none if d has no label): on the side of the point towards
# the middle of the x range, so that it stays inside the panel
labeltext <- function(d, xvar, yvar){
  if (is.null(d$label) || all(is.na(d$label))) return(NULL)
  i <- !is.na(d$label)
  ld <- data.frame(x=d[[xvar]][i], y=d[[yvar]][i], label=d$label[i])
  if (!is.null(d$panel)) ld$panel <- d$panel[i]
  right <- ld$x < mean(range(d[[xvar]], na.rm=TRUE))
  ggplot2::geom_text(data=ld, ggplot2::aes(x=x, y=y, label=label), inherit.aes=FALSE,
                     hjust=ifelse(right, -0.25, 1.25), size=3.5)
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
utils::globalVariables(c("xv", "res", "x", "y", "theo", "intercept", "slope", "time", "cov", "curve", "panel"))
