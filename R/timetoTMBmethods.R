#' timetoTMB Methods
#'
#' @param x,object A timetoTMB object, the result of \code{\link{timeto}}.
#' @param n Integer or NULL. The maximum number of times printed; NULL (default) for the option \code{gcTMB.print.max}, 100 unless set, e.g. \code{options(gcTMB.print.max = 200)}; \code{Inf} for all. When times are left out, a note at the end says how many, and how to show or extract all of them
#' @param ... Not used.
#'
#' @return Estimated times to reach the levels of the asymptote, and their covariance matrix
#'
#' @name timetoTMB-methods
#' @aliases print.timetoTMB coef.timetoTMB vcov.timetoTMB
NULL

#' @rdname timetoTMB-methods
#' @method print timetoTMB
#' @export
print.timetoTMB <- function(x, n=NULL, ...){
    ests <- x$estimates
    vc <- x$variance
    std <- sqrt(diag(vc))
    rdat <- data.frame(Estimate=ests, Std.Error=std)
    cat(if (identical(x$type, "absolute")) "Times to reach the responses:\n\n" else "Times to reach the percentages of the asymptote:\n\n")
    k <- printrows(n, nrow(rdat))
    print(rdat[seq_len(k), , drop=FALSE], digits=3)
    printomitted(k, nrow(rdat), "times", get="coef(x) and vcov(x) return all of them, and confint(x) their intervals")
    invisible(x)
}

#' @rdname timetoTMB-methods
#' @export
coef.timetoTMB <- function(object, ...){
    return(object$estimates)
}

#' @rdname timetoTMB-methods
#' @export
vcov.timetoTMB <- function(object, ...){
    return(object$variance)
}
