#' timetoTMB Methods
#'
#' @param x,object A timetoTMB object, the result of \code{\link{timeto}}.
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
print.timetoTMB <- function(x, ...){
    ests <- x$estimates
    vc <- x$variance
    std <- sqrt(diag(vc))
    rdat <- data.frame(Estimate=ests, Std.Error=std)
    cat(if (identical(x$type, "absolute")) "Times to reach the responses:\n\n" else "Times to reach the percentages of the asymptote:\n\n")
    print(rdat, digits=3)
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
