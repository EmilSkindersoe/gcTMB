# the number of rows the print methods show: n, or by default the option gcTMB.print.max (100 unless
# set), at most the total
printrows <- function(n, total){
  if (is.null(n)) n <- getOption("gcTMB.print.max", 100)
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 0) stop("n must be a non-negative number, or Inf for all rows.")
  as.integer(min(total, floor(n)))
}

# the footer of a print method that left rows out: how many, and how to see or extract them all; get
# names the functions that return all of them, extra adds a note on the rows left out
printomitted <- function(shown, total, what="rows", get=NULL, extra=NULL){
  if (shown >= total) return(invisible())
  cat("\n... ", total - shown, " more ", what, " not shown (", shown, " of ", total, ")", extra, ".\n",
      "Use print(x, n = Inf) to show all of them, or options(gcTMB.print.max = Inf) for all printing",
      if (!is.null(get)) paste(";\n", get, sep=""), ".\n", sep="")
  invisible()
}
