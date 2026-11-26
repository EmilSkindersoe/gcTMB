#' Viewer of observed growth curves
#'
#' A Shiny app to step through the individuals of a data set and see the observed curve of each, with
#' the other individuals in the background, filters by covariates, colour and line type by
#' covariates, and the growth rate between successive observations (\code{output="html"}); or two
#' ggplot2 graphs of all the individuals, of the curves and of their growth rates, with the curves
#' of the largest increase and decrease between successive observations in red (\code{output="jpeg"}).
#'
#' The individuals are the combinations of the \code{id} columns, e.g. \code{id=c("Plot", "Plant")}
#' for plants numbered within plots, and are labelled by their values joined by ":". Without
#' \code{id}, all observations are one curve. The rows are sorted by \code{x} within each individual,
#' and rows with a missing \code{x} or \code{y} are left out.
#'
#' Each covariate gets a filter: a choice of one value for a character, factor or logical covariate,
#' or a numeric one with at most 12 distinct values, and a range for the other numeric covariates.
#' The filters select observations, and the previous and next buttons (or the left and right arrow
#' keys) step through the individuals with observations left. The covariates can also set the colour
#' of the curves (a continuous scale for numeric covariates with more than 12 values) or their line
#' type (covariates with at most 6 values). The covariates of the selected individual, from its first
#' observation, are shown below the title of the graph. With "Same axes for all individuals" (on by
#' default), both graphs keep their axes as the individuals change: the graph of the curves the range
#' of \code{x} and \code{y} of all observations, and the graph of the growth rate the range of the
#' growth rates of all individuals. Moving the mouse over a graph shows the observation nearest to it:
#' its individual, \code{x} and \code{y} (the selected individual first, then the others shown), or on
#' the graph of the growth rate the interval and the rate.
#'
#' The growth rate is the difference of \code{y} between successive observations divided by the
#' difference of \code{x}, drawn at the middle of each interval; observations at the same \code{x} give no
#' rate.
#'
#' The app needs the package shiny. It is returned as a \code{shiny.appobj}, which starts when it is
#' printed, e.g. by calling \code{viewcurves()} at the console, or with \code{shiny::runApp()}.
#'
#' With \code{output="jpeg"}, the two graphs are drawn side by side, or written to the JPEG file
#' \code{file}. Without \code{colour} and \code{linetype}, the first two covariates set them: the one
#' with more distinct values the colour, the other the line type (if it has at most 6 values; the
#' first of them on a tie). A single covariate sets the colour. The \code{outliers} individuals
#' with the largest positive growth rate (the steepest increase between two successive observations)
#' and the \code{outliers} individuals with the most negative growth rate (the steepest decrease),
#' each individual counted at its own steepest interval, are drawn in red in the graph of the curves
#' and labelled at their last observation, the others fainter; they are reported by a message and
#' returned. A single wrong measurement shows up as an increase followed by a decrease, or the reverse.
#'
#' @param data A data frame.
#' @param x Character string. The column of the times (the x axis).
#' @param y Character string. The column of the response (the y axis).
#' @param id Character vector or NULL. The columns that together identify an individual; NULL (default)
#'   for a single curve of all observations.
#' @param covariates Character vector or NULL. Columns for the filters and the colour and line type of
#'   the curves.
#' @param colour,linetype Character string or NULL. The covariate that sets the colour or the line type;
#'   NULL (default) for none when the app starts (both can be changed in the app), and for the first
#'   two covariates with \code{output="jpeg"}.
#' @param output Character string. "html" (default) for the Shiny app; "jpeg" for the two graphs of
#'   all the individuals.
#' @param outliers Non-negative integer. With \code{output="jpeg"}, the number of individuals drawn in
#'   red for each direction, the largest increase and the largest decrease between successive
#'   observations, so up to twice as many individuals in all; 0 for none. Default 1. The same
#'   argument of \code{\link{plot.gcTMB}} counts the observations with the largest absolute residuals.
#' @param file Character string or NULL. With \code{output="jpeg"}, the JPEG file to write the two
#'   graphs to (12 by 5 inches at 150 dpi); NULL (default) draws them.
#'
#' @return With \code{output="html"}, a \code{shiny.appobj}. With \code{output="jpeg"}, invisibly, a
#'   list of the two ggplot graphs, \code{curves} and \code{rate}, and \code{outliers}, a data frame
#'   of the individuals drawn in red: \code{id}, \code{direction} ("increase" or "decrease"),
#'   \code{rate} (the growth rate) and \code{from} and \code{to} (the interval of \code{x}).
#'
#' @examples
#' if (interactive() && requireNamespace("shiny", quietly=TRUE)) {
#'   viewcurves(FabaHeights26, x="DAS", y="Height", id="Plant_ID",
#'              covariates=c("Genotype", "Treatment", "Soil"), colour="Genotype")
#' }
#' # the curves and growth rates of all plants, coloured by genotype, line type by treatment, with
#' # the plants of the steepest increase and decrease between two measurements in red
#' vc <- viewcurves(FabaHeights26, x="DAS", y="Height", id="Plant_ID",
#'                  covariates=c("Genotype", "Treatment"), output="jpeg")
#' vc$outliers
#' @export
viewcurves <- function(data, x, y, id=NULL, covariates=NULL, colour=NULL, linetype=NULL,
                       output=c("html", "jpeg"), outliers=1, file=NULL){
  output <- match.arg(output)
  d <- curvedata(data, x, y, id, covariates)
  ids <- unique(as.character(d$id))
  numx <- function(v) is.numeric(v) && length(unique(v)) > 12
  colchoices <- c("None", covariates)
  ltchoices <- c("None", Filter(function(v) !numx(d[[v]]) && length(unique(d[[v]])) <= 6, covariates))
  idname <- if (is.null(id)) "Individual" else paste(id, collapse=":")
  if (output == "jpeg"){
    if (!is.numeric(outliers) || length(outliers) != 1 || is.na(outliers) || outliers < 0 || outliers != round(outliers)){
      stop("outliers must be a non-negative integer.", call.=FALSE)
    }
    # the first two covariates: colour by the one with more values, line type by the other
    if (is.null(colour) && is.null(linetype) && length(covariates) > 0){
      cv <- covariates[seq_len(min(2, length(covariates)))]
      nv <- vapply(cv, function(v) length(unique(d[[v]][!is.na(d[[v]])])), numeric(1))
      ordv <- order(-nv)
      colour <- cv[ordv[1]]
      if (length(cv) == 2 && cv[ordv[2]] %in% ltchoices) linetype <- cv[ordv[2]]
    }
  }
  if (!is.null(colour) && !(colour %in% colchoices)) stop("colour must be one of the covariates.", call.=FALSE)
  if (!is.null(linetype) && !(linetype %in% ltchoices)) stop("linetype must be one of the covariates with at most 6 values.", call.=FALSE)
  if (output == "jpeg") return(staticcurves(d, x, y, idname, colour, linetype, numx, outliers, file))
  if (!requireNamespace("shiny", quietly=TRUE)) stop("viewcurves() needs the package shiny.", call.=FALSE)
  if (is.null(colour)) colour <- "None"
  if (is.null(linetype)) linetype <- "None"
  xr <- range(d$x)
  yr <- range(d$y)
  # the range of the growth rates of all individuals, for the fixed axis of the growth rate graph
  rr <- allrates(d)$y
  rr <- if (length(rr) > 0) range(rr) else NULL

  # a filter for each covariate: a choice of one value, or a range
  filterinput <- function(k){
    v <- d[[covariates[k]]]
    if (numx(v)){
      shiny::sliderInput(paste0("filter", k), covariates[k], min=min(v, na.rm=TRUE), max=max(v, na.rm=TRUE),
                         value=range(v, na.rm=TRUE))
    } else {
      shiny::selectInput(paste0("filter", k), covariates[k], c("All", sort(unique(as.character(v[!is.na(v)])))))
    }
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel(paste("Observed curves of", y)),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width=3,
        lapply(seq_along(covariates), filterinput),
        shiny::selectizeInput("ind", idname, choices=ids),
        shiny::fluidRow(
          shiny::column(6, shiny::actionButton("prev", "< Previous", width="100%")),
          shiny::column(6, shiny::actionButton("nxt", "Next >", width="100%"))
        ),
        shiny::tags$p(shiny::textOutput("position"), style="margin-top: 8px; color: grey;"),
        shiny::hr(),
        if (length(covariates) > 0) shiny::selectInput("colour", "Colour by", colchoices, selected=colour),
        if (length(ltchoices) > 1) shiny::selectInput("linetype", "Line type by", ltchoices, selected=linetype),
        shiny::checkboxInput("others", "Show the other individuals", TRUE),
        shiny::checkboxInput("rate", paste("Show the growth rate of", y), FALSE),
        shiny::checkboxInput("fixy", "Same axes for all individuals", TRUE),
        shiny::hr(),
        shiny::tableOutput("info")
      ),
      shiny::mainPanel(
        width=9,
        # the graphs with a tooltip of the observation nearest to the mouse
        shiny::div(style="position: relative;",
                   shiny::plotOutput("curve", height="450px", hover=shiny::hoverOpts("curvehover", delay=50, delayType="throttle")),
                   shiny::uiOutput("curvetip")),
        shiny::conditionalPanel("input.rate",
          shiny::div(style="position: relative;",
                     shiny::plotOutput("growth", height="300px", hover=shiny::hoverOpts("growthhover", delay=50, delayType="throttle")),
                     shiny::uiOutput("growthtip"))),
        shiny::div(shiny::tableOutput("obs"), style="overflow-x: auto;")
      )
    ),
    # the left and right arrow keys step through the individuals
    shiny::tags$script(shiny::HTML("
      document.addEventListener('keydown', function(e) {
        if (['INPUT', 'TEXTAREA'].includes(document.activeElement.tagName)) return;
        if (e.key === 'ArrowLeft')  document.getElementById('prev').click();
        if (e.key === 'ArrowRight') document.getElementById('nxt').click();
      });"))
  )

  server <- function(input, output, session){
    # the observations that pass the filters
    kept <- shiny::reactive({
      keep <- rep(TRUE, nrow(d))
      for (k in seq_along(covariates)){
        f <- input[[paste0("filter", k)]]
        if (is.null(f)) next
        v <- d[[covariates[k]]]
        if (numx(v)) keep <- keep & !is.na(v) & v >= f[1] & v <= f[2]
        else if (f != "All") keep <- keep & !is.na(v) & as.character(v) == f
      }
      d[keep, , drop=FALSE]
    })
    shown <- shiny::reactive(intersect(ids, as.character(kept()$id)))

    shiny::observeEvent(shown(), {
      sel <- if (isTRUE(input$ind %in% shown())) input$ind else shown()[1]
      shiny::updateSelectizeInput(session, "ind", choices=shown(), selected=sel)
    })
    step <- function(k){
      s <- shown()
      if (length(s) > 0) shiny::updateSelectizeInput(session, "ind", selected=stepid(s, input$ind, k))
    }
    shiny::observeEvent(input$prev, step(-1))
    shiny::observeEvent(input$nxt, step(1))

    output$position <- shiny::renderText({
      s <- shown()
      if (length(s) == 0) return("No observations pass the filters.")
      paste0(idname, " ", match(input$ind, s), " of ", length(s), " (arrow keys step through them)")
    })

    # the observations of the selected individual that pass the filters, with the mapped covariates
    mapped <- function(dd){
      cv <- if (is.null(input$colour)) "None" else input$colour
      lt <- if (is.null(input$linetype)) "None" else input$linetype
      if (cv != "None") dd$colour <- if (numx(d[[cv]])) dd[[cv]] else factor(dd[[cv]], levels=sort(unique(d[[cv]])))
      if (lt != "None") dd$linetype <- factor(dd[[lt]], levels=sort(unique(d[[lt]])))
      dd
    }
    ind <- shiny::reactive({
      shiny::req(input$ind)
      dd <- kept()
      dd[dd$id == input$ind, , drop=FALSE]
    })

    output$curve <- shiny::renderPlot({
      p1 <- mapped(ind())
      shiny::req(nrow(p1) > 0)
      bg <- mapped(kept())
      bg <- bg[bg$id != input$ind, , drop=FALSE]
      hascol <- !is.null(p1$colour)
      haslt <- !is.null(p1$linetype)
      aesl <- if (hascol && haslt) ggplot2::aes(x=x, y=y, group=id, colour=colour, linetype=linetype)
              else if (hascol) ggplot2::aes(x=x, y=y, group=id, colour=colour)
              else if (haslt) ggplot2::aes(x=x, y=y, group=id, linetype=linetype)
              else ggplot2::aes(x=x, y=y, group=id)
      g <- ggplot2::ggplot(mapping=aesl)
      if (input$others && nrow(bg) > 0){
        g <- g + if (hascol) ggplot2::geom_line(data=bg, alpha=0.3) else ggplot2::geom_line(data=bg, colour="grey80")
      }
      # the selected individual: with a colour covariate, outlined in black to stand out from the others
      if (hascol){
        g <- g + ggplot2::geom_line(data=p1, colour="black", linewidth=2.6, show.legend=FALSE) +
          ggplot2::geom_line(data=p1, linewidth=1.4) +
          ggplot2::geom_point(data=p1, ggplot2::aes(fill=colour), shape=21, colour="black", size=3, show.legend=FALSE)
        if (is.numeric(p1$colour)) g <- g + ggplot2::scale_fill_viridis_c(end=0.9)
        else g <- g + ggplot2::scale_fill_discrete(drop=FALSE)
      } else {
        g <- g + ggplot2::geom_line(data=p1, linewidth=1.3, colour="steelblue") +
          ggplot2::geom_point(data=p1, size=2.5, colour="steelblue")
      }
      if (hascol && is.numeric(p1$colour)) g <- g + ggplot2::scale_colour_viridis_c(end=0.9)
      if (hascol && is.factor(p1$colour)) g <- g + ggplot2::scale_colour_discrete(drop=FALSE)
      if (haslt) g <- g + ggplot2::scale_linetype_discrete(drop=FALSE)
      sub <- if (length(covariates) > 0) paste(paste(covariates, vapply(covariates, function(v) as.character(p1[[v]][1]), ""), sep=": "), collapse=", ")
      g <- g + ggplot2::labs(title=paste(idname, input$ind), subtitle=sub, x=x, y=y,
                             colour=if (hascol) input$colour, linetype=if (haslt) input$linetype) +
        ggplot2::theme_bw(base_size=14)
      if (input$fixy) g <- g + ggplot2::coord_cartesian(xlim=xr, ylim=yr)
      g
    })

    output$growth <- shiny::renderPlot({
      r <- growthrate(ind())
      shiny::req(nrow(r) > 0)
      g <- ggplot2::ggplot(r, ggplot2::aes(x=mid, y=rate)) +
        ggplot2::geom_hline(yintercept=0, colour="grey50") +
        ggplot2::geom_line(colour="darkorange", linewidth=1) +
        ggplot2::geom_point(colour="darkorange", size=2.5) +
        ggplot2::labs(x=paste(x, "(middle of interval)"), y=paste("Growth rate of", y, "per unit of", x)) +
        ggplot2::theme_bw(base_size=14)
      if (input$fixy) g <- g + ggplot2::coord_cartesian(xlim=xr, ylim=rr)
      g
    })

    # tooltips: the observation of the selected individual nearest to the mouse, else that of the other
    # individuals shown; on the growth rate graph, the interval of the selected individual
    tip <- function(hover, lines){
      if (is.null(hover) || length(lines) == 0) return(NULL)
      shiny::div(shiny::HTML(paste(lines, collapse="<br>")),
                 style=paste0("position: absolute; pointer-events: none; z-index: 100; left: ", hover$coords_css$x + 14,
                              "px; top: ", hover$coords_css$y + 14, "px; background: rgba(255, 255, 255, 0.92); ",
                              "border: 1px solid #999; border-radius: 4px; padding: 4px 8px; font-size: 13px; white-space: nowrap;"))
    }
    output$curvetip <- shiny::renderUI({
      h <- input$curvehover
      if (is.null(h)) return(NULL)
      near <- function(dd) if (nrow(dd) == 0) dd else shiny::nearPoints(dd, h, xvar="x", yvar="y", threshold=10, maxpoints=1)
      pt <- near(ind())
      if (nrow(pt) == 0 && isTRUE(input$others)){
        dd <- kept()
        pt <- near(dd[dd$id != input$ind, , drop=FALSE])
      }
      if (nrow(pt) == 0) return(NULL)
      tip(h, obstext(pt, idname, x, y))
    })
    output$growthtip <- shiny::renderUI({
      h <- input$growthhover
      if (is.null(h)) return(NULL)
      r <- growthrate(ind())
      if (nrow(r) == 0) return(NULL)
      pt <- shiny::nearPoints(r, h, xvar="mid", yvar="rate", threshold=10, maxpoints=1)
      if (nrow(pt) == 0) return(NULL)
      tip(h, ratetext(pt, x, y))
    })

    output$info <- shiny::renderTable({
      p1 <- ind()
      shiny::req(nrow(p1) > 0)
      n <- nrow(p1)
      data.frame(` `=c("Observations", "First", "Last", "Largest"),
                 Value=c(n, paste(format(p1$y[1]), "at", format(p1$x[1])), paste(format(p1$y[n]), "at", format(p1$x[n])),
                         format(max(p1$y))), check.names=FALSE)
    })

    output$obs <- shiny::renderTable({
      p1 <- ind()
      m <- rbind(format(p1$x, trim=TRUE), format(p1$y, trim=TRUE))
      rownames(m) <- c(x, y)
      m
    }, rownames=TRUE, colnames=FALSE)
  }

  shiny::shinyApp(ui, server)
}

# the lines of the tooltip of an observation (a row of the data of viewcurves): its individual, x and y
obstext <- function(pt, idname, xname, yname){
  c(paste0("<b>", idname, " ", htmlesc(pt$id[1]), "</b>"),
    paste0(htmlesc(xname), ": ", format(pt$x[1])), paste0(htmlesc(yname), ": ", format(pt$y[1])))
}

# the lines of the tooltip of a growth rate (a row of growthrate()): its interval and the rate
ratetext <- function(pt, xname, yname){
  c(paste0("<b>", htmlesc(xname), " ", format(pt$from[1]), " to ", format(pt$to[1]), "</b>"),
    paste0("Growth rate of ", htmlesc(yname), ": ", format(pt$rate[1], digits=4)))
}

htmlesc <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", as.character(s), fixed=TRUE), fixed=TRUE), fixed=TRUE)

# the data of viewcurves: columns x, y and id (the individual), the covariates under their own names,
# without the rows with a missing x or y, sorted by x within the individuals
curvedata <- function(data, x, y, id=NULL, covariates=NULL){
  if (!is.data.frame(data)) stop("data must be a data frame.", call.=FALSE)
  one <- function(v, what) if (!is.character(v) || length(v) != 1) stop(what, " must be a single column name.", call.=FALSE)
  one(x, "x")
  one(y, "y")
  if (!is.null(id) && !is.character(id)) stop("id must be a character vector of column names or NULL.", call.=FALSE)
  if (!is.null(covariates) && !is.character(covariates)) stop("covariates must be a character vector of column names or NULL.", call.=FALSE)
  miss <- setdiff(c(x, y, id, covariates), names(data))
  if (length(miss) > 0) stop("columns not in data: ", paste(miss, collapse=", "), ".", call.=FALSE)
  if (!is.numeric(data[[x]]) || !is.numeric(data[[y]])) stop("x and y must be numeric columns.", call.=FALSE)
  bad <- intersect(covariates, c("x", "y", "id", "colour", "linetype"))
  if (length(bad) > 0) stop("covariates cannot be named ", paste(bad, collapse=", "), "; rename the columns.", call.=FALSE)
  d <- data.frame(x=data[[x]], y=data[[y]],
                  id=if (is.null(id)) "all" else do.call(paste, c(lapply(id, function(v) as.character(data[[v]])), sep=":")))
  for (v in covariates) d[[v]] <- data[[v]]
  d <- d[!is.na(d$x) & !is.na(d$y), , drop=FALSE]
  if (nrow(d) == 0) stop("no rows with both x and y observed.", call.=FALSE)
  # individuals in the order of the id columns (factor levels or sorted values), then by x
  ord <- if (is.null(id)) order(d$x) else {
    keys <- lapply(id, function(v) data[[v]][!is.na(data[[x]]) & !is.na(data[[y]])])
    do.call(order, c(keys, list(d$x)))
  }
  d <- d[ord, , drop=FALSE]
  rownames(d) <- NULL
  d
}

# the individual k steps from cur in s, wrapping around at the ends; the first if cur is not in s
stepid <- function(s, cur, k){
  i <- match(cur, s)
  s[if (is.na(i)) 1 else (i - 1 + k) %% length(s) + 1]
}

# growth rate between successive observations of one individual, at the middle of each interval, and
# the end of the interval (to)
growthrate <- function(d){
  n <- nrow(d)
  if (n < 2) return(data.frame(mid=numeric(0), rate=numeric(0), from=numeric(0), to=numeric(0)))
  dx <- diff(d$x)
  r <- data.frame(mid=(d$x[-n] + d$x[-1])/2, rate=diff(d$y)/dx, from=d$x[-n], to=d$x[-1])
  r[dx > 0, , drop=FALSE]
}

# the growth rates of all individuals, in columns x (middle of the interval), y (rate) and id, with the
# colour and linetype columns of d (if any) of each individual's first observation
allrates <- function(d){
  rd <- do.call(rbind, lapply(split(d, factor(d$id, levels=unique(d$id))), function(di){
    r <- growthrate(di)
    if (nrow(r) == 0) return(NULL)
    out <- data.frame(x=r$mid, y=r$rate, id=di$id[1])
    if (!is.null(di$colour)) out$colour <- di$colour[1]
    if (!is.null(di$linetype)) out$linetype <- di$linetype[1]
    out
  }))
  if (is.null(rd)) rd <- data.frame(x=numeric(0), y=numeric(0), id=character(0))
  rownames(rd) <- NULL
  rd
}

# the individuals with the largest increase (largest positive growth rate) and the largest decrease
# (most negative growth rate) between successive observations, n of each, each at its own extreme
largestchanges <- function(d, n){
  none <- data.frame(id=character(0), direction=character(0), rate=numeric(0), from=numeric(0), to=numeric(0))
  if (n == 0) return(none)
  ch <- do.call(rbind, lapply(split(d, factor(d$id, levels=unique(d$id))), function(di){
    r <- growthrate(di)
    if (nrow(r) == 0) return(NULL)
    rbind(data.frame(id=di$id[1], direction="increase", r[which.max(r$rate), c("rate", "from", "to")]),
          data.frame(id=di$id[1], direction="decrease", r[which.min(r$rate), c("rate", "from", "to")]))
  }))
  if (is.null(ch)) return(none)
  inc <- ch[ch$direction == "increase" & ch$rate > 0, ]
  dec <- ch[ch$direction == "decrease" & ch$rate < 0, ]
  out <- rbind(utils::head(inc[order(-inc$rate), ], n), utils::head(dec[order(dec$rate), ], n))
  rownames(out) <- NULL
  out
}

# the two graphs of output="jpeg": the curves of all individuals, with the highlighted ones in red and
# labelled, and their growth rates (in columns x and y, so that they use the same mapping)
staticcurves <- function(d, xname, yname, idname, colour, linetype, numx, outliers, file){
  hascol <- !is.null(colour)
  haslt <- !is.null(linetype)
  if (hascol) d$colour <- if (numx(d[[colour]])) d[[colour]] else factor(d[[colour]], levels=sort(unique(d[[colour]])))
  if (haslt) d$linetype <- factor(d[[linetype]], levels=sort(unique(d[[linetype]])))
  aesl <- if (hascol && haslt) ggplot2::aes(x=x, y=y, group=id, colour=colour, linetype=linetype)
          else if (hascol) ggplot2::aes(x=x, y=y, group=id, colour=colour)
          else if (haslt) ggplot2::aes(x=x, y=y, group=id, linetype=linetype)
          else ggplot2::aes(x=x, y=y, group=id)
  lines <- function(alpha) if (hascol) ggplot2::geom_line(alpha=alpha) else ggplot2::geom_line(alpha=alpha, colour="grey40")
  styled <- function(g){
    if (hascol && is.numeric(d$colour)) g <- g + ggplot2::scale_colour_viridis_c(end=0.9)
    g + ggplot2::labs(colour=colour, linetype=linetype) + ggplot2::theme_bw()
  }

  hl <- largestchanges(d, outliers)
  hd <- d[d$id %in% hl$id, , drop=FALSE]
  gc <- ggplot2::ggplot(d, aesl) + lines(if (nrow(hd) > 0) 0.35 else 0.7)
  if (nrow(hd) > 0){
    last <- hd[!duplicated(hd$id, fromLast=TRUE), , drop=FALSE]
    gc <- gc + ggplot2::geom_line(data=hd, colour="red", linewidth=1.1, show.legend=FALSE) +
      ggplot2::geom_text(data=last, ggplot2::aes(label=id), colour="red", hjust=-0.15, size=3.5, show.legend=FALSE) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0.05, 0.12)))
  }
  gc <- styled(gc) + ggplot2::labs(title=paste("Observed curves of", yname), x=xname, y=yname)

  rd <- allrates(d)
  gr <- ggplot2::ggplot(rd, aesl) + ggplot2::geom_hline(yintercept=0, colour="grey50") + lines(0.5)
  gr <- styled(gr) + ggplot2::labs(title=paste("Growth rate of", yname), x=paste(xname, "(middle of interval)"),
                                   y=paste("Growth rate per unit of", xname))

  for (i in seq_len(nrow(hl))){
    message("Largest ", hl$direction[i], " between successive observations: ", idname, " ", hl$id[i], " (growth rate ",
            format(hl$rate[i], digits=3), " from ", xname, " ", format(hl$from[i]), " to ", format(hl$to[i]), ")")
  }
  plots <- list(curves=gc, rate=gr)
  if (!is.null(file)){
    grDevices::jpeg(file, width=12, height=5, units="in", res=150, quality=95)
    on.exit(grDevices::dev.off())
  }
  drawplots(plots, ncol=2)
  invisible(c(plots, list(outliers=hl)))
}

utils::globalVariables(c("mid", "rate", "id", "colour", "linetype", "label"))
