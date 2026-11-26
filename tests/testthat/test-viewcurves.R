test_that("curvedata builds the individuals, sorted by x, without missing x or y", {
  d <- data.frame(t=c(2, 1, 1, 2, 3), h=c(5, 3, 1, NA, 4), plot=c(1, 1, 2, 2, 2), plant=c("a", "a", "a", "a", "a"), g=c("u", "u", "v", "v", "v"))
  cd <- gcTMB:::curvedata(d, "t", "h", c("plot", "plant"), "g")
  expect_equal(cd$id, c("1:a", "1:a", "2:a", "2:a"))
  expect_equal(cd$x, c(1, 2, 1, 3))
  expect_equal(cd$y, c(3, 5, 1, 4))
  expect_equal(cd$g, c("u", "u", "v", "v"))
  expect_equal(unique(gcTMB:::curvedata(d, "t", "h", NULL, NULL)$id), "all")
  expect_error(gcTMB:::curvedata(d, "t", "nonsense", "plot"), "columns not in data: nonsense")
  expect_error(gcTMB:::curvedata(d, "t", "g"), "must be numeric")
  expect_error(gcTMB:::curvedata(d, c("t", "h"), "h"), "single column name")
})

test_that("stepping wraps around the ends", {
  s <- c("a", "b", "c")
  expect_equal(gcTMB:::stepid(s, "a", 1), "b")
  expect_equal(gcTMB:::stepid(s, "c", 1), "a")
  expect_equal(gcTMB:::stepid(s, "a", -1), "c")
  expect_equal(gcTMB:::stepid(s, "z", 1), "a")
})

test_that("the growth rate is the difference quotient at the middle of the intervals", {
  r <- gcTMB:::growthrate(data.frame(x=c(0, 2, 2, 6), y=c(1, 5, 6, 10)))
  expect_equal(r$mid, c(1, 4))
  expect_equal(r$rate, c(2, 1))
  expect_equal(nrow(gcTMB:::growthrate(data.frame(x=1, y=1))), 0)
})

test_that("viewcurves gives a shiny app whose filters and buttons step through the individuals", {
  skip_if_not_installed("shiny")
  app <- viewcurves(FabaHeights26, x="DAS", y="Height", id="Plant_ID",
                    covariates=c("Genotype", "Treatment", "DAS"), colour="Genotype", linetype="Treatment")
  expect_s3_class(app, "shiny.appobj")
  expect_error(viewcurves(FabaHeights26, "DAS", "Height", "Plant_ID", covariates="Genotype", colour="Soil"), "colour must be one of")
  expect_error(viewcurves(FabaHeights26, "DAS", "Height", "Plant_ID", covariates="DAS", linetype="DAS"), "at most 6 values")
  ids <- levels(FabaHeights26$Plant_ID)
  shiny::testServer(app, {
    session$setInputs(filter1="All", filter2="All", filter3=range(FabaHeights26$DAS), ind=ids[1],
                      colour="Genotype", linetype="Treatment", others=TRUE, rate=TRUE, fixy=TRUE)
    expect_equal(length(shown()), length(ids))
    expect_match(output$curve$src, "^data:image/png")
    expect_match(output$growth$src, "^data:image/png")
    g1 <- FabaHeights26$Genotype[1]
    session$setInputs(filter1=g1)
    expect_equal(length(shown()), length(unique(FabaHeights26$Plant_ID[FabaHeights26$Genotype == g1])))
    # a range filter of a numeric covariate keeps the observations within it
    session$setInputs(filter3=c(40, 60))
    expect_true(all(kept()$x >= 40 & kept()$x <= 60))
  })
  # the fixed axis of the growth rate graph: the range of the growth rates of all plants, from the
  # drop of 46_5 (132 cm at DAS 83 to 77 cm at DAS 90) to the rise of 39_4 (6.375 cm/day)
  rd <- gcTMB:::allrates(gcTMB:::curvedata(FabaHeights26, "DAS", "Height", "Plant_ID"))
  expect_equal(range(rd$y), c((77 - 132)/7, 6.375))
})

test_that("the growth rate keeps the interval it belongs to", {
  r <- gcTMB:::growthrate(data.frame(x=c(0, 2, 4, 8), y=c(0, 2, 8, 12)))
  expect_equal(r$from, c(0, 2, 4))
  expect_equal(r$to, c(2, 4, 8))
})

test_that("output='jpeg' maps the first two covariates and highlights the steepest increase and decrease", {
  # four individuals: a straight line, a jump upwards (b), a drop (c), a mild rise (d)
  d <- data.frame(t=rep(1:5, 4), id=rep(c("a", "b", "c", "d"), each=5),
                  g=rep(c("g1", "g2", "g3", "g1"), each=5), trt=rep(c("t1", "t2", "t1", "t2"), each=5))
  d$h <- c(1:5, c(1, 2, 3, 10, 17), c(1, 2, 3, 0, -3), c(1, 2, 3, 5, 7))
  pdf(NULL)
  on.exit(dev.off())
  expect_message(vc <- viewcurves(d, x="t", y="h", id="id", covariates=c("trt", "g"), output="jpeg"),
                 "Largest increase between successive observations: id b .growth rate 7 from t 3 to 4.")
  expect_equal(vc$outliers$id, c("b", "c"))
  expect_equal(vc$outliers$direction, c("increase", "decrease"))
  expect_equal(vc$outliers$rate, c(7, -3))
  expect_equal(vc$outliers$from, c(3, 3))
  expect_equal(vc$outliers$to, c(4, 4))
  expect_true(all(vapply(vc[c("curves", "rate")], inherits, logical(1), "ggplot")))
  # colour by g (three values), line type by trt (two values)
  expect_equal(vc$curves$labels$colour, "g")
  expect_equal(vc$curves$labels$linetype, "trt")
  expect_equal(vc$rate$labels$colour, "g")
  # the highlighted curves in red, labelled by their id
  red <- Filter(function(l) identical(l$aes_params$colour, "red"), vc$curves$layers)
  expect_length(red, 2)
  expect_setequal(red[[2]]$data$id, c("b", "c"))
  # reported in the console, not in the graph
  expect_null(vc$curves$labels$subtitle)
  # two of each direction, where there are growth rates of that sign (b and d increase most, only c decreases); none
  expect_equal(suppressMessages(viewcurves(d, "t", "h", "id", output="jpeg", outliers=2))$outliers$id, c("b", "d", "c"))
  v0 <- viewcurves(d, "t", "h", "id", output="jpeg", outliers=0)
  expect_equal(nrow(v0$outliers), 0)
  expect_length(Filter(function(l) identical(l$aes_params$colour, "red"), v0$curves$layers), 0)
  expect_null(v0$curves$labels$colour)
  # explicit colour and line type are kept
  vx <- suppressMessages(viewcurves(d, "t", "h", "id", covariates=c("trt", "g"), colour="trt", output="jpeg"))
  expect_equal(vx$curves$labels$colour, "trt")
  expect_null(vx$curves$labels$linetype)
  expect_error(viewcurves(d, "t", "h", "id", output="jpeg", outliers=-1), "non-negative integer")
  # written to a JPEG file
  f <- tempfile(fileext=".jpg")
  suppressMessages(viewcurves(d, "t", "h", "id", covariates="g", output="jpeg", file=f))
  expect_true(file.exists(f) && file.size(f) > 0)
})

test_that("in FabaHeights26 the steepest increase is plant 39_4 and the steepest decrease plant 46_5", {
  pdf(NULL)
  on.exit(dev.off())
  vc <- suppressMessages(viewcurves(FabaHeights26, "DAS", "Height", "Plant_ID", output="jpeg"))
  expect_equal(vc$outliers$id, c("39_4", "46_5"))
  # 46_5: 132 cm at DAS 83, 77 cm at DAS 90
  expect_equal(vc$outliers$rate[2], (77 - 132)/7)
  expect_equal(unlist(vc$outliers[2, c("from", "to")], use.names=FALSE), c(83, 90))
})

test_that("the tooltips show the individual, x and y, or the interval and the growth rate", {
  pt <- data.frame(x=83, y=132, id="46_5")
  expect_equal(gcTMB:::obstext(pt, "Plant_ID", "DAS", "Height"), c("<b>Plant_ID 46_5</b>", "DAS: 83", "Height: 132"))
  r <- gcTMB:::growthrate(data.frame(x=c(83, 90), y=c(132, 77)))
  expect_equal(gcTMB:::ratetext(r, "DAS", "Height"), c("<b>DAS 83 to 90</b>", "Growth rate of Height: -7.857"))
  # names are escaped for the HTML of the tooltip
  expect_equal(gcTMB:::obstext(data.frame(x=1, y=2, id="a<b"), "id", "x", "y & z")[c(1, 3)], c("<b>id a&lt;b</b>", "y &amp; z: 2"))
})
