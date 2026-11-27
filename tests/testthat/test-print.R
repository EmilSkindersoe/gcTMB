# the print methods show at most n rows (option gcTMB.print.max, 100 by default), with a note on the rest

test_that("print methods cut long output, say so, and show all with n = Inf or the option", {
  f <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25)
  ids <- unique(FabaHeights25$Plant_ID)
  mr <- maxrate(f, newdata=data.frame(Plant_ID=ids))
  ntot <- length(coef(mr))
  # 100 rows by default
  expect_equal(gcTMB:::printrows(NULL, 250), 100L)
  expect_equal(gcTMB:::printrows(NULL, 30), 30L)
  # the option; the note says how many are left out and how to get them
  old <- options(gcTMB.print.max=50)
  out <- capture.output(print(mr))
  options(old)
  expect_true(any(grepl(paste("... ", ntot - 50, " more estimates not shown (50 of ", ntot, ")", sep=""), out, fixed=TRUE)))
  expect_true(any(grepl("print(x, n = Inf)", out, fixed=TRUE)))
  expect_true(any(grepl("options(gcTMB.print.max = Inf)", out, fixed=TRUE)))
  expect_true(any(grepl("coef(x) and vcov(x) return all of them", out, fixed=TRUE)))
  # n; nothing is left out with Inf, and no note
  expect_false(any(grepl("not shown", capture.output(print(mr, n=Inf)))))
  expect_false(any(grepl("not shown", capture.output(print(mr)))))
  expect_true(any(grepl(paste(ntot - 10, "more estimates not shown"), capture.output(print(mr, n=10)))))
  expect_error(print(mr, n=-1), "non-negative")
  # the object itself keeps everything
  expect_length(coef(mr), ntot)

  # timeto
  tt <- timeto(f, 50, newdata=data.frame(Plant_ID=ids[1:5]))
  expect_true(any(grepl("2 more times not shown", capture.output(print(tt, n=3)))))
  # the coefficients of a fit and its summary
  f2 <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ factor(Treatment), data=FabaHeights25)
  ncf <- nrow(summary(f2)$coefficients)
  expect_true(any(grepl(paste(ncf - 2, "more coefficients not shown"), capture.output(print(f2, n=2)))))
  so <- capture.output(print(summary(f2), n=2))
  expect_true(any(grepl(paste(ncf - 2, "more coefficients not shown"), so)))
  expect_true(any(grepl("summary(fit)$coefficients", so, fixed=TRUE)))
  expect_false(any(grepl("not shown", capture.output(summary(f2)))))
  # partest: the note counts the significant comparisons left out
  pt <- partest(f2, newdata=data.frame(Treatment=sort(unique(FabaHeights25$Treatment))), parm="b1", contrasts="tukey")
  expect_gt(nrow(pt), 1)
  po <- capture.output(print(pt, n=1))
  nsig <- sum(as.data.frame(pt)$p.value[-1] < 1 - attr(pt, "level"))
  expect_true(any(grepl(paste(nrow(pt) - 1, " more comparisons not shown (1 of ", nrow(pt), "), ", nsig, " of them with p <", sep=""), po, fixed=TRUE)))
  expect_true(any(grepl("as.data.frame(x) returns all of them", po, fixed=TRUE)))
  expect_s3_class(as.data.frame(pt), "data.frame")
  expect_false(inherits(as.data.frame(pt), "partest.gcTMB"))
})
