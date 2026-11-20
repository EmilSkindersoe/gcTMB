# The data sets of the package from other fields: their structure, and the fits of their examples

test_that("the data sets from other packages have the documented structure", {
  expect_equal(dim(BerkeleyGrowth), c(4884, 5))
  expect_equal(names(BerkeleyGrowth), c("id", "sex", "age", "height", "weight"))
  expect_equal(nlevels(BerkeleyGrowth$id), 136)
  expect_equal(as.vector(table(BerkeleyGrowth$sex[!duplicated(BerkeleyGrowth$id)])), c(66, 70))
  expect_equal(levels(BerkeleyGrowth$sex), c("male", "female"))
  expect_equal(range(BerkeleyGrowth$age), c(0, 21))

  expect_equal(dim(PineHeightDiameter), c(1678, 6))
  expect_equal(names(PineHeightDiameter), c("plot", "dbh", "height", "ntrees", "dmean", "dvar"))
  expect_equal(nlevels(PineHeightDiameter$plot), 56)
  expect_false(anyNA(PineHeightDiameter))
  expect_true(all(PineHeightDiameter$height > 1.3))

  expect_equal(dim(SorghumMaize), c(236, 5))
  expect_equal(names(SorghumMaize), c("DOY", "Block", "Input", "Crop", "Biomass"))
  expect_equal(levels(SorghumMaize$Crop), c("Fiber sorghum", "Maize", "Sweet sorghum"))
  expect_equal(levels(SorghumMaize$Input), c("Low", "High"))
  expect_true(all(SorghumMaize$Biomass[SorghumMaize$DOY == 141] == 0))
})

test_that("the examples of the data sets converge", {
  fhd <- gcTMB(height ~ dbh, rform=b1 ~ 1 | plot, data=PineHeightDiameter, model="chapmanrichards", start=c(b6=1.3))
  expect_equal(fhd$opt$convergence, 0)
  # the curve starts at breast height
  expect_equal(unname(predict(fhd, newdata=data.frame(dbh=0), random=FALSE)), 1.3)

  fsm <- gcTMB(Biomass ~ DOY, fform=b1 ~ Crop * Input, data=SorghumMaize, model="beta", start=c(b4=141))
  expect_equal(fsm$opt$convergence, 0)
  expect_equal(unname(predict(fsm, newdata=data.frame(DOY=141, Crop="Maize", Input="High"))), 0)

  girls <- subset(BerkeleyGrowth, sex == "female" & age >= 8)
  fb <- gcTMB(height ~ age, data=girls, truncated=FALSE, fix=c(b6=FALSE))
  expect_equal(fb$opt$convergence, 0)
  # the growth spurt of the girls is around age 10-12, above a height of 100-140 cm
  b <- coef(fb)
  expect_true(b[["b2:(Intercept)"]] > 9 && b[["b2:(Intercept)"]] < 13)
  expect_true(b[["b6:(Intercept)"]] > 100 && b[["b6:(Intercept)"]] < 140)
})
