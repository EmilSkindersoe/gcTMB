# The data sets of the package from other fields: their structure, and the fits of their examples

test_that("the data sets from other packages have the documented structure", {
  expect_equal(dim(BerkeleyGrowth), c(4884, 5))
  expect_equal(names(BerkeleyGrowth), c("ID", "Sex", "Age", "Height", "Weight"))
  expect_equal(nlevels(BerkeleyGrowth$ID), 136)
  expect_equal(as.vector(table(BerkeleyGrowth$Sex[!duplicated(BerkeleyGrowth$ID)])), c(66, 70))
  expect_equal(levels(BerkeleyGrowth$Sex), c("male", "female"))
  expect_equal(range(BerkeleyGrowth$Age), c(0, 21))

  expect_equal(dim(PineHeightDiameter), c(1678, 6))
  expect_equal(names(PineHeightDiameter), c("Plot", "DBH", "Height", "Ntrees", "Dmean", "Dvar"))
  expect_equal(nlevels(PineHeightDiameter$Plot), 56)
  expect_false(anyNA(PineHeightDiameter))
  expect_true(all(PineHeightDiameter$Height > 1.3))

  expect_equal(dim(SorghumMaize), c(236, 5))
  expect_equal(names(SorghumMaize), c("DOY", "Block", "Input", "Crop", "Biomass"))
  expect_equal(levels(SorghumMaize$Crop), c("Fiber sorghum", "Maize", "Sweet sorghum"))
  expect_equal(levels(SorghumMaize$Input), c("Low", "High"))
  expect_true(all(SorghumMaize$Biomass[SorghumMaize$DOY == 141] == 0))
})

test_that("the examples of the data sets converge", {
  fhd <- gcTMB(Height ~ DBH, rform=b1 ~ 1 | Plot, data=PineHeightDiameter, model="chapmanrichards", start=c(b4=1.3))
  expect_equal(fhd$opt$convergence, 0)
  # the curve starts at breast height
  expect_equal(unname(predict(fhd, newdata=data.frame(DBH=0), random=FALSE)), 1.3)

  fsm <- gcTMB(Biomass ~ DOY, fform=b1 ~ Crop * Input, data=SorghumMaize, model="beta", start=c(b4=141))
  expect_equal(fsm$opt$convergence, 0)
  expect_equal(unname(predict(fsm, newdata=data.frame(DOY=141, Crop="Maize", Input="High"))), 0)

  girls <- subset(BerkeleyGrowth, Sex == "female" & Age >= 2)
  fb <- gcTMB(Height ~ Age, data=girls, model="preecebaines")
  expect_equal(fb$opt$convergence, 0)
  # the girls are about 160-170 cm as adults, with the middle of puberty around 11-13 years and the
  # peak height velocity a little earlier
  b <- coef(fb)
  expect_true(b[["b1:(Intercept)"]] > 160 && b[["b1:(Intercept)"]] < 170)
  expect_true(b[["b2:(Intercept)"]] > 11 && b[["b2:(Intercept)"]] < 13)
  mr <- coef(maxrate(fb))
  expect_true(mr[["tmax"]] > 9 && mr[["tmax"]] < b[["b2:(Intercept)"]])
})
