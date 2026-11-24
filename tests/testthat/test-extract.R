fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)

test_that("fixef, ranef, VarCorr and sigma agree with the fit", {
  fe <- fixef(fit)
  expect_equal(names(fe), c("b1:(Intercept)", "b2:(Intercept)", "b3:(Intercept)"))
  expect_equal(unname(fe), unname(coef(fit)[1:3]))
  re <- ranef(fit)$Plant_ID
  expect_equal(colnames(re), c("b1", "b2"))
  # the parameters of a plant are the fixed effects plus its random effects
  pp <- predict(fit, newdata=data.frame(Plant_ID=rownames(re)[1:3]), type="parameters")
  expect_equal(unname(pp[, 1:2]), unname(sweep(as.matrix(re[1:3, ]), 2, fe[1:2], "+")))
  vc <- VarCorr(fit)$Plant_ID
  expect_equal(unname(sqrt(diag(vc))), unname(exp(coef(fit)[c("log_sd:b1|Plant_ID", "log_sd:b2|Plant_ID")])))
  expect_equal(sigma(fit), unname(exp(coef(fit)["log_sigma"])))
  expect_equal(attr(VarCorr(fit), "sc"), sigma(fit))
  expect_length(corpars(fit), 0)
  expect_length(varpars(fit), 0)
})

test_that("fixef includes coefficients fixed by the user", {
  fr <- gcTMB(Height ~ DAS, data=FabaHeights25, start=c(b5=4))
  expect_equal(unname(fixef(fr)["b5:(Intercept)"]), 4)
  expect_false("b5:(Intercept)" %in% names(fixef(fit)))
})

test_that("criteria computes the information and goodness-of-fit criteria", {
  f0 <- gcTMB(Height ~ DAS, data=FabaHeights25)
  cr <- criteria(mixed=fit, fixed=f0)
  expect_equal(rownames(cr), c("mixed", "fixed"))
  expect_equal(cr$AIC, c(AIC(fit), AIC(f0)))
  expect_equal(cr$BIC, c(BIC(fit), BIC(f0)))
  r <- residuals(fit, type="response")
  expect_equal(cr$RMSE.c[1], sqrt(mean(r^2)))
  expect_equal(cr$ME.c[1], mean(r))
  # without random effects: conditional equals marginal, no conditional RTSE
  expect_equal(cr$R2.c[2], cr$R2.m[2])
  expect_true(is.na(cr$RTSE.c[2]))
  expect_equal(criteria(list(a=fit, b=f0)), cr, ignore_attr=TRUE)
  sel <- criteria(mixed=fit, fixed=f0, which=c("AIC", "R2", "RMSE.m"))
  expect_equal(names(sel), c("AIC", "R2.c", "R2.m", "RMSE.m"))
  expect_equal(sel, cr[names(sel)])
  expect_error(criteria(fit, which="nonsense"), "unknown criteria")
})

test_that("the generics work with and without nlme, for gcTMB and nlme fits", {
  skip_if_not_installed("nlme")
  fe <- gcTMB::fixef(fit)
  expect_equal(nlme::fixef(fit), fe)
  expect_equal(nlme::VarCorr(fit), gcTMB::VarCorr(fit))
  lf <- nlme::lme(Height ~ DAS, random=~1 | Plant_ID, data=FabaHeights25)
  expect_equal(gcTMB::fixef(lf), nlme::fixef(lf))
  expect_equal(gcTMB::ranef(lf), nlme::ranef(lf))
})
