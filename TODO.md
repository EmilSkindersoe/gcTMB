# gcTMB to-do list

## First priorities

- [ ] **Autocorrelation** of the residuals within individuals over time
  (e.g. AR(1) and continuous-time AR(1) for unequally spaced times), for
  repeated measurements of the same individual.
- [ ] **Variance functions** for heteroscedastic residuals (e.g. variance
  as a power of the fitted mean, or depending on time or a covariate,
  as in `nlme::varPower` / `varExp` / `varIdent`).

## Done in 0.1.0

- [x] Truncation of the curves at a threshold b4 (`truncated = TRUE`).
- [x] Log and logit links.
- [x] Hill and Weibull growth models.
- [x] Asymmetry parameter b5 for all models (Richards curve for the logistic model).
- [x] Truncation is the default, with b4 fixed at 0 (`fix = c(b4 = FALSE)` estimates it,
  `truncated = FALSE` turns it off).
- [x] `maxrate()`: time and value of the maximum growth rate and the fraction of the
  asymptote reached then, with Wald, profile and bootstrap intervals.
- [x] Normal and gamma CDF growth models, and a linear model b1*t + b2 for comparison.
