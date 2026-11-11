# gcTMB to-do list

## First priorities

- [ ] **Autocorrelation** of the residuals within individuals over time
  (e.g. AR(1) and continuous-time AR(1) for unequally spaced times), for
  repeated measurements of the same individual.
- [ ] **Variance functions** for heteroscedastic residuals (e.g. variance
  as a power of the fitted mean, or depending on time or a covariate,
  as in `nlme::varPower` / `varExp` / `varIdent`).

## Later

- [ ] More growth models: Hill, Weibull, and Richards (logistic with an
  asymmetry parameter). When added, generalise the model description in the
  vignette, README, DESCRIPTION and help pages.
