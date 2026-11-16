# gcTMB to-do list

## First priorities

- [ ] **Autocorrelation** of the residuals within individuals over time
  (e.g. AR(1) and continuous-time AR(1) for unequally spaced times), for
  repeated measurements of the same individual.
- [ ] **Variance functions** for heteroscedastic residuals (e.g. variance
  as a power of the fitted mean, or depending on time or a covariate,
  as in `nlme::varPower` / `varExp` / `varIdent`).
- [ ] **Test of equal maximum growth rates** between groups: extend `maxrate()`
  to an ANOVA-type test of whether the maximum growth rate (and its timing) is
  the same across groups, e.g. treatments in `newdata`. Both a Wald test by the
  delta method (the covariance of the rates of the groups) and a parametric
  bootstrap test.
- [ ] **Convert the nlme project to gcTMB** (waiting for the files):
  - the files go in `nlme-project/`, and the result in `nlme-project-gcTMB/`; both
    folders are ignored by git and by the package build;
  - the result is a self-contained folder to upload to a cloud service and run in
    batches on Linux, with parallel computing (e.g. `parallel::mclapply` or one
    job per batch);
  - review the nlme code and convert it to gcTMB;
  - assess the differences between the description in the manuscript and the
    actual implementation;
  - compare the computational time of nlme and gcTMB;
  - investigate whether the time-out after 120 seconds is needed, and drop the
    resampling of starting values, which gcTMB handles itself;
  - repeat the resampling until there is a true model and at least one
    competitor, instead of giving up after 10 attempts;
  - drop `abs(curve + eps)` for the simulated heights;
  - model type 3: fix the Richards asymmetry at v = 4 (the nlme scripts never set
    lv = log(4), so they simulate with v = 0.5);
  - add an extra oracle criterion, the marginal KL: the closed-form KL between
    normal distributions of one curve's observations at the design times, with
    mean and covariance (J D J' replaced by the spread of the curves, plus
    sigma^2 R for CAR-1) from the unscented transform (2q + 1 sigma points of
    the random effects), optionally Gauss-Hermite quadrature (3^q points) as a
    check; multiply by N. No simulation (too slow), and first-order
    linearisation was too inaccurate.
- [ ] **Simulation at new time points**: extend `simulate()` so that new data can
  be generated at custom time points, not only at the times in the data.
- [ ] **Profile `gcTMB()`** to find out what takes the most time in a fit
  (data setup, starting values, `MakeADFun`, optimisation, `sdreport`).

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
