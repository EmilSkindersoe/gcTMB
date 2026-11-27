# Changelog

## gcTMB 0.1.0.9027 (2026-11-27)

- New variance functions:
  - [`varExp()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md),
    g = exp(t v), as in nlme.
  - [`varPoly()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md),
    a polynomial standard deviation of a chosen `degree` (default 2), g
    = 1 + d1 (v - m) + … + dk (v - m)^k, centred at the mean m of v in
    the data. The likelihood is zero where the polynomial is 0 or below
    at an observation, and
    [`predict()`](https://rdrr.io/r/stats/predict.html) and
    [`simulate()`](https://rdrr.io/r/stats/simulate.html) give an error
    where it is 0 or below at new values.
  - [`varComb()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md),
    the product of the g of two or more variance functions, as in nlme,
    e.g. `varComb(varPower(), varExp())` for a standard deviation
    largest in the middle of growth. Its parameters are named `A.`,
    `B.`, and so on.
- `korf`: the coefficient `b2` has the log link by default. It grows
  with t^(b3 - 1), often to 1e8 or more, where the optimiser stopped
  with “singular convergence”; on the log scale the fits converge. The
  coefficients of `b2` are now on the log scale, and random effects of
  `b2` act on log `b2`; `link = c(b2 = "identity")` gives the old
  parameterisation.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html): a warning about
  non-finite simulated values suggests the log link of `b1` for count
  models: a simulated random asymptote below zero on the identity scale
  gives a negative mean.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html): the argument
  `random` is replaced by `re`, as `random` meant the opposite of its
  meaning in [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html). `re = "new"`
  (default, the old `random = TRUE`) draws new random effects,
  `re = "fitted"` (the old `random = FALSE`) keeps the fitted ones, and
  the new `re = "none"` simulates around the population curve. Scripts
  with `simulate(random = )` must be updated; they now give an error
  naming the replacement.
- [`predict()`](https://rdrr.io/r/stats/predict.html) without `newdata`
  failed for times stored as integers (e.g. `DAS` of the package
  datasets).

## gcTMB 0.1.0.9026 (2026-11-27)

### Starting values

- A coefficient of a covariate in `fform` without a starting value
  (`NA`) is now computed from starting values found separately for each
  combination of the covariates. Given values are kept, so `0` for the
  contrasts gives the previous pooled start.
- `start` can be a list with the components `fixef`, `VarCorr` and
  `sigma`, to start the random effects and the residual standard
  deviation as well.
- Preece-Baines models have a more robust start.
- [`corCAR1()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md):
  by default, rho starts where the correlation is 0.2 at the median
  spacing of the times (before: rho = 0.2).

### `update()` and `refit()`

- [`update()`](https://rdrr.io/r/stats/update.html) now only sets
  parameters of a fit: `fixef`, `VarCorr`, `sigma`, `corpar` and
  `varpar`. The argument `rho` is renamed `corpar`, and `rho` is no
  longer accepted.
- New
  [`refit()`](https://emilskindersoe.github.io/gcTMB/reference/refit.md)
  fits the model again with other components (e.g. `fform`, `rform`,
  `correlation`, `variance`, `model`, `family`) and reuses the other
  arguments. `startfrom = "fit"` starts from the estimates of the fit
  where they carry over; `startfrom = "data"` (default) uses the usual
  starting values.

### Output

- Print methods show at most 100 rows by default (`n` in
  [`print()`](https://rdrr.io/r/base/print.html), or
  `options(gcTMB.print.max = )`), with a note on what was left out and
  how to get it.
- Variance components, variance-function parameters and correlation
  parameters are printed with two decimals.
- `plot(type = "validation2")` has one row per grouping factor, in the
  order of `rform`, and one column per parameter.
- [`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
  warns about missing or unfitted levels of the grouping variables in
  `newdata`, as
  [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
  does.

## gcTMB 0.1.0.9025 (2026-11-26)

- [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md):
  rows of `newdata` with a fitted level of a grouping variable give the
  curve of that group, with its predicted random effects, as in
  [`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md);
  new argument `random`. Levels that are not in the fit give the
  population curve with a warning (before, the grouping column was
  silently ignored).
- [`confint()`](https://rdrr.io/r/stats/confint.html) for fits,
  [`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
  and
  [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
  gives Wald intervals by default; `method = "profile"` or `"boot"` for
  the others.
- `plot(type = "curves")` draws 95% confidence bands of the population
  curves.
- [`viewcurves()`](https://emilskindersoe.github.io/gcTMB/reference/viewcurves.md)
  shows numbers with the number of decimals of the data.

## gcTMB 0.1.0.9024 (2026-11-26)

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  `outliers` observations are drawn in red, with red labels, in the
  validation graphs.
- New website article “Exploring observed curves”, with the
  [`viewcurves()`](https://emilskindersoe.github.io/gcTMB/reference/viewcurves.md)
  app running in the browser.

## gcTMB 0.1.0.9023 (2026-11-26)

- [`viewcurves()`](https://emilskindersoe.github.io/gcTMB/reference/viewcurves.md):
  the argument `highlight` is renamed `outliers`, as in
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and the
  returned `$highlighted` is now `$outliers`.

## gcTMB 0.1.0.9022 (2026-11-26)

- New
  [`viewcurves()`](https://emilskindersoe.github.io/gcTMB/reference/viewcurves.md):
  a Shiny app to step through the observed curves of a data set, with
  filters, colours and line types by covariates, and growth rates; or
  plots of all individuals with the steepest increases and decreases
  marked (`output = "jpeg"`). shiny is in Suggests.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): `type` is
  now `"validation1"` (residuals, QQ plot and curves), `"validation2"`
  (QQ plots of the random effects) or `"curves"`. The observations with
  the largest absolute residuals (`outliers`, default 1) are labelled by
  their individual in models with random effects.

## gcTMB 0.1.0.9021 (2026-11-25)

- New variance functions for the `variance` argument of
  [`gcTMB()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md):
  [`varIdent()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md),
  [`varPower()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md)
  and
  [`varConstPower()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md),
  with the arguments and parameterisation of nlme, for the gaussian and
  lognormal families. The covariate can be a variable of the data,
  `fitted(.)` or the growth rate `rate(.)` of the fitted curve.
- New
  [`varpars()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md);
  `update(varpar = )`; print, summary,
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), residuals,
  simulation, prediction intervals and `timeto(residual = TRUE)` use the
  variance function.
- The default residuals are normalized residuals for the gaussian and
  lognormal families, and quantile residuals otherwise.
- New vignette “Variance functions”.

## gcTMB 0.1.0.9020 (2026-11-25)

- New spatial correlation structures
  [`corExp()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  and
  [`corGaus()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md),
  with an optional nugget and several position variables
  (`form = ~ x + y | group`) and euclidean, maximum or manhattan
  distance.
- [`corpars()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`update()`](https://rdrr.io/r/stats/update.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and normalized
  residuals support the new structures.
- `PineHeightDiameter` gains the tree positions `X` and `Y`.
- The vignette “Autocorrelation and simulation” is split into
  “Correlation structures” and “Simulation”.

## gcTMB 0.1.0.9019 (2026-11-25)

- [`summary()`](https://rdrr.io/r/base/summary.html) shows the
  unadjusted Wald test of each coefficient again;
  [`partest()`](https://emilskindersoe.github.io/gcTMB/reference/partest.md)
  gives tests of the differences between curves, adjusted for their
  number.

## gcTMB 0.1.0.9018 (2026-11-24)

- `plot(type = "validation")` adds the data with the fitted curves to
  the page; `type = "curves"` draws two panels: the individual curves
  with the data, and the population curves.

## gcTMB 0.1.0.9017 (2026-11-24)

- All column names of the package datasets start with a capital letter
  (e.g. `FabaHeights25$Height`, `BerkeleyGrowth$Age`). The values are
  unchanged; scripts using the old names must be updated.

## gcTMB 0.1.0.9016 (2026-11-24)

- `gcTMBcontrol(onerror = "warning")` returns a fit that fails with a
  warning instead of an error, e.g. in a loop over models; the error is
  kept in `fit$fiterror` and reported by print, summary and
  [`converged()`](https://emilskindersoe.github.io/gcTMB/reference/converged.md).
- [`summary()`](https://rdrr.io/r/base/summary.html) repeats the
  convergence warnings of the fit.
- `correlation` is now the 7th argument of
  [`gcTMB()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
  after `fix`, and `weights` comes after `b1type`. Calls that pass these
  arguments by position must be updated.

## gcTMB 0.1.0.9015 (2026-11-24)

- `pairwise()` is renamed
  [`partest()`](https://emilskindersoe.github.io/gcTMB/reference/partest.md).
  New argument `contrasts`: `"dunnett"` (default, each row of `newdata`
  against a reference) or `"tukey"` (all pairs), with
  `adjust = "single-step"`.
- [`summary()`](https://rdrr.io/r/base/summary.html) shows estimates and
  standard errors without tests, and no longer the Pearson dispersion
  statistic.

## gcTMB 0.1.0.9014 (2026-11-24)

- The shape `b3` of the Weibull, gamma and Korf models has the new
  default link `"logm1"`, `log(b3 - 1)`, which keeps it above 1. Its
  coefficients are now on that scale; `link = c(b3 = "identity")` gives
  the old parameterisation.

## gcTMB 0.1.0.9013 (2026-11-22)

- `chapmanrichards` is now `b4 + b1 (1 - exp(-b2 t))^b3` (Richards
  1959), with the rate `b2` and the shape `b3`.
- The `wykoff` model is replaced by `korf`,
  `b4 + b1 exp(b2/((1 - b3) t^(b3 - 1)))` (Kuzelka and Marusak 2015).
  Scripts fitting `model = "wykoff"` or `"chapmanrichards"` must be
  updated.
- The intercept `b4` of both forestry models is fixed by default at the
  breast height 1.3.

## gcTMB 0.1.0.9012 (2026-11-22)

- New datasets `FabaFlowers26` (flowers and pods per plant and date) and
  `FabaNodes26` (flowers and pods per node) from the 2026 field trial.

## gcTMB 0.1.0.9011 (2026-11-21)

- The additive parameter `b6` of 0.1.0.9010 is removed again; the models
  have the parameters `b1` to `b5`.
- New models: Wykoff (replaced by Korf in 0.1.0.9013), the beta growth
  function of Yin et al. (2003) as `"beta"` and `"beta_capped"`, and
  Preece-Baines model 1 of human growth, whose peak height velocity
  [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
  gives in closed form.
- `plot(type = "curves")` shows a single numeric covariate on a
  continuous colour scale.

## gcTMB 0.1.0.9010 (2026-11-20)

- New Chapman-Richards and beta growth models.
- New datasets `BerkeleyGrowth`, `PineHeightDiameter` and
  `SorghumMaize`.
- The x-axis need not be time; the residual plot is labelled with the x
  variable.

## gcTMB 0.1.0.9009 (2026-11-20)

- No user-visible changes from 0.1.0.9005 (2026-11-19) to 0.1.0.9009: R
  CMD check fixes on Linux, macOS, Windows and R 4.2, and checks on
  earlier versions of R.

## gcTMB 0.1.0.9004 (2026-11-19)

- Package website at <https://emilskindersoe.github.io/gcTMB/>.

## gcTMB 0.1.0.9003 (2026-11-18)

- `control = list(restart = TRUE)` optimises a fit once more from its
  estimates when it did not converge, and keeps the better optimum.
- [`converged()`](https://emilskindersoe.github.io/gcTMB/reference/converged.md)
  accepts a restarted fit that the restart improved by less than
  `lltol`.

## gcTMB 0.1.0.9002 (2026-11-18)

- Autocorrelated residuals with
  [`corAR1()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  and
  [`corCAR1()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md),
  whose parameter is called rho.
- Independent random effects with `||` in `rform`.
- Likelihood ratio tests are named `LR` in
  [`anova()`](https://rdrr.io/r/stats/anova.html).

## gcTMB 0.1.0.9001 (2026-11-16)

- Extractors,
  [`criteria()`](https://emilskindersoe.github.io/gcTMB/reference/criteria.md),
  [`update()`](https://rdrr.io/r/stats/update.html),
  [`converged()`](https://emilskindersoe.github.io/gcTMB/reference/converged.md),
  and simulation with new data.
- `control` option to skip the Hessian and the standard errors for
  faster fits.
- Plot types,
  [`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
  and
  [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md),
  `b1` as the maximum growth rate, and pairwise comparisons of curves.

## gcTMB 0.1.0.9000 (2026-11-16)

- First version of gcTMB, developed from drcTMB to fit growth curves
  only: sigmoid growth models with asymmetry and default truncation,
  [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md),
  the gaussian, lognormal, binomial, beta-binomial, beta, Poisson and
  negative binomial families, random effects, and bootstrap intervals.
