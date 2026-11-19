# gcTMB: Fitting Growth Curves Using the Template Model Builder

Fits non-linear growth curves over time by maximum likelihood using the
Template Model Builder ('TMB'), with covariates on the parameters and
random effects, which are integrated out by the Laplace approximation,
and autocorrelated residuals. The main function gcTMB() fits the
logistic (Richards), Gompertz, Hill, Weibull, normal and gamma growth
curves, each with an asymmetry parameter, and a straight line for
comparison. The curves are parameterised by the asymptote (or, with
b1type = "rate", the maximum growth rate), a time (of the maximum growth
rate or of half the asymptote), a scale, slope or shape parameter, a
threshold time before which a truncated curve is zero, and the
asymmetry, each with an identity, log, logit or inverse link, and any of
them can be fixed. The response can follow a Gaussian, lognormal,
Poisson, binomial, beta-binomial, beta or negative binomial (nbinom1,
nbinom2) distribution, with case weights. Fixed effects are formulas for
each parameter, e.g. b1 + b2 ~ treatment, and random effects are written
as in 'nlme', e.g. b1 + b2 ~ 1 \| individual: they can be correlated
across the parameters, with several terms, nested or crossed grouping
factors and random slopes. The Gaussian and lognormal residuals can be
autocorrelated within individuals, with the first-order autoregressive
structures corAR1() and continuous-time corCAR1() of 'nlme', whose
objects are accepted as well. Starting values are found by a
linearisation of the curves or a grid search (gcTMBstart()), and the
optimisation is controlled by gcTMBcontrol(). Inference covers Wald,
profile likelihood and parametric bootstrap confidence intervals
(confint(), profile(), bootstrap()), small-sample t and F inference for
Gaussian and lognormal models without random effects, likelihood ratio
tests with chi-square, F or parametric bootstrap p-values (anova()),
information criteria such as AIC, AICc, BIC and HQIC and goodness-of-fit
criteria (criteria()), and pairwise comparisons of a parameter between
curves with single-step multiplicity adjustment (pairwise()). The times
at which the curves reach a percentage of their asymptote or a given
response (timeto()), and the maximum growth rates, their timing and the
fraction of the asymptote reached then (maxrate()), are derived with
standard errors by automatic differentiation and the delta method, and
with profile likelihood and bootstrap intervals, for the population
curves, and for timeto() also for fitted individuals. Methods for gcTMB
fits include summary(), predict() with confidence and prediction
intervals of the curves or their parameters, fitted() and residuals()
(Pearson, randomised quantile and normalized), plot() of validation
plots and fitted curves, simulate() at the observed or new times and
individuals, update() to refit the model or set its parameters as the
true model of a simulation study, logLik(), coef(), vcov(), converged(),
and the extractors fixef(), ranef(), VarCorr(), sigma() and corpars(),
compatible with 'nlme' and 'lme4'. The data sets FabaHeights25 (heights
of faba bean plants under five treatments) and FabaCounts23 (pod counts
of six faba bean genotypes) illustrate the methods. Four vignettes cover
the growth models and their parameters, fixed and random effects, the
timing and maximum growth rates, and autocorrelation and simulation, and
an overview shows most of the features in one document.

## Growth models and fitting

[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
fits the growth curves: logistic (Richards), Gompertz, Hill, Weibull,
normal and gamma, each with an asymmetry b5, and a straight line. The
parameters are the asymptote b1 (or the maximum growth rate,
`b1type="rate"`), the time b2 of the maximum growth rate or of half the
asymptote, the scale, slope or shape b3, the threshold b4 of a truncated
curve and the asymmetry b5, each with a link and optionally fixed.
[`gcTMBstart`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md)
gives the starting values and
[`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md)
the settings of the optimisation and the inference.

## Effects, distributions and autocorrelation

Fixed effects of covariates on the parameters (`fform`) and random
effects (`rform`), correlated across parameters, with several terms,
nested or crossed factors and random slopes; the families gaussian,
lognormal, poisson, binomial, betabinomial, beta, nbinom1 and nbinom2;
and autocorrelated residuals within individuals,
[`corAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
and
[`corCAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md).

## Inference

- [`summary.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md),
  `coef`, `vcov`, `logLik`,
  [`converged`](https://emilskindersoe.github.io/gcTMB/reference/converged.md);

- confidence intervals:
  [`confint.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md)
  (profile likelihood, Wald, bootstrap),
  [`profile.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md),
  [`bootstrap`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md);

- tests and model choice:
  [`anova.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/anova.gcTMB.md)
  (likelihood ratio, F and parametric bootstrap tests),
  [`criteria`](https://emilskindersoe.github.io/gcTMB/reference/criteria.md)
  (information and goodness-of-fit criteria),
  [`pairwise`](https://emilskindersoe.github.io/gcTMB/reference/pairwise.md)
  (pairwise comparisons with single-step adjustment).

## Timing and growth rates

[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md),
the times to reach a percentage of the asymptote or a given response,
and
[`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md),
the maximum growth rate, its timing and the fraction of the asymptote
reached then, with standard errors and `confint` and `profile` methods.

## Predictions, diagnostics and simulation

[`predict.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
(curves and parameters, with confidence and prediction intervals),
`fitted`,
[`residuals.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
(Pearson, randomised quantile, normalized),
[`plot.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/plot.gcTMB.md)
(validation plots and curves),
[`simulate.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
(at the observed or new times and individuals), and
[`update.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/update.gcTMB.md)
(refit, or set the parameters of a true model for simulation).

## Extractors

[`fixef`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
[`ranef`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
[`VarCorr`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
`sigma`,
[`corpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
and
[`varpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
compatible with nlme and lme4 (see
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
Relation to nlme).

## Data and vignettes

The data sets
[`FabaHeights25`](https://emilskindersoe.github.io/gcTMB/reference/FabaHeights25.md)
and
[`FabaCounts23`](https://emilskindersoe.github.io/gcTMB/reference/FabaCounts23.md).
The vignettes:
[`vignette("gcTMB")`](https://emilskindersoe.github.io/gcTMB/articles/gcTMB.md),
an overview of most features;
[`vignette("growth-models")`](https://emilskindersoe.github.io/gcTMB/articles/growth-models.md),
[`vignette("fixed-random-effects")`](https://emilskindersoe.github.io/gcTMB/articles/fixed-random-effects.md),
[`vignette("timing-rates")`](https://emilskindersoe.github.io/gcTMB/articles/timing-rates.md)
and
[`vignette("autocorrelation-simulation")`](https://emilskindersoe.github.io/gcTMB/articles/autocorrelation-simulation.md).

## Author

**Maintainer**: Emil Skindersoe <e.skindersoe@gmail.com>
([ORCID](https://orcid.org/0009-0000-2172-1614))

Authors:

- Emil Skindersoe <e.skindersoe@gmail.com>
  ([ORCID](https://orcid.org/0009-0000-2172-1614))
