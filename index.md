# gcTMB

gcTMB fits non-linear growth curves over time by maximum likelihood with
the Template Model Builder ([TMB](https://github.com/kaskr/adcomp)). The
parameters of the curves can depend on covariates, have random effects
for repeated measurements of the same individuals, and the residuals can
be autocorrelated over time; the response can be continuous, a count or
a proportion. Derived quantities such as the time to reach 90% of the
asymptote or the maximum growth rate come with standard errors by
automatic differentiation.

## Features

- **Growth models**: general growth curves (logistic or Richards,
  Gompertz, Hill, Weibull, normal and gamma), each with an asymmetry
  parameter and optionally truncated at the start of growth, before
  which the response is zero (e.g. emergence, or the first pods); the
  beta growth function of determinate growth, e.g. of crops; for
  forestry, the Chapman-Richards and Korf curves with an intercept,
  e.g. height-diameter models of trees from breast height; for human
  growth, the Preece-Baines model of growth in height through childhood
  and puberty; and a straight line for comparison. The x-axis need not
  be time.
- **Interpretable parameters**: the asymptote (or the maximum growth
  rate, `b1type = "rate"`), the time of the maximum growth rate or of
  half the asymptote, a scale, slope or shape, the start of growth and
  the asymmetry; each with an identity, log, logit, inverse or
  log(b3 - 1) link, and any of them can be fixed.
- **Distributions**: gaussian, lognormal, poisson, binomial,
  beta-binomial, beta, and negative binomial (`nbinom1`, `nbinom2`).
- **Fixed and random effects**: formulas for each parameter,
  e.g. `fform = b1 + b2 ~ treatment` and `rform = b1 + b2 ~ 1 | plant`,
  as in [`nlme::nlme`](https://rdrr.io/pkg/nlme/man/nlme.html); random
  effects correlated across parameters, with several terms, nested or
  crossed grouping factors and random slopes, integrated out by the
  Laplace approximation.
- **Autocorrelation**: AR(1) and continuous-time AR(1) residuals within
  individuals, `correlation = corCAR1(form = ~ day | plant)`, with the
  correlation structures of nlme.
- **Inference**: profile likelihood, Wald and parametric bootstrap
  confidence intervals; likelihood ratio tests with chi-square, F or
  parametric bootstrap p-values; information criteria; tests of the
  differences between curves
  ([`partest()`](https://emilskindersoe.github.io/gcTMB/reference/partest.md)),
  against a reference or between all pairs, with simultaneous confidence
  intervals; small-sample t and F inference for gaussian models without
  random effects.
- **Derived quantities**: times to reach a percentage of the asymptote
  or a given response
  ([`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)),
  and the maximum growth rate, its timing and the fraction of the
  asymptote reached then
  ([`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)),
  with standard errors and confidence intervals.
- **Prediction and simulation**: predictions of the curves and their
  parameters with confidence and prediction intervals, validation plots,
  and simulation from the fitted model at new times and for new
  individuals, e.g. for simulation studies.

## Installation

gcTMB is installed from GitHub; it compiles its TMB model, so a C++
compiler is needed (Rtools on Windows).

``` r

pak::pak("EmilSkindersoe/gcTMB")
```

## Documentation

The website <https://emilskindersoe.github.io/gcTMB/> has all the help
pages and the vignettes as web pages.

The vignettes introduce the package with worked examples on faba bean
data:

1.  [Growth models and their
    parameters](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/growth-models.pdf):
    the growth models, the meaning of their parameters, inflection
    points and maximum growth rates, model comparison, confidence
    intervals, and simple treatment effects.
2.  [Fixed and random
    effects](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/fixed-random-effects.pdf):
    covariates and random effects, and likelihood ratio and F-tests of
    their inclusion, for heights and counts.
3.  [Timing and maximum growth
    rates](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/timing-rates.pdf):
    [`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
    and
    [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
    with treatment and genotype effects, the maximum growth rate as a
    parameter, and pairwise comparisons, for normal and count data.
4.  [Autocorrelation and
    simulation](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/autocorrelation-simulation.pdf):
    autocorrelated residuals within individuals, and simulation studies,
    including the recovery of the autocorrelation parameter.

[The
overview](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/gcTMB.pdf)
shows most of the features in one document. Each function has a help
page, and
[`?gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
describes the models in detail.

## Growth models

The main function is
[`gcTMB()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).
The response is a curve f(t) of the explanatory variable t, usually
time. The first curves below (logistic, Gompertz, Hill, Weibull, beta,
normal and gamma) are sigmoid growth curves, often used in agronomy and
ecology, e.g. for the growth of plants, crops and animals over a season.
The Chapman-Richards and Korf curves are height-diameter models, common
in forestry, with the diameter of the trees as t. The Preece-Baines
curve is a model of human growth in height through childhood and
puberty.

| `model` | f(t) | b2 | b3 | b4 |
|----|----|----|----|----|
| `"logistic"` | b1 / (1 + exp(-(t - b2)/b3)) | inflection point | scale | start of growth |
| `"gompertz"` | b1 exp(-exp(-(t - b2)/b3)) | inflection point | scale | start of growth |
| `"hill"` | b1 / (1 + (b2/t)^b3) | time of half the asymptote | slope on the log time axis | start of growth |
| `"weibull"` | b1 (1 - exp(-(b3 - 1)/b3 (t/b2)^b3)) | inflection point | shape, \> 1 | start of growth |
| `"beta"` | b1 (1 + (b3 - t)/(b3 - b2)) ((t - b4)/(b3 - b4))^((b3 - b4)/(b3 - b2)), back to 0 at 2 b3 - b2 | inflection point | end of growth | start of growth |
| `"beta_capped"` | the same, staying at b1 from b3 on | inflection point | end of growth | start of growth |
| `"normal"` | b1 Phi((t - b2)/b3) | inflection point | scale | start of growth |
| `"gamma"` | b1 Fgamma(t; b3, b2/(b3 - 1)) | inflection point | shape, \> 1 | start of growth |
| `"chapmanrichards"` | b4 + b1 (1 - exp(-b2 t))^b3 | rate | shape | intercept |
| `"korf"` | b4 + b1 exp(b2/((1 - b3) t^(b3 - 1))) | coefficient, \> 0 | shape, \> 1 | intercept |
| `"preecebaines"` | b1 - 2 (b1 - b3)/(exp(b4 (t - b2)) + exp(b5 (t - b2))) | middle of puberty | height at b2 | rate before puberty |
| `"linear"` | b1 t + b2 | intercept | not used | not used |

b1 is the asymptote approached at late times (the maximum, reached at
the end of growth b3, for the beta curves; the adult height for
Preece-Baines), and the slope of the linear model, which has no
asymptote. Phi is the standard normal distribution function, and
Fgamma(t; k, theta) the distribution function of the gamma distribution
with shape k and scale theta. The Hill, Weibull and gamma curves are for
positive times. The links of the parameters can be chosen,
e.g. `link = c(b1 = "log", b3 = "inverse")`, with the identity
(default), log, logit (b1) or inverse (the scale b3 of the logistic,
Gompertz and normal models, the rate 1/b3) links. The shape b3 of the
Weibull, gamma and Korf models, which must be larger than 1, has the
link log(b3 - 1) (`"logm1"`) by default; `link = c(b3 = "identity")`
gives the shape itself.

By default, the logistic, Gompertz, Hill, Weibull, normal and gamma
curves f(t) are truncated at the start of growth b4: shifted down to
zero at b4 and rescaled to keep the asymptote b1, b1 (f(t) -
f(b4))/(b1 - f(b4)) for t \> b4, so that the response is zero up to b4.
b4 is fixed at 0, so that the curve starts at zero at t = 0, unless
estimated with `fix = c(b4 = FALSE)`; `truncated = FALSE` gives the
curves without truncation.

Each of these six models has an asymmetry parameter b5 \> 0, fixed at 1
(the curves above) unless estimated with `fix = c(b5 = FALSE)`: the
logistic model becomes the Richards curve b1 / (1 + b5 exp(-(t -
b2)/b3))^(1/b5), the Hill model b1 / (1 + (b2/t)^(b3))b5, the Weibull,
normal and gamma models b1 (f(t)/b1)^b5, a power of the curve relative
to its asymptote, and the Gompertz model b1 (1 - (1 - exp(-exp(-(t -
b2)/b3)))^b5), a power of its survival function.

The other curves are never truncated, and only the general growth curves
have an asymmetry b5.

- **Determinate growth**: the beta growth function (Yin et al. 2003)
  grows from its start b4 to its maximum b1 at the end of growth b3.
  `"beta"` continues the curve of Yin et al. beyond b3, back to zero at
  2 b3 - b2, for responses that decline after their maximum;
  `"beta_capped"` stays at b1 from b3 on, the beta sigmoid growth
  function of Yin et al.
- **Forestry**: the Chapman-Richards curve of Richards (1959) and the
  Korf curve, as written by Kuželka and Marušák (2015), are
  height-diameter models of trees, with the diameter as t, in the
  parameters usual in forestry: the asymptote b1, the rate b2 and the
  shape b3 of Chapman-Richards (model M8 of Ciceu et al. 2020), and the
  asymptote b1, the coefficient b2 \> 0 and the shape b3 \> 1 of Korf.
  Their intercept b4 is fixed at the breast height of 1.3 m unless given
  in `start`, e.g. `start = c(b4 = 0)` for a curve from zero, or
  estimated with `fix = c(b4 = FALSE)`.
- **Human growth**: model 1 of Preece and Baines (1978) describes the
  growth in height through childhood and puberty, with the adult height
  b1, the age b2 in the middle of puberty, the height b3 at b2, and the
  rates b4 and b5 of growth before and during puberty, estimated by
  default.
  [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
  gives the peak height velocity of the pubertal growth spurt.

With `b1type = "rate"`, b1 is the maximum growth rate of the curve
instead of its asymptote, so that covariates and random effects on b1
act on the maximum growth rate.

Starting values come from a linearisation of the curve, e.g. log(y/(b1 -
y)) = (t - b2)/b3 for the logistic curve, a straight line in t (or in
log t for the Hill and Weibull curves) whose intercept and slope give b2
and b3, profiled over the asymptote b1, or from a grid search
([`gcTMBstart()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md)).

## Examples

### A single curve

``` r

library(gcTMB)
set.seed(1)
d <- data.frame(day = rep(seq(0, 24, 2), each = 3))
d$weight <- 50 / (1 + exp(-(d$day - 10) / 2.5)) + rnorm(nrow(d), 0, 1.5)
fit <- gcTMB(weight ~ day, data = d, model = "logistic")
fit
#> Call:
#> gcTMB(form = weight ~ day, data = d, model = "logistic")
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>         50.099           9.770           2.647  
#> 
#> Residual standard deviation: 1.426
```

[`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
gives the times at which the curve reaches a percentage of its
asymptote, with standard errors from automatic differentiation and
profile likelihood intervals by
[`confint()`](https://rdrr.io/r/stats/confint.html), and
[`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
the maximum growth rate and its timing:

``` r

tt <- timeto(fit, level = c(10, 50, 90))
confint(tt)
#>         2.5 %    97.5 %
#> 10%  4.192012  4.902495
#> 50%  9.664311 10.137903
#> 90% 15.044105 16.338969
maxrate(fit)
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>          Estimate Std.Error
#> tmax        9.770    0.1170
#> rate        4.850    0.1821
#> fraction    0.488    0.0022
```

### Treatment effects, random effects and autocorrelation

The data set `FabaHeights25` has the heights of 25 faba bean plants
under five treatments, measured on eight days after sowing. Treatment
effects on the asymptote and the time of the maximum growth rate,
correlated random effects of the plants on both, and continuous-time
AR(1) residuals within the plants:

``` r

FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fit1 <- gcTMB(Height ~ DAS, fform = b1 + b2 ~ trt, rform = b1 + b2 ~ 1 | Plant_ID,
              data = FabaHeights25, correlation = corCAR1(form = ~ DAS | Plant_ID))
summary(fit1)
#> Call:
#> gcTMB(form = Height ~ DAS, fform = b1 + b2 ~ trt, rform = b1 + 
#>     b2 ~ 1 | Plant_ID, data = FabaHeights25, correlation = corCAR1(form = ~DAS | 
#>     Plant_ID))
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 200, groups: Plant_ID 25
#> Log-likelihood: -480.8 (df = 16), AIC: 993.5, BIC: 1046
#> 
#> Coefficients:
#>                Estimate Std. Error
#> b1:(Intercept)  114.315      4.209
#> b1:trt0.3        -1.967      5.796
#> b1:trt0.5        12.109      5.843
#> b1:trt0.6        13.840      5.862
#> b1:trt0.8        11.427      5.863
#> b2:(Intercept)   30.099      1.157
#> b2:trt0.3        -1.480      1.609
#> b2:trt0.5         1.212      1.610
#> b2:trt0.6         1.884      1.616
#> b2:trt0.8         2.224      1.627
#> b3:(Intercept)   12.948      0.363
#> Tests of the differences between curves, with p-values and confidence intervals: partest().
#> 
#> Random effects:
#>  Groups   Name Std.Dev. Corr 
#>  Plant_ID b1   8.256         
#>           b2   1.814    -0.54
#> 
#> Residual standard deviation: 3.29
#> 
#> Correlation structure: continuous-time AR(1) in DAS within Plant_ID
#> rho: 0.9594
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 0.00024
```

Likelihood ratio tests of the treatment effects and of the
autocorrelation:

``` r

fit0 <- update(fit1, fform = NULL)
fit2 <- update(fit1, correlation = NULL)
anova(fit0, fit1)
#> Likelihood ratio tests of gcTMB models
#> 
#>      npar    AIC    BIC  logLik deviance     LR Df Pr(>LR)  
#> fit0    8 994.54 1020.9 -489.27   978.54                    
#> fit1   16 993.54 1046.3 -480.77   961.54 16.991  8  0.0302 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
anova(fit2, fit1)
#> Likelihood ratio tests of gcTMB models
#> 
#>      npar     AIC    BIC  logLik deviance     LR Df   Pr(>LR)    
#> fit2   15 1046.34 1095.8 -508.17  1016.34                        
#> fit1   16  993.54 1046.3 -480.77   961.54 54.799  1 1.335e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The curves of the plants and the population curves of the treatments:

``` r

plot(fit1, type = "curves", times = seq(0, 61, by = 0.5))
```

![Fitted growth curves of the faba bean plants and the population curves
of the treatments](reference/figures/README-curves-1.png)

The times to reach 90% of the asymptote for each treatment:

``` r

timeto(fit1, level = 90, newdata = data.frame(trt = levels(FabaHeights25$trt)), random = FALSE)
#> Times to reach the percentages of the asymptote:
#> 
#>             Estimate Std.Error
#> 90%:trt=0.1     59.9      1.36
#> 90%:trt=0.3     58.6      1.32
#> 90%:trt=0.5     61.0      1.37
#> 90%:trt=0.6     61.6      1.39
#> 90%:trt=0.8     61.9      1.40
```

### Counts

The data set `FabaCounts23` has the numbers of pods of 22 plants of six
genotypes. A Poisson model with an estimated start of growth, the time
of the first pods, and the maximum growth rate as the parameter b1, so
that the genotypes are compared by their maximum growth rates:

``` r

genotypes <- data.frame(Genotype = sort(unique(FabaCounts23$Genotype)))
fitp <- gcTMB(Pods ~ DAS, fform = b1 + b2 ~ Genotype, rform = b1 ~ 1 | Tube_ID, data = FabaCounts23,
              family = "poisson", fix = c(b4 = FALSE), b1type = "rate", link = c(b1 = "log"))
partest(fitp, newdata = genotypes)
#> Differences to the reference of b1 (maximum growth rate) on the scale of its log link:
#> 
#>                                    Estimate Std.Error z value p.value   lower
#> Genotype=Capri - Genotype=Aurora     -0.300     0.232  -1.296   0.542 -0.8801
#> Genotype=ILB938 - Genotype=Aurora    -0.457     0.231  -1.976   0.171 -1.0366
#> Genotype=Lynx - Genotype=Aurora       0.275     0.231   1.192   0.620 -0.3031
#> Genotype=Melodie - Genotype=Aurora    0.537     0.245   2.191   0.108 -0.0767
#> Genotype=Navara - Genotype=Aurora     0.177     0.232   0.762   0.900 -0.4041
#>                                    upper ratio ratio.lower ratio.upper
#> Genotype=Capri - Genotype=Aurora   0.280 0.741       0.415        1.32
#> Genotype=ILB938 - Genotype=Aurora  0.122 0.633       0.355        1.13
#> Genotype=Lynx - Genotype=Aurora    0.854 1.317       0.739        2.35
#> Genotype=Melodie - Genotype=Aurora 1.150 1.711       0.926        3.16
#> Genotype=Navara - Genotype=Aurora  0.758 1.193       0.668        2.13
#> 
#> Wald z tests; p-values and confidence limits adjusted by the single-step (Dunnett) method, 95% confidence limits, and ratios.
```

## Data sets

The package has eight data sets, from plant science, human growth,
forestry and crop science. Each is documented in its help page,
e.g. [`?FabaHeights26`](https://emilskindersoe.github.io/gcTMB/reference/FabaHeights26.md),
and used in the vignettes:

| Data set | Content | Response | Used in |
|----|----|----|----|
| `FabaHeights26` | Heights of 320 faba bean plants in a 2026 field trial: four genotypes on two soil types, irrigated or not, 14 days between 28 and 122 days after sowing | Height (cm) | [`vignette("gcTMB")`](https://emilskindersoe.github.io/gcTMB/articles/gcTMB.md) |
| `FabaFlowers26` | Flowers and pods of 96 plants of the same field trial: two genotypes on two soil types, irrigated or not, eight days between 59 and 84 days after sowing; new and running totals per plant | Counts | the example in [`?FabaFlowers26`](https://emilskindersoe.github.io/gcTMB/reference/FabaFlowers26.md) |
| `FabaNodes26` | The flowers and pods of `FabaFlowers26` per node: date and number of the flowers and of the pods of each of 879 nodes of the same field trial. Not a growth curve: for the pod set per node, to combine with `FabaHeights26` and `FabaFlowers26` | Pods per flower | the example in [`?FabaNodes26`](https://emilskindersoe.github.io/gcTMB/reference/FabaNodes26.md) |
| `FabaHeights25` | Heights of 25 faba bean plants, five at each of five treatment levels, eight days between 14 and 61 days after sowing (2025) | Height | [`vignette("growth-models")`](https://emilskindersoe.github.io/gcTMB/articles/growth-models.md), the examples above |
| `FabaCounts23` | Counts of flowers and pods of 22 faba bean plants of six genotypes, twelve days between 61 and 79 days after sowing (2023) | Counts | [`vignette("fixed-random-effects")`](https://emilskindersoe.github.io/gcTMB/articles/fixed-random-effects.md), the example above |
| `BerkeleyGrowth` | Heights and weights of 66 boys and 70 girls of the Berkeley Child Guidance Study, from birth to 21 years (from the package sitar) | Height, weight | [`vignette("timing-rates")`](https://emilskindersoe.github.io/gcTMB/articles/timing-rates.md), with the Preece-Baines curve |
| `PineHeightDiameter` | Heights and diameters at breast height of 1678 Scots pine trees on 56 sample plots in eastern Finland (from the package lmfor) | Height against diameter | [`vignette("fixed-random-effects")`](https://emilskindersoe.github.io/gcTMB/articles/fixed-random-effects.md), with the Chapman-Richards curve |
| `SorghumMaize` | Above-ground biomass of fiber sorghum, sweet sorghum and maize at two input levels in four blocks, central Greece (from the package nlraa) | Biomass | [`vignette("autocorrelation-simulation")`](https://emilskindersoe.github.io/gcTMB/articles/autocorrelation-simulation.md), with the beta curve |

The raw faba bean data and the scripts that prepare all eight data sets
are in the directory `data-raw` of the source package.

## Functions and methods

| Purpose | Functions |
|----|----|
| Fitting | [`gcTMB()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md), [`gcTMBstart()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md), [`gcTMBcontrol()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md), [`corAR1()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md), [`corCAR1()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md), [`update()`](https://rdrr.io/r/stats/update.html) |
| Results | [`summary()`](https://rdrr.io/r/base/summary.html), [`coef()`](https://rdrr.io/r/stats/coef.html), [`vcov()`](https://rdrr.io/r/stats/vcov.html), [`logLik()`](https://rdrr.io/r/stats/logLik.html), [`converged()`](https://emilskindersoe.github.io/gcTMB/reference/converged.md), [`fixef()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md), [`ranef()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md), [`VarCorr()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md), [`sigma()`](https://rdrr.io/r/stats/sigma.html), [`corpars()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md) |
| Confidence intervals | [`confint()`](https://rdrr.io/r/stats/confint.html), [`profile()`](https://rdrr.io/r/stats/profile.html), [`bootstrap()`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md) |
| Tests and model choice | [`anova()`](https://rdrr.io/r/stats/anova.html), [`criteria()`](https://emilskindersoe.github.io/gcTMB/reference/criteria.md), [`partest()`](https://emilskindersoe.github.io/gcTMB/reference/partest.md) |
| Derived quantities | [`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md), [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md) |
| Prediction and diagnostics | [`predict()`](https://rdrr.io/r/stats/predict.html), [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`plot()`](https://rdrr.io/r/graphics/plot.default.html) |
| Simulation | [`simulate()`](https://rdrr.io/r/stats/simulate.html), [`update()`](https://rdrr.io/r/stats/update.html) |

The extractors
[`fixef()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
[`ranef()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
and
[`VarCorr()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
work for gcTMB, nlme and lme4 fits alike, and the formulas of the
effects and the correlation structures follow nlme; gcTMB fits by
maximum likelihood, so its likelihoods compare with nlme fits with
`method = "ML"`. See
[`?gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
*Relation to nlme*.

## Acknowledgements

gcTMB started as a fork of
[drcTMB](https://github.com/daniel-gerhard/drcTMB) by Daniel Gerhard,
which fits dose-response curves with the Template Model Builder. It has
since been developed into a standalone package for growth curves, but
builds on the ideas and code of drcTMB.

## License

GPL (\>= 2).
