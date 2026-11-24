# Fit a growth curve

Fits a nonlinear growth curve over time by maximum likelihood with the
Template Model Builder (TMB), optionally with covariates on the
parameters and random effects, which are integrated out by the Laplace
approximation.

## Usage

``` r
gcTMB(
  form,
  fform = NULL,
  rform = NULL,
  data,
  start = NULL,
  fix = NULL,
  correlation = NULL,
  lower = -Inf,
  upper = Inf,
  family = "gaussian",
  model = "logistic",
  truncated = TRUE,
  link = NULL,
  b1type = c("asymptote", "rate"),
  weights = NULL,
  gsl = NULL,
  gsu = NULL,
  gsn = NULL,
  control = list()
)
```

## Arguments

- form:

  Formula, with the names of the response and the explanatory variable
  x, e.g. `weight ~ day`. x is usually time, but can be any variable
  along which the response grows towards an asymptote, e.g. the diameter
  of trees in a height-diameter model `height ~ dbh`.

- fform:

  Formula, list of formulas, or NULL. Fixed effects for the model
  parameters: a formula of the form b1 + b2 ~ Treatment, or a list of
  such formulas, e.g. list(b1 ~ Treatment, b2 ~ Treatment + Site).
  Parameters not listed have a single intercept (~ 1); NULL (default)
  gives every parameter a single intercept, one curve for all the data.
  A list of five right-hand side formulas, one for each of b1-b5 (X1 = ~
  1, ..., X5 = ~ 1), is also accepted.

- rform:

  Formula, list of formulas, or NULL. Random effects: a formula of the
  form b1 + b2 ~ 1 \| ID for random intercepts, or b1 ~ x \| ID for
  random intercepts and slopes in a covariate x, or a list of such
  formulas, e.g. list(b1 + b2 ~ 1 \| ID, b1 ~ 1 \| site). The random
  effects within a formula are correlated across the listed parameters
  and their intercepts and slopes; different formulas and the levels of
  nested factors (1 \| site/ID) are independent terms. With a double
  bar, all the random coefficients of the formula are independent,
  across the parameters as well as their intercepts and slopes: b1 + b2
  ~ 1 \|\| ID is list(b1 ~ 1 \| ID, b2 ~ 1 \| ID), a diagonal covariance
  matrix as random=pdDiag(b1 + b2 ~ 1) in nlme. A list of formulas
  combines the two into a block-diagonal covariance matrix: list(b1 + b2
  ~ 1 \| ID, b3 ~ 1 \| ID) gives correlated random effects of b1 and b2,
  and a random effect of b3 independent of both, as
  random=pdBlocked(list(pdSymm(b1 + b2 ~ 1), pdSymm(b3 ~ 1))) in nlme.
  The random slopes of a parameter usually go with a fixed effect of the
  covariate in fform. NULL (default) for no random effects.

- data:

  Data frame. Observations with missing values in the variables of the
  model or in the weights are left out.

- start:

  Numeric vector or NULL. Starting values for each model coefficient. NA
  entries are replaced by the starting values of
  [`gcTMBstart`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md),
  and NULL (default) uses them for all coefficients. Alternatively, a
  vector named by parameter, e.g. c(b1=10), gives the same starting
  value to all coefficients of each named parameter.

- fix:

  Logical vector or NULL. One entry for each model coefficient. If TRUE,
  a model coefficient is fixed to its value in start. Alternatively, a
  vector named by parameter, e.g. c(b1=TRUE), fixes all coefficients of
  the named parameters; their values must then be given in start, e.g.
  start=c(b1=10). NULL (default) fixes only the parameters that are
  fixed by default. The threshold (start of growth) b4 and the asymmetry
  b5 are fixed by default, at b4 = 0 (the curve starts at zero at t = 0)
  and b5 = 1 (the symmetric curve) unless given in start; so is the
  intercept b4 of the chapmanrichards and korf models, at the breast
  height of 1.3, while the rates b4 and b5 of the preecebaines model are
  estimated; fix=c(b4=FALSE) or fix=c(b5=FALSE), or an unnamed fix
  vector, estimates them. Without truncation, b4 is not used and always
  fixed, at 0 unless given in start.

- correlation:

  A correlation structure or NULL. The autocorrelation of the residuals
  within groups over time,
  [`corAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  or
  [`corCAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  (or the same structures of the package nlme), e.g.
  `corCAR1(form=~ time | ID)`, for the gaussian and lognormal (on the
  log scale) families. NULL (default) for independent residuals. See
  [`corStruct`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md);
  with a correlation structure, the weights are precision weights, with
  residual variance sigma^2/w.

- lower:

  Numeric, a single value or a vector. Lower bounds for the model
  coefficients b1-b5, a single value or a vector with the same length as
  start. Bounds of fixed coefficients are ignored; variance parameters
  are unbounded.

- upper:

  Numeric, a single value or a vector. Upper bounds for the model
  coefficients b1-b5, a single value or a vector with the same length as
  start.

- family:

  Character string. The distribution of the response, one of "gaussian",
  "lognormal", "binomial", "betabinomial", "beta", "poisson", "nbinom2",
  and "nbinom1". For the lognormal family, log(y) is normal with mean
  log f(t) and standard deviation sigma: the curve f is the median of a
  positive response, with a constant coefficient of variation; the mean
  is f(t) exp(sigma^2/2). The negative binomial families are for
  overdispersed counts, with variance mu (1 + mu/theta) for nbinom2 and
  mu (1 + phi) for nbinom1, and the beta-binomial family for
  overdispersed binomial counts, with variance n mu (1 - mu) (1 + (n -
  1)/(phi + 1)) for n trials, as in glmmTMB.

- model:

  Character string. The growth model: "logistic" (default), "gompertz",
  "hill", "weibull", "normal" or "gamma", each with an asymmetry
  parameter b5, "beta" or "beta_capped" for determinate growth,
  "chapmanrichards" or "korf" for height-diameter models in forestry,
  "preecebaines" for human growth in height, or "linear" for a straight
  line, see Details.

- truncated:

  Logical. If TRUE (default), the curve is truncated at a threshold time
  b4, by default fixed at 0: the response is zero up to b4 and then
  grows, see Details. FALSE for the curves without truncation. The beta,
  beta_capped, chapmanrichards, korf and preecebaines curves are never
  truncated.

- link:

  Named character vector or NULL. Links of the parameters b1-b5: NULL
  for the defaults, or links named by parameter, e.g.
  `link=c(b1="log", b3="inverse")`. The default is the identity link for
  all parameters, except the shape b3 of the weibull, gamma and korf
  models, which must be larger than 1 and has the link "logm1", log(b3 -
  1), by default (`link=c(b3="identity")` gives the shape itself). The
  links available for each parameter depend on its role in the model:

  - the asymptote b1: "identity", "log" or "logit" (for proportions);
    the slope b1 of the linear model: "identity" or "log";

  - the time b2 (inflection point, half the asymptote, middle of
    puberty), the rate b2 of chapmanrichards and the coefficient b2 of
    korf: "identity" or "log"; the intercept b2 of the linear model:
    "identity";

  - a scale b3 (logistic, gompertz, normal): "identity", "inverse", with
    which its coefficients are the growth rate 1/b3, or "log"; the shape
    b3 of weibull, gamma and korf: "logm1", "identity" or "log"; any
    other b3 (the slope of hill, the shape of chapmanrichards, the end
    of growth of beta, the height of preecebaines): "identity" or "log";

  - b4 (the threshold or start of growth, the intercept of
    chapmanrichards and korf, the rate of preecebaines) and the
    asymmetry b5: "identity" or "log".

  A link that a model does not allow gives an error naming the allowed
  ones. The coefficients, and start, fix, lower, upper and curves of
  [`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md),
  are on the scale of the links, see Details.

- b1type:

  Character string. The meaning of b1: "asymptote" (default), or "rate"
  for the maximum growth rate of the curve, see Details. With "rate",
  covariates and random effects on b1 are effects on the maximum growth
  rate, which
  [`anova.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/anova.gcTMB.md)
  and
  [`partest`](https://emilskindersoe.github.io/gcTMB/reference/partest.md)
  can test.

- weights:

  Numeric vector or NULL. Weights, one for each row of data; NULL
  (default) for a weight of 1 for every row. Unlike the `weights` of
  nlme, which are variance functions (`varFunc`), these are numeric case
  weights of the log-likelihood; with a correlation structure, they are
  precision weights, as `weights=varFixed(~ 1/w)` in nlme.

- gsl:

  Numeric vector or NULL. Lower bounds for a grid search, same length as
  start; NA for the default range. If gsl, gsu and gsn are all NULL
  (default), no grid search is performed and the starting values come
  from linear models (method "lm" of gcTMBstart); if any of them is
  given, NULL for the others means the default ranges and 10 grid
  points. See
  [`gcTMBstart`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md).

- gsu:

  Numeric vector or NULL. Upper bounds for a grid search, same length as
  start; NA for the default range, see gsl for NULL.

- gsn:

  Integer, a single number or a vector, or NULL. Number of grid points
  for each parameter, a single number or a vector with the same length
  as start. `gsn=10` alone runs a grid search over the default ranges;
  see gsl for NULL.

- control:

  List. Control parameters, the arguments of
  [`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md),
  e.g. `list(iter.max=2000)`. See
  [`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md)
  for details.

## Value

A gcTMB object, with methods for `print`, `summary`, `coef`, `vcov`,
`predict`, `confint`, `anova`, `plot` and others, see
[`gcTMB-methods`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md).

## Details

The logistic, gompertz, hill, weibull, normal and gamma curves are
sigmoid growth curves, often used in agronomy and ecology; they increase
from 0 towards the asymptote b1:

- logistic: f(t) = b1 / (1 + exp(-(t - b2)/b3)), symmetric around b2,
  where the curve reaches half of its asymptote;

- gompertz: f(t) = b1 exp(-exp(-(t - b2)/b3)); it reaches exp(-1), about
  37\\ asymptote at b2 and approaches the asymptote more slowly than it
  leaves zero;

- hill: f(t) = b1 / (1 + (b2/t)^b3), a logistic curve in log t
  (log-logistic), b1 / (1 + exp(-b3 (log t - log b2))); b2 is the time
  at which the curve reaches half of its asymptote, and b3 the slope on
  the log time axis;

- weibull: f(t) = b1 (1 - exp(-(b3 - 1)/b3 (t/b2)^b3)); b2 is the
  inflection point and b3 \> 1 the shape;

- normal: f(t) = b1 Phi((t - b2)/b3), with the standard normal
  distribution function Phi; symmetric around b2, where the curve
  reaches half of its asymptote;

- gamma: f(t) = b1 F(t), with F the cumulative distribution function
  (CDF) of the gamma distribution with shape b3 \> 1 and scale b2/(b3 -
  1), whose density has its mode at b2: b2 is the inflection point and
  b3 the shape;

- linear: f(t) = b1 t + b2, a straight line with slope b1 and intercept
  b2, for comparison with the growth curves. It has no asymptote, does
  not use b3-b5 and is not truncated;
  [`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
  gives the times at which it reaches given responses
  (`type="absolute"`), and
  [`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
  its constant rate b1.

Five more curves are never truncated: the beta curves of determinate
growth, the height-diameter models chapmanrichards and korf, common in
forestry, and the preecebaines model of human growth:

- chapmanrichards: f(t) = b4 + b1 (1 - exp(-b2 t))^b3, the
  Chapman-Richards curve of Richards (1959) in its usual form for
  height-diameter models (model M8 of Ciceu et al. 2020), zero above the
  intercept b4 up to t = 0: b1 is the asymptote, b2 the rate, b3 the
  shape and b4 the intercept, fixed by default at the breast height of
  1.3 m (`start=c(b4=0)` for a curve from zero, `fix=c(b4=FALSE)` to
  estimate it). For b3 \> 1 the inflection point is at log(b3)/b2, where
  the curve reaches (1 - 1/b3)^b3 of its asymptote; for b3 \<= 1 it
  grows fastest at t = 0;

- beta: the beta growth function of Yin et al. (2003), eqn 8 (with b4 =
  0 there), f(t) = b1 (1 + (b3 - t)/(b3 - b2)) ((t - b4)/(b3 -
  b4))^((b3 - b4)/(b3 - b2)) for b4 \< t \< 2 b3 - b2, and zero
  otherwise: b1 is the maximum w_max, b2 the time t_m of the maximum
  growth rate, b3 the end of growth t_e, where the curve reaches b1, and
  b4 the start of growth t_b, b4 \< b2 \< b3. After b3 the curve falls
  back to zero at 2 b3 - b2, for responses that decline after their
  maximum, e.g. biomass lost after maturity;

- beta_capped: the same curve up to b3, and b1 from b3 on, the beta
  sigmoid growth function of Yin et al. for determinate growth that
  stops at its maximum;

- korf: f(t) = b4 + b1 exp(b2/((1 - b3) t^(b3 - 1))) for t \> 0, and b4
  up to t = 0, the Korf curve: b1 is the asymptote above the intercept
  b4, b2 \> 0 a coefficient, b3 \> 1 the shape, and b4 the intercept as
  for chapmanrichards. The inflection point is at t = (b2/b3)^(1/(b3 -
  1)), where the curve reaches exp(-b3/(b3 - 1)) of its asymptote;

- preecebaines: model 1 of Preece and Baines (1978) for human growth in
  height, f(t) = b1 - 2 (b1 - b3)/(exp(b4 (t - b2)) + exp(b5 (t - b2))):
  the adult height b1, the age b2 in the middle of puberty, the height
  b3 at b2, and the rates b4 \< b5 of growth before and during puberty,
  which are estimated by default.
  [`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
  gives the peak of the pubertal growth spurt (the peak height
  velocity); `b1type="rate"` is not available.

**Relation to nlme.** Parts of the interface follow the package nlme
(Pinheiro and Bates, 2000), so that models can be moved between
[`nlme::nlme`](https://rdrr.io/pkg/nlme/man/nlme.html)/[`nlme::gnls`](https://rdrr.io/pkg/nlme/man/gnls.html)
and gcTMB:

- the fixed and random effects are formulas for the parameters, as the
  arguments `fixed` and `random` of
  [`nlme::nlme`](https://rdrr.io/pkg/nlme/man/nlme.html):
  `fform=b1 + b2 ~ trt` or `fform=list(b1 ~ trt, b2 ~ 1)`, and
  `rform=b1 + b2 ~ 1 | ID`. The random effects of a formula have an
  unstructured covariance matrix, as `pdSymm`; uncorrelated random
  effects, as `pdDiag`, are separate terms, e.g.
  `rform=list(b1 ~ 1 | ID, b2 ~ 1 | ID)` or `b1 ~ x || ID`. Nested
  factors (`1 | site/ID`) and the double bar use the syntax of lme4;

- the correlation structures
  [`corAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  and
  [`corCAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  have the arguments and the parameterisation of nlme, whose objects are
  accepted as well, and `residuals(fit, type="normalized")` gives the
  normalized residuals of nlme;

- [`fixef`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
  [`ranef`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  and
  [`VarCorr`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  are generics with the names of those of nlme, and work for gcTMB, nlme
  and lme4 fits whichever package is attached last;

- `predict` and `fitted` with `random=FALSE` correspond to `level=0` of
  nlme (the population curves), and with `random=TRUE` to the innermost
  level.

The differences: gcTMB fits by maximum likelihood (the Laplace
approximation with random effects), so its log-likelihoods, AIC and
likelihood ratio tests compare with nlme fits with `method="ML"`, not
with the default REML of
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) and
[`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html); the integration
over the random effects is exact for gaussian responses with random
effects that enter the curve linearly, e.g. random asymptotes, where
nlme's linearisation can stop slightly short of the maximum. `coef`
returns the estimated coefficients and variance parameters, not the
coefficients of each group as `coef` of nlme (use
`predict(fit, type="parameters")` for those); `weights` are numeric
weights, not variance functions; and the results of `ranef` and
`VarCorr` have the format of lme4 (a list of data frames, and of
covariance matrices).
[`confint.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md)
gives the intervals of `intervals` of nlme.

The asymmetry b5, the truncation at the threshold b4, the maximum growth
rates and inflection points, the rate 1/b3 and the small-sample
inference are described in the vignette,
[`vignette("gcTMB")`](https://emilskindersoe.github.io/gcTMB/articles/gcTMB.md);
the links in the argument `link`.

With `b1type="rate"`, b1 is the maximum growth rate of the curve, the
slope at its inflection point (see
[`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)),
instead of the asymptote; the asymptote is b1 divided by the maximum
growth rate of the curve with asymptote 1, which depends on b2-b5. The
model is the same as with the asymptote, unless b1 has random effects,
which are then random maximum growth rates. `predict(type="parameters")`
gives the asymptote as well. The maximum growth rate has a closed form
for all the models with b5 = 1, and for the logistic (Richards) and hill
curves also with an asymmetry b5; the gompertz, weibull, normal and
gamma models need b5 fixed at 1. The hill curve has a maximum growth
rate at a positive time only for b3 b5 \> 1. For the linear model, b1 is
the rate anyway.

## References

Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and
S-PLUS. Springer, New York.

Ciceu, A., Garcia-Duro, J., Seceleanu, I. and Badea, O. (2020). A
generalized nonlinear mixed-effects height-diameter model for Norway
spruce in mixed-uneven aged stands. Forest Ecology and Management 477,
118507.

Preece, M.A. and Baines, M.J. (1978). A new family of mathematical
models describing the human growth curve. Annals of Human Biology 5,
1-24.

Richards, F.J. (1959). A flexible growth function for empirical use.
Journal of Experimental Botany 10, 290-300.

Yin, X., Goudriaan, J., Lantinga, E.A., Vos, J. and Spiertz, H.J.
(2003). A flexible sigmoid function of determinate growth. Annals of
Botany 91, 361-371.

## See also

[`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
for the maximum growth rate and its timing,
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
for the times to reach a level of the asymptote,
[`gcTMBstart`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md)
for the starting values.

## Examples

``` r
# logistic growth of the faba bean plants
fit <- gcTMB(Height ~ DAS, data=FabaHeights25, model="logistic")
summary(fit)
#> Call:
#> gcTMB(form = Height ~ DAS, data = FabaHeights25, model = "logistic")
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 200
#> Log-likelihood: -699.2 (df = 4), AIC: 1406, BIC: 1420
#> 
#> Coefficients:
#>                Estimate Std. Error t value Pr(>|t|)    
#> b1:(Intercept) 120.7851     3.6524   33.07   <2e-16 ***
#> b2:(Intercept)  30.7300     0.6822   45.05   <2e-16 ***
#> b3:(Intercept)  12.8808     0.9769   13.19   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard deviation: 8.043 on 197 degrees of freedom
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 3.4e-08

# the other models, compared by AIC
fitg <- gcTMB(Height ~ DAS, data=FabaHeights25, model="gompertz")
fith <- gcTMB(Height ~ DAS, data=FabaHeights25, model="hill")
fitw <- gcTMB(Height ~ DAS, data=FabaHeights25, model="weibull")
fitn <- gcTMB(Height ~ DAS, data=FabaHeights25, model="normal")
fitgam <- gcTMB(Height ~ DAS, data=FabaHeights25, model="gamma")
fitlin <- gcTMB(Height ~ DAS, data=FabaHeights25, model="linear")
# the forestry models from zero instead of their default intercept, the breast height of 1.3 m
fitcr <- gcTMB(Height ~ DAS, data=FabaHeights25, model="chapmanrichards", start=c(b4=0))
fitb <- gcTMB(Height ~ DAS, data=FabaHeights25, model="beta_capped")
fitk <- gcTMB(Height ~ DAS, data=FabaHeights25, model="korf", start=c(b4=0))
AIC(fit, fitg, fith, fitw, fitn, fitgam, fitlin, fitcr, fitb, fitk)
#>        df      AIC
#> fit     4 1406.470
#> fitg    4 1407.658
#> fith    4 1409.892
#> fitw    4 1407.757
#> fitn    4 1406.412
#> fitgam  4 1409.685
#> fitlin  3 1443.592
#> fitcr   4 1410.100
#> fitb    4 1406.751
#> fitk    4 1415.341

# a height-diameter model from the breast height of 1.3 m, the default intercept b4
gcTMB(Height ~ DBH, data=PineHeightDiameter, model="chapmanrichards")
#> Call:
#> gcTMB(form = Height ~ DBH, data = PineHeightDiameter, model = "chapmanrichards")
#> 
#> Model:  chapmanrichards
#> Family: gaussian, identity link for all parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>       45.71868         0.01965         1.05978  
#> Fixed: b4:(Intercept) = 1.3
#> 
#> Residual standard deviation: 2.434

# the Richards curve: the logistic model with an estimated asymmetry b5
fitr5 <- gcTMB(Height ~ DAS, data=FabaHeights25, fix=c(b5=FALSE))
coef(fitr5)
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b5:(Intercept)      log_sigma 
#>    120.9525457     30.6657471     12.9844198      0.9728521      2.0772341 

# the rate 1/b3 instead of the scale b3
gcTMB(Height ~ DAS, data=FabaHeights25, link=c(b3="inverse"))
#> Call:
#> gcTMB(form = Height ~ DAS, data = FabaHeights25, link = c(b3 = "inverse"))
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, inverse link for b3, identity link for the other parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>      120.78512        30.72998         0.07764  
#> 
#> Residual standard deviation: 8.043

# treatment effects on the asymptote as ratios, with the log link
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fitl <- gcTMB(Height ~ DAS, fform=b1 ~ trt, data=FabaHeights25, link=c(b1="log"))
exp(coef(fitl)[2:5])
#> b1:trt0.3 b1:trt0.5 b1:trt0.6 b1:trt0.8 
#>  1.007483  1.077149  1.075851  1.050422 

# treatment effects on the asymptote and random effects of the plants
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fitre <- gcTMB(Height ~ DAS, fform=b1 ~ trt, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
summary(fitre)
#> Call:
#> gcTMB(form = Height ~ DAS, fform = b1 ~ trt, rform = b1 + b2 ~ 
#>     1 | Plant_ID, data = FabaHeights25)
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 200, groups: Plant_ID 25
#> Log-likelihood: -509.5 (df = 11), AIC: 1041, BIC: 1077
#> 
#> Coefficients:
#>                Estimate Std. Error z value Pr(>|z|)    
#> b1:(Intercept) 115.7419     4.0576  28.525   <2e-16 ***
#> b1:trt0.3       -3.1569     5.6806  -0.556   0.5784    
#> b1:trt0.5        9.8502     5.5933   1.761   0.0782 .  
#> b1:trt0.6       11.1356     5.6175   1.982   0.0474 *  
#> b1:trt0.8        8.9832     5.6351   1.594   0.1109    
#> b2:(Intercept)  30.8745     0.5402  57.151   <2e-16 ***
#> b3:(Intercept)  12.8441     0.2498  51.418   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Unadjusted Wald tests; tests of the differences between curves adjusted for their number: partest().
#> 
#> Random effects:
#>  Groups   Name Std.Dev. Corr 
#>  Plant_ID b1   8.842         
#>           b2   2.550    -0.33
#> 
#> Residual standard deviation: 2.075
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 5.2e-05

# treatment effects on the maximum growth rate, b1 with b1type = "rate", tested by a likelihood
# ratio test and compared with the first treatment
fitr <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25, b1type="rate")
fitr0 <- gcTMB(Height ~ DAS, fform=b2 ~ trt, data=FabaHeights25, b1type="rate")
anova(fitr0, fitr)
#> Likelihood ratio tests of gcTMB models
#> 
#>       npar    AIC    BIC  logLik deviance     LR Df  Pr(>LR)   
#> fitr0    8 1409.4 1435.8 -696.71   1393.4                      
#> fitr    12 1401.8 1441.4 -688.90   1377.8 15.619  4 0.003576 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
partest(fitr, newdata=data.frame(trt=levels(FabaHeights25$trt)))
#> Differences to the reference of b1 (maximum growth rate):
#> 
#>                   Estimate Std.Error t value p.value    lower upper
#> trt=0.3 - trt=0.1  -0.0149    0.0776  -0.193  0.9991 -0.20635 0.176
#> trt=0.5 - trt=0.1   0.2053    0.0822   2.498  0.0462  0.00247 0.408
#> trt=0.6 - trt=0.1   0.2180    0.0835   2.611  0.0344  0.01197 0.424
#> trt=0.8 - trt=0.1   0.1624    0.0839   1.937  0.1705 -0.04457 0.369
#> 
#> Wald t tests on 189 degrees of freedom; p-values and confidence limits adjusted by the single-step (Dunnett) method, 95% confidence limits.

# residuals of each plant correlated over time, continuous-time AR(1)
fitc <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25,
              correlation=corCAR1(form=~ DAS | Plant_ID))
corpars(fitc)
#>      rho 
#> 0.953607 

# pod counts with a poisson distribution
gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson")
#> Call:
#> gcTMB(form = Pods ~ DAS, data = FabaCounts23, family = "poisson")
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: poisson, identity link for all parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>         16.361          71.947           1.569  

# no pods before a threshold b4, after which the count grows
fitt <- gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE))
coef(fitt)
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b4:(Intercept) 
#>      17.690093      72.296815       2.117993      65.000000 
```
