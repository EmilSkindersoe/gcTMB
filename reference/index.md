# Package index

## Fitting growth curves

- [`gcTMB()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
  : Fit a growth curve
- [`gcTMBcontrol()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md)
  : Set control parameters for gcTMB
- [`gcTMBstart()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md)
  : Starting values for a gcTMB model
- [`corAR1()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  [`corCAR1()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  [`corExp()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  [`corGaus()`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  : Autocorrelation of the residuals within groups
- [`varIdent()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md)
  [`varPower()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md)
  [`varConstPower()`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md)
  : Variance functions for heteroscedastic residuals
- [`refit()`](https://emilskindersoe.github.io/gcTMB/reference/refit.md)
  : Refit a gcTMB model with other components
- [`update(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/update.gcTMB.md)
  : Set the parameters of a gcTMB fit
- [`converged()`](https://emilskindersoe.github.io/gcTMB/reference/converged.md)
  : Convergence of a gcTMB fit

## Methods and extractors

- [`print(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`summary(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`print(`*`<summary.gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`logLik(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`nobs(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`df.residual(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`fitted(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`residuals(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`predict(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`coef(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`vcov(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  [`simulate(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
  : gcTMB Methods
- [`fixef()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  [`ranef()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  [`VarCorr()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  [`sigma(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  [`corpars()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  [`varpars()`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  : Extract the components of a gcTMB fit
- [`plot(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/plot.gcTMB.md)
  : Plots of gcTMB models: validation and fitted curves
- [`criteria()`](https://emilskindersoe.github.io/gcTMB/reference/criteria.md)
  : Model selection criteria

## Exploring data

- [`viewcurves()`](https://emilskindersoe.github.io/gcTMB/reference/viewcurves.md)
  : Viewer of observed growth curves

## Inference

- [`anova(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/anova.gcTMB.md)
  : Likelihood ratio tests for gcTMB models
- [`confint(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md)
  [`confint(`*`<timetoTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md)
  : Confidence intervals for gcTMB models and times to reach a level of
  the asymptote
- [`profile(`*`<gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md)
  [`confint(`*`<profile.gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md)
  [`plot(`*`<profile.gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md)
  : Likelihood profiles of a gcTMB model
- [`bootstrap()`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md)
  [`print(`*`<bootstrap.gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md)
  [`confint(`*`<bootstrap.gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md)
  : Parametric bootstrap of a gcTMB model
- [`partest()`](https://emilskindersoe.github.io/gcTMB/reference/partest.md)
  [`print(`*`<partest.gcTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/partest.md)
  : Tests of the differences of a parameter between growth curves

## Timing and maximum growth rates

- [`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
  : Time to reach a level of the asymptote
- [`print(`*`<timetoTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/timetoTMB-methods.md)
  [`coef(`*`<timetoTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/timetoTMB-methods.md)
  [`vcov(`*`<timetoTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/timetoTMB-methods.md)
  : timetoTMB Methods
- [`profile(`*`<timetoTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/profile.timetoTMB.md)
  : Likelihood profiles of times to reach a level of the asymptote
- [`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
  : Maximum growth rate and its timing
- [`print(`*`<maxrateTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/maxrateTMB-methods.md)
  [`coef(`*`<maxrateTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/maxrateTMB-methods.md)
  [`vcov(`*`<maxrateTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/maxrateTMB-methods.md)
  : maxrateTMB Methods
- [`confint(`*`<maxrateTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/confint.maxrateTMB.md)
  : Confidence intervals for the maximum growth rate and its timing
- [`profile(`*`<maxrateTMB>`*`)`](https://emilskindersoe.github.io/gcTMB/reference/profile.maxrateTMB.md)
  : Profile likelihood of the maximum growth rate and its timing

## Data

- [`FabaHeights25`](https://emilskindersoe.github.io/gcTMB/reference/FabaHeights25.md)
  : Heights of faba bean plants
- [`FabaHeights26`](https://emilskindersoe.github.io/gcTMB/reference/FabaHeights26.md)
  : Heights of faba bean plants in a field trial
- [`FabaFlowers26`](https://emilskindersoe.github.io/gcTMB/reference/FabaFlowers26.md)
  : Flower and pod counts of faba bean plants in a field trial
- [`FabaNodes26`](https://emilskindersoe.github.io/gcTMB/reference/FabaNodes26.md)
  : Flowers and pods per node of faba bean plants in a field trial
- [`FabaCounts23`](https://emilskindersoe.github.io/gcTMB/reference/FabaCounts23.md)
  : Flower and pod counts of faba bean plants
- [`BerkeleyGrowth`](https://emilskindersoe.github.io/gcTMB/reference/BerkeleyGrowth.md)
  : Heights and weights of children in the Berkeley Child Guidance Study
- [`PineHeightDiameter`](https://emilskindersoe.github.io/gcTMB/reference/PineHeightDiameter.md)
  : Heights and diameters of Scots pine trees
- [`SorghumMaize`](https://emilskindersoe.github.io/gcTMB/reference/SorghumMaize.md)
  : Biomass of sorghum and maize

## Package

- [`gcTMB-package`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-package.md)
  : gcTMB: Fitting Growth Curves Using the Template Model Builder
