# Flowers and pods per node of faba bean plants in a field trial

This data set is not a growth curve and is not meant for the growth
models of the package: it describes the pod set (pods per flower) of
each node, e.g. for a binomial model as in the example. It belongs to
the same field trial (2026) as
[`FabaHeights26`](https://emilskindersoe.github.io/gcTMB/reference/FabaHeights26.md)
and
[`FabaFlowers26`](https://emilskindersoe.github.io/gcTMB/reference/FabaFlowers26.md),
with the same plants and `Plant_ID`, so that it can be merged with them,
e.g. to relate the pod set of a plant to its height or to the timing of
its flowering.

## Usage

``` r
FabaNodes26
```

## Format

A data frame with 879 rows and 11 variables:

- Genotype:

  genotype: Navarra or Stella

- Plot_ID:

  identifier of the plot

- Plant_ID:

  identifier of the plant, as in `FabaFlowers26` and `FabaHeights26`,
  e.g. "31_3"

- Soil:

  soil type: Loamy or Sandy

- Treatment:

  Irrigated or NotIrrigated

- Node:

  node number on the stem

- DAS_flowers:

  days after sowing of the count of the flowers of the node, `NA`
  without a flower record

- Flowers:

  number of flowers of the node, `NA` if not counted

- DAS_pods:

  days after sowing of the count of the pods of the node, `NA` without
  pods

- Pods:

  number of pods of the node

- Warning:

  1 if the record of the node cannot be right: more pods than flowers,
  or the pods counted on the same date as the flowers; 0 otherwise

Pods appear one to two weeks after flowering, so nodes that flowered in
the last week of the counts may not have had the time to set pods. The
raw data and the script that converts them, `DIGIFABA_field_26.R`, are
in the directory data-raw of the source package.

## Details

The flowers and pods of
[`FabaFlowers26`](https://emilskindersoe.github.io/gcTMB/reference/FabaFlowers26.md)
by node instead of by date: one row for each node of the 96 plants that
had flowers or pods, with the date and number of its flowers and the
date and number of its pods, for the pod set (pods per flower) of each
node. Each node was counted once, on the first date it had flowers and
the first date it had pods; nodes counted again on a later date are left
out of that date. A node that flowered but has no pod record set no pods
(`Pods` = 0). The number of flowers is missing (`NA`) where the flowers
of the node were not counted: a blank count, or pods on a node without a
flower record.

## Examples

``` r
# pod set per node by soil, irrigation and node position, from nodes with a reliable count
ok <- subset(FabaNodes26, !is.na(Flowers) & Warning == 0 & DAS_flowers <= 73)
fit <- glm(cbind(Pods, Flowers - Pods) ~ Soil * Treatment + Genotype + Node, family=binomial,
           data=ok)
summary(fit)
#> 
#> Call:
#> glm(formula = cbind(Pods, Flowers - Pods) ~ Soil * Treatment + 
#>     Genotype + Node, family = binomial, data = ok)
#> 
#> Coefficients:
#>                                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                      0.89433    0.17401   5.139 2.76e-07 ***
#> SoilSandy                        0.10949    0.09627   1.137    0.255    
#> TreatmentNotIrrigated           -0.48369    0.09041  -5.350 8.80e-08 ***
#> GenotypeStella                  -0.56841    0.07171  -7.926 2.26e-15 ***
#> Node                            -0.07932    0.01299  -6.106 1.02e-09 ***
#> SoilSandy:TreatmentNotIrrigated -0.63003    0.14946  -4.215 2.49e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1184.44  on 751  degrees of freedom
#> Residual deviance:  991.52  on 746  degrees of freedom
#> AIC: 2124
#> 
#> Number of Fisher Scoring iterations: 4
#> 
```
