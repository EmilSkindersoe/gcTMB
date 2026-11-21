# Heights and weights of children in the Berkeley Child Guidance Study

Longitudinal heights and weights of 66 boys and 70 girls from Berkeley,
California, born in 1928-29, measured from birth to 21 years: at birth,
0.085 years, every three months up to 2 years, every year up to 8 years
and every six months up to 21 years. Heights were measured from 1.75
years on; the heights of younger children are missing. The data were
published as an appendix to Tuddenham and Snyder (1954). The heights of
the first 18 years of most of the children are the `growth` data of the
package fda.

## Usage

``` r
BerkeleyGrowth
```

## Format

A data frame with 4884 rows and 5 variables:

- id:

  identifier of the child, a factor: 201-278 for the boys, 301-385 for
  the girls

- sex:

  male or female

- age:

  age in years

- height:

  height in cm

- weight:

  weight in kg

## Source

The data set `berkeley` of the package sitar (version 1.5.0, licence GPL
(\>= 2)), with the transcription errors corrected there, without the
other body measurements. The script that converts it is in the directory
data-raw of the source package.

## References

Tuddenham, R.D. and Snyder, M.M. (1954). Physical growth of California
boys and girls from birth to eighteen years. University of California
Publications in Child Development 1, 183-364.

Cole, T.J. (2025). sitar: Super Imposition by Translation and Rotation
Growth Curve Analysis. R package,
<https://CRAN.R-project.org/package=sitar>.

## Examples

``` r
# the growth of the girls from 2 years on: the Preece-Baines curve, with the adult height b1 and
# the middle of puberty b2, and the peak height velocity of the pubertal growth spurt
girls <- subset(BerkeleyGrowth, sex == "female" & age >= 2)
fit <- gcTMB(height ~ age, data=girls, model="preecebaines")
coef(fit)
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b4:(Intercept) b5:(Intercept) 
#>    166.3370146     11.9345925    154.0119514      0.1180604      0.8819423 
#>      log_sigma 
#>      1.7673499 
maxrate(fit)
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>          Estimate Std.Error
#> tmax       10.866   0.19993
#> rate        6.830   0.18381
#> fraction    0.883   0.00762
```
