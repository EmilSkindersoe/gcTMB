# Heights and diameters of Scots pine trees

Heights and diameters at breast height (1.3 m above the ground) of 1678
Scots pine trees on 56 fixed-area sample plots in Ilomantsi, eastern
Finland, collected by Timo Pukkala, for height-diameter models: the
height grows with the diameter towards an asymptote, from the breast
height of 1.3 m.

## Usage

``` r
PineHeightDiameter
```

## Format

A data frame with 1678 rows and 6 variables:

- plot:

  identifier of the sample plot, a factor

- dbh:

  diameter at breast height in cm

- height:

  tree height in m

- ntrees:

  number of trees on the plot

- dmean:

  mean diameter of the trees on the plot, in cm

- dvar:

  variance of the diameters of the trees on the plot

## Source

The data set `spati2` of the package lmfor (version 1.7, licence GPL-2),
with the variables renamed. The script that converts it is in the
directory data-raw of the source package.

## References

Pukkala, T. (1989). Prediction of tree diameter and height in a Scots
pine stand as a function of the spatial pattern of trees. Silva Fennica
23, 83-99. [doi:10.14214/sf.a15532](https://doi.org/10.14214/sf.a15532)

Mehtätalo, L. and Lappi, J. (2020). Biometry for Forestry and
Environmental Data: with Examples in R. Chapman and Hall/CRC, New York.
[doi:10.1201/9780429173462](https://doi.org/10.1201/9780429173462)

## Examples

``` r
# a Chapman-Richards curve from the breast height, the intercept b4 = 1.3, with a random asymptote
# of each plot
fit <- gcTMB(height ~ dbh, rform=b1 ~ 1 | plot, data=PineHeightDiameter,
             model="chapmanrichards", start=c(b4=1.3))
summary(fit)
#> Call:
#> gcTMB(form = height ~ dbh, rform = b1 ~ 1 | plot, data = PineHeightDiameter, 
#>     start = c(b4 = 1.3), model = "chapmanrichards")
#> 
#> Model:  chapmanrichards
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 1678, groups: plot 56
#> Log-likelihood: -2740 (df = 5), AIC: 5490, BIC: 5517
#> 
#> Coefficients:
#>                Estimate Std. Error
#> b1:(Intercept)   19.371      0.791
#> b2:(Intercept)   -1.454      0.634
#> b3:(Intercept)   17.180      1.074
#> Fixed: b4:(Intercept) = 1.3
#> 
#> Random effects:
#>  Groups Name Std.Dev.
#>  plot   b1   5.163   
#> 
#> Residual standard deviation: 1.141
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 1e-04
```
