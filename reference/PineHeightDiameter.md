# Heights and diameters of Scots pine trees

Heights and diameters at breast height (1.3 m above the ground) of 1678
Scots pine trees on 56 fixed-area sample plots in Ilomantsi, eastern
Finland, collected by Timo Pukkala, for height-diameter models: the
height grows with the diameter towards an asymptote, from the breast
height of 1.3 m. The positions of the trees within their plots allow
spatial correlation of the residuals between neighbouring trees, e.g.
`corExp(form=~ X + Y | Plot)`.

## Usage

``` r
PineHeightDiameter
```

## Format

A data frame with 1678 rows and 8 variables:

- Plot:

  identifier of the sample plot, a factor

- DBH:

  diameter at breast height in cm

- Height:

  tree height in m

- Ntrees:

  number of trees on the plot

- Dmean:

  mean diameter of the trees on the plot, in cm

- Dvar:

  variance of the diameters of the trees on the plot

- X, Y:

  position of the tree within its plot, in m

## Source

The data set `spati2` of the package lmfor (version 1.7, licence GPL-2),
with the variables renamed, and the positions of the same trees from the
data set `spati` of lmfor, converted from decimetres to metres. The
script that converts them is in the directory data-raw of the source
package.

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
fit <- gcTMB(Height ~ DBH, rform=b1 ~ 1 | Plot, data=PineHeightDiameter,
             model="chapmanrichards", start=c(b4=1.3))
summary(fit)
#> Call:
#> gcTMB(form = Height ~ DBH, rform = b1 ~ 1 | Plot, data = PineHeightDiameter, 
#>     start = c(b4 = 1.3), model = "chapmanrichards")
#> 
#> Model:  chapmanrichards
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 1678, groups: Plot 56
#> Log-likelihood: -2740 (df = 5), AIC: 5490, BIC: 5517
#> 
#> Coefficients:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> b1:(Intercept) 19.370992   0.790586   24.50   <2e-16 ***
#> b2:(Intercept)  0.058208   0.003636   16.01   <2e-16 ***
#> b3:(Intercept)  0.918873   0.029403   31.25   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Fixed: b4:(Intercept) = 1.3
#> 
#> Random effects:
#>  Groups Name Std.Dev.
#>  Plot   b1   5.16    
#> 
#> Residual standard deviation: 1.141
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 0.0019
```
