# Biomass of sorghum and maize

Above-ground biomass of fiber sorghum, sweet sorghum and maize in
central Greece, at two input levels (water and nitrogen) in four blocks,
on 11 days of the year from sowing, day 141, to day 303. Not every crop
was sampled on every day.

## Usage

``` r
SorghumMaize
```

## Format

A data frame with 236 rows and 5 variables:

- DOY:

  day of the year

- Block:

  block of the experimental design, 1-4

- Input:

  input level: Low or High

- Crop:

  Fiber sorghum, Maize or Sweet sorghum

- Biomass:

  biomass in Mg/ha

## Source

The data set `sm` of the package nlraa (version 1.9.10, licence GPL-3),
with the factors labelled. The script that converts it is in the
directory data-raw of the source package.

## References

Danalatos, N.G., Archontoulis, S.V. and Tsiboukas, K. (2009).
Comparative analysis of sorghum versus corn growing under optimum and
under water/nitrogen limited conditions in central Greece. In:
Proceedings of the 17th European Biomass Conference, Hamburg, 538-544.

Miguez, F. (2025). nlraa: Nonlinear Regression for Agricultural
Applications. R package, <https://CRAN.R-project.org/package=nlraa>.

## Examples

``` r
# growth from sowing, b4 = 141, to the end of growth b3, and the decline after it: the beta growth
# function, with a maximum of each crop and input level
fit <- gcTMB(Biomass ~ DOY, fform=b1 ~ Crop * Input, data=SorghumMaize, model="beta", start=c(b4=141))
summary(fit)
#> Call:
#> gcTMB(form = Biomass ~ DOY, fform = b1 ~ Crop * Input, data = SorghumMaize, 
#>     start = c(b4 = 141), model = "beta")
#> 
#> Model:  beta
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 236
#> Log-likelihood: -650.7 (df = 9), AIC: 1319, BIC: 1351
#> 
#> Coefficients:
#>                                Estimate Std. Error t value Pr(>|t|)    
#> b1:(Intercept)                  28.7264     1.0693  26.864  < 2e-16 ***
#> b1:CropMaize                   -11.5069     1.7359  -6.629 2.42e-10 ***
#> b1:CropSweet sorghum             2.0471     1.4496   1.412    0.159    
#> b1:InputHigh                     5.8764     1.4487   4.056 6.85e-05 ***
#> b1:CropMaize:InputHigh           0.8143     2.4385   0.334    0.739    
#> b1:CropSweet sorghum:InputHigh   0.9234     2.0462   0.451    0.652    
#> b2:(Intercept)                 228.0015     1.7396 131.066  < 2e-16 ***
#> b3:(Intercept)                 278.5168     1.1671 238.647  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Unadjusted Wald tests; tests of the differences between curves adjusted for their number: partest().
#> Fixed: b4:(Intercept) = 141
#> 
#> Residual standard deviation: 3.878 on 228 degrees of freedom
#> 
#> Optimiser: both X-convergence and relative convergence (5), maximum absolute gradient 1.4e-13
```
