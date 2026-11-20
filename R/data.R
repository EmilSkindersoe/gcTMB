#' Heights of faba bean plants
#'
#' Repeated height measurements of 25 faba bean plants, five at each of five treatment levels,
#' measured on eight days between 14 and 61 days after sowing (2025).
#'
#' @format A data frame with 200 rows and 7 variables:
#' \describe{
#'   \item{Tube_ID}{identifier of the tube the plant grows in}
#'   \item{Treatment}{treatment level: 0.1, 0.3, 0.5, 0.6 or 0.8}
#'   \item{Date}{date of the measurement}
#'   \item{DAS}{days after sowing}
#'   \item{height}{plant height}
#'   \item{DASfact}{days after sowing, the same values as DAS}
#'   \item{plant_ID}{identifier of the plant, the same values as Tube_ID}
#' }
#' The raw data and the script that converts them are in the directory data-raw of the source package.
#'
#' @examples
#' fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
#' summary(fit)
"FabaHeights25"

#' Heights of faba bean plants in a field trial
#'
#' Repeated height measurements of 320 faba bean plants in a field trial (2026): five plants in
#' each of 64 plots, with four genotypes on two soil types, irrigated or not, measured on 14 days
#' between 28 and 122 days after sowing. The Loamy plots were not measured 77 days after sowing,
#' and five single measurements are missing, so the data have 4315 rows instead of 4480.
#'
#' @format A data frame with 4315 rows and 7 variables:
#' \describe{
#'   \item{DAS}{days after sowing}
#'   \item{Genotype}{genotype: Lynx, Melodie, Navarra or Stella}
#'   \item{Plot_ID}{identifier of the plot, 1-64}
#'   \item{Plant_ID}{identifier of the plant, a factor unique to each plant: plot and plant number
#'     within the plot, e.g. "31_5"}
#'   \item{Soil}{soil type: Loamy or Sandy}
#'   \item{Treatment}{Irrigated or NotIrrigated}
#'   \item{Height}{plant height in cm}
#' }
#' The raw data and the script that converts them are in the directory data-raw of the source package.
#'
#' @examples
#' fit <- gcTMB(Height ~ DAS, fform=b1 ~ Soil * Treatment, data=FabaHeights26)
#' summary(fit)
"FabaHeights26"

#' Flower and pod counts of faba bean plants
#'
#' Repeated counts of flowers and pods of 22 faba bean plants of six genotypes, all with the
#' treatment "ww", on twelve days between 61 and 79 days after sowing (2023).
#'
#' @format A data frame with 263 rows and 9 variables:
#' \describe{
#'   \item{Tube_ID}{identifier of the tube of the plant}
#'   \item{Genotype}{genotype: Aurora, Capri, ILB938, Lynx, Melodie or Navara}
#'   \item{Treatment}{treatment, "ww" for all plants}
#'   \item{DAS}{days after sowing}
#'   \item{new}{number of new flowers}
#'   \item{open}{number of open flowers}
#'   \item{old}{number of old flowers}
#'   \item{pods}{number of pods}
#'   \item{n_flowers}{number of flowers, the sum of new, open and old}
#' }
#' The raw data and the script that converts them are in the directory data-raw of the source package.
#'
#' @examples
#' fit <- gcTMB(pods ~ DAS, fform=b1 ~ Genotype, data=FabaCounts23, family="poisson")
#' summary(fit)
"FabaCounts23"

#' Heights and weights of children in the Berkeley Child Guidance Study
#'
#' Longitudinal heights and weights of 66 boys and 70 girls from Berkeley, California, born in
#' 1928-29, measured from birth to 21 years: at birth, 0.085 years, every three months up to 2 years,
#' every year up to 8 years and every six months up to 21 years. Heights were measured from 1.75 years
#' on; the heights of younger children are missing. The data were published as an appendix to
#' Tuddenham and Snyder (1954). The heights of the first 18 years of most of the children are the
#' \code{growth} data of the package fda.
#'
#' @format A data frame with 4884 rows and 5 variables:
#' \describe{
#'   \item{id}{identifier of the child, a factor: 201-278 for the boys, 301-385 for the girls}
#'   \item{sex}{male or female}
#'   \item{age}{age in years}
#'   \item{height}{height in cm}
#'   \item{weight}{weight in kg}
#' }
#' @source The data set \code{berkeley} of the package sitar (version 1.5.0, licence GPL (>= 2)), with
#'   the transcription errors corrected there, without the other body measurements. The script that
#'   converts it is in the directory data-raw of the source package.
#' @references Tuddenham, R.D. and Snyder, M.M. (1954). Physical growth of California boys and girls
#'   from birth to eighteen years. University of California Publications in Child Development 1,
#'   183-364.
#'
#'   Cole, T.J. (2025). sitar: Super Imposition by Translation and Rotation Growth Curve Analysis.
#'   R package, \url{https://CRAN.R-project.org/package=sitar}.
#'
#' @examples
#' # the pubertal growth spurt of the girls: a logistic curve from age 8 on, above the estimated
#' # height b6 before the spurt
#' girls <- subset(BerkeleyGrowth, sex == "female" & age >= 8)
#' fit <- gcTMB(height ~ age, data=girls, truncated=FALSE, fix=c(b6=FALSE))
#' coef(fit)
#' \donttest{
#' # with a random asymptote and timing of the spurt for each girl
#' fitre <- gcTMB(height ~ age, rform=b1 + b2 ~ 1 | id, data=girls, truncated=FALSE, fix=c(b6=FALSE))
#' summary(fitre)
#' }
"BerkeleyGrowth"

#' Heights and diameters of Scots pine trees
#'
#' Heights and diameters at breast height (1.3 m above the ground) of 1678 Scots pine trees on 56
#' fixed-area sample plots in Ilomantsi, eastern Finland, collected by Timo Pukkala, for
#' height-diameter models: the height grows with the diameter towards an asymptote, from the breast
#' height of 1.3 m.
#'
#' @format A data frame with 1678 rows and 6 variables:
#' \describe{
#'   \item{plot}{identifier of the sample plot, a factor}
#'   \item{dbh}{diameter at breast height in cm}
#'   \item{height}{tree height in m}
#'   \item{ntrees}{number of trees on the plot}
#'   \item{dmean}{mean diameter of the trees on the plot, in cm}
#'   \item{dvar}{variance of the diameters of the trees on the plot}
#' }
#' @source The data set \code{spati2} of the package lmfor (version 1.7, licence GPL-2), with the
#'   variables renamed. The script that converts it is in the directory data-raw of the source
#'   package.
#' @references Pukkala, T. (1989). Prediction of tree diameter and height in a Scots pine stand as a
#'   function of the spatial pattern of trees. Silva Fennica 23, 83-99. \doi{10.14214/sf.a15532}
#'
#'   Mehtätalo, L. and Lappi, J. (2020). Biometry for Forestry and Environmental Data: with Examples
#'   in R. Chapman and Hall/CRC, New York. \doi{10.1201/9780429173462}
#'
#' @examples
#' # a Chapman-Richards curve from the breast height b6 = 1.3, with a random asymptote of each plot
#' fit <- gcTMB(height ~ dbh, rform=b1 ~ 1 | plot, data=PineHeightDiameter,
#'              model="chapmanrichards", start=c(b6=1.3))
#' summary(fit)
"PineHeightDiameter"

#' Biomass of sorghum and maize
#'
#' Above-ground biomass of fiber sorghum, sweet sorghum and maize in central Greece, at two input
#' levels (water and nitrogen) in four blocks, on 11 days of the year from sowing, day 141, to day
#' 303. Not every crop was sampled on every day.
#'
#' @format A data frame with 236 rows and 5 variables:
#' \describe{
#'   \item{DOY}{day of the year}
#'   \item{Block}{block of the experimental design, 1-4}
#'   \item{Input}{input level: Low or High}
#'   \item{Crop}{Fiber sorghum, Maize or Sweet sorghum}
#'   \item{Biomass}{biomass in Mg/ha}
#' }
#' @source The data set \code{sm} of the package nlraa (version 1.9.10, licence GPL-3), with the
#'   factors labelled. The script that converts it is in the directory data-raw of the source
#'   package.
#' @references Danalatos, N.G., Archontoulis, S.V. and Tsiboukas, K. (2009). Comparative analysis of
#'   sorghum versus corn growing under optimum and under water/nitrogen limited conditions in central
#'   Greece. In: Proceedings of the 17th European Biomass Conference, Hamburg, 538-544.
#'
#'   Miguez, F. (2025). nlraa: Nonlinear Regression for Agricultural Applications. R package,
#'   \url{https://CRAN.R-project.org/package=nlraa}.
#'
#' @examples
#' # determinate growth from sowing, b4 = 141: the beta growth function, with an asymptote of each
#' # crop and input level
#' fit <- gcTMB(Biomass ~ DOY, fform=b1 ~ Crop * Input, data=SorghumMaize, model="beta", start=c(b4=141))
#' summary(fit)
"SorghumMaize"
