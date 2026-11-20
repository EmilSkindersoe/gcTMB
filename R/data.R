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
