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
#'   \item{Height}{plant height}
#'   \item{DASfact}{days after sowing, the same values as DAS}
#'   \item{Plant_ID}{identifier of the plant, the same values as Tube_ID}
#' }
#' The raw data and the script that converts them are in the directory data-raw of the source package.
#'
#' @examples
#' fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
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

#' Flower and pod counts of faba bean plants in a field trial
#'
#' Repeated counts of flowers and pods of 96 faba bean plants in a field trial (2026): three plants
#' in each of 32 plots, with two genotypes on two soil types, irrigated or not, counted on eight
#' days between 59 and 84 days after sowing. The plants are plants 1-3 of the five of each plot in
#' \code{\link{FabaHeights26}}, with the same \code{Plant_ID}, so the two data sets can be merged.
#' The counts were recorded per node, each node only once: on the first date it had flowers (or
#' pods), as no more appear on a node afterwards. \code{Flowers} and \code{Pods} are therefore the
#' new flowers and pods since the previous date, summed over the nodes of each plant, and
#' \code{FlowersTotal} and \code{PodsTotal} the running totals, which follow a growth curve. A plant
#' without new flowers or pods on a date, including a plant that was no longer recorded, has 0 new
#' flowers or pods. Flowering had begun before the first count, so the flowers of the first date
#' include those of earlier days, and \code{FlowersTotal} does not start at 0.
#'
#' @format A data frame with 768 rows and 10 variables:
#' \describe{
#'   \item{DAS}{days after sowing}
#'   \item{Genotype}{genotype: Navarra or Stella}
#'   \item{Plot_ID}{identifier of the plot}
#'   \item{Plant_ID}{identifier of the plant, a factor unique to each plant: plot and plant number
#'     within the plot, e.g. "31_3", as in \code{FabaHeights26}}
#'   \item{Soil}{soil type: Loamy or Sandy}
#'   \item{Treatment}{Irrigated or NotIrrigated}
#'   \item{Flowers}{number of new flowers on the plant since the previous date}
#'   \item{Pods}{number of new pods on the plant since the previous date}
#'   \item{FlowersTotal}{number of flowers the plant has produced up to and including the date}
#'   \item{PodsTotal}{number of pods the plant has set up to and including the date}
#' }
#' The raw data and the script that converts them, \code{DIGIFABA_field_26.R}, are in the directory
#' data-raw of the source package.
#'
#' @examples
#' fit <- gcTMB(PodsTotal ~ DAS, fform=b1 + b2 ~ Soil * Treatment, rform=b1 ~ 1 | Plant_ID,
#'              data=FabaFlowers26)
#' summary(fit)
#' timeto(fit, level=50, newdata=expand.grid(Soil=c("Loamy", "Sandy"),
#'                                           Treatment=c("Irrigated", "NotIrrigated")), random=FALSE)
"FabaFlowers26"

#' Flowers and pods per node of faba bean plants in a field trial
#'
#' This data set is not a growth curve and is not meant for the growth models of the package: it
#' describes the pod set (pods per flower) of each node, e.g. for a binomial model as in the
#' example. It belongs to the same field trial (2026) as \code{\link{FabaHeights26}} and
#' \code{\link{FabaFlowers26}}, with the same plants and \code{Plant_ID}, so that it can be merged
#' with them, e.g. to relate the pod set of a plant to its height or to the timing of its
#' flowering.
#'
#' The flowers and pods of \code{\link{FabaFlowers26}} by node instead of by date: one row for each
#' node of the 96 plants that had flowers or pods, with the date and number of its flowers and the
#' date and number of its pods, for the pod set (pods per flower) of each node. Each node was
#' counted once, on the first date it had flowers and the first date it had pods; nodes counted
#' again on a later date are left out of that date. A node that flowered but has no pod record set
#' no pods (\code{Pods} = 0). The number of flowers is missing (\code{NA}) where the flowers of the
#' node were not counted: a blank count, or pods on a node without a flower record.
#'
#' @format A data frame with 879 rows and 11 variables:
#' \describe{
#'   \item{Genotype}{genotype: Navarra or Stella}
#'   \item{Plot_ID}{identifier of the plot}
#'   \item{Plant_ID}{identifier of the plant, as in \code{FabaFlowers26} and \code{FabaHeights26},
#'     e.g. "31_3"}
#'   \item{Soil}{soil type: Loamy or Sandy}
#'   \item{Treatment}{Irrigated or NotIrrigated}
#'   \item{Node}{node number on the stem}
#'   \item{DAS_flowers}{days after sowing of the count of the flowers of the node, \code{NA} without
#'     a flower record}
#'   \item{Flowers}{number of flowers of the node, \code{NA} if not counted}
#'   \item{DAS_pods}{days after sowing of the count of the pods of the node, \code{NA} without pods}
#'   \item{Pods}{number of pods of the node}
#'   \item{Warning}{1 if the record of the node cannot be right: more pods than flowers, or the pods
#'     counted on the same date as the flowers; 0 otherwise}
#' }
#' Pods appear one to two weeks after flowering, so nodes that flowered in the last week of the
#' counts may not have had the time to set pods. The raw data and the script that converts them,
#' \code{DIGIFABA_field_26.R}, are in the directory data-raw of the source package.
#'
#' @examples
#' # pod set per node by soil, irrigation and node position, from nodes with a reliable count
#' ok <- subset(FabaNodes26, !is.na(Flowers) & Warning == 0 & DAS_flowers <= 73)
#' fit <- glm(cbind(Pods, Flowers - Pods) ~ Soil * Treatment + Genotype + Node, family=binomial,
#'            data=ok)
#' summary(fit)
"FabaNodes26"

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
#'   \item{New}{number of new flowers}
#'   \item{Open}{number of open flowers}
#'   \item{Old}{number of old flowers}
#'   \item{Pods}{number of pods}
#'   \item{N_flowers}{number of flowers, the sum of New, Open and Old}
#' }
#' The raw data and the script that converts them are in the directory data-raw of the source package.
#'
#' @examples
#' fit <- gcTMB(Pods ~ DAS, fform=b1 ~ Genotype, data=FabaCounts23, family="poisson")
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
#'   \item{ID}{identifier of the child, a factor: 201-278 for the boys, 301-385 for the girls}
#'   \item{Sex}{male or female}
#'   \item{Age}{age in years}
#'   \item{Height}{height in cm}
#'   \item{Weight}{weight in kg}
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
#' # the growth of the girls from 2 years on: the Preece-Baines curve, with the adult height b1 and
#' # the middle of puberty b2, and the peak height velocity of the pubertal growth spurt
#' girls <- subset(BerkeleyGrowth, Sex == "female" & Age >= 2)
#' fit <- gcTMB(Height ~ Age, data=girls, model="preecebaines")
#' coef(fit)
#' maxrate(fit)
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
#'   \item{Plot}{identifier of the sample plot, a factor}
#'   \item{DBH}{diameter at breast height in cm}
#'   \item{Height}{tree height in m}
#'   \item{Ntrees}{number of trees on the plot}
#'   \item{Dmean}{mean diameter of the trees on the plot, in cm}
#'   \item{Dvar}{variance of the diameters of the trees on the plot}
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
#' # a Chapman-Richards curve from the breast height, the intercept b4 = 1.3, with a random asymptote
#' # of each plot
#' fit <- gcTMB(Height ~ DBH, rform=b1 ~ 1 | Plot, data=PineHeightDiameter,
#'              model="chapmanrichards", start=c(b4=1.3))
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
#' # growth from sowing, b4 = 141, to the end of growth b3, and the decline after it: the beta growth
#' # function, with a maximum of each crop and input level
#' fit <- gcTMB(Biomass ~ DOY, fform=b1 ~ Crop * Input, data=SorghumMaize, model="beta", start=c(b4=141))
#' summary(fit)
"SorghumMaize"
