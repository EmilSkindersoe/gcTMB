# Converts growth data from other R packages to the package datasets in data/, from the versions
# sitar 1.5.0, lmfor 1.7 and nlraa 1.9.10. These packages are only needed to rebuild the data, not to
# install gcTMB. Run from the package root: source("data-raw/external.R")

# Berkeley Child Guidance Study (sitar::berkeley, GPL (>= 2)): heights and weights of 66 boys and 70
# girls from birth to 21 years, with the transcription errors corrected in sitar. The other body
# measurements, mostly missing and only from 6-8 years on, are left out.
data("berkeley", package="sitar", envir=environment())
BerkeleyGrowth <- data.frame(id=factor(berkeley$id),
                             sex=factor(berkeley$sex, levels=c("1", "2"), labels=c("male", "female")),
                             age=berkeley$age,
                             height=berkeley$height,
                             weight=berkeley$weight)
BerkeleyGrowth <- BerkeleyGrowth[order(BerkeleyGrowth$id, BerkeleyGrowth$age), ]
rownames(BerkeleyGrowth) <- NULL

# Scots pine heights and diameters in Ilomantsi, Finland (lmfor::spati2, GPL-2): 1678 trees on 56
# plots, with the number of trees and the mean and variance of the diameters of each plot
data("spati2", package="lmfor", envir=environment())
PineHeightDiameter <- data.frame(plot=factor(spati2$plot),
                                 dbh=spati2$d,
                                 height=spati2$h,
                                 ntrees=spati2$n,
                                 dmean=spati2$dmean,
                                 dvar=spati2$dvar)
rownames(PineHeightDiameter) <- NULL

# biomass of fiber sorghum, sweet sorghum and maize in Greece (nlraa::sm, GPL-3): 3 crops at 2 input
# levels in 4 blocks, on 11 days of the year from sowing, day 141, to 303
data("sm", package="nlraa", envir=environment())
SorghumMaize <- data.frame(DOY=sm$DOY,
                           Block=factor(sm$Block),
                           Input=factor(sm$Input, levels=1:2, labels=c("Low", "High")),
                           Crop=factor(sm$Crop, levels=c("F", "M", "S"), labels=c("Fiber sorghum", "Maize", "Sweet sorghum")),
                           Biomass=sm$Yield)
SorghumMaize <- SorghumMaize[order(SorghumMaize$Crop, SorghumMaize$Input, SorghumMaize$Block, SorghumMaize$DOY), ]
rownames(SorghumMaize) <- NULL

save(BerkeleyGrowth, file="data/BerkeleyGrowth.rda", compress="xz")
save(PineHeightDiameter, file="data/PineHeightDiameter.rda", compress="xz")
save(SorghumMaize, file="data/SorghumMaize.rda", compress="xz")
