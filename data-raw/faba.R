# Converts the raw faba bean data in data-raw/ to the package datasets in data/.
# Run from the package root: source("data-raw/faba.R")

# plant heights: without the row numbers of the csv, and with the date as a Date
FabaHeights25 <- read.csv("data-raw/FabaHeights25.csv")
FabaHeights25$X <- NULL
FabaHeights25$Date <- as.Date(FabaHeights25$Date)

# flower and pod counts: without the row numbers and the duplicate column Genoype (identical to
# Genotype)
FabaCounts23 <- read.csv("data-raw/FabaCounts23.csv")
stopifnot(identical(FabaCounts23$Genotype, FabaCounts23$Genoype))
FabaCounts23$X <- NULL
FabaCounts23$Genoype <- NULL

save(FabaHeights25, file="data/FabaHeights25.rda", compress="xz")
save(FabaCounts23, file="data/FabaCounts23.rda", compress="xz")
