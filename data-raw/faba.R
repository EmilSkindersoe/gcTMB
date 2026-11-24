# Converts the raw faba bean data in data-raw/ to the package datasets in data/.
# Run from the package root: source("data-raw/faba.R")

# plant heights: without the row numbers of the csv, and with the date as a Date
FabaHeights25 <- read.csv("data-raw/FabaHeights25.csv")
FabaHeights25$X <- NULL
FabaHeights25$Date <- as.Date(FabaHeights25$Date)
# every column name starts with a capital letter, as in the other datasets
names(FabaHeights25)[names(FabaHeights25) == "height"] <- "Height"
names(FabaHeights25)[names(FabaHeights25) == "plant_ID"] <- "Plant_ID"

# flower and pod counts: without the row numbers and the duplicate column Genoype (identical to
# Genotype)
FabaCounts23 <- read.csv("data-raw/FabaCounts23.csv")
stopifnot(identical(FabaCounts23$Genotype, FabaCounts23$Genoype))
FabaCounts23$X <- NULL
FabaCounts23$Genoype <- NULL
# every column name starts with a capital letter, as in the other datasets
names(FabaCounts23) <- paste0(toupper(substring(names(FabaCounts23), 1, 1)), substring(names(FabaCounts23), 2))

save(FabaHeights25, file="data/FabaHeights25.rda", compress="xz")
save(FabaCounts23, file="data/FabaCounts23.rda", compress="xz")
