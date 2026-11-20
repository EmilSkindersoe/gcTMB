# Merges the two raw height files of the 2026 faba bean field trial into one clean table: the
# package dataset FabaHeights26 (data/FabaHeights26.rda) and a copy in data-raw/FabaHeights26.csv.
#   - data-raw/Height.xlsx      : all measurement dates except 6 July (13 dates)
#   - data-raw/Height 0607.xlsx : 6 July 2026 (DAS 83) only, without Genotype
# Run from the package root: source("data-raw/heights26.R")
# Every cleaning step ends in a stopifnot() check, so the script stops if the raw data are not
# as described in the comments.

library(readxl)

sowing <- as.Date("2026-04-14")  # DAS 0; every date in both files is sowing + DAS

# The two highlighted blocks in Height 0607.xlsx (plots 33-40 and 49-56) may have been entered
# the wrong way round: it was not written down which block was which, so the recorded order is
# a best guess. swap_highlighted = TRUE swaps their heights (plot 33 <-> 49, 34 <-> 50, ...,
# 40 <-> 56, same Plant_ID). Set it before sourcing to choose without editing this file:
#   swap_highlighted <- TRUE; source("data-raw/heights26.R")
# The script ends with a comparison of both versions against the neighbouring dates.
if (!exists("swap_highlighted")) swap_highlighted <- FALSE
stopifnot(isTRUE(swap_highlighted) || isFALSE(swap_highlighted))
swap_blocks <- function(d) {
  partner <- ifelse(d$Plot_ID %in% 33:40, d$Plot_ID + 16L,
                    ifelse(d$Plot_ID %in% 49:56, d$Plot_ID - 16L, d$Plot_ID))
  d$Height_cm <- d$Height_cm[match(paste(partner, d$Plant_ID), paste(d$Plot_ID, d$Plant_ID))]
  d
}

# ---- Height.xlsx ------------------------------------------------------------------------------
# Read everything as text: dates and heights come in mixed formats.
main <- as.data.frame(read_excel("data-raw/Height.xlsx", col_types = "text",
                                 .name_repair = "minimal"))
stopifnot(ncol(main) == 13,
          identical(names(main)[1:9], c("Date", "DAS", "Plot_ID", "Plant_ID", "Location",
                                        "Treatment", "Genotype", "Height_cm", "Notes")))
names(main)[10:13] <- c("X10", "X11", "X12", "X13")  # unnamed columns J-M
for (v in c("DAS", "Plot_ID", "Plant_ID")) main[[v]] <- as.integer(main[[v]])
# One Plant_ID is blank (DAS 56, plot 52); it sits between plants 3 and 5, so it is plant 4.
blank <- which(is.na(main$Plant_ID))
stopifnot(length(blank) == 1, main$DAS[blank] == 56, main$Plot_ID[blank] == 52,
          main$Plant_ID[blank + c(-1, 1)] == c(3, 5))
main$Plant_ID[blank] <- 4L
stopifnot(!anyNA(main[c("DAS", "Plot_ID", "Plant_ID")]))

# Dates. Text dates are d/m/Y, one with the typo "02026". Four dates were turned into Excel date
# numbers with day and month swapped (e.g. 12/5 became 5 December, serial 46361).
parse_date <- function(x) {
  out <- as.Date(rep(NA_character_, length(x)))
  serial <- grepl("^[0-9]+$", x)
  wrong <- as.Date(as.numeric(x[serial]), origin = "1899-12-30")
  out[serial] <- as.Date(format(wrong, "%Y-%d-%m"))  # swap day and month back
  out[!serial] <- as.Date(sub("/02026$", "/2026", x[!serial]), format = "%d/%m/%Y")
  out
}
main$Date <- parse_date(main$Date)
stopifnot(!anyNA(main$Date), all(main$Date == sowing + main$DAS))

# Column M only holds a note on the field layout of plot 19 (the data are in the right order).
stopifnot(identical(main$X13[!is.na(main$X13)],c("Note:", "Ploy 19", "wrong order", "1;2;5;4;3")))
main$X13 <- NULL

# Heights: "Dead" becomes NA; one value was typed with quotes ("49").
main$Height_cm[main$Height_cm %in% "Dead"] <- NA
main$Height_cm <- gsub('"', "", main$Height_cm)
stopifnot(!anyNA(suppressWarnings(as.numeric(na.omit(main$Height_cm)))))
main$Height_cm <- as.numeric(main$Height_cm)

# DAS 77 was not measured on the Loamy plots (33-64): their heights are empty, and these rows are
# left out. The unnamed columns J-L are not used.
is77L <- main$DAS == 77 & main$Plot_ID >= 33
stopifnot(sum(is77L) == 160, all(is.na(main$Height_cm[is77L])))
main <- main[!is77L, ]
key <- function(d) paste(d$Plot_ID, d$Plant_ID)

# ---- Height 0607.xlsx -------------------------------------------------------------------------
jul <- as.data.frame(read_excel("data-raw/Height 0607.xlsx",
                                col_types = c("date", "numeric", "numeric", "text", "numeric",
                                              "text", "text")))
stopifnot(nrow(jul) == 320, all(jul$DAS == 83), all(as.Date(jul$Date) == sowing + 83))
jul$Date <- as.Date(jul$Date)
for (v in c("DAS", "Plot_ID", "Plant_ID")) jul[[v]] <- as.integer(jul[[v]])
stopifnot(all(is.na(jul$Height_cm) | jul$Height_cm == "NA" | grepl("^[0-9.]+$", jul$Height_cm)))
jul$Height_cm[jul$Height_cm %in% "NA"] <- NA  # missing heights were typed as the text "NA"
jul$Height_cm <- as.numeric(jul$Height_cm)
# "Non irrigated" and "Non Irrigated" are the same treatment, called "NotIrrigated" in Height.xlsx.
stopifnot(all(jul$Treatment %in% c("Irrigated", "Non irrigated", "Non Irrigated")))
jul$Treatment[jul$Treatment != "Irrigated"] <- "NotIrrigated"

jul_recorded <- jul
if (swap_highlighted) jul <- swap_blocks(jul)
stopifnot(all.equal(swap_blocks(swap_blocks(jul_recorded)), jul_recorded))  # swapping twice = no swap

# Genotype, Location and Treatment are fixed per plot: take Genotype from Height.xlsx and check
# that the two files agree on the rest.
plots <- unique(main[c("Plot_ID", "Location", "Treatment", "Genotype")])
stopifnot(nrow(plots) == 64, !anyDuplicated(plots$Plot_ID))
m <- match(jul$Plot_ID, plots$Plot_ID)
stopifnot(!anyNA(m), all(jul$Location == plots$Location[m]),
          all(jul$Treatment == plots$Treatment[m]))
jul$Genotype <- plots$Genotype[m]

# ---- merge ------------------------------------------------------------------------------------
# Only the columns needed for the analysis; the soil type is called Location in the raw files.
# The date is left out: it is sowing + DAS (checked above).
cols <- c("DAS", "Genotype", "Plot_ID", "Plant_ID", "Location", "Treatment", "Height_cm")
FabaHeights26 <- rbind(main[cols], jul[cols])
names(FabaHeights26)[names(FabaHeights26) == "Location"] <- "Soil"
FabaHeights26 <- FabaHeights26[order(FabaHeights26$DAS, FabaHeights26$Plot_ID,
                                     FabaHeights26$Plant_ID), ]
rownames(FabaHeights26) <- NULL
n_per_DAS <- table(FabaHeights26$DAS)
stopifnot(nrow(FabaHeights26) == 14 * 64 * 5 - 160,
          all(n_per_DAS[names(n_per_DAS) != "77"] == 320), n_per_DAS[["77"]] == 160,
          !anyDuplicated(FabaHeights26[c("DAS", "Plot_ID", "Plant_ID")]))

# Plants without a height on a date are left out: blank or "NA" cells, "Dead" and "Missing plant".
missing <- is.na(FabaHeights26$Height_cm)
stopifnot(sum(missing) == 5) #We know of 5 missing plants, which comes from various sources such as missing or dead.
cat("Left out, no height recorded:\n")
print(FabaHeights26[missing, c("DAS", "Plot_ID", "Plant_ID")], row.names = FALSE)
FabaHeights26 <- FabaHeights26[!missing, ]
rownames(FabaHeights26) <- NULL
heights_long <- FabaHeights26  # with the raw Plant_ID (1-5), for the comparison at the end

# Plant_ID 1-5 is repeated in every plot; make it unique to the plant: plot_plant, e.g. "31_5".
FabaHeights26$Plant_ID <- interaction(FabaHeights26$Plot_ID, FabaHeights26$Plant_ID,
                                      sep = "_", lex.order = TRUE, drop = TRUE)
stopifnot(nlevels(FabaHeights26$Plant_ID) == 64 * 5,
          !anyDuplicated(FabaHeights26[c("DAS", "Plant_ID")]))
names(FabaHeights26)[names(FabaHeights26) == "Height_cm"] <- "Height"

write.csv(FabaHeights26, "data-raw/FabaHeights26.csv", row.names = FALSE)
save(FabaHeights26, file = "data/FabaHeights26.rda", compress = "xz")  # the package dataset

# Summary to check by eye: heights per date (rows) and soil x treatment (columns)
print(round(tapply(FabaHeights26$Height, list(FabaHeights26$DAS,
  paste(FabaHeights26$Soil, FabaHeights26$Treatment)), mean, na.rm = TRUE), 1))
cat("Plants measured per date:\n")
print(table(FabaHeights26$DAS))

# ---- which order of the highlighted blocks looks more realistic? -------------------------------
# Compares the DAS 83 heights with each plant's own heights at DAS 69 and DAS 90 (the Loamy plots
# have no DAS 77), for the recorded order and for the swapped one. A realistic DAS 83 value lies
# close to the straight line between DAS 69 and DAS 90 (small deviation, high correlation) and is
# not clearly below DAS 69 (plants do not shrink). The Loamy plots outside the blocks show what to
# expect from correctly entered data.
compare_orders <- function() {
  other <- heights_long[heights_long$DAS %in% c(69, 90), ]
  h <- function(das, d) other$Height_cm[other$DAS == das][
    match(paste(d$Plot_ID, d$Plant_ID), key(other[other$DAS == das, ]))]
  versions <- list(recorded = jul_recorded, swapped = swap_blocks(jul_recorded))
  out <- NULL
  for (v in names(versions)) {
    d <- versions[[v]]
    d <- d[d$Location == "Loamy", ]
    h69 <- h(69, d); h90 <- h(90, d)
    expected <- h69 + (83 - 69) / (90 - 69) * (h90 - h69)
    d$group <- ifelse(d$Plot_ID %in% 33:40, "plots 33-40",
                      ifelse(d$Plot_ID %in% 49:56, "plots 49-56", "other Loamy"))
    groups <- if (v == "recorded") unique(d$group) else c("plots 33-40", "plots 49-56")
    for (g in groups) {
      i <- d$group == g
      out <- rbind(out, data.frame(
        version = v, group = g, plants = sum(i),
        mean_abs_dev_cm = round(mean(abs(d$Height_cm[i] - expected[i]), na.rm = TRUE), 1),
        cor_expected = round(cor(d$Height_cm[i], expected[i], use = "complete.obs"), 2),
        below_DAS69_by_3cm = sum(d$Height_cm[i] < h69[i] - 3, na.rm = TRUE),
        above_DAS90_by_3cm = sum(d$Height_cm[i] > h90[i] + 3, na.rm = TRUE)))
    }
  }
  out$version[out$group == "other Loamy"] <- "(reference)"
  out[order(out$group, out$version), ]
}
cat("\nDAS 83: recorded vs swapped highlighted blocks (the csv uses",
    if (swap_highlighted) "the SWAPPED order)\n" else "the RECORDED order)\n")
print(compare_orders(), row.names = FALSE)
