# Converts the raw data of the 2026 faba bean field trial to two package datasets, each also
# written as a csv to data-raw/:
#   - FabaHeights26 (data/FabaHeights26.rda): the two raw height files merged into one table
#       - data-raw/Height.xlsx      : all measurement dates except 6 July (13 dates)
#       - data-raw/Height 0607.xlsx : 6 July 2026 (DAS 83) only, without Genotype
#   - FabaFlowers26 (data/FabaFlowers26.rda): new flowers and pods per plant and date, from
#       - data-raw/flowers.xlsx     : counts per node, 8 dates from DAS 59 to 84
# Run from the package root: source("data-raw/DIGIFABA_field_26.R")
# Every cleaning step ends in a stopifnot() check, so the script stops if the raw data are not
# as described in the comments.

library(readxl)

sowing <- as.Date("2026-04-14")  # DAS 0; every date in both files is sowing + DAS

# The two highlighted blocks in Height 0607.xlsx (plots 33-40 and 49-56) may have been entered
# the wrong way round: it was not written down which block was which, so the recorded order is
# a best guess. swap_highlighted = TRUE swaps their heights (plot 33 <-> 49, 34 <-> 50, ...,
# 40 <-> 56, same Plant_ID). Set it before sourcing to choose without editing this file:
#   swap_highlighted <- TRUE; source("data-raw/DIGIFABA_field_26.R")
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

# ==== FabaFlowers26: flowers.xlsx ==============================================================
# The raw file has one row per node with flowers or pods: two Node_ID columns, one for the node
# of the flowers (NumFlowers) and one for the node of the pods (NumPods), so a row can describe
# two different nodes. Each node is counted once, on the first date it has flowers (or pods): no
# more flowers or pods appear on a node afterwards, and those already there fall off. So the
# counts of a date are the new flowers and pods, on the nodes that started flowering (or setting
# pods) since the previous date. A plant without new flowers or pods on a date has one row with
# blank or "NA" counts, or no row at all. The counts are summed over the nodes of each plant, and
# every plant gets a row on every date, with 0 where nothing was recorded.
fl <- as.data.frame(suppressMessages(read_excel("data-raw/flowers.xlsx")))
stopifnot(identical(names(fl), c("Date", "DAS", "Location", "Treatment", "Plot_ID", "Genotype",
                                 "Plant_ID", "Node_ID...8", "NumFlowers", "Node_ID...10",
                                 "NumPods")))
names(fl)[c(8, 10)] <- c("Node_flowers", "Node_pods")
for (v in c("DAS", "Plot_ID", "Plant_ID")) fl[[v]] <- as.integer(fl[[v]])
stopifnot(!anyNA(fl[c("DAS", "Plot_ID", "Plant_ID", "Location", "Treatment", "Genotype")]),
          all(as.Date(fl$Date) == sowing + fl$DAS),  # the date is sowing + DAS, so it is left out
          all(fl$Plant_ID %in% 1:3))

# "Not irrigated" is called "NotIrrigated" in Height.xlsx and FabaHeights26.
stopifnot(all(fl$Treatment %in% c("Irrigated", "Not irrigated")))
fl$Treatment[fl$Treatment == "Not irrigated"] <- "NotIrrigated"

# NumFlowers is text because of "NA" cells; blank and "NA" both mean no flowers counted. 28 rows
# have a flower node but no count; they are counted as 0 flowers as well.
stopifnot(all(fl$NumFlowers %in% c(NA, "NA", as.character(1:20))))
fl$NumFlowers <- suppressWarnings(as.integer(fl$NumFlowers))
stopifnot(sum(is.na(fl$NumFlowers) & !(is.na(fl$Node_flowers) | fl$Node_flowers == "NA")) == 28,
          all(is.na(fl$NumPods) == is.na(fl$Node_pods)))
fl$Flowers_counted <- !is.na(fl$NumFlowers)  # for FabaNodes26, where a blank count is missing
fl$NumFlowers[is.na(fl$NumFlowers)] <- 0L
fl$NumPods <- as.integer(fl$NumPods)
fl$NumPods[is.na(fl$NumPods)] <- 0L

# Each node is on a plant only once per date. (Plot 36 on DAS 77 had plant 2 with pod nodes 13-16
# twice; the second set belongs to plant 3 and is corrected in flowers.xlsx.)
fl$Excel_row <- seq_len(nrow(fl)) + 1L  # row in flowers.xlsx, below the header
fl_pods <- fl[!is.na(fl$Node_pods), ]
fl_flowers <- fl[!(is.na(fl$Node_flowers) | fl$Node_flowers == "NA"), ]
stopifnot(!anyDuplicated(fl_pods[c("DAS", "Plot_ID", "Plant_ID", "Node_pods")]),
          !anyDuplicated(fl_flowers[c("DAS", "Plot_ID", "Plant_ID", "Node_flowers")]))

# A node should be counted on one date only, but some are counted on several dates. Only the
# first date is kept; the counts on the later dates are removed. A raw row can hold both a flower
# node and a pod node, so only the repeated flower (or pod) entry is removed, not the whole row.
# The repeats are printed with their rows in flowers.xlsx, so they can be checked against the
# field sheets.
repeated_nodes <- function(d, node, count) {
  d <- d[order(d$Plot_ID, d$Plant_ID, as.integer(d[[node]]), d$DAS), ]
  key <- paste(d$Plot_ID, d$Plant_ID, d[[node]])
  later <- duplicated(key)
  d <- d[key %in% key[later], ]
  data.frame(Excel_row = d$Excel_row, DAS = d$DAS, Plant = paste(d$Plot_ID, d$Plant_ID, sep = "_"),
             Node = d[[node]], Count = d[[count]],
             Action = ifelse(later[key %in% key[later]], "REMOVED", "kept"))
}
rep_flowers <- repeated_nodes(fl_flowers, "Node_flowers", "NumFlowers")
rep_pods <- repeated_nodes(fl_pods, "Node_pods", "NumPods")
drop_f <- fl$Excel_row %in% rep_flowers$Excel_row[rep_flowers$Action == "REMOVED"]
drop_p <- fl$Excel_row %in% rep_pods$Excel_row[rep_pods$Action == "REMOVED"]
fl$Node_flowers[drop_f] <- NA
fl$NumFlowers[drop_f] <- 0L
fl$Node_pods[drop_p] <- NA
fl$NumPods[drop_p] <- 0L
if (nrow(rep_flowers) + nrow(rep_pods) > 0) {
  line <- strrep("!", 70)
  cat("\n", line, "\n!! flowers.xlsx: NODES COUNTED ON MORE THAN ONE DATE",
      "\n!! only the first date is kept, the later counts are removed\n", line, "\n", sep = "")
  if (nrow(rep_flowers) > 0) {
    cat("\nFlower nodes:\n")
    print(rep_flowers, row.names = FALSE)
  }
  if (nrow(rep_pods) > 0) {
    cat("\nPod nodes:\n")
    print(rep_pods, row.names = FALSE)
  }
  cat(line, "\n\n", sep = "")
  warning("flowers.xlsx: removed ", sum(drop_f), " flower and ", sum(drop_p),
          " pod counts of nodes already counted on an earlier date; see the list printed above",
          call. = FALSE)
}

# Soil, treatment and genotype are fixed per plot and agree with Height.xlsx.
fl_plots <- unique(fl[c("Plot_ID", "Location", "Treatment", "Genotype")])
m <- match(fl_plots$Plot_ID, plots$Plot_ID)
stopifnot(nrow(fl_plots) == 32, !anyDuplicated(fl_plots$Plot_ID), !anyNA(m),
          all(fl_plots$Location == plots$Location[m]),
          all(fl_plots$Treatment == plots$Treatment[m]),
          all(fl_plots$Genotype == plots$Genotype[m]))

# Sum over the nodes of each plant. Plants stop being recorded once no new flowers or pods appear
# (only 4 plants on DAS 84), so a plant without a row on a date has 0 new flowers and 0 new pods.
fl_sums <- aggregate(cbind(Flowers = NumFlowers, Pods = NumPods) ~ DAS + Plot_ID + Plant_ID,
                     data = fl, FUN = sum)
fl_grid <- expand.grid(DAS = sort(unique(fl$DAS)), Plot_ID = fl_plots$Plot_ID, Plant_ID = 1:3)
FabaFlowers26 <- merge(fl_grid, fl_sums, all.x = TRUE)
FabaFlowers26$Flowers[is.na(FabaFlowers26$Flowers)] <- 0L
FabaFlowers26$Pods[is.na(FabaFlowers26$Pods)] <- 0L
FabaFlowers26 <- merge(FabaFlowers26, fl_plots)
stopifnot(nrow(FabaFlowers26) == 8 * 32 * 3,
          sum(FabaFlowers26$Flowers) == sum(fl$NumFlowers),
          sum(FabaFlowers26$Pods) == sum(fl$NumPods))

# Same columns as FabaHeights26; Plant_ID 1-3 is repeated in every plot, so make it unique to the
# plant: plot_plant, e.g. "31_3".
names(FabaFlowers26)[names(FabaFlowers26) == "Location"] <- "Soil"
FabaFlowers26 <- FabaFlowers26[order(FabaFlowers26$DAS, FabaFlowers26$Plot_ID,
                                     FabaFlowers26$Plant_ID),
                               c("DAS", "Genotype", "Plot_ID", "Plant_ID", "Soil", "Treatment",
                                 "Flowers", "Pods")]
rownames(FabaFlowers26) <- NULL
FabaFlowers26$Plant_ID <- interaction(FabaFlowers26$Plot_ID, FabaFlowers26$Plant_ID,
                                      sep = "_", lex.order = TRUE, drop = TRUE)
stopifnot(nlevels(FabaFlowers26$Plant_ID) == 32 * 3,
          !anyDuplicated(FabaFlowers26[c("DAS", "Plant_ID")]))
# The plants are the same as in FabaHeights26 (plants 1-3 of the five measured in each plot), so
# the two datasets can be merged by Plant_ID.
stopifnot(all(levels(FabaFlowers26$Plant_ID) %in% levels(FabaHeights26$Plant_ID)))

# Running totals per plant: all flowers (pods) the plant has produced up to and including the date.
# The rows are sorted by DAS, so a cumulative sum within each plant follows time.
FabaFlowers26$FlowersTotal <- ave(FabaFlowers26$Flowers, FabaFlowers26$Plant_ID, FUN = cumsum)
FabaFlowers26$PodsTotal <- ave(FabaFlowers26$Pods, FabaFlowers26$Plant_ID, FUN = cumsum)
last <- FabaFlowers26$DAS == max(FabaFlowers26$DAS)
stopifnot(!is.unsorted(FabaFlowers26$DAS),
          all.equal(FabaFlowers26$FlowersTotal[last],
                    as.vector(tapply(FabaFlowers26$Flowers, FabaFlowers26$Plant_ID, sum))),
          all.equal(FabaFlowers26$PodsTotal[last],
                    as.vector(tapply(FabaFlowers26$Pods, FabaFlowers26$Plant_ID, sum))))

write.csv(FabaFlowers26, "data-raw/FabaFlowers26.csv", row.names = FALSE)
save(FabaFlowers26, file = "data/FabaFlowers26.rda", compress = "xz")  # the package dataset

# Summary to check by eye: mean new flowers and pods per plant per date (rows) and soil x treatment
cat("\nFabaFlowers26, mean new flowers, new pods, total flowers and total pods per plant:\n")
for (y in c("Flowers", "Pods", "FlowersTotal", "PodsTotal"))
  print(round(tapply(FabaFlowers26[[y]], list(FabaFlowers26$DAS,
    paste(FabaFlowers26$Soil, FabaFlowers26$Treatment)), mean), 1))

# ==== FabaNodes26: one row per node of each plant ==============================================
# From the same cleaned rows as FabaFlowers26 (repeated nodes removed above), the date and number
# of the flowers of each node and the date and number of its pods, for the pod set per node.
#   - a node that flowered but has no pod row set no pods: Pods = 0, DAS_pods = NA
#   - a node with a flower node but a blank count, or with pods but no flower row: Flowers = NA
#   - Warning = 1 where the record cannot be right: more pods than flowers, or the pods on the
#     same date as the flowers; 0 otherwise
nd_f <- fl[!(is.na(fl$Node_flowers) | fl$Node_flowers == "NA"), ]
nd_f <- data.frame(Plot_ID = nd_f$Plot_ID, Plant_ID = nd_f$Plant_ID,
                   Node = as.integer(nd_f$Node_flowers), DAS_flowers = nd_f$DAS,
                   Flowers = ifelse(nd_f$Flowers_counted, nd_f$NumFlowers, NA_integer_))
nd_p <- fl[!is.na(fl$Node_pods), ]
nd_p <- data.frame(Plot_ID = nd_p$Plot_ID, Plant_ID = nd_p$Plant_ID,
                   Node = as.integer(nd_p$Node_pods), DAS_pods = nd_p$DAS, Pods = nd_p$NumPods)
stopifnot(!anyNA(nd_f$Node), !anyNA(nd_p$Node),
          !anyDuplicated(nd_f[c("Plot_ID", "Plant_ID", "Node")]),  # repeats were removed above
          !anyDuplicated(nd_p[c("Plot_ID", "Plant_ID", "Node")]))
FabaNodes26 <- merge(nd_f, nd_p, all = TRUE)
FabaNodes26$Pods[is.na(FabaNodes26$DAS_pods)] <- 0L
FabaNodes26$Warning <- as.integer(
  (!is.na(FabaNodes26$Flowers) & FabaNodes26$Pods > FabaNodes26$Flowers) |
  (!is.na(FabaNodes26$DAS_pods) & FabaNodes26$DAS_pods %in% FabaNodes26$DAS_flowers &
     FabaNodes26$DAS_pods == FabaNodes26$DAS_flowers))
FabaNodes26$Warning[is.na(FabaNodes26$Warning)] <- 0L
FabaNodes26 <- merge(FabaNodes26, fl_plots)
names(FabaNodes26)[names(FabaNodes26) == "Location"] <- "Soil"
FabaNodes26 <- FabaNodes26[order(FabaNodes26$Plot_ID, FabaNodes26$Plant_ID, FabaNodes26$Node),
                           c("Genotype", "Plot_ID", "Plant_ID", "Soil", "Treatment", "Node",
                             "DAS_flowers", "Flowers", "DAS_pods", "Pods", "Warning")]
rownames(FabaNodes26) <- NULL
FabaNodes26$Plant_ID <- factor(paste(FabaNodes26$Plot_ID, FabaNodes26$Plant_ID, sep = "_"),
                               levels = levels(FabaFlowers26$Plant_ID))
stopifnot(!anyNA(FabaNodes26$Plant_ID),
          sum(FabaNodes26$Flowers, na.rm = TRUE) == sum(fl$NumFlowers),
          sum(FabaNodes26$Pods) == sum(fl$NumPods),
          all(FabaNodes26$Warning %in% 0:1))

write.csv(FabaNodes26, "data-raw/FabaNodes26.csv", row.names = FALSE)
save(FabaNodes26, file = "data/FabaNodes26.rda", compress = "xz")  # the package dataset

cat("\nFabaNodes26:", nrow(FabaNodes26), "nodes;", sum(is.na(FabaNodes26$Flowers)),
    "without a flower count;", sum(FabaNodes26$Warning), "with Warning = 1:\n")
print(FabaNodes26[FabaNodes26$Warning == 1, c("Plant_ID", "Node", "DAS_flowers", "Flowers",
                                               "DAS_pods", "Pods")], row.names = FALSE)
