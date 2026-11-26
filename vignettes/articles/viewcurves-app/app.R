# The viewcurves() app of the website article "Exploring observed curves", run in the browser by
# Shinylive (webR). export.R puts the source of viewcurves() (R/viewcurves.R) and the data next to this
# file, since gcTMB, with its compiled TMB model, cannot be installed in the browser; viewcurves()
# itself needs only shiny and ggplot2.
library(shiny)
library(ggplot2)
source("viewcurves.R", local=TRUE)
load("FabaHeights26.rda")

viewcurves(FabaHeights26, x="DAS", y="Height", id="Plant_ID", covariates=c("Genotype", "Soil"),
           output="html")
