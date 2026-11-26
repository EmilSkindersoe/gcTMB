# Exports the viewcurves() app of the article "Exploring observed curves" with Shinylive into the
# built website, as <site>/viewcurves-app/, which the article embeds. Run from the package folder
# after building the site, as the pkgdown workflow does:
#   Rscript vignettes/articles/viewcurves-app/export.R docs
# The exported app (a few MB of Shinylive files) goes only into the website, not into the repository.
args <- commandArgs(trailingOnly=TRUE)
site <- if (length(args) > 0) args[1] else "docs"
if (!dir.exists(site)) stop("the site folder ", site, " does not exist; build the site first.")

app <- file.path(tempdir(), "viewcurves-app")
unlink(app, recursive=TRUE)
dir.create(app)
stopifnot(file.copy(c("vignettes/articles/viewcurves-app/app.R", "R/viewcurves.R", "data/FabaHeights26.rda"), app))
shinylive::export(app, file.path(site, "viewcurves-app"))
