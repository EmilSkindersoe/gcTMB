tmb_name <- "gcTMB_TMBExports"
tmb_flags <- commandArgs(trailingOnly = TRUE)

if(file.exists(paste0(tmb_name, ".cpp"))) {
  if(length(tmb_flags) == 0) tmb_flags <- ""
  # make only compares the .cpp with the object file, not the included .hpp and .h files with
  # the models: remove the old build, so that changes of the models are always compiled
  unlink(paste0(tmb_name, c(".o", .Platform$dynlib.ext)))
  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               safebounds = FALSE, safeunload = FALSE)
  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}

# cleanup done in ../Makevars[.win]
