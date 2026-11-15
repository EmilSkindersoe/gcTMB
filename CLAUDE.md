# Notes for Claude

- Don't `cat` whole large files. Read only the relevant part (Read with offset/limit, or grep), or write long command output to a file and show only the lines that matter.

- We have the following interpretation of all parameters throughout: b1 is the asymptotic value, b2 is the timing of maximum growth (timing of half saturation for Hill),
b3 is a scale parameter, b4 is a parameter related to the timing when the function can be strictly above 0 and b5 is a parameter for asymmetry.

- The vignette (`vignettes/gcTMB.rmd`) must always describe every growth function the package offers (formula, parameter meaning, inflection point, fraction of the asymptote reached at the inflection point, and maximum growth rate) and show each one fitted to the height data. When a growth model is added, update the vignette accordingly.