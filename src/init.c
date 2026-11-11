// Registration of native routines for the package library gcTMB. The TMB templates are compiled
// separately into gcTMB_TMBExports (see Makevars and TMB/compile.R), which registers its own
// routines through TMB_LIB_INIT. This library has no routines of its own; it is the package
// library R CMD check expects in src, and turns off dynamic symbol lookup.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_gcTMB(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
