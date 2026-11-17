#include "gcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// time at which a growth curve reaches a level of its asymptote, for use with
// MakeADFun(..., ADreport=TRUE): fn() returns the time, gr() its gradient w.r.t. (b1-b5, u1-u5, r)
template<class Type>
Type gcTime(objective_function<Type>* obj) {
  DATA_INTEGER(mod);
  DATA_INTEGER(trunc);    // 1: curve truncated at the threshold b4
  DATA_INTEGER(b1rate);   // 1: b1 is the maximum growth rate instead of the asymptote
  DATA_IVECTOR(lnk);      // links of b1-b5: 1 = identity, 2 = inverse, 3 = log, 4 = logit
  DATA_INTEGER(absolute); // 0: q is the fraction of the asymptote, 1: q is the response
  DATA_SCALAR(q);
  PARAMETER_VECTOR(b);    // b1-b5
  PARAMETER_VECTOR(u);    // random effects u1-u5, evaluated at zero
  PARAMETER(r);           // residual, evaluated at zero

  vector<Type> eta = b + u;
  Type f2 = linkinv(lnk(1), eta(1));
  Type f3 = linkinv(lnk(2), eta(2));
  Type f4 = linkinv(lnk(3), eta(3));
  Type f5 = linkinv(lnk(4), eta(4));
  Type f1 = gcasym(b1rate, mod, trunc, linkinv(lnk(0), eta(0)), f2, f3, f4, f5);
  // fraction of the asymptote at which the curve reaches the response q
  Type qf = absolute == 1 ? Type(q / f1) : q;
  Type time = gctime(mod, trunc, qf, f1, f2, f3, f4, f5, r);
  ADREPORT(time);
  return 0;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
