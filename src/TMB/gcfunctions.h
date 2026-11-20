// Helper functions shared by the gcTMB objective functions.
// Kept in a .h file (not .hpp) so TMBtools does not treat it as a model.

#ifndef GCFUNCTIONS_H
#define GCFUNCTIONS_H

// inverse link of a parameter from its linear predictor eta
// lnk: 1 = identity, 2 = inverse (the linear predictor of the scale b3 is then the rate 1/b3),
//      3 = log, 4 = logit
template<class Type>
Type linkinv(int lnk, Type eta){
  switch (lnk){
  case 2:
    return 1/eta;
  case 3:
    return exp(eta);
  case 4:
    return invlogit(eta);
  default:
    return eta;
  }
}

// The growth models are f(t) = f1 F(z(t)), with a distribution function F of a base distribution and a
// transformed time z(t), increasing in t:
//   mod 1 logistic: F logistic,          z = (t - f2) / f3
//   mod 2 gompertz: F gompertz,          z = (t - f2) / f3
//   mod 3 hill:     F logistic,          z = f3 (log t - log f2), f(t) = f1 / (1 + (f2 / t)^f3)
//   mod 4 weibull:  F complementary log-log, z = f3 (log t - log f2) + log((f3 - 1) / f3),
//                   f(t) = f1 (1 - exp(-(f3 - 1) / f3 (t / f2)^f3))
//   mod 5 normal:   F standard normal,   z = (t - f2) / f3
//   mod 6 gamma:    F gamma with shape f3, z = t (f3 - 1) / f2, the gamma distribution function with
//                   shape f3 and scale f2 / (f3 - 1), whose density has its mode at f2
//   mod 7 linear:   f(t) = f1 t, for comparison, plus the intercept f6 added by the templates; f2-f5 are
//                   not used, and it is not truncated
//   mod 8-9 chapmanrichards and beta: curves that start at f4, see gcstart below
// The templates add the additive parameter f6 to the curve, f6 + f(t) (the intercept of the linear model).
// f1 is the asymptote. f2 is the time of the maximum growth rate (the inflection point) for logistic,
// gompertz, weibull, normal and gamma, and the time of half the asymptote for hill. f3 is the scale on
// the time axis (logistic, gompertz, normal), positive, the slope of the hill
// curve on the log time axis, or the shape of the weibull and gamma curves, larger than 1. hill,
// weibull and gamma are zero up to t = 0.
// f5 is the asymmetry of the curve, with f5 = 1 the curves above; see gclogG below.
// The parameters are on their original scale, after the inverse links.

// base distribution of a model: 1 = logistic, F(z) = 1 / (1 + exp(-z)); 2 = gompertz,
// F(z) = exp(-exp(-z)); 3 = complementary log-log, F(z) = 1 - exp(-exp(z)); 4 = standard normal;
// 5 = gamma with shape sh
inline int gcbase(int mod){
  switch (mod){
  case 2:
    return 2;
  case 4:
    return 3;
  case 5:
    return 4;
  case 6:
    return 5;
  default:
    return 1;
  }
}

// models for positive times, which are zero up to t = 0
inline bool gclogtime(int mod){
  return mod == 3 || mod == 4 || mod == 6;
}

// transformed time z(t) of a model, and its inverse t(z)
template<class Type>
Type gcz(int mod, Type t, Type f2, Type f3){
  switch (mod){
  case 3:
    return f3 * (log(t) - log(f2));
  case 4:
    return f3 * (log(t) - log(f2)) + log((f3 - 1) / f3);
  case 6:
    return t * (f3 - 1) / f2;
  default:
    return (t - f2) / f3;
  }
}

template<class Type>
Type gczinv(int mod, Type z, Type f2, Type f3){
  switch (mod){
  case 3:
    return f2 * exp(z / f3);
  case 4:
    return f2 * exp((z - log((f3 - 1) / f3)) / f3);
  case 6:
    return z * f2 / (f3 - 1);
  default:
    return f2 + f3 * z;
  }
}

// log of the survival function S(z) = 1 - F(z) of a base distribution, for the truncated curves; sh is
// the shape of the gamma base
template<class Type>
Type gclogS(int base, Type z, Type sh){
  switch (base){
  case 2:
    // 1 - exp(-exp(-z))
    return logspace_sub(Type(0), Type(-exp(-z)));
  case 3:
    // exp(-exp(z))
    return -exp(z);
  case 4:
    return log(pnorm(Type(-z)));
  case 5:
    return logspace_sub(Type(0), Type(log(pgamma(z, sh))));
  default:
    // 1 / (1 + exp(z))
    return -logspace_add(Type(0), z);
  }
}

// inverse of the log survival function: z with log S(z) = ls
template<class Type>
Type gcSinvlog(int base, Type ls, Type sh){
  switch (base){
  case 2:
    return -log(-logspace_sub(Type(0), ls));
  case 3:
    return log(-ls);
  case 4:
    return -qnorm(Type(exp(ls)));
  case 5:
    return qgamma(Type(1 - exp(ls)), sh);
  default:
    return logspace_sub(Type(-ls), Type(0));
  }
}

// log of the distribution function F of a base distribution
template<class Type>
Type gclogF(int base, Type z, Type sh){
  switch (base){
  case 2:
    return -exp(-z);
  case 3:
    return logspace_sub(Type(0), Type(-exp(z)));
  case 4:
    return log(pnorm(z));
  case 5:
    return log(pgamma(z, sh));
  default:
    return -logspace_add(Type(0), Type(-z));
  }
}

// inverse of the log distribution function: z with log F(z) = lf (the gompertz base is only used
// through its survival function, see gclogG)
template<class Type>
Type gcFinvlog(int base, Type lf, Type sh){
  switch (base){
  case 2:
    return -log(-lf);
  case 3:
    return log(-logspace_sub(Type(0), lf));
  case 4:
    return qnorm(Type(exp(lf)));
  case 5:
    return qgamma(Type(exp(lf)), sh);
  default:
    return -logspace_sub(Type(-lf), Type(0));
  }
}

// Asymmetry f5 (b5 = 1: the symmetric curve of the base distribution). The curve is f1 G(z), with
//   logistic (Richards): G(z) = (1 + f5 exp(-z))^(-1/f5) = F(z - log f5)^(1/f5)
//   hill: G(z) = F(z)^f5, f(t) = f1 / (1 + (f2 / t)^f3)^f5
//   weibull, normal, gamma: G(z) = F(z)^f5
//   gompertz: G(z) = 1 - (1 - F(z))^f5, a power of the survival function, as a power of the gompertz
//     F is only a shift of z
// so that log G = p log F(z - s) (power of F) or log(1 - G) = p log S(z) (power of S). f3 is passed
// on as the shape of the gamma base.
inline bool gcpowS(int mod){
  return mod == 2;
}

template<class Type>
Type gcpow(int mod, Type f5){
  return mod == 1 ? Type(1 / f5) : f5;
}

template<class Type>
Type gcshift(int mod, Type f5){
  return mod == 1 ? Type(log(f5)) : Type(0);
}

// log G(z) and log(1 - G(z)) of the asymmetric curve
template<class Type>
Type gclogG(int mod, Type z, Type f3, Type f5){
  int base = gcbase(mod);
  if (gcpowS(mod)) return logspace_sub(Type(0), Type(f5 * gclogS(base, z, f3)));
  return gcpow(mod, f5) * gclogF(base, Type(z - gcshift(mod, f5)), f3);
}

template<class Type>
Type gclogSG(int mod, Type z, Type f3, Type f5){
  int base = gcbase(mod);
  if (gcpowS(mod)) return f5 * gclogS(base, z, f3);
  return logspace_sub(Type(0), Type(gcpow(mod, f5) * gclogF(base, Type(z - gcshift(mod, f5)), f3)));
}

// z with log(1 - G(z)) = ls, and z with G(z) = q
template<class Type>
Type gcSGinvlog(int mod, Type ls, Type f3, Type f5){
  int base = gcbase(mod);
  if (gcpowS(mod)) return gcSinvlog(base, Type(ls / f5), f3);
  return gcshift(mod, f5) + gcFinvlog(base, Type(logspace_sub(Type(0), ls) / gcpow(mod, f5)), f3);
}

template<class Type>
Type gcGinv(int mod, Type q, Type f3, Type f5){
  int base = gcbase(mod);
  if (gcpowS(mod)) return gcSinvlog(base, Type(log(1 - q) / f5), f3);
  return gcshift(mod, f5) + gcFinvlog(base, Type(log(q) / gcpow(mod, f5)), f3);
}

// log(1 - G) at the threshold f4 of a truncated curve; for the models on positive times, a threshold
// up to zero does not truncate the curve (1 - G = 1). Both branches of a conditional expression are
// evaluated, the other one at a valid time.
template<class Type>
Type gclogS4(int mod, Type f2, Type f3, Type f4, Type f5){
  if (!gclogtime(mod)) return gclogSG(mod, gcz(mod, f4, f2, f3), f3, f5);
  Type f4p = CppAD::CondExpGt(f4, Type(0), f4, Type(1));
  return CppAD::CondExpGt(f4, Type(0), gclogSG(mod, gcz(mod, f4p, f2, f3), f3, f5), Type(0));
}

// Models that start at the threshold f4, where the curve itself is zero (it is zero up to f4 as well),
// so that a truncation at f4 changes nothing (they are always fitted as truncated curves); they do
// not use f5:
//   mod 8 chapmanrichards: f(t) = f1 (1 - exp(-(t - f4) / f3))^c, c = exp((f2 - f4) / f3), with the
//     inflection point f2 for c > 1 (f2 > f4), the scale f3 and the start f4
//   mod 9 beta (Yin et al. 2003): f(t) = f1 (1 + (te - t) / (te - f2)) ((t - f4) / (te - f4))^((te - f4) / (te - f2))
//     for f4 < t < te, and f1 from the end of growth te = f2 + f3 on: the time f2 of the maximum
//     growth rate, the duration f3 of the growth after it, and the start f4. With u = (t - f4) / D,
//     D = te - f4, and d = f3 / D: G(u) = (1 + (1 - u) / d) u^(1 / d), 0 < u < 1.
inline bool gcstart(int mod){
  return mod >= 8;
}

// exponent c of the chapmanrichards curve
template<class Type>
Type gccr(Type f2, Type f3, Type f4){
  return exp((f2 - f4) / f3);
}

// fraction G of the asymptote reached at time t by a curve starting at f4, for t > f4; t is clamped
// to just after the start (and to the end of growth of the beta curve), so that both branches of the
// conditional expressions in gcmean are finite
template<class Type>
Type gcGstart(int mod, Type t, Type f2, Type f3, Type f4){
  Type eps = Type(1e-10);
  if (mod == 9){
    Type D = f2 + f3 - f4;
    Type d = f3 / D;
    Type u = (t - f4) / D;
    u = CppAD::CondExpGt(u, eps, u, eps);
    u = CppAD::CondExpLt(u, Type(1), u, Type(1));
    return (1 + (1 - u) / d) * exp(log(u) / d);
  }
  Type z = (t - f4) / f3;
  z = CppAD::CondExpGt(z, eps, z, eps);
  return exp(gccr(f2, f3, f4) * logspace_sub(Type(0), Type(-z)));
}

// time at which a curve starting at f4 reaches the fraction q of its asymptote. The beta curve has no
// closed form: Newton steps from the root found in double precision by bisection, whose derivatives
// are those of the implicit function (the starting point is a constant of the tape, which is only
// evaluated at the parameters it was taped at, as in gcTime)
template<class Type>
Type gctimestart(int mod, Type q, Type f2, Type f3, Type f4){
  if (mod != 9){
    Type z = -logspace_sub(Type(0), Type(log(q) / gccr(f2, f3, f4)));
    return f4 + f3 * z;
  }
  Type D = f2 + f3 - f4;
  Type d = f3 / D;
  double dd = asDouble(d), qd = asDouble(q);
  double lo = 0, hi = 1;
  for (int k = 0; k < 100; k++){
    double m = (lo + hi) / 2;
    double g = (1 + (1 - m) / dd) * pow(m, 1 / dd);
    if (g < qd) lo = m; else hi = m;
  }
  Type u = Type((lo + hi) / 2);
  for (int k = 0; k < 3; k++){
    Type g = (1 + (1 - u) / d) * exp(log(u) / d);
    Type dg = exp((1 / d - 1) * log(u)) / d * (1 - u) * (1 + 1 / d);
    u = u - (g - q) / dg;
  }
  return f4 + D * u;
}

// maximum growth rate of a curve starting at f4 with asymptote 1: at the inflection point f2 of the
// chapmanrichards curve, (1 - 1/c)^(c - 1) / f3 (guarded for curves without an inflection point after
// the start); at f2 for the beta curve, (1 + 1/d) (1 - d)^(1/d - 1) / D
template<class Type>
Type gcrate1start(int mod, Type f2, Type f3, Type f4){
  if (mod == 9){
    Type D = f2 + f3 - f4;
    Type d = f3 / D;
    return (1 + 1 / d) * exp((1 / d - 1) * log(1 - d)) / D;
  }
  Type c = gccr(f2, f3, f4);
  Type cc = CppAD::CondExpGt(c, Type(1), c, Type(2));
  return CppAD::CondExpGt(c, Type(1), Type(exp((cc - 1) * log(1 - 1 / cc)) / f3), Type(1 / f3));
}

// growth curve f(t) = f1 G(z(t)), with the asymmetry f5. trunc = 1: G truncated at the threshold f4,
// as a distribution function truncated from below: (G(z) - G(z4)) / (1 - G(z4)) = 1 - S(z) / S(z4)
// for t > f4 and 0 for t <= f4, z4 = z(f4), S = 1 - G; computed by the survival function for
// thresholds far beyond f2. The time t is data, so the models on positive times can branch on t <= 0.
// The linear model is f1 t (plus f6).
template<class Type>
Type gcmean(int mod, int trunc, Type t, Type f1, Type f2, Type f3, Type f4, Type f5){
  if (mod == 7) return f1 * t;
  if (gcstart(mod)) return f1 * CppAD::CondExpGt(t, f4, gcGstart(mod, t, f2, f3, f4), Type(0));
  if (gclogtime(mod) && asDouble(t) <= 0) return Type(0) * f1;
  Type z = gcz(mod, t, f2, f3);
  if (trunc == 1){
    // both branches are evaluated: below the threshold, the ratio of the survival functions is finite
    Type Ft = 1 - exp(gclogSG(mod, z, f3, f5) - gclogS4(mod, f2, f3, f4, f5));
    return f1 * CppAD::CondExpGt(t, f4, Ft, Type(0));
  }
  return f1 * exp(gclogG(mod, z, f3, f5));
}

// maximum growth rate of the curve with asymptote f1 = 1, for b1 as the maximum growth rate
// (b1rate = 1): the asymptote is then b1 / gcrate1. Closed forms at the inflection point tm of the
// curve without truncation, for the asymmetry f5 of the logistic (Richards) and hill curves, and
// f5 = 1 for the others (checked in R). Truncated at f4 (trunc = 1): the curve without truncation
// scaled by 1 / S(z4) beyond f4; if the inflection point is not after the threshold, the curve grows
// fastest just after it, with the slope at f4 by a central difference (evaluated at a valid time,
// as both branches of the conditional expression are evaluated). The linear model has rate b1.
template<class Type>
Type gcrate1(int mod, int trunc, Type f2, Type f3, Type f4, Type f5){
  if (mod == 7) return Type(1);
  if (gcstart(mod)) return gcrate1start(mod, f2, f3, f4);
  Type tm = f2;
  Type r0;
  switch (mod){
  case 1:
    r0 = pow(Type(1 + f5), Type(-(1 + f5) / f5)) / f3;
    break;
  case 2:
    r0 = exp(Type(-1)) / f3;
    break;
  case 3: {
    // v = (f2 / tm)^f3 at the inflection point
    Type v = (f3 + 1) / (f3 * f5 - 1);
    tm = f2 * pow(v, Type(-1 / f3));
    r0 = f3 * f5 * v * pow(Type(1 + v), Type(-f5 - 1)) / tm;
    break;
  }
  case 4: {
    Type cc = (f3 - 1) / f3;
    r0 = f3 * cc * exp(-cc) / f2;
    break;
  }
  case 5:
    // standard normal density at 0, 1 / sqrt(2 pi)
    r0 = Type(0.3989422804014327) / f3;
    break;
  default:
    // gamma with shape f3 and scale f2 / (f3 - 1): the density of z = t (f3 - 1) / f2 at its mode
    r0 = dgamma(Type(f3 - 1), f3, Type(1), false) * (f3 - 1) / f2;
  }
  if (trunc != 1) return r0;
  Type h = Type(1e-4);
  Type t4 = gclogtime(mod) ? Type(CppAD::CondExpGt(f4, Type(2 * h), f4, Type(1))) : f4;
  Type g4 = (exp(gclogG(mod, gcz(mod, Type(t4 + h), f2, f3), f3, f5)) -
             exp(gclogG(mod, gcz(mod, Type(t4 - h), f2, f3), f3, f5))) / (2 * h);
  return CppAD::CondExpGt(tm, f4, r0, g4) / exp(gclogS4(mod, f2, f3, f4, f5));
}

// asymptote of a curve from f1, which is the asymptote (b1rate = 0) or the maximum growth rate
// (b1rate = 1)
template<class Type>
Type gcasym(int b1rate, int mod, int trunc, Type f1, Type f2, Type f3, Type f4, Type f5){
  if (b1rate != 1) return f1;
  return f1 / gcrate1(mod, trunc, f2, f3, f4, f5);
}

// time at which the curve reaches the fraction q of the asymptote f1, the inverse of gcmean:
// t(Ginv(q)); r: shift on the response scale (residual). Truncated at f4 (trunc = 1):
// 1 - S(z) / S(z4) = q gives t(Sinv((1 - q) S(z4))), on the log scale of S. The linear model has no
// asymptote: q f1 is the response, as for the absolute times, reached at (q f1 - r - f2) / f1.
template<class Type>
Type gctime(int mod, int trunc, Type q, Type f1, Type f2, Type f3, Type f4, Type f5, Type r){
  Type qr = q - r / f1;
  if (mod == 7) return qr;
  if (gcstart(mod)) return gctimestart(mod, qr, f2, f3, f4);
  Type z;
  if (trunc == 1){
    z = gcSGinvlog(mod, Type(log(1 - qr) + gclogS4(mod, f2, f3, f4, f5)), f3, f5);
  } else {
    z = gcGinv(mod, qr, f3, f5);
  }
  return gczinv(mod, z, f2, f3);
}

// mean f clamped into the open unit interval, for the beta and beta-binomial log-likelihoods
template<class Type>
Type unitclamp(Type f){
  Type lo = Type(1e-12), hi = Type(1 - 1e-12);
  return CppAD::CondExpLt(f, lo, lo, CppAD::CondExpGt(f, hi, hi, f));
}

// -Inf for a mean f outside the open unit interval, else ll; ll is computed at the clamped mean,
// as the beta functions are finite for most negative arguments (lgamma), and the derivatives of
// both branches of a conditional expression are evaluated
template<class Type>
Type inunit(Type f, Type ll){
  Type ninf = Type(-std::numeric_limits<double>::infinity());
  return CppAD::CondExpLe(f, Type(0), ninf, CppAD::CondExpGe(f, Type(1), ninf, ll));
}

// log-likelihood contribution of a single observation
// fam: 1 = gaussian (phi = sd), 2 = poisson, 3 = binomial (bn = number of trials),
//      4 = beta (phi = precision), 5 = nbinom2 (variance f (1 + f/phi), phi = size),
//      6 = nbinom1 (variance f (1 + phi)), 7 = betabinomial (bn = number of trials, phi = precision),
//      8 = lognormal (median f, phi = sd of log(y))
// A mean outside the range of the family gives -Inf (beta, betabinomial) or NaN (log of a negative mean
// or probability), which the optimisers treat as invalid.
template<class Type>
Type loglik0(int fam, Type y, Type f, Type bn, Type phi){
  switch (fam){
  case 1:
    return dnorm(y, f, phi, true);
  case 2:
    return dpois(y, f, true);
  case 3:
    return dbinom(y, bn, f, true);
  case 4: {
    Type fc = unitclamp(f);
    return inunit(f, Type(lgamma(phi) - lgamma(fc*phi) - lgamma((1-fc)*phi) + (fc*phi-1)*log(y) + ((1-fc)*phi - 1)*log(1 - y)));
  }
  case 5:
    // on the log scale of the mean and of the variance minus the mean, for numerical stability
    return dnbinom_robust(y, log(f), Type(2*log(f) - log(phi)), 1);
  case 6:
    return dnbinom_robust(y, log(f), Type(log(f) + log(phi)), 1);
  case 7: {
    // binomial with a beta distributed probability, with mean f and precision phi
    Type fc = unitclamp(f);
    return inunit(f, Type(lgamma(bn + 1) - lgamma(y + 1) - lgamma(bn - y + 1)
      + lgamma(y + fc*phi) + lgamma(bn - y + (1-fc)*phi) - lgamma(bn + phi)
      - lgamma(fc*phi) - lgamma((1-fc)*phi) + lgamma(phi)));
  }
  case 8:
    // lognormal with median f and standard deviation phi of log(y); density of y, with the Jacobian
    return dnorm(log(y), log(f), phi, true) - log(y);
  default:
    return 0;
  }
}

// log-likelihood, with a mean of exactly zero for the count families (poisson, binomial, nbinom2,
// nbinom1, betabinomial), as a truncated curve before its threshold: a zero count has probability
// one and a positive count zero. The log-likelihood of the other branch is computed at a valid mean,
// as the derivatives of both branches of a conditional expression are evaluated.
template<class Type>
Type loglik(int fam, Type y, Type f, Type bn, Type phi){
  if (fam == 2 || fam == 3 || fam == 5 || fam == 6 || fam == 7){
    Type ninf = Type(-std::numeric_limits<double>::infinity());
    Type ll = loglik0(fam, y, Type(CppAD::CondExpEq(f, Type(0), Type(0.5), f)), bn, phi);
    return CppAD::CondExpEq(f, Type(0), Type(CppAD::CondExpEq(y, Type(0), Type(0), ninf)), ll);
  }
  return loglik0(fam, y, f, bn, phi);
}

// simulate a single observation with mean f, same parameterisation as loglik; a mean of zero gives zero
template<class Type>
Type simobs(int fam, Type f, Type bn, Type phi){
  if (fam != 1 && asDouble(f) == 0) return Type(0);
  switch (fam){
  case 1:
    return rnorm(f, phi);
  case 2:
    return rpois(f);
  case 3:
    return rbinom(bn, f);
  case 4:
    return rbeta(f*phi, (1-f)*phi);
  case 5:
    return rnbinom2(f, Type(f*(1 + f/phi)));
  case 6:
    return rnbinom2(f, Type(f*(1 + phi)));
  case 7:
    return rbinom(bn, rbeta(f*phi, (1-f)*phi));
  case 8:
    return exp(rnorm(log(f), phi));
  default:
    return 0;
  }
}

#endif
