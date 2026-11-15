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
//   mod 7 linear:   f(t) = f1 t + f2, for comparison; f3-f5 are not used, and it is not truncated
// f1 is the asymptote. f2 is the time of the maximum growth rate (the inflection point) for logistic,
// gompertz, weibull, normal and gamma, and the time of half the asymptote for hill. f3 is the scale on
// the time axis (logistic, gompertz, normal), positive for an increasing curve, the slope of the hill
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

// growth curve f(t) = f1 G(z(t)), with the asymmetry f5. trunc = 1: G truncated at the threshold f4,
// as a distribution function truncated from below: (G(z) - G(z4)) / (1 - G(z4)) = 1 - S(z) / S(z4)
// for t > f4 and 0 for t <= f4, z4 = z(f4), S = 1 - G; computed by the survival function for
// thresholds far beyond f2. The time t is data, so the models on positive times can branch on t <= 0.
// The linear model is f1 t + f2.
template<class Type>
Type gcmean(int mod, int trunc, Type t, Type f1, Type f2, Type f3, Type f4, Type f5){
  if (mod == 7) return f1 * t + f2;
  if (gclogtime(mod) && asDouble(t) <= 0) return Type(0) * f1;
  Type z = gcz(mod, t, f2, f3);
  if (trunc == 1){
    // both branches are evaluated: below the threshold, the ratio of the survival functions is finite
    Type Ft = 1 - exp(gclogSG(mod, z, f3, f5) - gclogS4(mod, f2, f3, f4, f5));
    return f1 * CppAD::CondExpGt(t, f4, Ft, Type(0));
  }
  return f1 * exp(gclogG(mod, z, f3, f5));
}

// time at which the curve reaches the fraction q of the asymptote f1, the inverse of gcmean:
// t(Ginv(q)); r: shift on the response scale (residual). Truncated at f4 (trunc = 1):
// 1 - S(z) / S(z4) = q gives t(Sinv((1 - q) S(z4))), on the log scale of S. The linear model has no
// asymptote: q f1 is the response, as for the absolute times, reached at (q f1 - r - f2) / f1.
template<class Type>
Type gctime(int mod, int trunc, Type q, Type f1, Type f2, Type f3, Type f4, Type f5, Type r){
  Type qr = q - r / f1;
  if (mod == 7) return qr - f2 / f1;
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
