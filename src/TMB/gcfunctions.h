// Helper functions shared by the gcTMB objective functions.
// Kept in a .h file (not .hpp) so TMBtools does not treat it as a model.

#ifndef GCFUNCTIONS_H
#define GCFUNCTIONS_H

// inverse link of a parameter from its linear predictor eta
// lnk: 1 = identity, 2 = inverse (the linear predictor of the scale b3 is then the rate 1/b3)
template<class Type>
Type linkinv(int lnk, Type eta){
  switch (lnk){
  case 2:
    return 1/eta;
  default:
    return eta;
  }
}

// growth curve: f(t) = f1 * F((t - f2) / f3), with F increasing from 0 to 1
// mod: 1 = logistic, F(z) = 1 / (1 + exp(-z)); 2 = gompertz, F(z) = exp(-exp(-z))
// f1 is the asymptote, f2 the time of the maximum growth rate (the inflection point of the curve), and
// f3 the scale on the time axis, positive for an increasing curve. The parameters are on their
// original scale, after the inverse links. f4 and f5 are reserved for future models and not used.
template<class Type>
Type gcmean(int mod, Type t, Type f1, Type f2, Type f3, Type f4, Type f5){
  Type z = (t - f2) / f3;
  Type F = 0;
  switch (mod){
  case 1:
    // 1 / (1 + exp(-z)), on the log scale for numerical stability
    F = exp(-logspace_add(Type(0), Type(-z)));
    break;
  case 2:
    F = exp(-exp(-z));
    break;
  }
  return f1 * F;
}

// inverse of F at level q
template<class Type>
Type gcFinv(int mod, Type q){
  switch (mod){
  case 1:
    return log(q / (1 - q));
  case 2:
    return -log(-log(q));
  default:
    return 0;
  }
}

// time at which the curve reaches the fraction q of the asymptote f1, the inverse of gcmean:
// f2 + f3 Finv(q); r: shift on the response scale (residual)
template<class Type>
Type gctime(int mod, Type q, Type f1, Type f2, Type f3, Type f4, Type f5, Type r){
  Type qr = q - r / f1;
  return f2 + f3 * gcFinv(mod, qr);
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
Type loglik(int fam, Type y, Type f, Type bn, Type phi){
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

// simulate a single observation with mean f, same parameterisation as loglik
template<class Type>
Type simobs(int fam, Type f, Type bn, Type phi){
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
