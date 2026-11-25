#include "gcfunctions.h"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// correlation of residuals dt apart in time (CAR(1), 0 < rho < 1), or lag positions apart (AR(1), -1 < rho < 1,
// by repeated products, as rho may be negative)
template<class Type>
Type corlagcor(int cortype, Type rho, Type dt, int lag){
  if (cortype == 2) return exp(dt * log(rho));
  Type cr = 1;
  for (int k = 0; k < lag; k++) cr *= rho;
  return cr;
}

template<class Type>
Type gcTMB(objective_function<Type>* obj) {
  DATA_INTEGER(fam);
  DATA_INTEGER(mod);
  DATA_INTEGER(trunc);     // 1: curve truncated at the threshold b4
  DATA_INTEGER(b1rate);    // 1: b1 is the maximum growth rate instead of the asymptote
  DATA_IVECTOR(lnk);       // links of b1-b5: 1 = identity, 2 = inverse, 3 = log, 4 = logit
  DATA_IVECTOR(tnp);       // random effect terms: number of random coefficients per term
  DATA_IVECTOR(tnl);       // number of levels per term
  DATA_IVECTOR(tpar);      // parameter indices (0-4) of the random coefficients of all terms, concatenated
  DATA_IVECTOR(tzb);       // first column in Z (and Zp) of the design block of each random coefficient
  DATA_INTEGER(simre);     // simulate: 1 = new random effects, 0 = keep the fitted ones
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_VECTOR(w);
  DATA_VECTOR(bn);
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);
  DATA_MATRIX(X3);
  DATA_MATRIX(X4);
  DATA_MATRIX(X5);
  DATA_SPARSE_MATRIX(Z);   // design blocks of all terms, side by side: a column for each level, with the
                           // covariate of a random slope (1 for intercepts); zero columns without random effects
  DATA_VECTOR(xp);         // predictions (predict with se.fit): times, length zero otherwise
  DATA_MATRIX(Xp1);
  DATA_MATRIX(Xp2);
  DATA_MATRIX(Xp3);
  DATA_MATRIX(Xp4);
  DATA_MATRIX(Xp5);
  DATA_SPARSE_MATRIX(Zp);  // random effects of the fitted levels; zero rows for population predictions
  DATA_INTEGER(cortype);   // autocorrelation of the residuals within groups: 0 = none, 1 = AR(1), 2 = CAR(1)
  DATA_IVECTOR(corprev);   // previous observation (zero-based) of the same group in time, -1 for the first; length n
  DATA_VECTOR(cordt);      // CAR(1): time since the previous observation of the group; length n
  DATA_IVECTOR(corlag);    // AR(1): lag (in positions) since the previous observation of the group; length n
  DATA_IVECTOR(corord);    // observations in an order with each after its previous one, for simulation;
                           // spatial structures: the observations of each group consecutive
  DATA_IVECTOR(corgsize);  // spatial structures (3 = exponential, 4 = Gaussian): number of observations per
                           // group, in the order of corord; length zero otherwise
  DATA_VECTOR(cordist);    // spatial structures: distance matrix of each group, by column, concatenated
  PARAMETER_VECTOR(b1);
  PARAMETER_VECTOR(b2);
  PARAMETER_VECTOR(b3);
  PARAMETER_VECTOR(b4);
  PARAMETER_VECTOR(b5);
  PARAMETER_VECTOR(u);     // length zero without random effects
  PARAMETER_VECTOR(theta); // length zero without random effects
  PARAMETER(log_sigma);    // log sd (gaussian; of log(y) for lognormal), precision (beta, betabinomial), size (nbinom2) or dispersion (nbinom1);
                           // mapped (fixed) for poisson and binomial
  PARAMETER_VECTOR(cor_par); // autocorrelation on the unconstrained scale of nlme: log((1 + rho)/(1 - rho)) for AR(1),
                           // log(rho/(1 - rho)) for CAR(1), log(range) and logit(nugget) for the spatial
                           // structures (without nugget: length 1); mapped (fixed) without autocorrelation

  Type phi = exp(log_sigma);
  // correlated residuals (gaussian, and log(y) for lognormal): stationary, with the marginal standard
  // deviation sigma/sqrt(w) and the correlation rho^dt between observations dt apart in time within a group
  // (AR(1), CAR(1)), or (1 - nugget) exp(-d/range) (exponential) or (1 - nugget) exp(-(d/range)^2)
  // (Gaussian) between observations a distance d apart
  bool correlated = (cortype > 0) && (fam == 1 || fam == 8);
  bool spatial = (cortype >= 3);
  Type rho = 0, range = 1, nugget = 0;
  if (cortype == 1) rho = (exp(cor_par(0)) - Type(1)) / (exp(cor_par(0)) + Type(1));
  if (cortype == 2) rho = Type(1) / (Type(1) + exp(-cor_par(0)));
  if (spatial){
    range = exp(cor_par(0));
    if (cor_par.size() > 1) nugget = Type(1) / (Type(1) + exp(-cor_par(1)));
  }
  if (correlated && !spatial) ADREPORT(rho);
  if (fam == 1 || fam == 8){
    Type sigma = phi;
    ADREPORT(sigma);
  }
  if (fam == 4){
    ADREPORT(phi);
  }

  Type nll = 0;
  int n = y.size();

  // random effects: column k of Zu is added to the linear predictor of b(k+1)
  // term t has tnl(t) levels and tnp(t) random coefficients (intercepts or slopes) of the parameters
  // tpar, with design blocks at columns tzb of Z, correlated within a term and independent between terms
  matrix<Type> Zu(n, 5);
  Zu.setZero();
  // simulated random effects (the fitted ones for simre = 0), in the order of u, and means
  vector<Type> usim(u.size());
  vector<Type> fsim(n);
  int nt = tnp.size();
  if (nt > 0){
    int ncorr = 0;
    for (int t = 0; t < nt; t++) ncorr += tnp(t)*tnp(t);
    vector<Type> corr(ncorr);

    int uoff = 0, toff = 0, poff = 0, coff = 0;
    for (int t = 0; t < nt; t++){
      int np = tnp(t);
      int nl = tnl(t);
      vector<Type> ut = u.segment(uoff, nl*np);
      vector<int> dim(2);
      dim << nl, np;
      array<Type> U(ut, dim);

      vector<Type> logsd = theta.segment(toff, np);
      vector<Type> corr_transf = theta.segment(toff + np, np*(np-1)/2);
      vector<Type> sd = exp(logsd);
      density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
      density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);

      for (int j = 0; j < nl; j++){
        nll += scnldens(U.transpose().col(j));
      }

      // new random effects for each level
      SIMULATE {
        if (simre == 1){
          vector<Type> uj(np);
          for (int j = 0; j < nl; j++){
            scnldens.simulate(uj);
            for (int k = 0; k < np; k++) U(j, k) = uj(k);
          }
        }
        for (int k = 0; k < np; k++){
          for (int j = 0; j < nl; j++) usim(uoff + j + k*nl) = U(j, k);
        }
      }

      matrix<Type> ct = nldens.cov();
      for (int l = 0; l < np; l++){
        for (int k = 0; k < np; k++) corr(coff + k + l*np) = ct(k, l);
      }

      for (int k = 0; k < np; k++){
        vector<Type> Zuk = Z.middleCols(tzb(poff + k), nl) * U.col(k).matrix();
        Zu.col(tpar(poff + k)) += Zuk.matrix();
      }

      uoff += nl*np;
      toff += np + np*(np-1)/2;
      poff += np;
      coff += np*np;
    }
    ADREPORT(corr);
  }

  vector<Type> Xb1 = X1 * b1;
  vector<Type> Xb2 = X2 * b2;
  vector<Type> Xb3 = X3 * b3;
  vector<Type> Xb4 = X4 * b4;
  vector<Type> Xb5 = X5 * b5;

  Type f;
  vector<Type> fv(n);
  for (int i = 0; i < n; i++){
    Type f2 = linkinv(lnk(1), Type(Xb2(i) + Zu(i, 1)));
    Type f3 = linkinv(lnk(2), Type(Xb3(i) + Zu(i, 2)));
    Type f4 = linkinv(lnk(3), Type(Xb4(i) + Zu(i, 3)));
    Type f5 = linkinv(lnk(4), Type(Xb5(i) + Zu(i, 4)));
    Type f1 = gcasym(b1rate, mod, trunc, linkinv(lnk(0), Type(Xb1(i) + Zu(i, 0))), f2, f3, f4, f5);
    f = gcmean(mod, trunc, x(i), f1, f2, f3, f4, f5);
    fv(i) = f;
    if (!correlated){
      nll += -w(i) * loglik(fam, y(i), f, bn(i), phi);
      SIMULATE {
        y(i) = simobs(fam, f, bn(i), phi);
      }
    }
    SIMULATE {
      fsim(i) = f;
    }
  }

  if (correlated && spatial){
    // standardised residuals r as below, multivariate normal within each group with the correlation
    // matrix of the distances between its observations; observations with weight zero are left out
    vector<Type> r(n);
    for (int i = 0; i < n; i++){
      if (w(i) <= 0){
        r(i) = 0;
      } else {
        Type e = (fam == 8) ? Type(log(y(i)) - log(fv(i))) : Type(y(i) - fv(i));
        r(i) = e * sqrt(w(i)) / phi;
        // the Jacobian of the standardisation, and of the log for lognormal
        nll -= Type(0.5) * log(w(i)) - log(phi);
        if (fam == 8) nll += log(y(i));
      }
    }
    vector<Type> rs(n);
    SIMULATE {
      for (int i = 0; i < n; i++) rs(i) = rnorm(Type(0), Type(1));
    }
    int off = 0, doff = 0;
    for (int g = 0; g < corgsize.size(); g++){
      int m = corgsize(g);
      matrix<Type> R(m, m);
      vector<Type> rg(m);
      for (int a = 0; a < m; a++){
        rg(a) = r(corord(off + a));
        for (int b = 0; b < m; b++){
          if (a == b){
            R(a, b) = 1;
          } else {
            Type d = cordist(doff + a + b*m) / range;
            Type k = (cortype == 3) ? Type(exp(-d)) : Type(exp(-d*d));
            R(a, b) = (Type(1) - nugget) * k;
          }
        }
      }
      density::MVNORM_t<Type> mvn(R);
      nll += mvn(rg);
      SIMULATE {
        vector<Type> sg = mvn.simulate();
        for (int a = 0; a < m; a++) rs(corord(off + a)) = sg(a);
      }
      off += m;
      doff += m*m;
    }
    SIMULATE {
      for (int i = 0; i < n; i++){
        Type wi = (w(i) > 0) ? w(i) : Type(1);
        Type e = rs(i) * phi / sqrt(wi);
        y(i) = (fam == 8) ? Type(fv(i) * exp(e)) : Type(fv(i) + e);
      }
    }
  }

  if (correlated && !spatial){
    // standardised residuals r = (g(y) - g(f)) sqrt(w)/sigma, g the identity (gaussian) or log
    // (lognormal), with precision weights w; within a group in time, r given the previous residual is
    // normal with mean cr r_prev and variance 1 - cr^2, cr = rho^dt (the Markov property of AR(1)).
    // Observations with weight zero are not in the chains and do not contribute.
    // all residuals first: the previous observation of a group may come later in the data
    vector<Type> r(n);
    for (int i = 0; i < n; i++){
      if (w(i) <= 0){
        r(i) = 0;
      } else {
        Type e = (fam == 8) ? Type(log(y(i)) - log(fv(i))) : Type(y(i) - fv(i));
        r(i) = e * sqrt(w(i)) / phi;
      }
    }
    for (int i = 0; i < n; i++){
      if (w(i) <= 0) continue;
      int p = corprev(i);
      Type ld;
      if (p < 0){
        ld = dnorm(r(i), Type(0), Type(1), true);
      } else {
        Type cr = corlagcor(cortype, rho, cordt(i), corlag(i));
        ld = dnorm(r(i), cr * r(p), sqrt(Type(1) - cr*cr), true);
      }
      // density of y: the Jacobian of the standardisation, and of the log for lognormal
      ld += Type(0.5) * log(w(i)) - log(phi);
      if (fam == 8) ld -= log(y(i));
      nll -= ld;
    }
    SIMULATE {
      vector<Type> rs(n);
      rs.setZero();
      for (int k = 0; k < n; k++){
        int i = corord(k);
        int p = corprev(i);
        if (p < 0 || w(i) <= 0){
          rs(i) = rnorm(Type(0), Type(1));
        } else {
          Type cr = corlagcor(cortype, rho, cordt(i), corlag(i));
          rs(i) = cr * rs(p) + sqrt(Type(1) - cr*cr) * rnorm(Type(0), Type(1));
        }
        Type wi = (w(i) > 0) ? w(i) : Type(1);
        Type e = rs(i) * phi / sqrt(wi);
        y(i) = (fam == 8) ? Type(fv(i) * exp(e)) : Type(fv(i) + e);
      }
    }
  }
  SIMULATE {
    REPORT(y);
    REPORT(usim);
    REPORT(fsim);
  }

  // mean at the times xp, for standard errors of the predictions by the delta method
  int npred = xp.size();
  if (npred > 0){
    matrix<Type> Zup(npred, 5);
    Zup.setZero();
    int uoff = 0, poff = 0;
    for (int t = 0; t < nt; t++){
      int np = tnp(t);
      int nl = tnl(t);
      for (int k = 0; k < np; k++){
        vector<Type> uk = u.segment(uoff + k*nl, nl);
        vector<Type> Zuk = Zp.middleCols(tzb(poff + k), nl) * uk.matrix();
        Zup.col(tpar(poff + k)) += Zuk.matrix();
      }
      uoff += nl*np;
      poff += np;
    }
    vector<Type> Xpb1 = Xp1 * b1;
    vector<Type> Xpb2 = Xp2 * b2;
    vector<Type> Xpb3 = Xp3 * b3;
    vector<Type> Xpb4 = Xp4 * b4;
    vector<Type> Xpb5 = Xp5 * b5;
    vector<Type> mu(npred);
    for (int i = 0; i < npred; i++){
      Type f2 = linkinv(lnk(1), Type(Xpb2(i) + Zup(i, 1)));
      Type f3 = linkinv(lnk(2), Type(Xpb3(i) + Zup(i, 2)));
      Type f4 = linkinv(lnk(3), Type(Xpb4(i) + Zup(i, 3)));
      Type f5 = linkinv(lnk(4), Type(Xpb5(i) + Zup(i, 4)));
      Type f1 = gcasym(b1rate, mod, trunc, linkinv(lnk(0), Type(Xpb1(i) + Zup(i, 0))), f2, f3, f4, f5);
      mu(i) = gcmean(mod, trunc, xp(i), f1, f2, f3, f4, f5);
    }
    ADREPORT(mu);
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
