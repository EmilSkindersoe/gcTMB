# Properties of the growth models. The position of a model in the list is its index (mod) in the TMB
# templates, which define the mean function (gcmean) and the inverse of F (gcFinv) of each model; F
# and Finv below are the same functions in R. A new model needs an entry here, at the same index as in
# the templates.
#
# All models are of the form f(t) = b1 F((t - b2) / b3), with F increasing from 0 to 1.
# b1-b5: role of each parameter, NA if the model does not use it (it is then always fixed)
#   b1: "asymptote", the response approached at late times
#   b2: "inflection", the time of the maximum growth rate
#   b3: "scale" on the time axis, positive for an increasing curve
#   b4, b5: reserved for parameters of future models
# F(z): increasing from 0 to 1, with its maximum slope at z = 0
# Finv(q): inverse of F, which linearises the model for the starting values: Finv(y / b1) = (t - b2) / b3
gcmodels <- list(
  logistic=list(b1="asymptote", b2="inflection", b3="scale", b4=NA, b5=NA,
                F=function(z) plogis(z),
                Finv=function(q) qlogis(q)),
  gompertz=list(b1="asymptote", b2="inflection", b3="scale", b4=NA, b5=NA,
                F=function(z) exp(-exp(-z)),
                Finv=function(q) -log(-log(q)))
)

# links of the parameters b1-b5 of a model m (see gcmodels): the possible links of each parameter, the
# first one is the default. The scale b3 may have the inverse link, so that its linear predictor is
# the rate 1/b3.
modellinks <- function(m){
  list(b1="identity",
       b2="identity",
       b3=if (m$b3 %in% "scale") c("identity", "inverse") else "identity",
       b4="identity",
       b5="identity")
}

# names of the links, in the order of their index in the TMB templates
linknames <- c("identity", "inverse")

# index of each link, as in the TMB templates: 1 = identity, 2 = inverse
linkindex <- function(link){
  match(link, linknames)
}

# parameters b1-b5 from their linear predictors eta (a vector, or a matrix with a column for each
# parameter), by the inverse links lnk (see linkindex)
linkpars <- function(lnk, eta){
  if (is.matrix(eta)) return(do.call(cbind, lapply(1:5, function(k) linkfun(lnk[k], eta[, k], inverse=TRUE))))
  sapply(1:5, function(k) linkfun(lnk[k], eta[k], inverse=TRUE))
}

# link function (inverse=FALSE) or inverse link of a parameter, by its index (see linkindex); the
# inverse link is its own inverse
linkfun <- function(lnk, x, inverse=FALSE){
  switch(lnk, x, 1/x)
}

# derivative of the inverse link lnk at the linear predictor eta
dlinkinv <- function(lnk, eta){
  switch(lnk, rep(1, length(eta)), -1/eta^2)
}

# properties of a model, by name or index
gcmodel <- function(mod){
  gcmodels[[mod]]
}

# index of the model, as in the TMB template
modelindex <- function(model){
  which(model == names(gcmodels))
}

# parameters not used by the model, which cannot be estimated
unusedpars <- function(mod){
  m <- gcmodel(mod)
  which(is.na(c(m$b1, m$b2, m$b3, m$b4, m$b5)))
}

# growth curve, same as gcmean in the TMB template: f1 F((t - f2) / f3)
gcmeanR <- function(mod, t, f1, f2, f3, f4, f5){
  gcmodel(mod)$F((t - f2) / f3) * f1
}

# time at which a curve with parameters b2 and b3 reaches the fraction q of its asymptote, as gctime
# in the TMB template (residual zero)
gctimeR <- function(mod, q, b2, b3){
  b2 + b3 * gcmodel(mod)$Finv(q)
}
