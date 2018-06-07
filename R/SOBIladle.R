SOBIladle <- function(X, tau = 1:12, l = 20, sim = "geom", n.boot = 200,
                      ncomp = ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X) - 1),
                      maxiter = 1000, eps = 1e-06, ...) {
  data.name <- deparse(substitute(X))
  method <- "SOBI"
  
  if (length(tau) == 1) tau <- 1:tau
  
  n <- nrow(X)
  p <- ncol(X)
  MEAN <- colMeans(X)
  X.c <- sweep(X, 2, MEAN, "-")
  COV <- crossprod(X.c)/(n - 1)
  COV.EVD <- eigen(COV, symmetric = TRUE)
  COV.sqrt.i <- COV.EVD$vectors %*% tcrossprod(diag(COV.EVD$values^(-0.5)),
                                               COV.EVD$vectors)

  Z <- tcrossprod(X.c, COV.sqrt.i)

  M_array <- array(0, dim = c(p, p, length(tau)))
  for(t in 1:length(tau)) {
    M_array[, , t] <- crossprod(Z[1:(n - tau[t]), ], Z[(tau[t] + 1):n, ])/(n - tau[t])
    M_array[, , t] <- (M_array[, , t] + t(M_array[, , t]))/2
  }

  EV.Mdata <- MSobi(X, k_set = tau)
  frjddata <- frjd2(EV.Mdata, eps = eps, maxiter = maxiter)
  Dfrjddata <- diag(apply(frjddata$D^2, 1:2, sum))
  EVdata <- frjddata$V[, order(Dfrjddata, decreasing = TRUE)]

  RES <- tsboot(X, SOBIbootLADLE, R = n.boot, sim = sim, l = l, EVdata = EVdata, tau = tau, rank = ncomp, ...)

  fis <- RES$t
  fn0 <- c(0, colMeans(fis))
  fn <- fn0/(1 + sum(fn0))
  phin <- sort(Dfrjddata, decreasing = TRUE)[1:(ncomp + 1)]/(1 + sum(sort(Dfrjddata, decreasing = TRUE)[1:(ncomp + 1)]))
  gn <- fn + phin
  est.k <- which.min(gn) - 1
  W <- crossprod(EVdata, COV.sqrt.i)
  S <- ts(tcrossprod(X.c, W))
  if(!is.null(attributes(X)$tsp)) {
    attributes(S)$tsp <-  attributes(X)$tsp
  }
  colnames(S) <- paste0("Series", 1:p)
  

  RES <- list(method = method, k = est.k, fn = fn, phin = phin, data.name = data.name,
       gn = gn, lambda = sort(Dfrjddata, decreasing = TRUE)[1:(ncomp + 1)],
       W = W, S = S, MU = MEAN, sim = sim, tau = tau)
  
  class(RES) <- "ladle"
  RES
  
}