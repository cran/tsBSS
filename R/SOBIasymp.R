SOBIasymp <- function (X, k, tau = 1:12, eps = 1e-06, maxiter = 200 ) {
  DNAME <- deparse(substitute(X))
  # method <- match.arg(method, c("satterthwaite", "integration", "saddlepoint"))
  p <- ncol(X)
  n <- nrow(X)
  
  if (length(tau) == 1) tau <- 1:tau
  ntau <- length(tau)
  
  MU <- colMeans(X)
  X.C <- sweep(X, 2, MU, "-")
  COV <- crossprod(X.C)/n
  EVD.COV <- eigen(COV, symmetric = TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% tcrossprod(diag((1/EVD.COV$values)^0.5), EVD.COV$vectors)
  X.C.W <- tcrossprod(X.C, COV.inv.sqrt)
  
  R <- array(0, dim = c(p, p, ntau))
  n <- nrow(X)
  for (i in 1:ntau) {
    Yt <- X.C.W[1:(n - tau[i]), ]
    Yti <- X.C.W[(1 + tau[i]):n, ]
    Ri <- crossprod(Yt, Yti)/nrow(Yt)
    R[, , i] <- (Ri + t(Ri))/2
  }
  
  JD <- frjd(R, eps = eps, maxiter = maxiter)
  W <- crossprod(JD$V, COV.inv.sqrt)
  sumEVs <- apply(JD$D^2,1:2,sum)
  
  ORDER <- order(diag(sumEVs), decreasing = TRUE)
  D <- diag(sumEVs)[ORDER]
  W <- W[ORDER, ]
  Z <- tcrossprod(X.C, W)
  
  D2 <- sumEVs[ORDER, ORDER]
  Tk <- n * sum((D2[(k + 1):p, (k + 1):p]))
  
  Z <- ts(Z)
  if(!is.null(attributes(X)$tsp)) {
    attributes(Z)$tsp <-  attributes(X)$tsp
  }
  colnames(Z) <- paste0("Series", 1:p)
  
  names(Tk) <- "T"
  
  PARAMETER <- ntau/2 * (p - k) *((p - k) + 1)
  names(PARAMETER) <- c("df")
  PVAL <- 1 - pchisq(Tk, PARAMETER)
  METHOD <- c("SOBI test for white noise processes")
  ALTERNATIVE <- paste0("there are fewer than ", p - k, 
                        " white noise components")
  RES <-  list(statistic = Tk, p.value = PVAL, parameter = PARAMETER, 
               method = METHOD, data.name = DNAME, alternative = ALTERNATIVE, 
               k = k, W = W, S = Z, D = D, MU = MU)
  class(RES) <- c("ictest", "htest")
  RES
}
