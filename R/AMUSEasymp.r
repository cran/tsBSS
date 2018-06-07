AMUSEasymp <- function (X, k, tau = 1) {
  DNAME <- deparse(substitute(X))
  # method <- match.arg(method, c("satterthwaite", "integration", "saddlepoint"))
  p <- ncol(X)
  n <- nrow(X)
  
  MU <- colMeans(X)
  X.C <- sweep(X, 2, MU, "-")
  COV <- crossprod(X.C)/n
  EVD.COV <- eigen(COV, symmetric = TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% tcrossprod(diag((1/EVD.COV$values)^0.5), EVD.COV$vectors)
  X.C.W <- tcrossprod(X.C, COV.inv.sqrt)
  
  Yt <- X.C.W[1:(n - tau), ]
  Yti <- X.C.W[(1 + tau):n, ]
  R <- crossprod(Yt, Yti)/nrow(Yt)
  R <- (R + t(R))/2
  
  
  EVD <- eigen(R, symmetric = TRUE)
  W <- crossprod(EVD$vectors, COV.inv.sqrt)
  
  ORDER <- order(EVD$values^2, decreasing = TRUE)
  D <- EVD$values[ORDER]
  W <- W[ORDER, ]
  Z <- tcrossprod(X.C, W)
  
  D2 <- D^2
  Tk <- n * sum((D2[(k + 1):p]))
  
  Z <- ts(Z)
  if(!is.null(attributes(X)$tsp)) {
    attributes(Z)$tsp <-  attributes(X)$tsp
  }
  colnames(Z) <- paste0("Series", 1:p)
  
  names(Tk) <- "T"
  PARAMETER <- 0.5 * (p - k) *((p - k) + 1)
  names(PARAMETER) <- c("df")
  PVAL <- 1 - pchisq(Tk, PARAMETER)
  METHOD <- c("SOBI test for white noise processes")
  ALTERNATIVE <- paste0("there are fewer than ", p - k, 
                        " white noise components")
  RES <-  list(statistic = Tk, p.value = PVAL, parameter = PARAMETER, 
               method = METHOD, data.name = DNAME, alternative = ALTERNATIVE, 
               k = k, W = W, S = Z, D = D, MU = MU, tau = tau)
  class(RES) <- c("ictest", "htest")
  RES
}
