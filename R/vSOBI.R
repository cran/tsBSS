# Method vSOBI
vSOBI <- function(X, ...) UseMethod("vSOBI")

# main function for vSOBI
vSOBI.default <- function (X, k = 1:12, eps = 1e-06, maxiter = 1000, G = "pow", 
                           ordered = FALSE, acfk = NULL, original = TRUE, ...){
  G <- match.arg(G, c("pow", "lcosh"))
  MEAN <- colMeans(X)
  COV <- cov(X)
  EVD <- eigen(COV, symmetric = TRUE)
  COV.sqrt.i <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
  X.C <- sweep(X, 2, MEAN, "-")
  Y <- tcrossprod(X.C, COV.sqrt.i)
  p <- ncol(X)
  U <- diag(p) #Initial value for the orthogonal matrix U
  crit <- Inf
  iter <- 0
  if(length(which(k < 1) != 0)) stop("only non-zero lags allowed")
  nk <- length(k)
  Tk <- array(NA, dim = c(p, p, nk))
  if (G == "pow") {TIKuse <- TIKc} else {TIKuse <- TIKlcc}
  while(crit > eps) {
      for (i in 1:nk){
        Tk[ , , i] <- TIKuse(Y, U, k = k[i], method = 3)
      }
    TU <- apply(Tk, c(1, 2), sum)
    EVDt <- eigen(tcrossprod(TU), symmetric = TRUE)
    COVt.sqrt.i <- EVDt$vectors %*% tcrossprod(diag(EVDt$values^(-0.5)), EVDt$vectors)
    U.new <- COVt.sqrt.i %*% TU #Updated U
    crit <- sqrt(sum((abs(U.new) - abs(U))^2)) #Comparing the current and the new matrix U
    iter <- iter + 1
    if(iter > maxiter) stop("maxiter reached without convergence")
    U <- U.new
  } #While the criterion value is below tolerance value.
  W <- crossprod(U, COV.sqrt.i) #Unmixing matrix
  S <- tcrossprod(X.C, W)
  if (ordered == TRUE) { #Ordering by volatility
    if (is.null(acfk) == TRUE) { acfk <- k }
    ord <- ordf(S, acfk, p, ...)
    if (original == TRUE) {
      S <- ord$S # Original independent components
    } else {
      S <- ord$RS # Residuals based on ARMA fit, if applicable; otherwise otiginal IC's
    }
  }
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, k = k, S = S)
  if (ordered == TRUE) {
    RES$fits <- ord$fits
    RES$armaeff <- ord$armaeff
    RES$linTS <- ord$linTS
    RES$linP <- ord$linP
    RES$volTS <- ord$volTS
    RES$volP <- ord$volP
  }
  class(RES) <- c("bssvol", "bss")
  RES
}


vSOBI.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- vSOBI.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}
