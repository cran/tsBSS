#Principal Volatility Components
# Modified PVC from Miettinen et al. (2017)

PVC <- function(X, k = 1:12, ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, ...) {
  nk <- length(k)
  MEAN <- colMeans(X)
  COV <- cov(X)
  EVD <- eigen(COV, symmetric = TRUE)
  COV.sqrt.i <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
  X.C <- sweep(X, 2, MEAN, "-")
  Y <- tcrossprod(X.C, COV.sqrt.i)
  p <- ncol(X)
  R <- PVCkc(Y, k)
  U <- eigen(R, symmetric = TRUE)$vectors
  W <- crossprod(U, COV.sqrt.i)
  W <- diag(sign(rowMeans(W))) %*% W
  S <- tcrossprod(X.C, W)
  if (ordered == TRUE) { #Ordering by volatility
    if (is.null(acfk) == TRUE) { acfk <- k }
    ord <- ordf(S, acfk, p, W, alpha, ...)
    W <- ord$W
    if (original == TRUE) {
      S <- ord$S # Original independent components
    } else {
      S <- ord$RS # Residuals based on ARMA fit, if applicable; otherwise original IC's
      Sraw <- ord$S
      Sraw <- ts(Sraw, names = paste("Series", 1:p))
      if (is.ts(X)) attr(Sraw, "tsp") <- attr(X, "tsp")
    }
  }
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, k = k, S = S, MU = MEAN)
  if (ordered == TRUE) {
    if (original == FALSE) {
      RES$Sraw <- Sraw
    }
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

