# Method gFOBI
gFOBI <- function(X, ...) UseMethod("gFOBI")

# main function for gFOBI
gFOBI.default <- function (X, k = 0:12, eps = 1e-06, maxiter = 100, method = "frjd",
                           na.action = na.fail, weight = NULL,
                           ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, ...) 
{
  nk <- length(k)
  method <- match.arg(method, c("rjd", "frjd"))
  MEAN <- colMeans(X)
  COV <- cov(X)
  EVD <- eigen(COV, symmetric = TRUE)
  COV.sqrt.i <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
  X.C <- sweep(X, 2, MEAN, "-") #Centered
  Y <- tcrossprod(X.C, COV.sqrt.i) #Scaled
  p <- ncol(X)
  R <- array(0, dim = c(p, p, nk))
  n <- nrow(X)
  for (i in 1:nk) {
    Yt <- Y[1:(n - k[i]), ]
    Yti <- Y[(1 + k[i]):n, ]
    r <- sqrt(rowSums(Yt^2))
    Yu <- r*Yti
    Ri <- crossprod(Yu)/nrow(Yt)
    R[, , i] <- Ri
  }
  JD <- switch(method, frjd = {
    frjd(R, eps = eps, maxiter = maxiter, na.action = na.action, weight = weight)$V
  }, rjd = {
    rjd(R, eps = eps, maxiter = maxiter, na.action = na.action)$V
  })
  W <- crossprod(JD, COV.sqrt.i)
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

gFOBI.ts <- function(X, ...)
{
  x <- as.matrix(X)
  RES <- gFOBI.default(x,...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

