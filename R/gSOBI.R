# gSOBI algorithm (combination of SOBI and vSOBI)
gSOBI <- function (X, k1 = 1:12, k2 = 1:3, b = 0.9, eps = 1e-06, maxiter = 1000, 
                   ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, ...){
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
  K1 <- length(k1)
  K2 <- length(k2)
  Tk1 <- array(NA, dim = c(p, p, K1))
  Tk2 <- array(NA, dim = c(p, p, K2))
  while(crit > eps) {
      for (i in 1:K1) {
        Tk1[ , , i] <- TIK1c(Y, U, k = k1[i])
      }
      for (i in 1:K2) {
        Tk2[ , , i] <- TIKc(Y, U, k = k2[i], method = 3)
      }
    TU <- b*apply(Tk1, c(1, 2), sum) + (1 - b)*apply(Tk2, c(1, 2), sum)
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
    if (is.null(acfk) == TRUE) { acfk <- 1:max(k1, k2) }
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
  if (is.ts(X)) attr(S, "tsp") <- attr(X, "tsp")
  RES <- list(W = W, k1 = k1, k2 = k2, S = S, MU = MEAN)
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

#Function to order the components (calculates the volatilities of the components or their residuals)
# Used by the functions gFOBI, gJADE, gSOBI, vSOBI, PVC and FixNA
ordf <- function(S, acfk, p, W, alpha, ...) {
  lblin1 <- lbtest(S, acfk, "linear") #Linear autocorrelations exist?
  lineff <- lblin1$p_val
  armaeff <- rep(0, p)
  fits <- vector("list", p)
  S2 <- S
  coef <- numeric(p)
  for (i in 1:p) {
    if(lineff[i] < alpha) {
      fits[[i]] <- forecast::auto.arima(S[, i], stationary = T, seasonal = F, ...)
      S2[, i] <- fits[[i]]$residuals #Replaces series with residuals (if autocorrelation exists)
      armaeff[i] <- 1 #Logical vector: 1 if series is replaced with residuals; 0 if not
    }
  }
  vol <- lbtest(S2, acfk, "squared")
  ord <- vol$TS
  fits2 <- vector("list", p)
  j <- 1
  ordered <- order(ord, decreasing = T)
  for(i in ordered) {
    fits2[[j]] <- fits[[i]]
    j = j + 1
  }
  list(S = S[, ordered],
       RS = S2[, ordered], #Residuals if ARMA effects; otherwise original independent component
       W = W[ordered, ],
       fits = fits2,
       volTS = vol$TS[ordered], volP = vol$p_val[ordered],
       armaeff = armaeff[ordered], linTS = lblin1$TS[ordered], linP = lineff[ordered])
}


`print.bssvol` <- function(x, ...) {
  print.listof(x[(names(x) != "S") & (names(x) != "Sraw") & (names(x) != "MU") & (names(x) != "fits") & (names(x) != "armaeff") & (names(x) != "linTS")
                 & (names(x) != "linP") & (names(x) != "volTS") & (names(x) != "volP")], ...)
}
