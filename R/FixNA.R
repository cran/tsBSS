# Method FixNA (and FixNA2)
FixNA <- function(X,...) UseMethod("FixNA")

# main function for FixNA and FixNA2)
FixNA.default <- function (X, k = 1:12, eps = 1e-06, maxiter = 1000, G = "pow", method = "FixNA", ...){
  method <- match.arg(method, c("FixNA", "FixNA2"))
  G <- match.arg(G, c("pow", "lcosh"))
  if (method == "FixNA") met = 1 else met = 2 #FixNA = 1, FixNA2 = 2
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
  if (G == "pow") {TIKuse <- TIK} else {TIKuse <- TIKlc}
  while(crit > eps) {
    for (i in 1:nk) {
      Tk[ , , i] <- TIKuse(Y, U, k = k[i], method = met)
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
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, k = k, S = S)
  class(RES) <- "bss"
  RES
}


FixNA.ts <- function(X, ...)
{
  x <- as.matrix(X)
  RES <- FixNA.default(x,...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}