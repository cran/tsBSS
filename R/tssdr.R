# Method tssdr (time series supervised dimesion reduction)
tssdr <- function (y, X, algorithm = "TSIR", k = 1:12, H = 10, weight = 0.5, method = "frjd",
                   eps = 1e-06, maxiter = 1000, ...) {
  #weight only applies if algorithm = "TSSH"
  if (length(y) != nrow(X)) stop("y and X have different lengths!")
  algorithm <- match.arg(algorithm, c("TSIR", "TSAVE", "TSSH"))
  method <- match.arg(method, c("rjd", "frjd"))
  if((algorithm == "TSSH") & length(H) == 1) {
    H <- rep(H, 2)
    warning("H should be a 2-vector for TSSH. Using the given H for both parts.")
  }
  if((algorithm != "TSSH") & length(H) == 2) {
    stop('H should be a scalar for TSIR and TSAVE!')
  }
  MEAN <- colMeans(X)
  p <- ncol(X)
  COV <- cov(X)
  EVD <- eigen(COV, symmetric = TRUE)
  COV.sqrt.i <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
  X.C <- sweep(X, 2, MEAN, "-")
  XCS <- tcrossprod(X.C, COV.sqrt.i)
  nk <- length(k)
  if (length(H) == 1) { #For TSIR and TSAVE
    slices <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/H))),
                          include.lowest = TRUE, labels = FALSE))
  } else { #For TSSH
    slices <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/H[1]))),
                            include.lowest = TRUE, labels = FALSE))
    slices2 <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/H[2]))),
                            include.lowest = TRUE, labels = FALSE))
  }
  R <- array(0, dim = c(p, p, nk))
  switch(algorithm,
         TSIR  = {
           for (i in 1:length(k)) {
             R[, , i] <- TSIRc(XCS, slices = slices, k = k[i], h = H)
           }
         },
         TSAVE  = {
           for (i in 1:length(k)) {
             R[, , i] <- TSAVEc(XCS, slices = slices, k = k[i], h = H)
           }
         },
         TSSH  = {
           for (i in 1:length(k)) {
             R[, , i] <- (1 - weight)*TSIRc(XCS, slices = slices, k = k[i], h = H[1]) +
               weight*TSAVEc(XCS, slices = slices2, k = k[i], h = H[2])
           }
         }
  )
  # Joint diagonalization
  JD <- switch(method, frjd = {
    frjd(R, eps = eps, maxiter = maxiter, ...)
  }, rjd = {
    rjd(R, eps = eps, maxiter = maxiter, ...)
  })
  sumdg <- diag(apply(JD$D, 1:2, sum))
  ord <- order(sumdg, decreasing = TRUE)
  P <- diag(p)
  P <- P[ord, ]
  D <- JD$D
  DTable <- matrix(0, ncol = p, nrow = nk)
  for (j in 1:length(k)) {
    D[ , , j] <- P %*% tcrossprod(D[ , , j], P) #Diagonal elements are now in order
    DTable[j, ] <- diag(D[ , , j])
  }
  colnames(DTable) <- paste("Dir.", (1:p), sep = "")
  rownames(DTable) <- paste("Lag ", k, sep = "")
  V <- JD$V %*% t(P)
  W <- crossprod(V, COV.sqrt.i) #p*p matrix
  S <- tcrossprod(X.C, W) #Values for the all possible directions
  S <- ts(cbind(y, S), names = c("y", paste("Series", 1:p))) # Response included
  if (is.ts(X)) attr(S, "tsp") <- attr(X, "tsp")
  RES <- list(W = W, k = k, S = S, L = DTable/sum(DTable), H = H,
              yname = deparse(substitute(y)),
              Xname = deparse(substitute(X)),
              algorithm = algorithm)
  class(RES) <- "tssdr"
  RES
}

# Printing method for objects of class "tssdr"
`print.tssdr` <- function(x, digits = 3, ...) {
  print.listof(x[(names(x) == "W") | (names(x) == "k") | (names(x) == "L")], digits = digits, ...)
}

# Extracting the estimated directions:
# Components method for objects of class "tssdr"
`components.tssdr` <- function(x, ...) x$S[, -1]

# Plotting method for objects of class "tssdr" (R's basic time series plot)
`plot.tssdr` <- function(x, main = "The response and the directions", ...) {
  plot.ts(x$S, main = main, ...)
}
