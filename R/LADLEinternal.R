# Functions needed for AMUSE and SOBI ladles

# function for frjd without maxiter check
frjd2 <- function (X, maxiter = 100, eps = 1e-06) {
    dim.X <- dim(X)
        p <- dim.X[1]
        K <- dim.X[3]
        Xt <- aperm(X, c(1, 3, 2))
        X <- matrix(Xt, ncol = p)
    dim.X <- dim(X)

    p <- dim.X[2]
    K <- dim.X[1]/p
    weight <- rep(1, K)
    res <- .C("rjd2c", as.double(as.vector(X)), as.integer(c(K,
        p, maxiter)), as.double(as.vector(weight)), as.double(eps),
        res = double(p^2 + 1), PACKAGE = "tsBSS")$res
    iter <- res[p^2 + 1]
    V <- matrix(res[1:p^2], p, p)
    D <- X
    for (k in 1:K) {
        D[((k - 1) * p + 1):(k * p), ] <- crossprod(V, crossprod(X[((k -
            1) * p + 1):(k * p), ], V))
    }
    D <- array(t(D), c(p, p, K))
    list(V = V, D = D, iter = iter)
}

# function to measure the bootstrap variation of the eigenvectors
fi <- function (EVboot, EVdata, rank) {
    fni <- numeric(rank)
    for (ii in 1:rank) {
        fni[ii] <- det(crossprod(EVdata[, 1:ii], EVboot[, 1:ii]))
    }
    1 - abs(fni)
}

# Mfunction for AMUSE
MAmuse <- function (x, k) {
    n <- nrow(x)
    x.c <- sweep(x, 2, colMeans(x), "-")
    COV <- crossprod(x.c)/(n - 1)
    COV.EVD <- eigen(COV, symmetric = TRUE)
    COV.sqrt.i <- COV.EVD$vectors %*% tcrossprod(diag(COV.EVD$values^(-0.5)),
        COV.EVD$vectors)
    Z <- tcrossprod(x.c, COV.sqrt.i)
    M <- crossprod(Z[1:(n - k), ], Z[(k + 1):n, ])/(n -
        k)
    M.sym <- (M + t(M))/2
    M.sym
    crossprod(M.sym)
}

# Bootstrapping function for AMUSE
AMUSEbootLADLE <- function (X, EVdata, tau, rank) {
    Mboot <- MAmuse(X,k = tau)
    EVboot <- eigen(Mboot, symmetric = TRUE)$vectors
    fi(EVboot, EVdata, rank)
}

# Mfunction for SOBI
MSobi <- function (x, k_set) {
  n <- nrow(x)
  p <- ncol(x)
  x.c <- sweep(x, 2, colMeans(x), "-")
  COV <- crossprod(x.c)/(n - 1)
  COV.EVD <- eigen(COV, symmetric = TRUE)
  COV.sqrt.i <- COV.EVD$vectors %*% tcrossprod(diag(COV.EVD$values^(-0.5)),
                                               COV.EVD$vectors)
  Z <- tcrossprod(x.c, COV.sqrt.i)

  M_array <- array(0, dim = c(p, p, length(k_set)))
  for(t in 1:length(k_set)){
    M_array[, , t] <- crossprod(Z[1:(n - k_set[t]), ], Z[(k_set[t] + 1):n, ])/(n - k_set[t])
    M_array[, , t] <- (M_array[, , t] + t(M_array[, , t]))/2
  }

  M_array
}

# Bootstrapping function for SOBI
SOBIbootLADLE <- function (X, EVdata, tau, rank) {
  Mboot <- MSobi(X, k_set = tau)
  frjdboot <- frjd2(Mboot)
  Dfrjd <- diag(apply(frjdboot$D^2, 1:2, sum))
  EVboot <- frjdboot$V[, order(Dfrjd, decreasing = TRUE)]
  fi(EVboot, EVdata, rank)
}
