AMUSEladle <- function(X, tau = 1, l = 20, sim = "geom", n.boot = 200,
                       ncomp = ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X) - 1), ...) {
    data.name <- deparse(substitute(X))
    method <- "AMUSE"

    p <- ncol(X)
    n <- nrow(X)
    MEAN <- colMeans(X)
    X.c <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(X.c)/(n - 1)
    COV.EVD <- eigen(COV, symmetric = TRUE)
    COV.sqrt.i <- COV.EVD$vectors %*% tcrossprod(diag(COV.EVD$values^(-0.5)),
        COV.EVD$vectors)

    Z <- tcrossprod(X.c, COV.sqrt.i)
    Mdata <- crossprod(Z[1:(n - tau), ], Z[(tau + 1):n, ])/(n - tau)
    Mdata.sym <- (Mdata + t(Mdata))/2

    Mdata.sym <- crossprod(Mdata.sym)

    EV.Mdata <- eigen(Mdata.sym, symmetric = TRUE)
    EVdata <- EV.Mdata$vectors


    RES <- tsboot(X, AMUSEbootLADLE, R = n.boot, sim = sim, l = l, EVdata = EVdata,
                  tau = tau, rank = ncomp, ...)
    fis <- RES$t
    fn0 <- c(0, colMeans(fis))
    fn <- fn0/(1 + sum(fn0))
    phin <- EV.Mdata$values[1:(ncomp + 1)]/(1 + sum(EV.Mdata$values[1:(ncomp + 1)]))
    gn <- fn + phin
    est.k <- which.min(gn) - 1

    W <- crossprod(EVdata, COV.sqrt.i)
    S <- ts(tcrossprod(X.c, W))
    if(!is.null(attributes(X)$tsp)) {
      attributes(S)$tsp <-  attributes(X)$tsp
    }
    colnames(S) <- paste0("Series", 1:p)
    
    
    RES <- list(method = method, k = est.k, fn = fn, phin = phin, data.name = data.name,
                gn = gn, lambda = sort(EV.Mdata$values, decreasing = TRUE)[1:(ncomp + 1)],
                W = W, S = S, MU = MEAN, sim = sim, lag = lag)
    
    class(RES) <- "ladle"
    RES
    }