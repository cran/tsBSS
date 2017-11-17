# Ljung-Box-type test
# Option "linear": Checking if linear autocorrelations are present
#  - whether ARMA-components are present or not
# Option "squared": Checking if second order autocorrelations are present
# - whether volatility is constant or not

lbtest <- function(X, k, type = "squared") {
  type <- match.arg(type, c("linear", "squared"))
  lbc <- switch(type, 
                linear = lblinMc,
                squared = lbsqMc)
  TS <- as.vector(lbc(X, k))
  K <- length(k)
  p_val <- 1 - pchisq(TS, K)
  RES <- list(TS = TS, p_val = p_val, Xname = deparse(substitute(X)),
              varnames = colnames(X), k = k, K = K, type = type)
  class(RES) <- "lbtest"
  RES
}


`print.lbtest` <- function(x, digits = 3, ...) {
  cat("\n")
  cat(strwrap(paste("Serial autocorrelation test for", x$Xname), prefix = "\t"), sep = "\n")
  cat("\n")
  cat(paste("Testing for", x$type), "autocorrelations based on lags ")
  cat(x$k, "\n")
  cat(paste("Based on a chi squared test with", x$K), "degrees of freedom \n\n")
  cat("The test statistic and the corresponding p-value for each series: \n")
  restab <- as.data.frame(cbind(x$varnames, round(x$TS, digits), format.pval(x$p_val, digits = digits)))
  colnames(restab) <- c("Series", "Statistic", "p-value")
  print(restab, row.names = FALSE, digits = 4)
  cat("\n")
  cat("For each series the alternative hypothesis is:", paste("serial", x$type), "correlation exists \n")
}
