lblinMc <- function(X, k)
{
  RES <- .Call( "lblinM", X, k, PACKAGE = "tsBSS")$RES
  RES
}

lbsqMc <- function(X, k)
{
  RES <- .Call("lbsqM", X, k, PACKAGE = "tsBSS")$RES
  RES
}

