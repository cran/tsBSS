PVCkc <- function(Y, k)
{
  RES <- .Call( "PVCk", Y, k, PACKAGE = "tsBSS")$R
  RES
}
