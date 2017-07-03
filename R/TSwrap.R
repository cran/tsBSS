TSIRc <- function(X, slices = slices, k = 1, h = 10)
    {
    RES <- .Call( "TSIR", X, slices, k, h, PACKAGE = "tsBSS")$RES
    RES
    }

TSAVEc <- function(X, slices = slices, k = 1, h = 10)
{
  RES <- .Call( "TSAVE", X, slices, k, h, PACKAGE = "tsBSS")$RES
  RES
}
