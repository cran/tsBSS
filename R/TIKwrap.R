TIK <- function(Y, U = U, k = 1, method = 3)
    {
    RES <- .Call( "TIK", Y, U, k, method, PACKAGE = "tsBSS")
    RES$Tik
    }


