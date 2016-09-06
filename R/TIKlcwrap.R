TIKlc <- function(Y, U = U, k = 1, method = 3)
    {
    RES <- .Call( "TIKlc", Y, U, k, method, PACKAGE = "tsBSS")
    RES$Tik
    }


