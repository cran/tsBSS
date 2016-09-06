#include <RcppArmadillo.h>

RcppExport SEXP CCK(SEXP Y, SEXP k);
RcppExport SEXP TIK(SEXP Y, SEXP U, SEXP k, SEXP method);
RcppExport SEXP TIKlc(SEXP Y, SEXP U, SEXP k, SEXP method);
