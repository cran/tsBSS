#include <RcppArmadillo.h>

RcppExport SEXP CCK(SEXP Y, SEXP k);
RcppExport SEXP TIK(SEXP Y, SEXP U, SEXP k, SEXP method);
RcppExport SEXP TIKlc(SEXP Y, SEXP U, SEXP k, SEXP method);
RcppExport SEXP TSAVE(SEXP X, SEXP slices, SEXP h, SEXP k);
RcppExport SEXP TSIR(SEXP X, SEXP slices, SEXP h, SEXP k);
