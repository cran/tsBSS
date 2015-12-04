#include "tsBSS.h"
#include <Rcpp.h>
//#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP CCK(SEXP Y, SEXP k)
  {
    mat X = as<arma::mat>(Y);
    vec tau = as<vec>(k);

    int p = X.n_cols;
    int n = X.n_rows;
    int nk = tau.n_elem;

    mat Ip = eye(p,p);
    mat CK;
    mat C23;

    cube CCK(p, p, p*p*nk);
    cube CC1(p, p, p*p*nk);
    cube CC23(p, p, p*p*nk);

    int ind = 0;

    for( int kk = 0; kk < nk; kk = kk + 1 )
      {
      int K = tau(kk);
      int nK = n-K;
      mat LAMBDA_K = trans(X.rows(0,nK-1))* X.rows(K,n-1) / nK;

      cube Hkt(p,p,nK);

      for( int hh=0; hh<nK; hh=hh+1)
        {
          Hkt.slice(hh) = trans(X.row(hh+K))*X.row(hh);
        }

      for( int ii=0; ii<p; ii=ii+1)
        {
          for( int jj=0; jj<p; jj=jj+1)
            {
            mat Eij = zeros(p, p);
            Eij(ii,jj) = 1;
            mat Eij2 = Eij + Eij.t();

            mat BK = zeros(p,p);

                for (int bb=0; bb < nK; bb=bb+1)
                    {
                    BK = BK + trans(Hkt.slice(bb)) * Eij * Hkt.slice(bb);
                    }


                BK = BK/nK;

                if(ii==jj)
                {
                C23 = trans(LAMBDA_K) * Eij2 * LAMBDA_K + Ip;
                CK = BK - C23;
                }
                else
                {
                C23 = trans(LAMBDA_K) * Eij2 * LAMBDA_K;
                CK = BK - C23;
                }

                CCK.slice(ind) = CK;
                CC1.slice(ind) = BK;
                CC23.slice(ind) = C23;
                ind = ind+1;
            }
        }
      }

    return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("k") = tau,
                            Rcpp::Named("CCK") = CCK
                            );
  }

SEXP TIK(SEXP Y, SEXP U, SEXP k)
  {
    mat X = as<arma::mat>(Y);
    mat O = as<arma::mat>(U);
    int tau = as<int>(k);

    int p = X.n_cols;
    int n = X.n_rows;

    mat Xt = X.rows(0,n-1-tau);
    mat Xttau = X.rows(tau,n-1);

    mat Tik = zeros(p,p);

    for( int ii = 0; ii < p; ii = ii + 1)
      {
        vec oy = Xt * O.col(ii);
        vec oytau = Xttau * O.col(ii);
        double Tik1 = mean(pow(oy % oytau,2));
        mat Tik2 = mean((Xt % repmat(2 * (oy % oytau) % oytau,1, p)), 0);
        mat Tik3 = mean((Xttau % repmat(2 * (oy % oytau) % oy,1, p)), 0);
        Tik.col(ii) = trans((Tik1-1) * (Tik2 + Tik3));
      }
      return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("k") = tau,
                            Rcpp::Named("Tik") = Tik
                            );
  }
