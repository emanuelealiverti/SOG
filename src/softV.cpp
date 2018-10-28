#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
/*
Last part is I( |x| - \lambda )_+ ). Trick with sign
*/

// [[Rcpp::export]]
arma::vec softV(arma::vec X, double lambda){

	arma::vec temp = sign(X) % ( abs(X) - lambda) % ( ( sign( abs(X) - lambda ) + 1 ) / 2 );
	arma::vec out = normalise(temp);

  return out;
}
