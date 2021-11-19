#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix GaussianKernel(NumericVector x1, NumericVector x2, List kernelPar) {
	double l = as<double>(kernelPar["l"]);
	double sigmaF = as<double>(kernelPar["sigmaF"]);
	int nx1 = x1.size();
	int nx2 = x2.size();
	NumericMatrix out(nx1, nx2);
	
	for (int i=0; i < nx1; i++) {
		for (int j=0; j < nx2; j++) {
			out(i, j) = pow(sigmaF, 2.0)*exp(-1.0/2.0*pow(x1[i]-x2[j],2.0)/pow(l, 2.0));
		}
	}
	return(wrap(out));
}

