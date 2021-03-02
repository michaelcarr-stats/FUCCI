
#include <cmath>
#include<Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix InitDomain_y(int &nrow, int &ncol, double &delta) {
    NumericMatrix domain_y(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            domain_y(i, j) = ((double) nrow - i - 1) * sqrt(3.0) / 2 * delta;
        }
    }
    return domain_y;
}

