#include<Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix InitDomain_x(int &nrow, int &ncol, double &delta) {
    NumericMatrix domain_x(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            if (((nrow - i) % 2) == 0) {
                domain_x(i, j) = j * delta;
            }
            else {
                domain_x(i, j) = (j + 0.5) * delta;
            }
        }
    }
    return domain_x;
}