#include<Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
IntegerMatrix InitDomain(NumericMatrix &domain_x, int &nrow, int &ncol, double &Xmax, double &InitialL, double &InitialDensity, NumericMatrix &percentIS) {
    
    //initialisation
    IntegerMatrix domain(nrow, ncol); 
    
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            
            //Initial Seeding LHS & RHS
                //check if node is within the initial seed domain
                if ((domain_x(i, j) >= 0 && domain_x(i, j) <= InitialL/2) || (domain_x(i, j) >= Xmax - InitialL/2 && domain_x(i, j) <= Xmax)) {
                    if (R::runif(0,1) < InitialDensity) {
                        double r = R::runif(0,1);
                        if (r <= percentIS(0,0)) {
                            domain(i, j) = 1;
                        }
                        else if (r > percentIS(0,0) && r <= (percentIS(0,0) + percentIS(1,0))) {
                            domain(i, j) = 2;
                        }
                        else {
                            domain(i, j) = 3;
                        }
                    } else {
                        domain(i, j) = 0;
                    } 
                } else {
                    domain(i, j) = 0;
                }
        }
    }

    return domain;
}
