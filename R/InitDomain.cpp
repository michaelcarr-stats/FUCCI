#include<Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
IntegerMatrix InitDomain(NumericMatrix &domain_x, int &nrow, int &ncol, double &Xmax, double &InitialL, double &InitialDensity, NumericMatrix &percentIS,  NumericMatrix &countIS, bool &SetCells) {
    
    //initialisation
    IntegerMatrix domain(nrow, ncol); 
    
    if (SetCells == false) { // use percentages to fill domain
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
    } else { //use initial cell counts
        int row;
        int col;
        int count = 0;
        
        IntegerVector Index = Rcpp::sample(nrow*ncol,nrow*ncol,false); //randomly sample indexs
        for (int i = 0; i < nrow*ncol; i++) {
            row = (Index(i) - 1) % nrow; //compute row index
            col = floor((Index(i) - 1)/nrow); //compute column index
            if ((domain_x(row, col) >= 0 && domain_x(row, col) <= InitialL/2) || (domain_x(row, col) >= Xmax - InitialL/2 && domain_x(row, col) <= Xmax)) {
                count += 1;
                if (count <= countIS[0]) { //ensure red cell limit not reached
                    domain(row, col) = 1; //assign red agent to location
                } else if (count <= countIS[0] + countIS[1]) { //ensure yellow cell limit not reached
                    domain(row, col) = 2; //assign yellow agent to location
                } else if (count <= countIS[0] + countIS[1] + countIS[2]) { //ensure green cell limit not reached
                    domain(row, col) = 3; //assign green agent to location
                }
            }
        }
    }
    
    return domain;
}
