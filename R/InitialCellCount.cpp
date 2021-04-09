#include<Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
IntegerMatrix InitialCellCount(IntegerMatrix &domain, NumericMatrix &domain_x, NumericMatrix &domain_y, int &nrow, int &ncol, double &Xmax, double &Ymax) {
  
  //initialise variables
  IntegerMatrix CellCounts(3,1);
  int Nred = 0;
  int Nyellow = 0;
  int Ngreen = 0;
  
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      if (domain_x(i,j) < Xmax && domain_y(i,j) < Ymax) { //ensure cells inside domain
        if (domain(i,j) == 1){ 
          Nred += 1;
        } else if (domain(i,j) == 2) {
          Nyellow += 1;
        } else if (domain(i,j) == 3) {
          Ngreen += 1;
        }
      } 
    }
  }
  
  CellCounts(0,0) = Nred;
  CellCounts(1,0) = Nyellow;
  CellCounts(2,0) = Ngreen;
  
  return CellCounts;
  
}
