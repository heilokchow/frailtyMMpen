#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector FaccCal(NumericVector y, NumericVector x, NumericVector Target, NumericVector A, NumericVector C) {
  
  Rcpp::Dimension d = y.attr("dim");
  int a = d[0], b = d[1];
  double temp(0.0);
  
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < b; j++) {
      temp = 0.0;
      for (int k = 0; k < a; k++) {
        for (int h = 0; h < b; h++) {
          temp += A[k] * (y[k+h*a] >= y[i+j*a]) * x[k+h*a] / C[k];
        }
      }
      Target[i+j*a] = temp;
    }
  }
  
  return Target;
}
