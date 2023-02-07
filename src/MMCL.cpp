#include <RcppGSL.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <string>
#include <algorithm>

using namespace Rcpp;
using namespace std;

struct Reorder {
  double a;
  double b;
  int id;
};

// [[Rcpp::export]]
List MMCL(NumericVector y, NumericVector X, NumericVector d, NumericVector coef, NumericVector lambda,
         double tht, int frailty, int penalty, double tune, int a, int b, int p) {

  int N = a*b;
  vector<Reorder> LaR(N);

  for (int i = 0; i < N; i++) {
    LaR[i].a = lambda[i];
    LaR[i].b = y[i];
    LaR[i].id = i;
  }

  sort(LaR.begin(), LaR.end(), [](const Reorder& z1, const Reorder& z2){return(z1.b < z2.b);});
  vector<double> cumLam(N);
  NumericVector La(N);
  NumericVector La1(N);
  vector<int> idd(N);

  for (int i = 0; i < N; i++) {
    cumLam[i] = LaR[i].a;
  }
  partial_sum(cumLam.begin(), cumLam.end(), cumLam.begin(), plus<double>());
  for (int i = 0; i < N; i++) {
    La[LaR[i].id] = cumLam[i];
    La1[i] = cumLam[i];
    idd[i] = LaR[i].id;
  }
  
//   vy = as.vector(y)
//   vd = as.vector(d)
//   
//   
//   La = (cumsum(lambda[order(vy)]))[rank(vy)]
//   La = matrix(La, a, b)
//   
//   YpreExp = matrix(0, a, b)
//   for (i in 1:a) {
//     YpreExp[i,] = exp(X[i,,] %*% coef)
//   }
//   
//   A = rowSums(La*YpreExp)
//     B = apply((lambda*YpreExp)^d, 1, prod)
//     D = rowSums(d) 
//     
//     int0 <- vector("numeric", length = a)
//     int1 <- vector("numeric", length = a)
//     
//     if (frailty == "LogN" || frailty == "InvGauss" || frailty == "PVF") {
//       for (i in 1:a) {  
//         int0[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
//                             i = i, est.tht = est.tht, A = A, B = B, D = D, frailty = frailty, power = power, mode = 0)$value
//       }
//       
//       if (any(int0 == 0)) {
//         return(list(coef = rep(1/p, p), est.tht = est.tht, lambda = rep(1/N, N)))
//       }
//       
//       for (i in 1:a) {  
//         int1[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
//                             i = i, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i], frailty = frailty, power = power, mode = 1)$value
//       }
//     }
//     
//     if (frailty == "Gamma") {
//       C2 = 1/est.tht + A
//       A2 = 1/est.tht + D
//       int1 = A2/C2
//     }
//     
// # Update lambda Variables
//     ME = matrix(int1, a, b)
//       E_0 = as.vector(ME*YpreExp)
//       SUM_0 = cumsum((E_0[order(vy)])[seq(N,1,-1)])
//       SUM_0 = (SUM_0[seq(N,1,-1)])[rank(vy)] 
//     
//     lambda = vd/SUM_0
//     
// # Update coefficients
//     AVE_X = apply(abs(X), c(1,2), sum)
//       for(k in 1:p) {
//         E_1 = as.vector(ME*X[,,k]*YpreExp)
//         E_2 = as.vector(ME*AVE_X*abs(X[,,k])*YpreExp)
//         
//         SUM_1 = cumsum((E_1[order(vy)])[seq(N, 1, -1)])
//         SUM_1 = (SUM_1[seq(N, 1, -1)])[rank(vy)] 
//         SUM_2 = cumsum((E_2[order(vy)])[seq(N, 1, -1)])
//         SUM_2 = (SUM_2[seq(N, 1, -1)])[rank(vy)] 
//         
//         DE_1 = sum(d*X[,,k]) - sum(vd*SUM_1/SUM_0) 
//         DE_2 = -sum(vd*SUM_2/SUM_0) 
//         
//         if (!is.null(penalty)) {
//           if (penalty == "LASSO") {
//             DE_1 = DE_1 - N*sign(coef[k])*tune
//             DE_2 = DE_2 - N*tune/abs(coef[k])
//           }
//           if (penalty == "MCP") {
//             DE_1 = DE_1 - N*sign(coef[k])*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)
//             DE_2 = DE_2 - N*(tune - abs(coef[k])/3)*(abs(coef[k]) <= 3*tune)/abs(coef[k])
//           }
//           if (penalty == "SCAD") {
//             DE_1 = DE_1 - N*sign(coef[k])*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)
//             DE_2 = DE_2 - N*(tune*(abs(coef[k]) <= tune) + max(0,3.7*tune - abs(coef[k]))*(abs(coef[k]) > tune)/2.7)/abs(coef[k])
//           }
//         }
//         
//         coef[k] = coef[k] - DE_1/DE_2
//       }
//       
  List ret = List::create(La, La1, y, idd);
  return ret;
}

