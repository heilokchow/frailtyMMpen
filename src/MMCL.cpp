#include <RcppGSL.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <string>
#include <algorithm>

using namespace Rcpp;

struct Reorder {
  double a;
  double b;
  int id;
};

struct intParams {
  double a;
  double b;
  double d;
  double s;
  double por;
};

double logN1int(double x, void * params) {
  intParams k = *(intParams *) params;
  double ret = 1/(x*std::sqrt(2*M_PI*k.s))*std::exp(-std::pow(std::log(x), 2.0)/(2*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x);
  return ret;
}

double logN2int(double x, void * params) {
  intParams k = *(intParams *) params;
  double ret = 1/(x*std::sqrt(2*M_PI*k.s))*std::exp(-std::pow(std::log(x), 2.0)/(2*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)*x/k.por;
  return ret;
}

// [[Rcpp::export]]
List MMCL(NumericVector y, NumericVector X, NumericVector d, NumericVector coef, NumericVector lambda,
         double tht, int frailty, int penalty, double tune, int a, int b, int p) {

  int N = a*b;
  std::vector<Reorder> LaR(N);
  std::vector<double> cumLam(N);
  NumericVector La(N);
  NumericVector YpreExp(N, 0.0);
  intParams kint;
  
  NumericVector AA(N);
  NumericVector BB(N);
  
  NumericVector A(a, 0.0);
  NumericVector B(a, 1.0);
  NumericVector D(a, 0.0);
  
  NumericVector int1(a, 0.0);
  NumericVector int2(a, 0.0);

  for (int i = 0; i < N; i++) {
    LaR[i].a = lambda[i];
    LaR[i].b = y[i];
    LaR[i].id = i;
  }

  std::sort(LaR.begin(), LaR.end(), [](const Reorder& z1, const Reorder& z2){return(z1.b < z2.b);});

  for (int i = 0; i < N; i++) {
    cumLam[i] = LaR[i].a;
  }
  
  std::partial_sum(cumLam.begin(), cumLam.end(), cumLam.begin(), std::plus<double>());
  
  for (int i = 0; i < N; i++) {
    La[LaR[i].id] = cumLam[i];
  }
  
  for (int i = 0; i < p; i++) {
    YpreExp += X[Range(i*N, (i+1)*N-1)] * coef[i];
  }

  YpreExp = exp(YpreExp);
  
  AA = La * YpreExp;
  BB = lambda * YpreExp;
  
  for (int j = 0; j < b; j++) {
    for (int i = 0; i < a; i++) {
      A[i] += AA[j*a + i];
      if (d[j*a + i] == 1) {
        B[i] *= BB[j*a + i];
        D[i]++;
      }
    }
  }
  
  if (frailty == 1) {
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2;
    double result(0.0), error;
    
    F1.function = &logN1int;
    F2.function = &logN2int;
    F1.params = &kint;
    F2.params = &kint;
    kint.s = tht;
    
    for (int i = 0; i < a; i++) {
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      kint.por = result;
      int1[i] = result;
      
      gsl_integration_qagiu (&F2, 0, 0, 1e-7, 1000, w, &result, &error);
      
      int2[i] = result;
    }
    
    gsl_integration_workspace_free (w);
  }
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
// if (frailty == "LogN" || frailty == "InvGauss") {
//   int2 <- vector("numeric", length = a)
//   for (i in 1:a) {  
//     int2[i] = integrate(int_tao, lower = 0, upper = 20, stop.on.error = FALSE,
//                         i = i, est.tht = est.tht, A = A, B = B, D = D, tao0 = int0[i], frailty = frailty, mode = 2)$value
//   }
//   est.tht = sum(int2)/a
// }
// 
// if (est.tht < 0) {
//   est.tht = 1
// } 
// 
// if (frailty == "Gamma") {
//   Q01 = a*(digamma(1/est.tht)+log(est.tht)-1)/(est.tht^2) + sum(A2/C2-digamma(A2)+log(C2))/(est.tht^2) 
//   Q02 = a*(3-2*digamma(1/est.tht)-2*log(est.tht))/(est.tht^3)+2*sum(digamma(A2)-log(C2)-A2/C2)/(est.tht^3) - a*trigamma(1/est.tht)/(est.tht^4)
//   
//   est.tht1 = est.tht - Q01/Q02
//   if(est.tht1 > 0) {
//     est.tht = est.tht1 
//   }
// }

  List ret = List::create(La, int1, int2, AA, BB, A, B, D);
  return ret;
}

