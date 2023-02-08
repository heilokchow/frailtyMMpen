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

double logN3int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::pow(std::log(x), 2.0)/(x*std::sqrt(2*M_PI*k.s))*std::exp(-std::pow(std::log(x), 2.0)/(2*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)/k.por;
  return ret;
  
}

double InvG1int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = 1/(std::sqrt(2*M_PI*k.s)*std::pow(x, 1.5))*std::exp(-std::pow(x-1, 2.0)/(2*x*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x);
  return ret;
  
}

double InvG2int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = 1/(std::sqrt(2*M_PI*k.s)*std::pow(x, 1.5))*std::exp(-std::pow(x-1, 2.0)/(2*x*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)*x/k.por;
  return ret;
  
}

double InvG3int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::pow(x-1, 2.0)/x/(std::sqrt(2*M_PI*k.s)*std::pow(x, 1.5))*std::exp(-std::pow(x-1, 2.0)/(2*x*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)/k.por;
  return ret;
  
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void computeLAM(NumericVector& LAM, const NumericVector& lambda, const NumericVector& y, const int& N, int mode) {
  
  std::vector<Reorder> LaR(N);
  std::vector<double> cumLam(N);
  
  for (int i = 0; i < N; i++) {
    LaR[i].a = lambda[i];
    LaR[i].b = y[i];
    LaR[i].id = i;
  }
  
  if (mode == 0) {
    std::sort(LaR.begin(), LaR.end(), [](const Reorder& z1, const Reorder& z2){return(z1.b < z2.b);});
  }
  
  if (mode == 1) {
    std::sort(LaR.begin(), LaR.end(), [](const Reorder& z1, const Reorder& z2){return(z1.b > z2.b);});
  }
  
  for (int i = 0; i < N; i++) {
    cumLam[i] = LaR[i].a;
  }
  
  std::partial_sum(cumLam.begin(), cumLam.end(), cumLam.begin(), std::plus<double>());
  
  for (int i = 0; i < N; i++) {
    LAM[LaR[i].id] = cumLam[i];
  }
  
}

// [[Rcpp::export]]
List MMCL(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0,
         const double& tht0, int frailty, int penalty, double tune, int a, int b, int p) {

  int N = a*b;

  NumericVector coef = clone(coef0);
  NumericVector lambda = clone(lambda0);
  double tht = tht0;
  
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
  NumericVector int3(a, 0.0);

  NumericVector ME(N, 0.0);
  NumericVector E0(N, 0.0);
  NumericVector SUM0(N, 0.0);
  
  NumericVector AVEX(N, 0.0);
  NumericVector E1(N, 0.0);
  NumericVector D11(N, 0.0);
  NumericVector SUM1(N, 0.0);
  NumericVector E2(N, 0.0);
  NumericVector SUM2(N, 0.0);
  NumericVector D22(N, 0.0);
  
  double D1(0.0), D2(0.0);
  
  computeLAM(La, lambda, y, N, 0);
  
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
  
  if (frailty == 1 || frailty == 2) {
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2, F3;
    double result(0.0), error;
    
    if (frailty == 1) {
      F1.function = &logN1int;
      F2.function = &logN2int;
      F3.function = &logN3int;
    }
    
    if (frailty == 2) {
      F1.function = &InvG1int;
      F2.function = &InvG2int;
      F3.function = &InvG3int;
    }
    
    F1.params = &kint;
    F2.params = &kint;
    F3.params = &kint;
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
      
      gsl_integration_qagiu (&F3, 0, 0, 1e-7, 1000, w, &result, &error);
      
      int3[i] = result;
    }
    
    gsl_integration_workspace_free (w);
    
    tht = std::accumulate(int3.begin(), int3.end(), 0.0) / a;
    
  }

  // Update lambda Variables
  
  for (int i = 0; i < b; i++) {
    ME[Range(i*a, (i+1)*a-1)] = int2;
  }
  
  E0 = ME * YpreExp;
  computeLAM(SUM0, E0, y, N, 1);
  lambda = d / SUM0;
  
  // Update coefficients
  
  for (int i = 0; i < p; i++) {
    AVEX += abs(X[Range(i*N, (i+1)*N-1)]);
  } 
  
  for (int i = 0; i < p; i++) {
    
    E1 = ME * X[Range(i*N, (i+1)*N-1)] * YpreExp;
    E2 = ME * AVEX *X[Range(i*N, (i+1)*N-1)] * YpreExp;
    
    computeLAM(SUM1, E1, y, N, 1);
    computeLAM(SUM2, E2, y, N, 1);
    
    D11 = d * X[Range(i*N, (i+1)*N-1)] - d * SUM1 / SUM0;
    D22 = - d * SUM2 / SUM0;
    
    D1 = std::accumulate(D11.begin(), D11.end(), 0.0);
    D2 = std::accumulate(D22.begin(), D22.end(), 0.0);
    
    if (penalty != 0) {
      
      if (penalty == 1) {
        D1 -= N*sgn(coef[i])*tune;
        D2 -= N*tune/std::abs(coef[i]);
      }
      
      if (penalty == 2) {
        D1 -= N*sgn(coef[i])*(tune - std::abs(coef[i])/3)*(std::abs(coef[i]) <= 3*tune);
        D2 -= N*(tune - std::abs(coef[i])/3)*(std::abs(coef[i]) <= 3*tune)/std::abs(coef[i]);
      }
      
      if (penalty == 3) {
        D1 -= N*sgn(coef[i])*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, 3.7*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/2.7);
        D2 -= N*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, 3.7*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/2.7)/std::abs(coef[i]);
      }
    }
    
    coef[i] -= D1/D2;
    
  }

  List ret = List::create(_["coef"] = coef, _["est.tht"] = tht, _["lambda"] = lambda, AVEX, SUM1, D2, ME, SUM0, int1, int2);
  return ret;
}

