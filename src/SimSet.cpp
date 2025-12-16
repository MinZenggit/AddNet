#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>

using namespace Rcpp;

#define PI 3.14159265

double fs(double t, const double& shift1, int i, const int &n);
double fr(double t, const double& shift2, int i, const int &n);
double fg(double t, const int& n);

std::mt19937 gen(123);

// [[Rcpp::export]]
List SimSetC(int n, double shift1, double shift2, NumericVector Zij) {
  
  if (Rf_isNull(Zij.attr("dim"))) {
    throw std::runtime_error("'x' does not have 'dim' attibute.");
  }
  Rcpp::Dimension d = Zij.attr("dim");
  if (d.size() != 3) {
    throw std::runtime_error("'x' must have 3 dimensions.");
  }
  
  std::size_t p = d[2];
  
  if (d[0] != n || d[1] != n)
    return 0;
  
  double maxit(10.0), temp(0.0), tp(0.0), rej(0.0);
  std::vector<int> se;
  std::vector<int> re;
  std::vector<double> te;
  
  int kp(0);
  std::poisson_distribution<int> Pdis(maxit);
  std::uniform_real_distribution<double> Udis(0.0,1.0);
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++ ) {
      if (i != j) {
        kp = Pdis(gen);
        
        if (kp != 0) {
          for (int z = 0; z < kp; z++) {
            tp = Udis(gen);
            temp = 2/sqrt(n) + fs(tp, shift1, i, n) + fr(tp, shift2, j, n);
            for (int di = 0; di < p; di++)
              temp += fg(tp, n) * Zij[i-1+n*(j-1)+n*n*di];
            
            rej = Udis(gen);
            if (rej < temp / maxit) {
              se.push_back(i);
              re.push_back(j);
              te.push_back(tp);
            }
          }
        }
        
      }
    }
  }
  
  return List::create(se, re, te);
}


// SIM1.6
double fs(double t, const double& shift1, int i, const int &n) {
  if (i <= n / 2) return (0.5 + 0.5*sin(2*PI*t))/sqrt(n);
  return -(0.5 + 0.5*sin(2*PI*t))/sqrt(n);
}

double fr(double t, const double& shift2, int i, const int &n) {
  if (i <= n / 2) return (0.5 + 0.5*cos(2*PI*t))/sqrt(n);
  return -(0.5 + 0.5*cos(2*PI*t))/sqrt(n);
}

double fg(double t, const int& n) {
  return 0.25/sqrt(n);
}