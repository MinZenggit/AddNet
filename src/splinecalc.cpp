#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void splinecalc(NumericVector Patemp, NumericVector N_tall, NumericVector N_tallM, NumericVector N_tallMC, NumericVector b, int n,
                       int nb, int z1, int z2, int z3, double h1, double t, NumericVector t_sep_t) {
  int z11 = z1 - 1;
  int z22 = z2 - 1;
  int z33 = z3 - 1;
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double diff(0.0);
  double ktemp(0.0);
  int n1 = n*(n-1);
  int n2 = n*n;

  for (int j = 0; j < 100; j++) {
    diff = (t - t_sep_t[j])/h1;
    ktemp = inv_sqrt_2pi / h1 * std::exp(-0.5 * diff * diff);
    N_tall[j*n1 + z33] += ktemp;
    N_tallMC[j*n1 + z33] += ktemp * ktemp;
    N_tallM[j*n2 + z22*n + z11] += ktemp * 0.28209/h1;
    // N_tallM[j*n2 + z22*n + z11] += ktemp * ktemp;
    for (int i = 0; i < nb; i++) {
      b[j*nb + i] += Patemp[i] * ktemp;
    }
  }
}

