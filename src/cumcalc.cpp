#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void cumcalc(NumericVector Patemp, NumericVector N_tC, NumericVector B, int n,
                int nb, int z3, double t, NumericVector t_sep_t) {
  int z33 = z3 - 1;
  double ktemp(0.0);
  int n1 = n*(n-1);
//对于每一个固定的时间节点j，如果这个事件发生的事件小于j，那么参数的积分就累加上
  for (int j = 0; j < 100; j++) {
    if (t < t_sep_t[j]) {
      for (int i = 0; i < nb; i++) { //nb = (2*n-1)+p
        B[j*nb + i] += Patemp[i];
      }
    }// z3 = hash[(z1-1)*n+z2]
    N_tC[j*n1 + z33] += 1; //这里可能要放到if条件句内部才对
  }
}

