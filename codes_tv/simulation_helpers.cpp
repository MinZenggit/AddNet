// 文件名: simulation_helpers.cpp (已更正)

#include <Rcpp.h>
#include <cmath>
#include <vector>

// 定义常量
#define PI 3.14159265

// 声明辅助函数
double fs(double t, int i, int n);
double fr(double t, int i, int n);
double fg(double t, int n);

// [[Rcpp::plugins(cpp11)]]

// 核心事件生成函数，用于单个时间段
// [[Rcpp::export]]
Rcpp::List generate_events_for_period_cpp(int n, double t_start, double t_end, 
                                          const Rcpp::NumericVector& Z_period, 
                                          double max_lambda, int p) {
  
  Rcpp::Dimension d = Z_period.attr("dim");
  if (d[0] != n || d[1] != n || d[2] != p) {
    Rcpp::stop("Dimension of Z_period does not match n and p.");
  }
  
  std::vector<int> se;
  std::vector<int> re;
  std::vector<double> te;
  
  double intensity_base, intensity_total;
  
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= n; ++j) {
      if (i == j) continue;
      
      // --- 更正开始 ---
      // 使用 R 的 C API 直接调用 rpois，更稳定高效
      // R::rpois 返回 double，需要转换为 int
      int num_potential_events = static_cast<int>(R::rpois(max_lambda * (t_end - t_start)));
      // --- 更正结束 ---
      
      if (num_potential_events > 0) {
        // 在循环内生成每个事件的时间和拒绝值
        for (int k = 0; k < num_potential_events; ++k) {
          // --- 更正开始 ---
          // 使用 R 的 C API 直接调用 runif
          double t_current = R::runif(t_start, t_end);
          double rejection_value = R::runif(0.0, max_lambda);
          // --- 更正结束 ---
          
          intensity_base = 2.0 / std::sqrt(n) + fs(t_current, i, n) + fr(t_current, j, n);
          
          double z_effect = 0.0;
          if (p > 0) {
            int z_idx = (i-1) + n * (j-1);
            z_effect = fg(t_current, n) * Z_period[z_idx];
          }
          
          intensity_total = intensity_base + z_effect;
          
          if (rejection_value < intensity_total) {
            se.push_back(i);
            re.push_back(j);
            te.push_back(t_current);
          }
        }
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("sender") = se,
    Rcpp::Named("receiver") = re,
    Rcpp::Named("time") = te
  );
}

// 强度函数的组成部分 (无变化)
double fs(double t, int i, int n) {
  if (i <= n / 2) return (0.5 + 0.5 * std::sin(2 * PI * t)) / std::sqrt(n);
  return -(0.5 + 0.5 * std::sin(2 * PI * t)) / std::sqrt(n);
}

double fr(double t, int i, int n) {
  if (i <= n / 2) return (0.5 + 0.5 * std::cos(2 * PI * t)) / std::sqrt(n);
  return -(0.5 + 0.5 * std::cos(2 * PI * t)) / std::sqrt(n);
}

double fg(double t, int n) {
  return 0.5 * std::sin(2 * PI * t) / std::sqrt(n); 
}
