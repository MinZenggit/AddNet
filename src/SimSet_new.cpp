#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <cstdio>

using namespace Rcpp;

// --- 全局常量定义 ---
#ifndef PI
#define PI 3.141592653589793
#endif

const double EFFECT_SIZE_MULTIPLIER = 10.0;

// 函数声明
double bs_new(double t, int n, double C_lambda, double Csparse);
double fs_new(double t, int i, int n, double Csparse);
double fr_new(double t, int j, int n, double Csparse);
double fg_new(double t, int n, double Csparse);

// 全局随机数生成器
std::mt19937 gen(123);

// [[Rcpp::export]]
List SimSetC(int n, double C_lambda, double Csparse, NumericVector Zij) {
  
  if (Rf_isNull(Zij.attr("dim"))) {
    throw std::runtime_error("'Zij' does not have 'dim' attribute.");
  }
  Rcpp::Dimension d = Zij.attr("dim");
  if (d.size() != 3) {
    throw std::runtime_error("'Zij' must have 3 dimensions.");
  }
  
  std::size_t p = d[2];
  
  if (d[0] != n || d[1] != n)
    Rcpp::stop("Dimensions of 'Zij' do not match 'n'.");
  
  double maxit(10.0), temp(0.0), tp(0.0), rej(0.0); 
  std::vector<int> se;
  std::vector<int> re;
  std::vector<double> te;
  
  int kp(0);
  std::poisson_distribution<int> Pdis(maxit);
  std::uniform_real_distribution<double> Udis(0.0,1.0);
  
  bool negative_intensity_detected = false;
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++ ) {
      if (i != j) {
        kp = Pdis(gen);
        
        if (kp != 0) {
          for (int z = 0; z < kp; z++) {
            tp = Udis(gen);
            
            double temp = bs_new(tp, n, C_lambda, Csparse) + 
              fs_new(tp, i, n, Csparse) + 
              fr_new(tp, j, n, Csparse);
            
            for (int di = 0; di < p; di++)
              temp += fg_new(tp, n, Csparse) * Zij[i-1+n*(j-1)+n*n*di];
            
            if (temp < 0) {
              if (!negative_intensity_detected) {
                char buffer[256];
                snprintf(buffer, sizeof(buffer),
                         "Warning: Negative intensity lambda_ij(t) = %.4f detected for pair (i=%d, j=%d) at t=%.4f. "
                         "The simulation will proceed by setting it to 0, but results may be biased. "
                         "Consider increasing C_lambda.",
                         temp, i, j, tp);
                Rcpp::warning(buffer);
                negative_intensity_detected = true;
              }
            }
            
            temp = std::max(0.0, temp);
            
            double rej = Udis(gen);
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
  
  return List::create(_["se"] = se, _["re"] = re, _["te"] = te);
}

// --- 强度函数 ---

// [[Rcpp::export]]
double bs_new(double t, int n, double C_lambda, double Csparse) {
  double n_double = static_cast<double>(n);
  double scale_factor = std::pow(n_double, -2.0 * Csparse);
  return C_lambda * scale_factor;
}

// [[Rcpp::export]]
double fs_new(double t, int i, int n, double Csparse) {
  // *** MODIFIED LOGIC ***
  // 节点效应被归一化，其波动范围不再与 n 成正比
  double n_double = static_cast<double>(n);
  double i_double = static_cast<double>(i);
  double scale_factor = std::pow(n_double, -2.0 * Csparse);
  
  // 将 (i - (n+1)/2) 归一化到 [-1, 1] 区间
  double normalized_effect = (i_double - (n_double + 1.0) / 2.0) / (n_double / 2.0);
  
  double effect = EFFECT_SIZE_MULTIPLIER * normalized_effect * std::sin(2.0 * PI * t);
  return effect * scale_factor;
}

// [[Rcpp::export]]
double fr_new(double t, int j, int n, double Csparse) {
  // *** MODIFIED LOGIC ***
  // 节点效应被归一化，其波动范围不再与 n 成正比
  double n_double = static_cast<double>(n);
  double j_double = static_cast<double>(j);
  double scale_factor = std::pow(n_double, -2.0 * Csparse);
  
  // 将 (j - (n+1)/2) 归一化到 [-1, 1] 区间
  double normalized_effect = (j_double - (n_double + 1.0) / 2.0) / (n_double / 2.0);
  
  double effect = EFFECT_SIZE_MULTIPLIER * normalized_effect * std::cos(2.0 * PI * t);
  return effect * scale_factor;
}

// [[Rcpp::export]]
double fg_new(double t, int n, double Csparse) {
  double n_double = static_cast<double>(n);
  double scale_factor = std::pow(n_double, -2.0 * Csparse);
  return EFFECT_SIZE_MULTIPLIER * std::sin(2.0 * PI * t) * scale_factor;
}


// --- 对应的积分函数 ---

// [[Rcpp::export]]
NumericVector bsI_cpp(NumericVector t, int n, double C_lambda, double Csparse) {
  double const_val = bs_new(0.0, n, C_lambda, Csparse);
  return t * const_val;
}

// [[Rcpp::export]]
NumericVector fsI_cpp(NumericVector t, int i, int n, double Csparse) {
  // *** MODIFIED LOGIC ***
  // 对应新的 fs_new 的积分
  double n_double = static_cast<double>(n);
  double i_double = static_cast<double>(i);
  double scale_factor = std::pow(n_double, -2.0 * Csparse);
  
  // 归一化后的效应部分
  double normalized_effect = (i_double - (n_double + 1.0) / 2.0) / (n_double / 2.0);
  
  double const_part = EFFECT_SIZE_MULTIPLIER * normalized_effect * scale_factor;
  return const_part * (1.0 - Rcpp::cos(2.0 * PI * t)) / (2.0 * PI);
}

// [[Rcpp::export]]
NumericVector frI_cpp(NumericVector t, int j, int n, double Csparse) {
  // *** MODIFIED LOGIC ***
  // 对应新的 fr_new 的积分
  double n_double = static_cast<double>(n);
  double j_double = static_cast<double>(j);
  double scale_factor = std::pow(n_double, -2.0 * Csparse);
  
  // 归一化后的效应部分
  double normalized_effect = (j_double - (n_double + 1.0) / 2.0) / (n_double / 2.0);
  
  double const_part = EFFECT_SIZE_MULTIPLIER * normalized_effect * scale_factor;
  return const_part * Rcpp::sin(2.0 * PI * t) / (2.0 * PI);
}

// [[Rcpp::export]]
NumericVector fgI_cpp(NumericVector t, int n, double Csparse) {
  double n_double = static_cast<double>(n);
  double scale_factor = std::pow(n_double, -2.0 * Csparse);
  double const_part = EFFECT_SIZE_MULTIPLIER * scale_factor;
  return const_part * (1.0 - Rcpp::cos(2.0 * PI * t)) / (2.0 * PI);
}
