// =============================================================================
// SIMULATION SETTING: Dominant Node Effects (Corrected Implementation)
//
// - alpha_1, beta_1: O(log(n)/n) and POSITIVE
// - alpha_i, beta_j (for i,j > 1): O(log(n)/n^2) and NEGATIVE
// - lambda_0 (baseline), gamma (covariate): O(log(n)/n^2)
//
// This corrected structure ensures positivity more robustly.
// ===============
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <numeric>

using namespace Rcpp; // <--- 在这里添加这一行

// 定义PI常量
#ifndef PI
#define PI 3.141592653589793
#endif

// 全局随机数生成器 (使用固定种子以保证可复现性)
std::mt19937 gen(123);

// 函数声明 (与原始文件保持一致)
double bs(double t, int n, double Clambda);
double fs(double t, int i, int n);
double fr(double t, int j, int n);
double fg(double t, int n);

// =============================================================================
// 核心模拟函数 (SimSetC)
// 签名不变，但内部逻辑为确保鲁棒性而优化
// =============================================================================
// [[Rcpp::export]]
List SimSetC(int n, double Clambda, double shift, Rcpp::NumericVector Zij) {
  
  // --- 输入验证 ---
  if (Rf_isNull(Zij.attr("dim"))) Rcpp::stop("'Zij' does not have a 'dim' attribute.");
  Rcpp::Dimension d = Zij.attr("dim");
  if (d.size() != 3) Rcpp::stop("'Zij' must have 3 dimensions.");
  if (d[0] != n || d[1] != n) Rcpp::stop("Dimensions of 'Zij' do not match 'n'.");
  
  int p = d[2];
  // 注意: 'shift' 参数保留以匹配签名，但在此模型中未使用。
  
  // --- 动态计算强度上界 (maxit) ---
  // 为了让 thinning 算法正确工作，maxit 必须大于等于任何时刻任何点对的最大强度
  double max_bs = bs(0.0, n, Clambda);
  double max_fs = fs(0.5, 1, n); // fs max at t=0.5 for i=1
  double max_fr = fr(0.0, 1, n); // fr max at t=0 for j=1
  
  double max_abs_fg = 0.0;
  if (p > 0) max_abs_fg = std::abs(fg(0.25, n));
  
  double max_abs_Zij = 0.0;
  if (p > 0) {
    for (int k=0; k<Zij.size(); ++k) {
      if (std::abs(Zij[k]) > max_abs_Zij) max_abs_Zij = std::abs(Zij[k]);
    }
  }
  
  double maxit = max_bs + max_fs + max_fr + max_abs_fg * max_abs_Zij * p;
  maxit = std::ceil(maxit * 1.1); // 增加10%的安全余量并向上取整
  
  // --- Thinning 算法 ---
  double temp(0.0), tp(0.0), rej(0.0);
  std::vector<int> se, re;
  std::vector<double> te;
  
  int kp(0);
  std::poisson_distribution<int> Pdis(maxit);
  std::uniform_real_distribution<double> Udis(0.0, 1.0);
  
  bool negative_intensity_detected = false;
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      if (i == j) continue;
      
      kp = Pdis(gen);
      
      if (kp > 0) {
        for (int z = 0; z < kp; z++) {
          tp = Udis(gen);
          
          // 计算总强度
          temp = bs(tp, n, Clambda) + fs(tp, i, n) + fr(tp, j, n);
          for (int di = 0; di < p; di++) {
            temp += fg(tp, n) * Zij[i-1 + n*(j-1) + n*n*di];
          }
          
          if (temp < 0) {
            if (!negative_intensity_detected) {
              Rcpp::warning("Calculated intensity is negative. Consider increasing 'Clambda'. Simulation proceeds by setting negative values to 0. (This warning appears only once).");
              negative_intensity_detected = true;
            }
            temp = 0.0; // 截断负值
          }
          
          if (temp > maxit) {
            Rcpp::Rcout << "Calculated intensity " << temp << " exceeds maxit " << maxit << std::endl;
            Rcpp::stop("Intensity exceeds maxit. Increase safety margin or check calculations.");
          }
          
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
  
  return List::create(_["se"] = se, _["re"] = re, _["te"] = te);
}


// =============================================================================
// 新设计的强度函数 (函数签名不变)
// =============================================================================

// --- 内部可调常量 ---
// 您可以修改这些值来改变 alpha, beta, gamma 的相对强度
const double CALPHA_INTERNAL = 3.0;
const double CBETA_INTERNAL = 3.0;
const double CGAMMA_INTERNAL = 5.0;

// --- Baseline: lambda_0(t) ---
// 数量级: O(log(n)/n^2), 正
// 由 Clambda 控制其大小以保证总强度为正
// [[Rcpp::export]]
double bs(double t, int n, double Clambda) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  return Clambda * std::log(n_val) / (n_val * n_val);
}

// --- Sender Effect: alpha_i(t) ---
// 约束: alpha_1 > 0, O(log(n)/n); alpha_i < 0, O(log(n)/n^2) for i>1; sum(alpha_i) = 0
// [[Rcpp::export]]
double fs(double t, int i, int n) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  double time_component = 0.5 + 0.5 * sin(2.0 * PI * t); // 范围 [0, 1]
  
  double alpha_1_val = CALPHA_INTERNAL * std::log(n_val) / n_val * time_component;
  
  if (i == 1) {
    return alpha_1_val; // 正, O(log(n)/n)
  } else {
    return -alpha_1_val / (n_val - 1.0); // 负, O(log(n)/n^2)
  }
}

// --- Receiver Effect: beta_j(t) ---
// 约束: beta_1 > 0, O(log(n)/n); beta_j < 0, O(log(n)/n^2) for j>1; sum(beta_j) = 0
// [[Rcpp::export]]
double fr(double t, int j, int n) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  double time_component = 0.5 + 0.5 * cos(2.0 * PI * t); // 范围 [0, 1]
  
  double beta_1_val = CBETA_INTERNAL * std::log(n_val) / n_val * time_component;
  
  if (j == 1) {
    return beta_1_val; // 正, O(log(n)/n)
  } else {
    return -beta_1_val / (n_val - 1.0); // 负, O(log(n)/n^2)
  }
}

// --- Covariate Effect: gamma(t) ---
// 数量级: O(log(n)/n^2)
// [[Rcpp::export]]
double fg(double t, int n) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  double scale_factor = std::log(n_val) / (n_val * n_val);
  return CGAMMA_INTERNAL * scale_factor * sin(2.0 * PI * t);
}


// =============================================================================
// 对应的积分函数 (从 0 到 t)
// =============================================================================

// [[Rcpp::export]]
NumericVector bsI_cpp(NumericVector t, int n, double Clambda) {
  double const_val = bs(0.0, n, Clambda); // bs is constant in t
  return t * const_val;
}

// [[Rcpp::export]]
NumericVector fsI_cpp(NumericVector t, int i, int n) {
  if (n <= 1) return t * 0.0;
  double n_val = static_cast<double>(n);
  
  NumericVector time_integral = 0.5 * t - Rcpp::cos(2.0 * PI * t) / (4.0 * PI) + 1.0 / (4.0 * PI);
  double alpha_1_scale = CALPHA_INTERNAL * std::log(n_val) / n_val;
  NumericVector alpha_1_integral = alpha_1_scale * time_integral;
  
  if (i == 1) {
    return alpha_1_integral;
  } else {
    return -alpha_1_integral / (n_val - 1.0);
  }
}

// [[Rcpp::export]]
NumericVector frI_cpp(NumericVector t, int j, int n) {
  if (n <= 1) return t * 0.0;
  double n_val = static_cast<double>(n);
  
  NumericVector time_integral = 0.5 * t + Rcpp::sin(2.0 * PI * t) / (4.0 * PI);
  double beta_1_scale = CBETA_INTERNAL * std::log(n_val) / n_val;
  NumericVector beta_1_integral = beta_1_scale * time_integral;
  
  if (j == 1) {
    return beta_1_integral;
  } else {
    return -beta_1_integral / (n_val - 1.0);
  }
}

// [[Rcpp::export]]
NumericVector fgI_cpp(NumericVector t, int n) {
  if (n <= 1) return t * 0.0;
  double n_val = static_cast<double>(n);
  double scale_factor = std::log(n_val) / (n_val * n_val);
  double coefficient = CGAMMA_INTERNAL * scale_factor;
  
  return coefficient * (1.0 - Rcpp::cos(2.0 * PI * t)) / (2.0 * PI);
}
