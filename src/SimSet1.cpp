#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <numeric>

using namespace Rcpp;

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
// =============================================================================
// [[Rcpp::export]]
List SimSetC(int n, double Clambda, double shift, Rcpp::NumericVector Zij) {
  
  // --- 输入验证 ---
  if (Rf_isNull(Zij.attr("dim"))) stop("'Zij' does not have a 'dim' attribute.");
  Dimension d = Zij.attr("dim");
  if (d.size() != 3) stop("'Zij' must have 3 dimensions.");
  if (d[0] != n || d[1] != n) stop("Dimensions of 'Zij' do not match 'n'.");
  
  int p = d[2];
  // 注意: 'shift' 参数保留以匹配签名，但在此模型中未使用。
  
  // --- 动态计算强度上界 (maxit) ---
  // maxit 必须大于等于任何时刻任何点对的最大强度
  double max_bs = bs(0.0, n, Clambda);
  double max_fs = fs(0.5, 1, n); // fs max at t=0.5 for i in positive group
  double max_fr = fr(0.0, 1, n); // fr max at t=0 for j in positive group
  
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
          
          temp = bs(tp, n, Clambda) + fs(tp, i, n) + fr(tp, j, n);
          for (int di = 0; di < p; di++) {
            temp += fg(tp, n) * Zij[i-1 + n*(j-1) + n*n*di];
          }
          
          if (temp < 0) {
            if (!negative_intensity_detected) {
              warning("Calculated intensity is negative. CRITICAL: You MUST increase 'Clambda' for valid results. Simulation proceeds by setting negative values to 0. (This warning appears only once).");
              negative_intensity_detected = true;
            }
            temp = 0.0;
          }
          
          if (temp > maxit) {
            Rcout << "Calculated intensity " << temp << " exceeds maxit " << maxit << std::endl;
            stop("Intensity exceeds maxit. Increase safety margin or check calculations.");
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
// 新设计的强度函数 (阶数 O(1/sqrt(n)))
// =============================================================================

// --- 内部可调常量 ---
const double CALPHA_INTERNAL = 1.0;
const double CBETA_INTERNAL = 1.0;
const double CGAMMA_INTERNAL = 0.05;

// --- Baseline: lambda_0(t) ---
// 阶数: O(1/sqrt(n)), 正
// [[Rcpp::export]]
double bs(double t, int n, double Clambda) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  return Clambda / std::sqrt(n_val);
}

// --- Sender Effect: alpha_i(t) ---
// 阶数: O(1/sqrt(n)), sum(alpha_i) = 0
// [[Rcpp::export]]
double fs(double t, int i, int n) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = CALPHA_INTERNAL / std::sqrt(n_val);
  double time_component = 0.5 + 0.5 * sin(2.0 * PI * t);
  double base_effect = scale_factor * time_component;
  
  if (i <= n_pos) {
    return base_effect; // 正效应
  } else {
    return -base_effect * (n_pos_d / n_neg_d); // 负效应，以保证总和为0
  }
}

// --- Receiver Effect: beta_j(t) ---
// 阶数: O(1/sqrt(n)), sum(beta_j) = 0
// [[Rcpp::export]]
double fr(double t, int j, int n) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = CBETA_INTERNAL / std::sqrt(n_val);
  double time_component = 0.5 + 0.5 * cos(2.0 * PI * t);
  double base_effect = scale_factor * time_component;
  
  if (j <= n_pos) {
    return base_effect; // 正效应
  } else {
    return -base_effect * (n_pos_d / n_neg_d); // 负效应，以保证总和为0
  }
}

// --- Covariate Effect: gamma(t) ---
// 阶数: O(1/sqrt(n))
// [[Rcpp::export]]
double fg(double t, int n) {
  if (n <= 1) return 0.0;
  double n_val = static_cast<double>(n);
  return CGAMMA_INTERNAL / std::sqrt(n_val) * sin(2.0 * PI * t);
}


// =============================================================================
// 对应的积分函数 (从 0 到 t)
// =============================================================================

// [[Rcpp::export]]
NumericVector bsI_cpp(NumericVector t, int n, double Clambda) {
  double const_val = bs(0.0, n, Clambda);
  return t * const_val;
}

// [[Rcpp::export]]
NumericVector fsI_cpp(NumericVector t, int i, int n) {
  if (n <= 1) return t * 0.0;
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = CALPHA_INTERNAL / std::sqrt(n_val);
  NumericVector time_integral = 0.5 * t - Rcpp::cos(2.0 * PI * t) / (4.0 * PI) + 1.0 / (4.0 * PI);
  NumericVector base_integral = scale_factor * time_integral;
  
  if (i <= n_pos) {
    return base_integral;
  } else {
    return -base_integral * (n_pos_d / n_neg_d);
  }
}

// [[Rcpp::export]]
NumericVector frI_cpp(NumericVector t, int j, int n) {
  if (n <= 1) return t * 0.0;
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = CBETA_INTERNAL / std::sqrt(n_val);
  NumericVector time_integral = 0.5 * t + Rcpp::sin(2.0 * PI * t) / (4.0 * PI);
  NumericVector base_integral = scale_factor * time_integral;
  
  if (j <= n_pos) {
    return base_integral;
  } else {
    return -base_integral * (n_pos_d / n_neg_d);
  }
}

// [[Rcpp::export]]
NumericVector fgI_cpp(NumericVector t, int n) {
  if (n <= 1) return t * 0.0;
  double n_val = static_cast<double>(n);
  double coefficient = CGAMMA_INTERNAL / std::sqrt(n_val);
  
  return coefficient * (1.0 - Rcpp::cos(2.0 * PI * t)) / (2.0 * PI);
}
