// =============================================================================
// SIMULATION SETTING: All terms are of order O(log(n)/n^2)
//
// This file defines the data generating process where the baseline rate,
// sender effects, receiver effects, and covariate effects are all scaled
// by log(n)/n^2. This setting is designed for testing scenarios with
// relatively weak signals.
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <cstdio> // *** MODIFIED: Added for snprintf (used in warning) ***

using namespace Rcpp;

// 定义PI常量
#ifndef PI
#define PI 3.141592653589793
#endif

// 函数声明
double bs(double t, const int& n, double Clambda);
double fs(double t, int i, const int &n);
double fr(double t, int j, const int &n);
double fg(double t, const int& n);

// 全局随机数生成器 (使用固定种子以保证可复现性)
std::mt19937 gen(123);

// [[Rcpp::export]]
List SimSetC(int n, double Clambda, double shift, Rcpp::NumericVector Zij) {
  
  // --- Input Validation ---
  if (Rf_isNull(Zij.attr("dim"))) {
    throw std::runtime_error("'Zij' does not have a 'dim' attribute.");
  }
  Rcpp::Dimension d = Zij.attr("dim");
  if (d.size() != 3) {
    throw std::runtime_error("'Zij' must have 3 dimensions.");
  }
  
  std::size_t p = d[2];
  
  if (d[0] != n || d[1] != n)
    Rcpp::stop("Dimensions of 'Zij' do not match 'n'.");
  
  // --- Simulation via Thinning Algorithm ---
  double maxit(1.0), temp(0.0), tp(0.0), rej(0.0);
  std::vector<int> se;
  std::vector<int> re;
  std::vector<double> te;
  
  int kp(0);
  std::poisson_distribution<int> Pdis(maxit);
  std::uniform_real_distribution<double> Udis(0.0,1.0);
  
  // *** MODIFIED: Added flag for one-time warning ***
  bool negative_intensity_detected = false;
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++ ) {
      if (i == j) continue; // Skip self-interactions
      
      kp = Pdis(gen);
      
      if (kp > 0) {
        for (int z = 0; z < kp; z++) {
          tp = Udis(gen); // Propose an event time
          
          // Calculate total intensity at time tp
          temp = bs(tp, n, Clambda) + fs(tp, i, n) + fr(tp, j, n);
          for (int di = 0; di < p; di++) {
            temp += fg(tp, n) * Zij[i-1 + n*(j-1) + n*n*di];
          }
          
          // *** MODIFIED: Add check and warning for negative intensity ***
          if (temp < 0) {
            if (!negative_intensity_detected) {
              char buffer[256];
              snprintf(buffer, sizeof(buffer),
                       "Warning: Calculated intensity is negative: %.4e for pair (i=%d, j=%d) at t=%.4f. "
                       "This can happen if Clambda is too small or covariate effects are strongly negative. "
                       "Simulation proceeds by setting negative values to 0. (This warning appears only once).",
                       temp, i, j, tp);
              Rcpp::warning(buffer);
              negative_intensity_detected = true;
            }
          }
          
          // Ensure intensity is non-negative (as a safeguard)
          temp = std::max(0.0, temp);
          
          // Acceptance-rejection step
          if (temp > maxit) {
            Rcpp::stop("Calculated intensity exceeds maxit. This should not happen in this setting.");
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


// --- Setting: Half-Positive/Half-Negative, All terms O(log(n)/n^2) ---

// [[Rcpp::export]]
double bs(double t, const int &n, double Clambda) {
  double n_val = static_cast<double>(n);
  return Clambda * std::log(n_val) / (n_val * n_val);
}

// [[Rcpp::export]]
double fs(double t, int i, const int &n) {
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = std::log(n_val) / (n_val * n_val);
  double time_component = 0.5 + 0.5 * sin(2.0 * PI * t);
  double base_effect = scale_factor * time_component;
  
  if (i <= n_pos) {
    return base_effect;
  } else {
    return -base_effect * (n_pos_d / n_neg_d);
  }
}

// [[Rcpp::export]]
double fr(double t, int j, const int &n) {
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = std::log(n_val) / (n_val * n_val);
  double time_component = 0.5 + 0.5 * cos(2.0 * PI * t);
  double base_effect = scale_factor * time_component;
  
  if (j <= n_pos) {
    return base_effect;
  } else {
    return -base_effect * (n_pos_d / n_neg_d);
  }
}

// [[Rcpp::export]]
// *** MODIFIED: Changed fg to a sin form ***
double fg(double t, const int& n) {
  double n_val = static_cast<double>(n);
  double scale_factor = std::log(n_val) / (n_val * n_val);
  // Covariate effect gamma: O(log(n)/n^2) with time-varying component
  return 5.0 * scale_factor * sin(2.0 * PI * t);
}

// --- Corresponding Integral Functions ---

// [[Rcpp::export]]
NumericVector bsI_cpp(NumericVector t, int n, double Clambda) {
  double const_val = bs(0.0, n, Clambda);
  return t * const_val;
}

// [[Rcpp::export]]
NumericVector fsI_cpp(NumericVector t, int i, int n) {
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = std::log(n_val) / (n_val * n_val);
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
  double n_val = static_cast<double>(n);
  int n_pos = static_cast<int>(std::floor(n_val / 2.0));
  double n_pos_d = static_cast<double>(n_pos);
  double n_neg_d = n_val - n_pos_d;
  
  double scale_factor = std::log(n_val) / (n_val * n_val);
  NumericVector time_integral = 0.5 * t + Rcpp::sin(2.0 * PI * t) / (4.0 * PI);
  NumericVector base_integral = scale_factor * time_integral;
  
  if (j <= n_pos) {
    return base_integral;
  } else {
    return -base_integral * (n_pos_d / n_neg_d);
  }
}

// [[Rcpp::export]]
// *** MODIFIED: Updated fgI_cpp to be the integral of the new fg ***
NumericVector fgI_cpp(NumericVector t, int n) {
  double n_val = static_cast<double>(n);
  double scale_factor = std::log(n_val) / (n_val * n_val);
  double coefficient = 5.0 * scale_factor;
  
  // The function to integrate is f(t) = C * sin(2*pi*t)
  // The integral F(t) = C * [-cos(2*pi*t) / (2*pi)]
  // Integral from 0 to t is F(t) - F(0)
  // F(t) - F(0) = C * [(-cos(2*pi*t)/2pi) - (-cos(0)/2pi)]
  //             = C * [1 - cos(2*pi*t)] / (2*pi)
  
  return coefficient * (1.0 - Rcpp::cos(2.0 * PI * t)) / (2.0 * PI);
}
