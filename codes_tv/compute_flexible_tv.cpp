// compute_flexible_tv.cpp
#include <RcppArmadillo.h>
#include <set>
#include <map>
#include <tuple>
#include <algorithm> // For std::upper_bound

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 辅助结构体 (无变化)
struct ParamSet {
  arma::mat gamma_hat, alpha_hat, beta_hat, alpha_var, beta_var;
  arma::cube gamma_var;
  arma::vec lambda0_hat, lambda0_var;
  
  ParamSet(int n, int p, int K) :
    gamma_hat(p, K, fill::zeros), alpha_hat(n, K, fill::zeros),
    beta_hat(n, K, fill::zeros), alpha_var(n, K, fill::zeros),
    beta_var(n, K, fill::zeros), gamma_var(p, p, K, fill::zeros),
    lambda0_hat(K, fill::zeros), lambda0_var(K, fill::zeros) {}
  
  List to_R_list() {
    return List::create(
      Named("gamma_hat") = gamma_hat, Named("gamma_var") = gamma_var,
      Named("alpha_hat") = alpha_hat, Named("alpha_var") = alpha_var,
      Named("beta_hat") = beta_hat,   Named("beta_var") = beta_var,
      Named("lambda0_hat") = lambda0_hat, Named("lambda0_var") = lambda0_var
    );
  }
};

// 核心计算函数 (无变化)
void process_k_slice_sparse(
    int k, int n, int p, double delta_t,
    const std::map<std::pair<int, int>, double>& events_k,
    const std::map<std::pair<int, int>, arma::vec>& z_map,
    const arma::mat& diff2_z_core, const arma::mat& inv_XtX,
    ParamSet& params
) {
  // --- 这部分代码与原始文件完全相同 ---
  // --- 点估计部分 ---
  arma::vec diff2_response_vec(((n-1)*(n-1)), fill::zeros);
  if (p > 0) {
    for(auto const& [coords, count] : events_k) {
      int r = coords.first; int c = coords.second;
      if (r > 0 && c > 0)   diff2_response_vec[(r-1) + (c-1)*(n-1)] += count;
      if (r < n-1 && c > 0) diff2_response_vec[r + (c-1)*(n-1)] -= count;
      if (r > 0 && c < n-1) diff2_response_vec[(r-1) + c*(n-1)] -= count;
      if (r < n-1 && c < n-1) diff2_response_vec[r + c*(n-1)] += count;
    }
    params.gamma_hat.col(k) = (diff2_z_core * diff2_response_vec) / delta_t;
  }
  arma::vec row_sums(n, fill::zeros), col_sums(n, fill::zeros);
  double total_events = 0;
  for(auto const& [coords, count] : events_k) {
    row_sums(coords.first) += count; col_sums(coords.second) += count; total_events += count;
  }
  arma::vec diff_alpha(n-1, fill::zeros); arma::vec diff_beta(n-1, fill::zeros);
  if (p > 0) {
    arma::mat z_diff_alpha_gamma(n-1, n, fill::zeros); arma::mat z_diff_beta_gamma(n, n-1, fill::zeros);
    for(auto const& [coords, z_vec] : z_map) {
      int r = coords.first; int c = coords.second;
      double z_gamma_val = dot(z_vec, params.gamma_hat.col(k));
      if (r > 0)   z_diff_alpha_gamma(r-1, c) += z_gamma_val;
      if (r < n-1) z_diff_alpha_gamma(r, c)   -= z_gamma_val;
      if (c > 0)   z_diff_beta_gamma(r, c-1)  += z_gamma_val;
      if (c < n-1) z_diff_beta_gamma(r, c)    -= z_gamma_val;
    }
    diff_alpha = (row_sums.subvec(1, n-1) - row_sums.subvec(0, n-2)) / (n * delta_t) - mean(z_diff_alpha_gamma, 1);
    diff_beta = (col_sums.subvec(1, n-1) - col_sums.subvec(0, n-2)) / (n * delta_t) - mean(z_diff_beta_gamma, 0).t();
  } else {
    diff_alpha = (row_sums.subvec(1, n-1) - row_sums.subvec(0, n-2)) / (n * delta_t);
    diff_beta = (col_sums.subvec(1, n-1) - col_sums.subvec(0, n-2)) / (n * delta_t);
  }
  arma::vec raw_alpha_hat = cumsum(join_cols(zeros<vec>(1), diff_alpha));
  params.alpha_hat.col(k) = raw_alpha_hat - mean(raw_alpha_hat);
  arma::vec raw_beta_hat = cumsum(join_cols(zeros<vec>(1), diff_beta));
  params.beta_hat.col(k) = raw_beta_hat - mean(raw_beta_hat);
  double Z_gamma_sum = 0;
  if (p > 0) { for(auto const& [coords, z_vec] : z_map) { Z_gamma_sum += dot(z_vec, params.gamma_hat.col(k)); } }
  double total_residual = total_events / delta_t - (n * accu(params.alpha_hat.col(k)) + n * accu(params.beta_hat.col(k)) + Z_gamma_sum);
  params.lambda0_hat(k) = total_residual / (n * n);
  
  // --- 方差计算部分 (无变化) ---
  arma::mat lambda_hat_k(n, n);
  lambda_hat_k.fill(params.lambda0_hat(k));
  lambda_hat_k.each_col() += params.alpha_hat.col(k);
  lambda_hat_k.each_row() += params.beta_hat.col(k).t();
  if (p > 0) {
    for(auto const& [coords, z_vec] : z_map) {
      lambda_hat_k(coords.first, coords.second) += dot(z_vec, params.gamma_hat.col(k));
    }
  }
  params.lambda0_var(k) = accu(lambda_hat_k) / (pow(n, 4.0) * delta_t);
  arma::vec row_sums_lambda = sum(lambda_hat_k, 1);
  arma::vec col_sums_lambda = sum(lambda_hat_k, 0).t();
  double common_denom = n * n * delta_t;
  arma::vec var_diff_alpha(n-1); for(int l=0; l < n-1; ++l) var_diff_alpha(l) = (row_sums_lambda(l) + row_sums_lambda(l+1)) / common_denom;
  arma::vec cov_diff_alpha(n-2); for(int l=0; l < n-2; ++l) cov_diff_alpha(l) = -row_sums_lambda(l+1) / common_denom;
  for(int i=0; i < n; ++i) { double v=0; for(int l=0; l<n-1; ++l) { double c=(l+1.0)/n-(l>=i?1:0); v+=c*c*var_diff_alpha(l); } for(int l=0; l<n-2; ++l) { double c1=(l+1.0)/n-(l>=i?1:0); double c2=(l+2.0)/n-((l+1)>=i?1:0); v+=2*c1*c2*cov_diff_alpha(l); } params.alpha_var(i,k)=v; }
  arma::vec var_diff_beta(n-1); for(int l=0; l < n-1; ++l) var_diff_beta(l) = (col_sums_lambda(l) + col_sums_lambda(l+1)) / common_denom;
  arma::vec cov_diff_beta(n-2); for(int l=0; l < n-2; ++l) cov_diff_beta(l) = -col_sums_lambda(l+1) / common_denom;
  for(int j=0; j < n; ++j) { double v=0; for(int l=0; l<n-1; ++l) { double c=(l+1.0)/n-(l>=j?1:0); v+=c*c*var_diff_beta(l); } for(int l=0; l<n-2; ++l) { double c1=(l+1.0)/n-(l>=j?1:0); double c2=(l+2.0)/n-((l+1)>=j?1:0); v+=2*c1*c2*cov_diff_beta(l); } params.beta_var(j,k)=v; }
  
  if (p > 0) {
    arma::mat Omega_k(p, p, fill::zeros);
    auto get_z = [&](int r, int c) -> arma::vec {
      if (r < 0 || r >= n || c < 0 || c >= n) return arma::zeros<vec>(p);
      auto it = z_map.find({r, c});
      return (it != z_map.end()) ? it->second : arma::zeros<vec>(p);
    };
    auto get_diff2_z = [&](int r, int c) -> arma::vec {
      if (r < 0 || r >= n-1 || c < 0 || c >= n-1) return arma::zeros<vec>(p);
      return get_z(r+1, c+1) - get_z(r, c+1) - get_z(r+1, c) + get_z(r, c);
    };
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        arma::vec x_ab     = get_diff2_z(a, b);
        arma::vec x_am1b   = get_diff2_z(a - 1, b);
        arma::vec x_abm1   = get_diff2_z(a, b - 1);
        arma::vec x_am1bm1 = get_diff2_z(a - 1, b - 1);
        arma::vec tilde_z_vec = x_ab - x_am1b - x_abm1 + x_am1bm1;
        Omega_k += tilde_z_vec * tilde_z_vec.t() * lambda_hat_k(a, b);
      }
    }
    params.gamma_var.slice(k) = (1.0 / delta_t) * inv_XtX * Omega_k * inv_XtX;
  }
}

// 主导出函数 (已修改以处理时变协变量)
// [[Rcpp::export]]
List compute_results_sparse_cpp_tv(
    int n, int p, int K,
    Rcpp::DataFrame events_df,
    Rcpp::Nullable<Rcpp::List> z_list_nullable,
    Rcpp::Nullable<Rcpp::NumericVector> z_times_nullable
) {
  // --- 1. 事件数据预处理 (无变化) ---
  auto group_events = [&](Rcpp::DataFrame df) {
    std::vector<std::map<std::pair<int, int>, double>> grouped(K);
    IntegerVector event_k = df["k"]; IntegerVector event_s = df["sender"];
    IntegerVector event_r = df["receiver"]; NumericVector event_c = df["count"];
    for(int i=0; i<event_k.size(); ++i) {
      grouped[event_k[i]][{event_s[i]-1, event_r[i]-1}] = event_c[i];
    }
    return grouped;
  };
  auto grouped_truncate = group_events(events_df);
  
  // --- 2. 初始化结果和协变量相关变量 (已修改) ---
  ParamSet params_truncate(n, p, K);
  double delta_t_truncate = 1.0 / static_cast<double>(K);
  
  // 协变量相关对象
  std::map<std::pair<int, int>, arma::vec> z_map;
  arma::mat diff2_z_core, inv_XtX;
  Rcpp::List z_dfs;
  Rcpp::NumericVector z_times;
  
  if (p > 0) {
    if (!z_list_nullable.isNotNull() || !z_times_nullable.isNotNull()) {
      stop("z_list and z_times must be provided when p > 0.");
    }
    z_dfs = Rcpp::as<Rcpp::List>(z_list_nullable.get());
    z_times = Rcpp::as<Rcpp::NumericVector>(z_times_nullable.get());
  }
  
  int active_z_idx = -1; // 跟踪当前活动的协变量集索引
  
  // --- 3. 主循环 (已修改) ---
  for (int k = 0; k < K; ++k) {
    if (p > 0) {
      // 确定当前时间片k应该使用哪个协变量集
      double current_t = (static_cast<double>(k) + 0.5) / K;
      
      // std::upper_bound 找到第一个大于 current_t 的时间点
      // 它返回一个迭代器，distance算出这个迭代器是第几个元素 (0-based)
      // 这正好是我们需要的协变量集的索引
      auto it = std::upper_bound(z_times.begin(), z_times.end(), current_t);
      int required_z_idx = std::distance(z_times.begin(), it);
      
      // **关键优化**: 仅当协变量集发生变化时才更新矩阵
      if (required_z_idx != active_z_idx) {
        active_z_idx = required_z_idx;
        
        // 清空并重新填充 z_map
        z_map.clear();
        Rcpp::DataFrame current_z_df = Rcpp::as<Rcpp::DataFrame>(z_dfs[active_z_idx]);
        IntegerVector z_senders = current_z_df["sender"];
        IntegerVector z_receivers = current_z_df["receiver"];
        arma::mat z_values = as<arma::mat>(current_z_df["values"]);
        for(int i=0; i<z_senders.size(); ++i) {
          z_map[{z_senders[i]-1, z_receivers[i]-1}] = z_values.row(i).t();
        }
        
        // 重新计算依赖于 z_map 的矩阵
        arma::mat diff2_z_mat((n-1)*(n-1), p, fill::zeros);
        for(auto const& [coords, z_vec] : z_map) {
          int r = coords.first; int c = coords.second;
          if (r > 0 && c > 0)   diff2_z_mat.row((r-1) + (c-1)*(n-1)) += z_vec.t();
          if (r < n-1 && c > 0) diff2_z_mat.row(r + (c-1)*(n-1)) -= z_vec.t();
          if (r > 0 && c < n-1) diff2_z_mat.row((r-1) + c*(n-1)) -= z_vec.t();
          if (r < n-1 && c < n-1) diff2_z_mat.row(r + c*(n-1)) += z_vec.t();
        }
        arma::mat XtX = diff2_z_mat.t() * diff2_z_mat;
        if (rcond(XtX) < 1e-12) {
          std::string msg = "Covariate matrix is singular for time interval " + std::to_string(active_z_idx + 1);
          stop(msg);
        }
        inv_XtX = arma::inv_sympd(XtX);
        diff2_z_core = inv_XtX * diff2_z_mat.t();
      }
    }
    
    // 使用当前（可能已更新）的协变量矩阵调用核心计算函数
    process_k_slice_sparse(k, n, p, delta_t_truncate, grouped_truncate[k], z_map, diff2_z_core, inv_XtX, params_truncate);
  }
  
  // --- 4. 返回结果 (已简化) ---
  return params_truncate.to_R_list();
}
