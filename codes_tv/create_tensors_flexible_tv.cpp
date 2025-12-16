#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <map>
#include <utility>
#include <tuple>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List create_event_tensors_sparse_cpp(
    const Rcpp::IntegerVector& senders, 
    const Rcpp::IntegerVector& receivers, 
    const Rcpp::NumericVector& event_times, 
    int n, 
    int K,
    double bandwidth,
    bool create_kernel = false) {
  
  std::map<std::tuple<int, int, int>, double> aggregator_truncate;
  std::map<std::tuple<int, int, int>, double> aggregator_kernel;
  
  int n_events = event_times.size();
  
  double norm_const = 0.0;
  std::vector<double> time_points;
  double prefactor = 0.0;
  
  if (create_kernel) {
    norm_const = 1.0 / (bandwidth * std::sqrt(2.0 * M_PI));
    time_points.resize(K);
    for (int k = 0; k < K; ++k) {
      time_points[k] = (static_cast<double>(k) + 0.5) / static_cast<double>(K);
    }
    prefactor = -0.5 / (bandwidth * bandwidth);
  }
  
  for (int i = 0; i < n_events; ++i) {
    int sender_idx = senders[i] - 1;
    int receiver_idx = receivers[i] - 1;
    double t_event = event_times[i];
    
    if (sender_idx < 0 || sender_idx >= n || receiver_idx < 0 || receiver_idx >= n) {
      continue;
    }
    
    // ************************** BUG FIX STARTS HERE **************************
    // Part 1: 更新截断法聚合器 (使用与基准代码一致的 ceil 逻辑)
    int k_idx_truncate = static_cast<int>(std::ceil(t_event * K));
    if (k_idx_truncate > K) k_idx_truncate = K;
    if (k_idx_truncate < 1) k_idx_truncate = 1;
    k_idx_truncate -= 1; // 将 [1, K] 范围映射到 [0, K-1] 索引
    // ************************** BUG FIX ENDS HERE ****************************
    
    aggregator_truncate[{k_idx_truncate, sender_idx, receiver_idx}]++;
    
    if (create_kernel) {
      for (int k_idx_kernel = 0; k_idx_kernel < K; ++k_idx_kernel) {
        double time_diff = t_event - time_points[k_idx_kernel];
        double weight = norm_const * std::exp(prefactor * time_diff * time_diff);
        if (weight > 1e-9) {
          aggregator_kernel[{k_idx_kernel, sender_idx, receiver_idx}] += weight;
        }
      }
    }
  }
  
  auto map_to_df = [](const std::map<std::tuple<int, int, int>, double>& aggregator) {
    std::vector<int> k_col, sender_col, receiver_col;
    std::vector<double> count_col;
    
    for (const auto& pair : aggregator) {
      k_col.push_back(std::get<0>(pair.first));
      sender_col.push_back(std::get<1>(pair.first) + 1);
      receiver_col.push_back(std::get<2>(pair.first) + 1);
      count_col.push_back(pair.second);
    }
    
    return Rcpp::DataFrame::create(
      Rcpp::Named("k") = k_col,
      Rcpp::Named("sender") = sender_col,
      Rcpp::Named("receiver") = receiver_col,
      Rcpp::Named("count") = count_col
    );
  };
  
  Rcpp::List final_list;
  final_list["truncate"] = map_to_df(aggregator_truncate);
  
  if (create_kernel) {
    final_list["kernel"] = map_to_df(aggregator_kernel);
  } else {
    final_list["kernel"] = R_NilValue;
  }
  
  return final_list;
}
