#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <map>
#include <utility>

// [[Rcpp::depends(RcppArmadillo)]]

// Helper to convert a map to a List for R
Rcpp::List map_to_list(const std::map<std::pair<int, int>, double>& event_map) {
  std::vector<int> senders, receivers;
  std::vector<double> counts;
  
  for (const auto& pair : event_map) {
    senders.push_back(pair.first.first + 1); // Convert back to 1-based index for R
    receivers.push_back(pair.first.second + 1);
    counts.push_back(pair.second);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("sender") = Rcpp::wrap(senders),
    Rcpp::Named("receiver") = Rcpp::wrap(receivers),
    Rcpp::Named("count") = Rcpp::wrap(counts)
  );
}


// [[Rcpp::export]]
Rcpp::List create_event_tensors_flexible_cpp(
    const Rcpp::IntegerVector& senders, 
    const Rcpp::IntegerVector& receivers, 
    const Rcpp::NumericVector& event_times, 
    int n, 
    int K,
    double bandwidth,
    bool create_kernel = false) {
  
  // Use maps to aggregate counts for sparse events
  std::vector<std::map<std::pair<int, int>, double>> aggregator_truncate(K);
  std::vector<std::map<std::pair<int, int>, double>> aggregator_kernel;
  
  if (create_kernel) {
    aggregator_kernel.resize(K);
  }
  
  int n_events = event_times.size();
  
  // Pre-calculate for kernel method
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
    int sender_idx = senders[i] - 1; // Use 0-based index internally
    int receiver_idx = receivers[i] - 1;
    double t_event = event_times[i];
    
    if (sender_idx < 0 || sender_idx >= n || receiver_idx < 0 || receiver_idx >= n) {
      continue;
    }
    
    // Part 1: Update Truncate Tensor
    int k_idx_truncate = static_cast<int>(std::floor(t_event * K));
    if (k_idx_truncate >= K) k_idx_truncate = K - 1;
    if (k_idx_truncate < 0) k_idx_truncate = 0;
    
    aggregator_truncate[k_idx_truncate][{sender_idx, receiver_idx}]++;
    
    // Part 2: Update Kernel Tensor
    if (create_kernel) {
      for (int k_idx_kernel = 0; k_idx_kernel < K; ++k_idx_kernel) {
        double time_diff = t_event - time_points[k_idx_kernel];
        double weight = norm_const * std::exp(prefactor * time_diff * time_diff);
        if (weight > 1e-9) { // Optimization: only add if weight is non-trivial
          aggregator_kernel[k_idx_kernel][{sender_idx, receiver_idx}] += weight;
        }
      }
    }
  }
  
  // Convert maps to lists of lists for R
  Rcpp::List final_truncate(K);
  for(int k = 0; k < K; ++k) {
    final_truncate[k] = map_to_list(aggregator_truncate[k]);
  }
  
  Rcpp::List final_list;
  final_list["truncate"] = final_truncate;
  
  if (create_kernel) {
    Rcpp::List final_kernel(K);
    for(int k = 0; k < K; ++k) {
      final_kernel[k] = map_to_list(aggregator_kernel[k]);
    }
    final_list["kernel"] = final_kernel;
  } else {
    final_list["kernel"] = R_NilValue;
  }
  
  return final_list;
}


// [[Rcpp::export]]
Rcpp::List create_z_matrices_cpp(
    const Rcpp::IntegerVector& senders,
    const Rcpp::IntegerVector& receivers,
    const arma::mat& z_values,
    int n,
    int p) {
  
  if (p == 0) {
    return Rcpp::List::create(
      Rcpp::Named("inv_XtX") = arma::mat(0, 0),
      Rcpp::Named("diff2_z_core") = arma::mat(0, (n-1)*(n-1))
    );
  }
  
  arma::mat diff2_z_mat((n - 1) * (n - 1), p, arma::fill::zeros);
  
  for (int i = 0; i < senders.size(); ++i) {
    int u = senders[i] - 1;
    int v = receivers[i] - 1;
    arma::rowvec z_row = z_values.row(i);
    
    // Contribution to diff2_z(u, v)
    if (u < n - 1 && v < n - 1) {
      diff2_z_mat.row(u + v * (n - 1)) += z_row;
    }
    // Contribution to diff2_z(u+1, v)
    if (u < n - 1 && v < n - 1 && u + 1 < n - 1) {
      diff2_z_mat.row((u + 1) + v * (n - 1)) -= z_row;
    }
    // Contribution to diff2_z(u, v+1)
    if (u < n - 1 && v < n - 1 && v + 1 < n - 1) {
      diff2_z_mat.row(u + (v + 1) * (n - 1)) -= z_row;
    }
    // Contribution to diff2_z(u+1, v+1)
    if (u < n - 1 && v < n - 1 && u + 1 < n - 1 && v + 1 < n - 1) {
      diff2_z_mat.row((u + 1) + (v + 1) * (n - 1)) += z_row;
    }
  }
  
  arma::mat XtX = diff2_z_mat.t() * diff2_z_mat;
  if (arma::rcond(XtX) < std::numeric_limits<double>::epsilon()) {
    Rcpp::stop("Covariate matrix is singular.");
  }
  arma::mat inv_XtX = arma::inv_sympd(XtX);
  arma::mat diff2_z_core = inv_XtX * diff2_z_mat.t();
  
  return Rcpp::List::create(
    Rcpp::Named("inv_XtX") = inv_XtX,
    Rcpp::Named("diff2_z_core") = diff2_z_core
  );
}
