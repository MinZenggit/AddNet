#include <Rcpp.h>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
Rcpp::List create_event_tensors_flexible_cpp(
    const Rcpp::IntegerVector& senders, 
    const Rcpp::IntegerVector& receivers, 
    const Rcpp::NumericVector& event_times, 
    int n, 
    int K,
    double bandwidth,
    bool create_kernel = false) { // New control parameter, defaults to false
  
  // Truncate tensor is always created
  Rcpp::IntegerVector counts_truncate(n * n * K, 0);
  
  // Kernel tensor is created only if requested
  Rcpp::NumericVector counts_kernel;
  if (create_kernel) {
    counts_kernel = Rcpp::NumericVector(n * n * K, 0.0);
  }
  
  int n_events = event_times.size();
  
  // Pre-calculate for kernel method only if needed
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
  
  // --- Single loop over all events ---
  for (int i = 0; i < n_events; ++i) {
    int sender_idx = senders[i] - 1;
    int receiver_idx = receivers[i] - 1;
    double t_event = event_times[i];
    
    if (sender_idx < 0 || sender_idx >= n || receiver_idx < 0 || receiver_idx >= n) {
      continue;
    }
    
    // --- Part 1: Update Truncate Tensor (Always runs) ---
    int k_idx_truncate = static_cast<int>(std::ceil(t_event * K));
    if (k_idx_truncate > K) k_idx_truncate = K;
    if (k_idx_truncate < 1) k_idx_truncate = 1;
    k_idx_truncate -= 1;
    
    int flat_index_truncate = sender_idx + receiver_idx * n + k_idx_truncate * n * n;
    counts_truncate[flat_index_truncate]++;
    
    // --- Part 2: Update Kernel Tensor (Only if requested) ---
    if (create_kernel) {
      for (int k_idx_kernel = 0; k_idx_kernel < K; ++k_idx_kernel) {
        double time_diff = t_event - time_points[k_idx_kernel];
        double weight = norm_const * std::exp(prefactor * time_diff * time_diff);
        int flat_index_kernel = sender_idx + receiver_idx * n + k_idx_kernel * n * n;
        counts_kernel[flat_index_kernel] += weight;
      }
    }
  }
  
  counts_truncate.attr("dim") = Rcpp::Dimension(n, n, K);
  
  if (create_kernel) {
    counts_kernel.attr("dim") = Rcpp::Dimension(n, n, K);
    return Rcpp::List::create(
      Rcpp::Named("truncate") = counts_truncate,
      Rcpp::Named("kernel") = counts_kernel
    );
  } else {
    // If kernel was not created, return a list with a NULL kernel element
    return Rcpp::List::create(
      Rcpp::Named("truncate") = counts_truncate,
      Rcpp::Named("kernel") = R_NilValue // This becomes NULL in R
    );
  }
}
