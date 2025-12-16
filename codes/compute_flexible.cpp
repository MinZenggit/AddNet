#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Helper struct to hold all parameter estimates for one method
struct ParamSet {
  arma::mat gamma_hat;
  arma::cube gamma_var;
  arma::mat alpha_hat;
  arma::mat beta_hat;
  arma::vec lambda0_hat;
  arma::mat alpha_var;
  arma::mat beta_var;
  arma::vec lambda0_var;
  
  ParamSet(int n, int p, int K) :
    gamma_hat(p, K, fill::zeros),
    gamma_var(p, p, K, fill::zeros),
    alpha_hat(n, K, fill::zeros),
    beta_hat(n, K, fill::zeros),
    lambda0_hat(K, fill::zeros),
    alpha_var(n, K, fill::zeros),
    beta_var(n, K, fill::zeros),
    lambda0_var(K, fill::zeros) {}
  
  List to_R_list() {
    return List::create(
      Named("gamma_hat") = gamma_hat, Named("gamma_var") = gamma_var,
      Named("alpha_hat") = alpha_hat, Named("alpha_var") = alpha_var,
      Named("beta_hat") = beta_hat,   Named("beta_var") = beta_var,
      Named("lambda0_hat") = lambda0_hat, Named("lambda0_var") = lambda0_var
    );
  }
};

// Helper function to process one slice (k) for one method
void process_k_slice(
    int k, int n, int p, double delta_t,
    const arma::mat& current_events,
    const arma::cube& z_constant,
    const arma::mat& diff2_z_core,
    const arma::mat& inv_XtX,
    const arma::cube& diff2_z_cube,
    const arma::cube& z_diff_alpha,
    const arma::cube& z_diff_beta,
    ParamSet& params // Pass by reference to modify
) {
  // --- Step 1: Estimate gamma_hat ---
  arma::mat diff2_response = (current_events.submat(1, 1, n-1, n-1) - 
    current_events.submat(0, 1, n-2, n-1) - 
    current_events.submat(1, 0, n-1, n-2) + 
    current_events.submat(0, 0, n-2, n-2));
  if (p > 0) {
    params.gamma_hat.col(k) = (diff2_z_core * vectorise(diff2_response)) / delta_t;
  }
  
  // --- Step 2: Estimate alpha_hat and beta_hat ---
  // A. Alpha
  arma::mat event_diff_alpha = (current_events.rows(1, n-1) - current_events.rows(0, n-2)) / delta_t;
  arma::mat correction_alpha(n-1, n, fill::zeros);
  if (p > 0) {
    for (int i = 0; i < n-1; ++i) {
      for (int j = 0; j < n; ++j) {
        correction_alpha(i, j) = dot(arma::vec(z_diff_alpha.tube(i, j)), params.gamma_hat.col(k));
      }
    }
  }
  arma::vec diff_alpha = mean(event_diff_alpha - correction_alpha, 1);
  arma::vec raw_alpha_hat = cumsum(join_cols(zeros<vec>(1), diff_alpha));
  params.alpha_hat.col(k) = raw_alpha_hat - mean(raw_alpha_hat);
  
  // B. Beta
  arma::mat event_diff_beta = (current_events.cols(1, n-1) - current_events.cols(0, n-2)) / delta_t;
  arma::mat correction_beta(n, n-1, fill::zeros);
  if (p > 0) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n-1; ++j) {
        correction_beta(i, j) = dot(arma::vec(z_diff_beta.tube(i, j)), params.gamma_hat.col(k));
      }
    }
  }
  arma::vec diff_beta = mean(event_diff_beta - correction_beta, 0).t();
  arma::vec raw_beta_hat = cumsum(join_cols(zeros<vec>(1), diff_beta));
  params.beta_hat.col(k) = raw_beta_hat - mean(raw_beta_hat);
  
  // --- Step 3 & 4: Estimate lambda0_hat and full intensity lambda_hat_k ---
  arma::mat Z_gamma_k(n, n, fill::zeros);
  if (p > 0) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        Z_gamma_k(i, j) = dot(arma::vec(z_constant.tube(i, j)), params.gamma_hat.col(k));
      }
    }
  }
  arma::mat lambda_hat_k = Z_gamma_k;
  lambda_hat_k.each_col() += params.alpha_hat.col(k);
  lambda_hat_k.each_row() += params.beta_hat.col(k).t();
  
  double total_residual = accu(current_events) / delta_t - accu(lambda_hat_k);
  params.lambda0_hat(k) = total_residual / (n * n);
  lambda_hat_k.fill(params.lambda0_hat(k));
  lambda_hat_k += Z_gamma_k;
  lambda_hat_k.each_col() += params.alpha_hat.col(k);
  lambda_hat_k.each_row() += params.beta_hat.col(k).t();
  
  // --- Step 5: Estimate Variances ---
  // A. Lambda0 Var
  params.lambda0_var(k) = accu(lambda_hat_k) / (pow(n, 4.0) * delta_t);
  
  // B. Alpha Var
  arma::vec row_sums = sum(lambda_hat_k, 1);
  arma::vec var_diff_alpha(n-1);
  arma::vec cov_diff_alpha(n-2);
  double common_denom = n * n * delta_t;
  for(int l=0; l < n-1; ++l) var_diff_alpha(l) = (row_sums(l) + row_sums(l+1)) / common_denom;
  for(int l=0; l < n-2; ++l) cov_diff_alpha(l) = -row_sums(l+1) / common_denom;
  for(int i=0; i < n; ++i) {
    double current_var = 0.0;
    for(int l=0; l < n-1; ++l) {
      double c_il = (double)(l + 1.0) / n - (l >= i ? 1.0 : 0.0);
      current_var += c_il * c_il * var_diff_alpha(l);
    }
    for(int l=0; l < n-2; ++l) {
      double c_il = (double)(l + 1.0) / n - (l >= i ? 1.0 : 0.0);
      double c_il_plus_1 = (double)(l + 2.0) / n - ((l+1) >= i ? 1.0 : 0.0);
      current_var += 2.0 * c_il * c_il_plus_1 * cov_diff_alpha(l);
    }
    params.alpha_var(i, k) = current_var;
  }
  
  // C. Beta Var
  arma::vec col_sums = sum(lambda_hat_k, 0).t();
  arma::vec var_diff_beta(n-1);
  arma::vec cov_diff_beta(n-2);
  for(int l=0; l < n-1; ++l) var_diff_beta(l) = (col_sums(l) + col_sums(l+1)) / common_denom;
  for(int l=0; l < n-2; ++l) cov_diff_beta(l) = -col_sums(l+1) / common_denom;
  for(int j=0; j < n; ++j) {
    double current_var = 0.0;
    for(int l=0; l < n-1; ++l) {
      double c_jl = (double)(l + 1.0) / n - (l >= j ? 1.0 : 0.0);
      current_var += c_jl * c_jl * var_diff_beta(l);
    }
    for(int l=0; l < n-2; ++l) {
      double c_jl = (double)(l + 1.0) / n - (l >= j ? 1.0 : 0.0);
      double c_jl_plus_1 = (double)(l + 2.0) / n - ((l+1) >= j ? 1.0 : 0.0);
      current_var += 2.0 * c_jl * c_jl_plus_1 * cov_diff_beta(l);
    }
    params.beta_var(j, k) = current_var;
  }
  
  // D. Gamma Var
  if (p > 0) {
    arma::mat Omega_k(p, p, fill::zeros);
    arma::vec z_zeros = zeros<vec>(p);
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        arma::vec v_ab     = (a < n-1 && b < n-1) ? arma::vec(diff2_z_cube.tube(a, b))     : z_zeros;
        arma::vec v_am1b   = (a > 0   && b < n-1) ? arma::vec(diff2_z_cube.tube(a-1, b))   : z_zeros;
        arma::vec v_abm1   = (a < n-1 && b > 0)   ? arma::vec(diff2_z_cube.tube(a, b-1))   : z_zeros;
        arma::vec v_am1bm1 = (a > 0   && b > 0)   ? arma::vec(diff2_z_cube.tube(a-1, b-1)) : z_zeros;
        arma::vec tilde_z_vec = v_ab - v_am1b - v_abm1 + v_am1bm1;
        Omega_k += tilde_z_vec * tilde_z_vec.t() * lambda_hat_k(a, b);
      }
    }
    params.gamma_var.slice(k) = (1.0 / delta_t) * inv_XtX * Omega_k * inv_XtX;
  }
}

// [[Rcpp::export]]
List compute_results_flexible_cpp(int n, int p, 
                                  const arma::cube& event_counts_truncate, 
                                  Rcpp::Nullable<Rcpp::NumericVector> event_counts_kernel_nullable,
                                  const arma::cube& z_constant, 
                                  const arma::mat& diff2_z_core,
                                  const arma::mat& inv_XtX,
                                  const arma::cube& diff2_z_cube) {
  int K = event_counts_truncate.n_slices;
  
  // Determine if we need to compute kernel results
  bool compute_kernel = event_counts_kernel_nullable.isNotNull();
  
  // Initialize result containers
  ParamSet params_truncate(n, p, K);
  ParamSet params_kernel(n, p, K); // Created but might not be filled
  
  // Precompute shared z differences
  arma::cube z_diff_alpha, z_diff_beta;
  if (p > 0) {
    z_diff_alpha = z_constant.rows(1, n-1) - z_constant.rows(0, n-2);
    z_diff_beta = z_constant.cols(1, n-1) - z_constant.cols(0, n-2);
  }
  
  // Effective delta_t for each method
  double delta_t_truncate = 1.0 / static_cast<double>(K);
  double delta_t_kernel = 1.0;
  
  // Convert kernel counts to arma::cube if it exists
  arma::cube event_counts_kernel;
  if (compute_kernel) {
    Rcpp::NumericVector kernel_vec = event_counts_kernel_nullable.get();
    event_counts_kernel = arma::cube(kernel_vec.begin(), n, n, K, false);
  }
  
  // --- Single loop over K, processing methods conditionally ---
  for (int k = 0; k < K; ++k) {
    // Process slice k for truncate method (always)
    process_k_slice(k, n, p, delta_t_truncate, event_counts_truncate.slice(k),
                    z_constant, diff2_z_core, inv_XtX, diff2_z_cube,
                    z_diff_alpha, z_diff_beta, params_truncate);
    
    // Process slice k for kernel method (conditionally)
    if (compute_kernel) {
      process_k_slice(k, n, p, delta_t_kernel, event_counts_kernel.slice(k),
                      z_constant, diff2_z_core, inv_XtX, diff2_z_cube,
                      z_diff_alpha, z_diff_beta, params_kernel);
    }
  }
  
  if (compute_kernel) {
    return List::create(
      Named("truncate") = params_truncate.to_R_list(),
      Named("kernel") = params_kernel.to_R_list()
    );
  } else {
    return List::create(
      Named("truncate") = params_truncate.to_R_list(),
      Named("kernel") = R_NilValue
    );
  }
}
