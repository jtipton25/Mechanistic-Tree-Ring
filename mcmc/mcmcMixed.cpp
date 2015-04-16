#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, myFunctions, bayesTreeRing)]]
#include "mcmcMix.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List makeMCMC(mat y, mat Temp, mat P, vec species, List params, List process,
              bool sim = false){
  
  // calibration model priors
  double mu_mu_beta_0 = as<double>(params["mu_mu_beta_0"]);
  double s2_mu_beta_0 = as<double>(params["s2_mu_beta_0"]);
  double alpha_s2_beta_0 = as<double>(params["alpha_s2_beta_0"]);
  double beta_s2_beta_0 = as<double>(params["beta_s2_beta_0"]);
  double mu_mu_beta_1 = as<double>(params["mu_mu_beta_1"]);
  double s2_mu_beta_1 = as<double>(params["s2_mu_beta_1"]);
  double alpha_s2_beta_1 = as<double>(params["alpha_s2_beta_1"]);
  double beta_s2_beta_1 = as<double>(params["beta_s2_beta_1"]);
  double alpha_s2 = as<double>(params["alpha_s2"]);
  double beta_s2 = as<double>(params["beta_s2"]);
  
  // VS growth model priors
  double alpha_Temp_min = as<double>(params["alpha_Temp_min"]);
  double beta_Temp_min = as<double>(params["beta_Temp_min"]);
  double Temp_min_lower = as<double>(params["Temp_min_lower"]);
  double Temp_min_upper = as<double>(params["Temp_min_upper"]);
  double alpha_Temp_max = as<double>(params["alpha_Temp_max"]);
  double beta_Temp_max = as<double>(params["beta_Temp_max"]);
  double Temp_max_lower = as<double>(params["Temp_max_lower"]);
  double Temp_max_upper = as<double>(params["Temp_max_upper"]);
  double alpha_P_min = as<double>(params["alpha_P_min"]);
  double beta_P_min = as<double>(params["beta_P_min"]);
  double P_min_lower = as<double>(params["P_min_lower"]);
  double P_min_upper = as<double>(params["P_min_upper"]);
  double alpha_P_max = as<double>(params["alpha_P_max"]);
  double beta_P_max = as<double>(params["beta_P_max"]);
  double P_max_lower = as<double>(params["P_max_lower"]);
  double P_max_upper = as<double>(params["P_max_upper"]);
  double Temp_min_tune_tmp = as<double>(params["Temp_min_tune"]);
  double Temp_max_tune_tmp = as<double>(params["Temp_max_tune"]);
  double P_min_tune_tmp = as<double>(params["P_min_tune"]);
  double P_max_tune_tmp = as<double>(params["P_max_tune"]);
 
  // probit growth model priors
  double mu_mu_gamma_T = as<double>(params["mu_mu_gamma_T"]);
  double s2_mu_gamma_T = as<double>(params["s2_mu_gamma_T"]);
  double alpha_s2_gamma_T = as<double>(params["alpha_s2_gamma_T"]);
  double beta_s2_gamma_T = as<double>(params["beta_s2_gamma_T"]);
  //  IG prior for xi2_T
  //  double alpha_xi2_T = as<double>(params["alpha_xi2_T"]);
  //  double beta_xi2_T = as<double>(params["beta_xi2_T"]);
  // logNormal prior for xi2_T
  double mu_mu_xi2_T = as<double>(params["mu_mu_xi2_T"]);
  double s2_mu_xi2_T = as<double>(params["s2_mu_xi2_T"]);
  double alpha_s2_xi2_T = as<double>(params["alpha_s2_xi2_T"]);
  double beta_s2_xi2_T = as<double>(params["beta_s2_xi2_T"]);
  double mu_mu_gamma_P = as<double>(params["mu_mu_gamma_P"]);
  double s2_mu_gamma_P = as<double>(params["s2_mu_gamma_P"]);
  double alpha_s2_gamma_P = as<double>(params["alpha_s2_gamma_P"]);
  double beta_s2_gamma_P = as<double>(params["beta_s2_gamma_P"]);
  //  IG prior for xi2_P
  //  double alpha_xi2_P = as<double>(params["alpha_xi2_P"]);
  //  double beta_xi2_P = as<double>(params["beta_xi2_P"]);
  // logNormal prior for xi2_P
  double mu_mu_xi2_P = as<double>(params["mu_mu_xi2_P"]);
  double s2_mu_xi2_P = as<double>(params["s2_mu_xi2_P"]);
  double alpha_s2_xi2_P = as<double>(params["alpha_s2_xi2_P"]);
  double beta_s2_xi2_P = as<double>(params["beta_s2_xi2_P"]);
  double gamma_T_tune_tmp = as<double>(params["gamma_T_tune"]);
  double xi2_T_tune_tmp = as<double>(params["xi2_T_tune"]);
  double gamma_P_tune_tmp = as<double>(params["gamma_P_tune"]);
  double xi2_P_tune_tmp = as<double>(params["xi2_P_tune"]);
 
  int n_mcmc = as<int>(params["n_mcmc"]);
  int n_thin = as<int>(params["n_thin"]);
  int num_species = as<int>(params["num_species"]);
  double psi = as<double>(params["psi"]);
 
  vec day_len = as<vec>(params["day_len"]);
  double W_T_tune = as<double>(params["W_T_tune"]);
  double W_P_tune = as<double>(params["W_P_tune"]);
  uvec T_obs_idx = as<uvec>(params["T_obs_idx"]);
  uvec P_obs_idx = as<uvec>(params["P_obs_idx"]);
 
  //
  // Load empirical Bayes output
  //
 
  vec phi_vec = as<vec>(process["phi_vec"]);
  mat Sigma_inv = as<mat>(process["Sigma_inv"]);
  double trend_0 = as<double>(process["trend_0"]);
  double trend_1 = as<double>(process["trend_1"]);
 
  //
  // Initialize variables
  //
 
  int t = y.n_rows;
  int p = y.n_cols;
  mat H(t, p);
  int N_obs = 0;
  for(int i = 0; i < p; i++){
    for(int k = 0; k < t; k++){
      if(y(k, i) == 0){
  H(k, i) = 0;
      } else {
  H(k, i) = 1;
  N_obs += 1;
      }
    }
  }
  
  vec H_col_sums(p);
  for(int i = 0; i < p; i++){
    H_col_sums(i) = sum(H.col(i));
  }
  vec species_m1 = species - 1;
  vec N_obs_species(num_species, fill::zeros);
  for(int i = 0; i < p; i++){
    int species_idx = species_m1(i);
    N_obs_species(species_idx) += H_col_sums(i); 
  }
 
  int N_T = T_obs_idx.n_elem; 
  int N_P = P_obs_idx.n_elem;
  mat one_mat(t, num_species, fill::ones);
  mat zero_mat(t, num_species, fill::zeros);
  double log_psi = log(psi);
  double log_one_minus_psi = log(1 - psi);
 
  vec J_num_species(num_species, fill::ones);
  vec J(24);
  for(int i = 0; i < 24; i++){
    if(i < 12){
      J(i) = 1;
    } else {
      J(i) = 0;
    }
  }
  vec Jtrend_0 = J * trend_0;
  vec Jtrend_1 = J * trend_1;
 
  vec mu_T = rowMeans(Temp);
  vec s_T = rowSds(Temp);
  vec mu_P = rowMeans(P);
  vec s_P = rowSds(P);
  mat mu_T_mat(12, t);
  mat s_T_mat(12, t);
  mat mu_P_mat(12, t);
  mat s_P_mat(12, t);
  for(int i = 0; i < t; i++){
    mu_T_mat.col(i) = mu_T;
    s_T_mat.col(i) = s_T;
    mu_P_mat.col(i) = mu_P;
    s_P_mat.col(i) = s_P;
  }
 
  //
  // initialize MCMC variables
  //
 
  double mu_beta_0 = R::rnorm(mu_mu_beta_0, sqrt(s2_mu_beta_0));
  double mu_beta_0_tilde = R::rnorm(mu_mu_beta_0, sqrt(s2_mu_beta_0));
  double s2_beta_0 = 1.0 / R::rgamma(alpha_s2_beta_0, 1.0 / beta_s2_beta_0);
  double s2_beta_0_tilde = 1.0 / R::rgamma(alpha_s2_beta_0, 1.0 / beta_s2_beta_0);
  double mu_beta_1 = R::rnorm(mu_mu_beta_1, sqrt(s2_mu_beta_1));
  double mu_beta_1_tilde = R::rnorm(mu_mu_beta_1, sqrt(s2_mu_beta_1));
  double s2_beta_1 = 1.0 / R::rgamma(alpha_s2_beta_1, 1.0 / beta_s2_beta_1);
  double s2_beta_1_tilde = 1.0 / R::rgamma(alpha_s2_beta_1, 1.0 / beta_s2_beta_1);
  double mu_gamma_T = R::rnorm(mu_mu_gamma_T, sqrt(s2_mu_gamma_T));
  double s2_gamma_T = 1.0;
//  / R::rgamma(alpha_s2_gamma_T, 1.0 / beta_s2_gamma_T);
  double mu_xi2_T = R::rnorm(mu_mu_xi2_T, sqrt(s2_mu_xi2_T));
  double s2_xi2_T = 1.0;
//  / R::rgamma(alpha_s2_xi2_T, 1.0 / beta_s2_xi2_T);
  double mu_gamma_P = R::rnorm(mu_mu_gamma_P, sqrt(s2_mu_gamma_P));
  double s2_gamma_P = 1.0;
//  / R::rgamma(alpha_s2_gamma_P, 1.0 / beta_s2_gamma_P);
  double mu_xi2_P = R::rnorm(mu_mu_xi2_P, sqrt(s2_mu_xi2_P));
  double s2_xi2_P = 1.0;
//  / R::rgamma(alpha_s2_xi2_P, 1.0 / beta_s2_xi2_P);
 
  vec beta_0(num_species);
  vec beta_1(num_species);
  vec s2(num_species); 
  vec beta_0_tilde(num_species);
  vec beta_1_tilde(num_species);
  vec s2_tilde(num_species);
  mat x(t, num_species);
  vec Temp_min(num_species);
  vec Temp_max(num_species);
  vec P_min(num_species);
  vec P_max(num_species);
  vec gamma_T(num_species);
  vec xi2_T(num_species);
  vec gamma_P(num_species);
  vec xi2_P(num_species);
 
  for(int i = 0; i < num_species; i++){
    for(int l = 0; l < t; l++){
      x(l, i) = R::rbinom(1, psi);
    }
    beta_0(i) = R::rnorm(mu_beta_0, sqrt(s2_beta_0));
    beta_1(i) = R::rnorm(mu_beta_1, sqrt(s2_beta_1));
    s2(i) = 1.0 / R::rgamma(alpha_s2, 1.0 / beta_s2);
    beta_0_tilde(i) = R::rnorm(mu_beta_0, sqrt(s2_beta_0));
    beta_1_tilde(i) = R::rnorm(mu_beta_1, sqrt(s2_beta_1));
    s2_tilde(i) = 1.0 / R::rgamma(alpha_s2, 1.0 / beta_s2);
    Temp_min(i) = R::rbeta(alpha_Temp_min, beta_Temp_min) * 
      (Temp_min_upper - Temp_min_lower) + Temp_min_lower;
    Temp_max(i) = R::rbeta(alpha_Temp_max, beta_Temp_max) * 
      (Temp_max_upper - Temp_max_lower) + Temp_max_lower;
    P_min(i) = R::rbeta(alpha_P_min, beta_P_max) * 
      (P_min_upper - P_min_lower) + P_min_lower;
    P_max(i) = R::rbeta(alpha_P_max, beta_P_max) * 
      (P_max_upper - P_max_lower) + P_max_lower;
    gamma_T(i) = R::rnorm(mu_gamma_T, sqrt(s2_gamma_T));
    xi2_T(i) = R::rlnorm(mu_xi2_T, sqrt(s2_xi2_T));
    if(xi2_T(i) > 50){
    	xi2_T(i) = 50;
    }
    gamma_P(i) = R::rnorm(mu_gamma_P, sqrt(s2_gamma_P));
    xi2_P(i) = R::rlnorm(mu_xi2_P, sqrt(s2_xi2_P));
    if(xi2_P(i) > 50){
    	xi2_P(i) = 50;
    }
  }
  vec xi_T = sqrt(xi2_T);
  vec xi_P = sqrt(xi2_P);
 
  // adaptive tuning
  vec Temp_min_tune = Temp_min_tune_tmp * J_num_species;
  vec Temp_max_tune = Temp_max_tune_tmp * J_num_species;
  vec P_min_tune = P_min_tune_tmp * J_num_species;
  vec P_max_tune = P_max_tune_tmp * J_num_species;
  vec gamma_T_tune = gamma_T_tune_tmp * J_num_species;
  vec xi2_T_tune = xi2_T_tune_tmp * J_num_species;
  vec gamma_P_tune = gamma_P_tune_tmp * J_num_species;
  vec xi2_P_tune = xi2_P_tune_tmp * J_num_species;
 
  mat Sigma_s_T(12, 12, fill::zeros);
  mat Sigma_s_P(12, 12, fill::zeros);
  Sigma_s_T.diag() = pow(s_T, 2);
  Sigma_s_P.diag() = pow(s_P, 2);
  mat W_T = trans(mvrnormArma(t, mu_T, Sigma_s_T));
  W_T.cols(T_obs_idx - 1) = Temp;
  mat W_T_tilde = (W_T - mu_T_mat) / s_T_mat;
 
  mat W_P = trans(mvrnormArma(t, mu_P, Sigma_s_P));
  W_P.cols(P_obs_idx - 1) = P;
  mat W_P_tilde = (W_P - mu_P_mat) / s_P_mat;
  mat W = rbindARMA(W_T, W_P);
  mat W_tilde = rbindARMA(W_T_tilde, W_P_tilde);
 
  mat zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, 
                         Temp_min, Temp_max, P_min, P_max, species_m1, x);
  //
  // Initialize save variables
  //
 
  int n_burn = n_mcmc / 2;
  int n_save = (n_mcmc - n_burn) / n_thin;
  mat beta_0_save(num_species, n_save, fill::zeros);
  vec mu_beta_0_save(n_save, fill::zeros);
  vec s2_beta_0_save(n_save, fill::zeros);
  mat beta_0_tilde_save(num_species, n_save, fill::zeros);
  vec mu_beta_0_tilde_save(n_save, fill::zeros);
  vec s2_beta_0_tilde_save(n_save, fill::zeros);
  mat beta_1_save(num_species, n_save, fill::zeros);
  vec mu_beta_1_save(n_save, fill::zeros);
  vec s2_beta_1_save(n_save, fill::zeros);
  mat beta_1_tilde_save(num_species, n_save, fill::zeros);
  vec mu_beta_1_tilde_save(n_save, fill::zeros);
  vec s2_beta_1_tilde_save(n_save, fill::zeros);
  mat s2_save(num_species, n_save, fill::zeros);
  mat s2_tilde_save(num_species, n_save, fill::zeros);
  mat Temp_min_save(num_species, n_save, fill::zeros);
  mat Temp_max_save(num_species, n_save, fill::zeros);
  mat P_min_save(num_species, n_save, fill::zeros);
  mat P_max_save(num_species, n_save, fill::zeros);
  mat gamma_T_save(num_species, n_save, fill::zeros);
  vec mu_gamma_T_save(n_save, fill::zeros);
  vec s2_gamma_T_save(n_save, fill::zeros);
  mat xi2_T_save(num_species, n_save, fill::zeros);
  vec mu_xi2_T_save(n_save, fill::zeros);
  vec s2_xi2_T_save(n_save, fill::zeros);
  mat gamma_P_save(num_species, n_save, fill::zeros);
  vec mu_gamma_P_save(n_save, fill::zeros);
  vec s2_gamma_P_save(n_save, fill::zeros);
  mat xi2_P_save(num_species, n_save, fill::zeros);
  vec mu_xi2_P_save(n_save, fill::zeros);
  vec s2_xi2_P_save(n_save, fill::zeros);
  cube W_save(24, t, n_save, fill::zeros);
  cube x_save(t, num_species, n_save, fill::zeros);
  vec Temp_min_accept(num_species, fill::zeros);
  vec Temp_max_accept(num_species, fill::zeros);
  vec P_min_accept(num_species, fill::zeros);
  vec P_max_accept(num_species, fill::zeros);
  vec gamma_T_accept(num_species, fill::zeros);
  vec xi2_T_accept(num_species, fill::zeros);
  vec gamma_P_accept(num_species, fill::zeros);
  vec xi2_P_accept(num_species, fill::zeros);
  double W_T_accept = 0;
  double W_P_accept = 0;
  vec Temp_min_accept_tmp(num_species, fill::zeros);
  vec Temp_max_accept_tmp(num_species, fill::zeros);
  vec P_min_accept_tmp(num_species, fill::zeros);
  vec P_max_accept_tmp(num_species, fill::zeros);
  vec gamma_T_accept_tmp(num_species, fill::zeros);
  vec xi2_T_accept_tmp(num_species, fill::zeros);
  vec gamma_P_accept_tmp(num_species, fill::zeros);
  vec xi2_P_accept_tmp(num_species, fill::zeros);
  double W_T_accept_tmp = 0;
  double W_P_accept_tmp = 0;
 
  //
  // start MCMC chain
  //
 
  Rcout << "\n\n" << "Starting Hierarchical Mixed model fit with random growth model, will run for " << n_mcmc << 
    " iterations\n\n";
  for(int k = 0; k < n_mcmc; k++){
    if((k + 1) % 100 == 0){
      Rcout << " " << k + 1;
    }
   
    Rcpp::checkUserInterrupt();
    
    //
    // sample W_T
    //
    vec like_T = makeLikelihoodMix(y, H, beta_0, beta_1, beta_0_tilde, 
                                   beta_1_tilde, zeta, t, p, species_m1, s2,
                                   s2_tilde, x);
    // sample for t = t-N_T to 2
    for(int i = t - N_T - 1; i > 0; i--){
      vec W_tilde_star_T = W_tilde.col(i);
      vec W_star_vec_T = W.col(i);
      for(int l = 0; l < 12; l++){
  W_tilde_star_T(l)  =  W_tilde_star_T(l) + R::rnorm(0, W_T_tune);
  W_star_vec_T(l) = W_tilde_star_T(l) * s_T(l) + mu_T(l);
      }
      rowvec x_ind_T = x.row(i);
      vec zeta_star_T = makeZetaIndividualMix(p, day_len, W_star_vec_T, gamma_T, 
                xi_T, gamma_P, xi_P, Temp_min, 
                Temp_max, P_min, P_max, 
                species_m1, x_ind_T);
      vec tmp1_star_T(24);
      vec tmp2_star_T(24);
      vec tmp1_T(24);
      vec tmp2_T(24);
      if (i == t - N_T - 1){
  tmp1_star_T = (W_tilde_star_T - Jtrend_0 - Jtrend_1) - phi_vec % 
    (W_tilde.col(i - 1) - Jtrend_0 - Jtrend_1);        
  tmp2_star_T = (W_tilde.col(i + 1) - Jtrend_0 - Jtrend_1) - 
    phi_vec % (W_tilde_star_T - Jtrend_0 - Jtrend_1);
  tmp1_T = (W_tilde.col(i) - Jtrend_0 - Jtrend_1) - phi_vec % 
    (W_tilde.col(i - 1) - Jtrend_0 - Jtrend_1);
  tmp2_T = (W_tilde.col(i + 1) - Jtrend_0 - Jtrend_1) - phi_vec % 
    (W_tilde.col(i) - Jtrend_0 - Jtrend_1);
      } else {
  tmp1_star_T = (W_tilde_star_T - Jtrend_0) - phi_vec % 
    (W_tilde.col(i - 1) - Jtrend_0);
  tmp2_star_T = (W_tilde.col(i + 1) - Jtrend_0) - phi_vec % 
    (W_tilde_star_T - Jtrend_0);
  tmp1_T = (W_tilde.col(i) - Jtrend_0) - phi_vec % 
    (W_tilde.col(i - 1) - Jtrend_0);
  tmp2_T = (W_tilde.col(i + 1) - Jtrend_0) - phi_vec % 
    (W_tilde.col(i) - Jtrend_0);
      }
      double mh1_T = makeLikelihoodIndividualMix(y, H, beta_0, beta_1, 
                                                 beta_0_tilde, beta_1_tilde, 
                                                 zeta_star_T, i, p, species_m1,
                                                 s2, s2_tilde, x) - 0.5 * 
  as_scalar(tmp1_star_T.t() * Sigma_inv * tmp1_star_T) - 
  0.5 * 
  as_scalar(tmp2_star_T.t() * Sigma_inv * tmp2_star_T);
      double mh2_T = like_T(i) - 0.5 * 
  as_scalar(tmp1_T.t() * Sigma_inv * tmp1_T) - 
  0.5 * as_scalar(tmp2_T.t() * Sigma_inv * tmp2_T);
      double mh_T = exp(mh1_T - mh2_T);
      if(mh_T > R::runif(0, 1)){
  W_tilde.col(i) = W_tilde_star_T;
  W.col(i) = W_star_vec_T;
  if(k > n_burn){
    W_T_accept += 1.0  / (n_burn * (t - N_T));
  }
  W_T_accept_tmp += 1.0  / (50 * (t - N_T));
      }
    }
    // sample for t = 1
    vec W_tilde_star_T = W_tilde.col(0);
    vec W_star_vec_T = W.col(0);
    for(int l = 0; l < 12; l++){
      W_tilde_star_T(l) = W_tilde_star_T(l) + R::rnorm(0, W_T_tune);
      W_star_vec_T(l) = (W_tilde_star_T(l) + trend_0)  * s_T(l) + mu_T(l);
    }
    rowvec x_ind_T = x.row(0);
    vec zeta_star_T = makeZetaIndividualMix(p, day_len, W_star_vec_T, gamma_T,
              xi_T, gamma_P, xi_P, Temp_min, 
              Temp_max, P_min, P_max, species_m1,
              x_ind_T);
    vec tmp_star_T = (W_tilde.col(1) - Jtrend_0) - 
      phi_vec % (W_tilde_star_T - Jtrend_0);
    int idx_val_T = 0;
    double mh1_T = makeLikelihoodIndividualMix(y, H, beta_0, beta_1, 
                                               beta_0_tilde, beta_1_tilde, 
                                               zeta_star_T, idx_val_T, p, 
                                               species_m1, s2, s2_tilde, x) -
      0.5 * as_scalar(tmp_star_T.t() * Sigma_inv * tmp_star_T);
    vec tmp_T = (W_tilde.col(1) - Jtrend_0) - phi_vec % 
      (W_tilde.col(0) - Jtrend_0);
    double mh2_T = - 0.5 * as_scalar(tmp_T.t() * Sigma_inv * tmp_T) + like_T(0);
    double mh_T = exp(mh1_T - mh2_T);
    if(mh_T > R::runif(0, 1)){
      W_tilde.col(0) = W_tilde_star_T;
      W.col(0) = W_star_vec_T;
      if(k > n_burn){
  W_T_accept += 1.0  / (n_burn * (t - N_T));  
      }
      W_T_accept_tmp += 1.0  / (50 * (t - N_T));  
    } 
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min, 
                       Temp_max, P_min, P_max, species_m1, x);
    
    //
    // sample W_P
    //
    
    vec like_P = makeLikelihoodMix(y, H, beta_0, beta_1, beta_0_tilde, 
                                   beta_1_tilde, zeta, t, p, species_m1, s2, 
                                   s2_tilde, x);
    // sample for t = t-N_T to 2
    for(int i = t - N_T - 1; i > 0; i--){
      vec W_tilde_star = W_tilde.col(i);
      vec W_star_vec = W.col(i);
      for(int l = 12; l < 24; l++){
  W_tilde_star(l)  =  W_tilde_star(l) + R::rnorm(0, W_P_tune);
  W_star_vec(l) = W_tilde_star(l) * s_P(l - 12) + mu_P(l - 12);
      }
      rowvec x_ind_P = x.row(i);
      vec zeta_star_P = makeZetaIndividualMix(p, day_len, W_star_vec, gamma_T,
                xi_T, gamma_P, xi_P, Temp_min, 
                Temp_max, P_min, P_max, 
                species_m1, x_ind_P);
      //     vec tmp1_star(24);
      //     vec tmp2_star(24);
      vec tmp1_star = W_tilde_star - phi_vec % W_tilde.col(i - 1);
      vec tmp2_star = W_tilde.col(i + 1) - phi_vec % W_tilde_star;
      double mh1 = makeLikelihoodIndividualMix(y, H, beta_0, beta_1, 
                                               beta_0_tilde, beta_1_tilde,
                                               zeta_star_P, i, p, species_m1,
                                               s2, s2_tilde, x)- 0.5 * 
  as_scalar(tmp1_star.t() * Sigma_inv * tmp1_star) - 
  0.5 * as_scalar(tmp2_star.t() * Sigma_inv * tmp2_star);
      vec tmp1 = W_tilde.col(i) - phi_vec % W_tilde.col(i - 1);
      vec tmp2 = W_tilde.col(i + 1) - phi_vec % W_tilde.col(i);
      double mh2 = like_P(i) - 0.5 * as_scalar(tmp1.t() * Sigma_inv * tmp1) - 
  0.5 * as_scalar(tmp2.t() * Sigma_inv * tmp2);
      double mh = exp(mh1 - mh2);
      if(mh > R::runif(0, 1)){
  W_tilde.col(i) = W_tilde_star;
  W.col(i) = W_star_vec;
  if(k > n_burn){
    W_P_accept += 1.0  / (n_burn * (t - N_P));
  }
  W_P_accept_tmp += 1.0  / (50 * (t - N_P));
      }
    }
    // sample for t = 1
    vec W_tilde_star = W_tilde.col(0);
    vec W_star_vec = W.col(0);
    for(int l = 12; l < 24; l++){
      W_tilde_star(l) = W_tilde_star(l) + R::rnorm(0, W_P_tune);
      W_star_vec(l) = W_tilde_star(l)  * s_P(l - 12) + mu_P(l - 12);
    }
    rowvec x_ind_P = x.row(0);
    vec zeta_star_P = makeZetaIndividualMix(p, day_len, W_star_vec, gamma_T, 
              xi_T, gamma_P, xi_P, Temp_min, 
              Temp_max, P_min, P_max, species_m1, 
              x_ind_P);
    vec tmp_star = W_tilde.col(1) - phi_vec % W_tilde_star;
    int idx_val = 0;
    double mh1 = makeLikelihoodIndividualMix(y, H, beta_0, beta_1, beta_0_tilde,
                                             beta_1_tilde, zeta_star_P, 
                                             idx_val, p, species_m1, s2, 
                                             s2_tilde, x) - 
      0.5 * as_scalar(tmp_star.t() * Sigma_inv * tmp_star);
    vec tmp = W_tilde.col(1) - phi_vec % W_tilde.col(0);
    double mh2 = like_P(0) - 0.5 * as_scalar(tmp.t() * Sigma_inv * tmp);
    double mh = exp(mh1 - mh2);
    if(mh > R::runif(0, 1)){
      W_tilde.col(0) = W_tilde_star;
      W.col(0) = W_star_vec;
      if(k > n_burn){
  W_P_accept += 1.0  / (n_burn * (t - N_P));
      }
      W_P_accept_tmp += 1.0  / (50 * (t - N_P));
    } 
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min, 
                       Temp_max, P_min, P_max, species_m1, x);
    
    //
    // sample beta_0
    //                        
    mat A0(num_species, num_species, fill::zeros);
    vec A_vec0 = makeBeta0VarMix(H, s2, species_m1, t, p, num_species, s2_beta_0, x);
    vec b0 = makeBeta0MeanMix(y, H, beta_1, zeta, species_m1, t, p, num_species, 
            s2, mu_beta_0, s2_beta_0, x);
    A0.diag() = A_vec0;
    beta_0 = rMVNArma(A0, b0);
   
    // sample mu_beta_0
    double a_mu_beta_0 = num_species / s2_beta_0 + 1.0 / s2_mu_beta_0;
    double b_mu_beta_0 = sum(beta_0) / s2_beta_0 + mu_mu_beta_0 / s2_mu_beta_0;
    mu_beta_0 = rMVNArmaScalar(a_mu_beta_0, b_mu_beta_0);
   
    // sample s2_beta_0
    s2_beta_0 = 1.0 / R::rgamma(alpha_s2_beta_0 + 0.5 * num_species, 1.0 /
        (beta_s2_beta_0 + 0.5 *
         sum(pow(beta_0 - mu_beta_0, 2))));
   
    //
    // sample beta_1
    //
   
    mat A1(num_species, num_species, fill::zeros);
    vec A_vec1 = makeBeta1VarMix(H, zeta, species_m1, t, p, num_species, s2, s2_beta_1, x);
    vec b1 = makeBeta1MeanMix(y, H, beta_0, zeta, species_m1, t, p, 
            num_species, s2, mu_beta_1, s2_beta_1, x);    
    A1.diag() = A_vec1;
    beta_1 = rMVNArma(A1, b1);
   
    // sample mu_beta_1
    double a_mu_beta_1 = num_species / s2_beta_1 + 1.0 / s2_mu_beta_1;
    double b_mu_beta_1 = sum(beta_1) / s2_beta_1 + mu_mu_beta_1 / s2_mu_beta_1;
    mu_beta_1 = rMVNArmaScalar(a_mu_beta_1, b_mu_beta_1);
   
    // sample s2_beta_1
    s2_beta_1 = 1.0 / R::rgamma(alpha_s2_beta_1 + 0.5 * num_species, 1.0 /
        (beta_s2_beta_1 + 0.5 *
         sum(pow(beta_1 - mu_beta_1, 2))));  
   
    //
    // sample beta_0_tilde
    //
   
    mat A0_tilde(num_species, num_species, fill::zeros);
    vec A_vec0_tilde = makeBeta0VarTildeMix(H, s2_tilde, species_m1, t, p, 
              num_species, s2_beta_0_tilde, x);
    vec b0_tilde = makeBeta0MeanTildeMix(y, H, beta_1_tilde, zeta, species_m1, t, p,
           num_species, s2_tilde, mu_beta_0_tilde, s2_beta_0_tilde, x);
    A0_tilde.diag() = A_vec0_tilde;
    beta_0_tilde = rMVNArma(A0_tilde, b0_tilde);
   
    // sample mu_beta_0_tilde
    double a_mu_beta_0_tilde = num_species / s2_beta_0_tilde + 1.0 / s2_mu_beta_0;
    double b_mu_beta_0_tilde = sum(beta_0_tilde) / s2_beta_0_tilde + mu_mu_beta_0 / s2_mu_beta_0;
    mu_beta_0_tilde = rMVNArmaScalar(a_mu_beta_0_tilde, b_mu_beta_0_tilde);
   
    // sample s2_beta_0_tilde
    s2_beta_0_tilde = 1.0 / R::rgamma(alpha_s2_beta_0 + 0.5 * num_species, 1.0 /
              (beta_s2_beta_0 + 0.5 *
               sum(pow(beta_0_tilde - mu_beta_0_tilde, 2))));
   
    //
    // sample beta_1_tilde
    //
   
    mat A1_tilde(num_species, num_species, fill::zeros);
    vec A_vec1_tilde = makeBeta1VarTildeMix(H, zeta, species_m1, t, p, num_species, 
              s2_tilde, s2_beta_1_tilde, x);
    vec b1_tilde = makeBeta1MeanTildeMix(y, H, beta_0_tilde, zeta, species_m1, t, p,
           num_species, s2_tilde, mu_beta_1_tilde, 
           s2_beta_1_tilde, x);          
    A1_tilde.diag() = A_vec1_tilde;
    beta_1_tilde = rMVNArma(A1_tilde, b1_tilde);
   
    // sample mu_beta_1_tilde
    double a_mu_beta_1_tilde = num_species / s2_beta_1_tilde + 
      1.0 / s2_mu_beta_1;
    double b_mu_beta_1_tilde = sum(beta_1_tilde) / s2_beta_1_tilde + 
      mu_mu_beta_1 / s2_mu_beta_1;
    mu_beta_1_tilde = rMVNArmaScalar(a_mu_beta_1_tilde, b_mu_beta_1_tilde);
   
    // sample s2_beta_1
    s2_beta_1_tilde = 1.0 / R::rgamma(alpha_s2_beta_1 + 0.5 * num_species, 1.0 /
              (beta_s2_beta_1 + 0.5 * 
               sum(pow(beta_1_tilde - mu_beta_1_tilde, 2))));  
   
    //
    // sample VS-Lite paramters
    //
   
    // sample Temp_min
    vec Temp_min_star(num_species);
    for(int i = 0; i < num_species; i++){
      Temp_min_star(i) = R::rnorm(Temp_min(i), Temp_min_tune(i));
      if(Temp_min_star(i) < Temp_min_lower){
  Temp_min_star(i) = Temp_min_lower;
      }
      if(Temp_min_star(i) > Temp_min_upper){
  Temp_min_star(i) = Temp_min_upper;
      }
    }
    mat zeta_star = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P,       
        Temp_min_star, Temp_max, P_min, P_max, 
        species_m1, x);
    vec like_star = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde,
              zeta_star, species_m1, 
              t, p, num_species, s2_tilde, x);
    vec like = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde, 
               zeta, species_m1, t, p, 
               num_species, s2_tilde, x);
    for(int i = 0; i < num_species; i++){
      double kernel_star = (Temp_min_star(i) - Temp_min_lower) / 
  (Temp_min_upper - Temp_min_lower);
      double kernel = (Temp_min(i) - Temp_min_lower) / 
  (Temp_min_upper - Temp_min_lower);
      double mh1 = like_star(i) + // conditional likelihood under proposed value
  // density of proposal under prior
  (alpha_Temp_min - 1.0) * log(kernel_star) + 
  (beta_Temp_min - 1.0) * log(1.0 - kernel_star);             
      double mh2 = like(i) + // conditional likelihood under current value
  // density of current under prior
  (alpha_Temp_min - 1.0) * log(kernel) + (beta_Temp_min - 1.0) * 
  log(1.0 - kernel);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
  Temp_min(i) = Temp_min_star(i);
  if(k > n_burn){
    Temp_min_accept(i) += 1.0 / n_burn;
  }
  Temp_min_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min,
           Temp_max, P_min, P_max, species_m1, x);
   
    // sample Temp_max
    vec Temp_max_star(num_species);
    for(int i = 0; i < num_species; i++){
      Temp_max_star(i) = R::rnorm(Temp_max(i), Temp_max_tune(i));
      if(Temp_max_star(i) < Temp_max_lower){
        Temp_max_star(i) = Temp_max_lower;
      }
      if(Temp_max_star(i) > Temp_max_upper){
        Temp_max_star(i) = Temp_max_upper;
      }
    }
    zeta_star = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P,
          Temp_min, Temp_max_star, P_min, P_max, species_m1, x);
    like_star = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde,
                zeta_star, species_m1, 
                t, p, num_species, s2_tilde, x);
    like = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde, 
           zeta, species_m1, t, p, 
           num_species, s2_tilde, x);
    for(int i = 0; i < num_species; i++){
      double kernel_star = (Temp_max_star(i) - Temp_max_lower) / 
  (Temp_max_upper - Temp_max_lower);
      double kernel = (Temp_max(i) - Temp_max_lower) / 
  (Temp_max_upper - Temp_max_lower);
      double mh1 = like_star(i) + // conditional likelihood under proposed value
  // density of proposal under prior
  (alpha_Temp_max - 1.0) * log(kernel_star) + 
  (beta_Temp_max - 1.0) * log(1.0 - kernel_star);             
      double mh2 = like(i) + // conditional likelihood under proposed value
  // density of proposal under prior
  (alpha_Temp_max - 1.0) * log(kernel) + (beta_Temp_max - 1.0) *
  log(1.0 - kernel);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
  Temp_max(i) = Temp_max_star(i);
  if(k > n_burn){
    Temp_max_accept(i) += 1.0 / n_burn;
  }
  Temp_max_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min,
           Temp_max, P_min, P_max, species_m1, x);
   
    // sample P_min
    vec P_min_star(num_species);
    for(int i = 0; i < num_species; i++){
      P_min_star(i) = R::rnorm(P_min(i), P_min_tune(i));
      if(P_min_star(i) < P_min_lower){
  P_min_star(i) = P_min_lower;
      }
      if(P_min_star(i) > P_min_upper){
  P_min_star(i) = P_min_upper;
      }
    }
   
    zeta_star = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P,
          Temp_min, Temp_max, P_min_star, P_max, species_m1, 
          x);
    like_star = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde,
                zeta_star, species_m1, 
                t, p, num_species, s2_tilde, x);
    like = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde, 
           zeta, species_m1, t, p, 
           num_species, s2_tilde, x);
    for(int i = 0; i < num_species; i++){
      double kernel_star = (P_min_star(i) - P_min_lower) / 
  (P_min_upper - P_min_lower);
      double kernel = (P_min(i) - P_min_lower) / (P_min_upper - P_min_lower);
      double mh1 = like_star(i) + // conditional likelihood under proposal
  // density of proposal under prior
  (alpha_P_min - 1.0) * log(kernel_star) + (beta_P_min - 1.0) * 
  log(1.0 - kernel_star);             
      double mh2 = like(i) + // conditional likelihood under current value
  // density of current value under prior
  (alpha_P_min - 1.0) * log(kernel) + (beta_P_min - 1.0) * 
  log(1.0 - kernel);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
  P_min(i) = P_min_star(i);
  if(k > n_burn){
    P_min_accept(i) += 1.0 / n_burn;
  }
  P_min_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min,
           Temp_max, P_min, P_max, species_m1, x);
   
    // sample P_max
    vec P_max_star(num_species);
    for(int i = 0; i < num_species; i++){
      P_max_star(i) = R::rnorm(P_max(i), P_max_tune(i));
      if(P_max_star(i) < P_max_lower){
  P_max_star(i) = P_max_lower;
      }
      if(P_max_star(i) > P_max_upper){
  P_max_star(i) = P_max_upper;
      }
    }
   
    zeta_star = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, 
          Temp_min, Temp_max, P_min, P_max_star, species_m1, x);
    like_star = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde,
                zeta_star, species_m1, 
                t, p, num_species, s2_tilde, x);
    like = makeLikelihoodSpeciesTildeMix(y, H, beta_0_tilde, beta_1_tilde, 
           zeta, species_m1, t, p, 
           num_species, s2_tilde, x);
    for(int i = 0; i < num_species; i++){
      double kernel_star = (P_max_star(i) - P_max_lower) / 
  (P_max_upper - P_max_lower);
      double kernel = (P_max(i) - P_max_lower) / (P_max_upper - P_max_lower);
      double mh1 = like_star(i) + // conditional likelihood under proposal
  // density of proposal under prior
  (alpha_P_max - 1.0) * log(kernel_star) + (beta_P_max - 1.0) * 
  log(1.0 - kernel_star);             
      double mh2 = like(i) + // conditional likelihood under current value
  // density of current value under prior
  (alpha_P_max - 1.0) * log(kernel) + (beta_P_max - 1.0) * 
  log(1.0 - kernel);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
  P_max(i) = P_max_star(i);
  if(k > n_burn){
    P_max_accept(i) += 1.0 / n_burn;
  }
  P_max_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min,
           Temp_max, P_min, P_max, species_m1, x);
   
    //
    // sample probit paramters
    //
   
    // sample gamma_T
    vec gamma_T_star(num_species); 
    for(int i = 0; i < num_species; i++){
      gamma_T_star(i) = R::rnorm(gamma_T(i), gamma_T_tune(i));
    }
    zeta_star = makeZetaMix(t, p, day_len, W, gamma_T_star, xi_T, gamma_P, xi_P, 
          Temp_min, Temp_max, P_min, P_max, species_m1, 
          x);
    like_star = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta_star, 
           species_m1, t, p, num_species, s2, x);
    like = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta, species_m1, t, p, 
            num_species, s2, x);
    for(int i = 0; i < num_species; i++){
      double mh1 = like_star(i) + // conditional likelihood under proposal
  // density of proposal under prior 
  - 0.5 / s2_gamma_T * pow(gamma_T_star(i) - mu_gamma_T, 2);
      double mh2 = like(i) + // conditional likelihood under current value
  // density of current value under prior
  - 0.5 / s2_gamma_T * pow(gamma_T(i) - mu_gamma_T, 2);
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
  gamma_T(i) = gamma_T_star(i);
  if(k > n_burn){
    gamma_T_accept(i) += 1.0 / n_burn;
  }
  gamma_T_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min, 
           Temp_max, P_min, P_max, species_m1, x);
   
    // sample mu_gamma_T
    double a_mu_gamma_T = num_species / s2_gamma_T + 1.0 / s2_mu_gamma_T;
    double b_mu_gamma_T = sum(gamma_T) / s2_gamma_T + mu_mu_gamma_T / s2_mu_gamma_T;
    mu_gamma_T = rMVNArmaScalar(a_mu_gamma_T, b_mu_gamma_T);
   
    // sample s2_gamma_T
    s2_gamma_T = 1.0 / R::rgamma(alpha_s2_gamma_T + 0.5 * num_species, 1.0 / 
         (beta_s2_gamma_T + 0.5 * sum(pow(gamma_T - mu_gamma_T, 2))));
   
    // sample xi2_T
    vec xi2_T_star(num_species);
    vec xi_T_star = xi_T;
    for(int i = 0; i < num_species; i++){
      xi2_T_star(i) = R::rnorm(xi2_T(i), xi2_T_tune(i));
      if(xi2_T_star(i) > 0){
  xi_T_star(i) = sqrt(xi2_T_star(i));
      }
    }
    zeta_star = makeZetaMix(t, p, day_len, W, gamma_T, xi_T_star, gamma_P, xi_P, Temp_min, 
          Temp_max, P_min, P_max, species_m1, x);
    // growth model likelihood by species for proposed value
    like_star = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta_star, 
           species_m1, t, p, num_species, s2, x);
    // growth model likelihood by species for current value
    like = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta, species_m1, t, p, 
            num_species, s2, x);
    for(int i = 0; i < num_species; i++){
      if(xi2_T_star(i) > 0){
  double mh1 = like_star(i) + // conditional likelihood under proposal
    //   // likelihood of proposal under IG prior
    //   - (alpha_xi2_T + 1.0) * log(xi2_T_star(i)) - beta_xi2_T / xi2_T_star(i);
    // likelihood of proposal under Lognormal prior
    - log(xi2_T_star(i)) - 0.5 * pow(log(xi2_T_star(i)) - mu_xi2_T, 2) / s2_xi2_T;
  double mh2 = like(i) + // conditional likelihood under current value
    //   // likelihood of current value under IG prior
    //   - (alpha_xi2_T + 1.0) * log(xi2_T(i)) - beta_xi2_T / xi2_T(i);
    // likelihood of current value under Lognormal prior
    - log(xi2_T(i)) - 0.5 * pow(log(xi2_T(i)) - mu_xi2_T, 2) / s2_xi2_T;
  double mh = exp(mh1 - mh2);   
  if(mh > R::runif(0, 1)){
    xi2_T(i) = xi2_T_star(i);
    xi_T(i) = xi_T_star(i);
    if(k > n_burn){
      xi2_T_accept(i) += 1.0 / n_burn;
    }
    xi2_T_accept_tmp(i) += 1.0 / 50;
  }
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min, 
           Temp_max, P_min, P_max, species_m1, x);
     
    // sample mu_xi2_T
    double a_mu_xi2_T = num_species / s2_xi2_T + 1.0 / s2_mu_xi2_T;
    double b_mu_xi2_T = sum(log(xi2_T)) / s2_xi2_T + mu_mu_xi2_T / s2_mu_xi2_T;
    mu_xi2_T = rMVNArmaScalar(a_mu_xi2_T, b_mu_xi2_T);
   
    // sample s2_xi2_T
    s2_xi2_T = 1.0 / R::rgamma(alpha_s2_xi2_T + 0.5 * num_species, 1.0 / 
             (beta_s2_xi2_T + 0.5 * sum(pow(log(xi2_T) - mu_xi2_T, 2))));
                              
    // sample gamma_P
    vec gamma_P_star = gamma_P; 
    for(int i = 0; i < num_species; i++){
      gamma_P_star(i) = R::rnorm(gamma_P(i), gamma_P_tune(i));
    }
    zeta_star = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P_star, xi_P, 
          Temp_min, Temp_max, P_min, P_max, species_m1, 
          x);
    // growth model likelihood by species for proposed value
    like_star = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta_star, 
           species_m1, t, p, num_species, s2, x);
    // growth model likelihood by species for current value
    like = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta, species_m1, t, p, 
            num_species, s2, x);
    for(int i = 0; i < num_species; i++){
      double mh1 = like_star(i) + // conditional likelihood under proposal
  // likelihood of proposal under prior
  - 0.5 / s2_gamma_P * pow(gamma_P_star(i) - mu_gamma_P, 2);
      double mh2 = like(i) + // conditional likelihood under proposal
  // likelihood of proposal under prior
  - 0.5 / s2_gamma_P * pow(gamma_P(i) - mu_gamma_P, 2);
      double mh = exp(mh1 -mh2);    
      if(mh > R::runif(0, 1)){
  gamma_P(i) = gamma_P_star(i);
  if(k > n_burn){
    gamma_P_accept(i) += 1.0 / n_burn;
  }
  gamma_P_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min,
           Temp_max, P_min, P_max, species_m1, x);
     
    // sample mu_gamma_P
    double a_mu_gamma_P = num_species / s2_gamma_P + 1.0 / s2_mu_gamma_P;
    double b_mu_gamma_P = sum(gamma_P) / s2_gamma_P +
      mu_mu_gamma_P / s2_mu_gamma_P;
    mu_gamma_P = rMVNArmaScalar(a_mu_gamma_P, b_mu_gamma_P);
     
    // sample s2_gamma_P
    s2_gamma_P = 1.0 / R::rgamma(alpha_s2_gamma_P + 0.5 * num_species, 1.0 / 
         (beta_s2_gamma_P + 0.5 * sum(pow(gamma_P - mu_gamma_P, 2))));
     
    // sample xi2_P
    vec xi2_P_star(num_species);
    vec xi_P_star = xi_P;
    for(int i = 0; i < num_species; i++){
      xi2_P_star(i) = R::rnorm(xi2_P(i), xi2_P_tune(i));
      if(xi2_P_star(i) > 0){
  xi_P_star(i) = sqrt(xi2_P_star(i));
      }
    }
    zeta_star = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P_star, Temp_min, 
          Temp_max, P_min, P_max, species_m1, x);
    // growth model likelihood by species for proposed value
    like_star = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta_star, 
           species_m1, t, p, num_species, s2, x);
    // growth model likelihood by species for current value
    like = makeLikelihoodSpeciesMix(y, H, beta_0, beta_1, zeta, species_m1, t, p, 
            num_species, s2, x);
    for(int i = 0; i < num_species; i++){
      if(xi2_P_star(i) > 0){
  double mh1 = like_star(i) + // conditional likelihood under proposal
    //   // likelihood of proposal under IG prior
    //   - (alpha_xi2_P + 1.0) * log(xi2_P_star(i)) - beta_xi2_P / xi2_P_star(i);
    // likelihood of proposal under Lognormal prior
    - log(xi2_P_star(i)) - 0.5 * pow(log(xi2_P_star(i)) - mu_xi2_P, 2) / s2_xi2_P;
  double mh2 = like(i) + // conditional likelihood under current value
    //   // likelihood of current value under IG prior
    //   - (alpha_xi2_P + 1.0) * log(xi2_P(i)) - beta_xi2_P / xi2_P(i);
    // likelihood of current value under Lognormal prior
    - log(xi2_P(i)) - 0.5 * pow(log(xi2_P(i)) - mu_xi2_P, 2) / s2_xi2_P;
  double mh = exp(mh1 - mh2);   
  if(mh > R::runif(0, 1)){
    xi2_P(i) = xi2_P_star(i);
    xi_P(i) = xi_P_star(i);
    if(k > n_burn){
      xi2_P_accept(i) += 1.0 / n_burn;
    }
    xi2_P_accept_tmp(i) += 1.0 / 50;
  }
      }
    }
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min, 
           Temp_max, P_min, P_max, species_m1, x);
     
    // sample mu_xi2_P
    double a_mu_xi2_P = num_species / s2_xi2_P + 1.0 / s2_mu_xi2_P;
    double b_mu_xi2_P = sum(log(xi2_P)) / s2_xi2_P + mu_mu_xi2_P / s2_mu_xi2_P;
    mu_xi2_P = rMVNArmaScalar(a_mu_xi2_P, b_mu_xi2_P);
   
    // sample s2_xi2_P
    s2_xi2_P = 1.0 / R::rgamma(alpha_s2_xi2_P + 0.5 * num_species, 1.0 / 
             (beta_s2_xi2_P + 0.5 * sum(pow(log(xi2_P) - mu_xi2_P, 2))));
     
    //
    // sample s2
    //         
     
    vec s2_alpha_vec = makeS2AlphaMix(H, species_m1, t, p, num_species, x);
    vec s2_beta_vec = makeS2BetaMix(y, H, beta_0, beta_1, zeta, species_m1, t, p, 
            num_species, x);
    for(int i = 0; i < num_species; i++){
      s2(i) = 1.0 / R::rgamma(alpha_s2 + 0.5 * s2_alpha_vec(i), 
            1.0 / (beta_s2 + 0.5 * s2_beta_vec(i)));
    }
       
    //
    // sample s2_tilde
    //
     
    vec s2_alpha_tilde_vec = makeS2AlphaTildeMix(H, species_m1, t, p, num_species, x);
    vec s2_beta_tilde_vec = makeS2BetaTildeMix(y, H, beta_0_tilde, beta_1_tilde, zeta,
                 species_m1, t, p, num_species, x);
    for(int i = 0; i < num_species; i++){
      s2_tilde(i) = 1.0 / R::rgamma(alpha_s2 + 0.5 * s2_alpha_tilde_vec(i), 
            1.0 / (beta_s2 + 0.5 * s2_beta_tilde_vec(i)));
    }          
     
    //
    // sample model choice matrix x
    //
     
    mat zetaPro = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P,        
            Temp_min, Temp_max, P_min, P_max, species_m1,
            one_mat);
    mat zetaVS = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, 
           Temp_min, Temp_max, P_min, P_max, species_m1,
           zero_mat);
    mat likePsiPro = makeLikelihoodPsi(y, H, beta_0, beta_1, zetaPro, 
               species_m1, num_species, t, p, s2);
    mat likePsiVS = makeLikelihoodPsi(y, H, beta_0_tilde, beta_1_tilde, zetaVS, 
              species_m1, num_species, t, p, s2_tilde);
    for(int i = 0; i < t; i++){
      for(int l = 0; l < num_species; l++){
  double psi_tilde = 1.0 / (1.0 + exp((likePsiVS(i, l) + log_one_minus_psi) - 
              (likePsiPro(i, l) + log_psi)));
  x(i, l) = R::rbinom(1, psi_tilde);
      }
    }
     
    zeta = makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min, 
           Temp_max, P_min, P_max, species_m1, x);
     
    //
    // Adjust tuning parameters
    //
     
    if((k + 1) % 50 == 0){ // start counting at 0
 		  // Adaptive tuning for Temp_min
    	updateTuningVec(k, Temp_min_accept_tmp, Temp_min_tune);
  	  // Adaptive tuning for Temp_max
  	  updateTuningVec(k, Temp_max_accept_tmp, Temp_max_tune);
  	  // Adaptive tuning for P_min
  	  updateTuningVec(k, P_min_accept_tmp, P_min_tune);
  	  // Adaptive tuning for P_max
  	  updateTuningVec(k, P_max_accept_tmp, P_max_tune);
  	  // Adaptive tuning for gamma_T
    	updateTuningVec(k, gamma_T_accept_tmp, gamma_T_tune);
  	  // Adaptive tuning for xi2_T
  	  updateTuningVec(k, xi2_T_accept_tmp, xi2_T_tune);
  	  // Adaptive tuning for gamma_P
  	  updateTuningVec(k, gamma_P_accept_tmp, gamma_P_tune);
  	  // Adaptive tuning for xi2_P
  	  updateTuningVec(k, xi2_P_accept_tmp, xi2_P_tune);
      // Adaptive tuning for W_T
      updateTuning(k, W_T_accept_tmp, W_T_tune);
      // Adaptive tuning for W_P
      updateTuning(k, W_P_accept_tmp, W_P_tune);


    }
      
    //
    // Save Variables
    //
     
    if(k > n_burn){
      if((k + 1) % n_thin == 0){
  int k_tmp = (k + 1 - n_burn) / n_thin - 1;
  beta_0_save.col(k_tmp) = beta_0;
  mu_beta_0_save(k_tmp) = mu_beta_0;
  s2_beta_0_save(k_tmp) = s2_beta_0;
  beta_1_save.col(k_tmp) = beta_1;
  mu_beta_1_save(k_tmp) = mu_beta_1;
  s2_beta_1_save(k_tmp) = s2_beta_1;
  beta_0_tilde_save.col(k_tmp) = beta_0_tilde;
  mu_beta_0_tilde_save(k_tmp) = mu_beta_0_tilde;
  s2_beta_0_tilde_save(k_tmp) = s2_beta_0_tilde;
  beta_1_tilde_save.col(k_tmp) = beta_1_tilde;
  mu_beta_1_tilde_save(k_tmp) = mu_beta_1_tilde;
  s2_beta_1_tilde_save(k_tmp) = s2_beta_1_tilde;
  s2_save.col(k_tmp) = s2;
  s2_tilde_save.col(k_tmp) = s2_tilde;
  Temp_min_save.col(k_tmp) = Temp_min;
  Temp_max_save.col(k_tmp) = Temp_max;
  P_min_save.col(k_tmp) = P_min;
  P_max_save.col(k_tmp) = P_max;
  gamma_T_save.col(k_tmp) = gamma_T;
  mu_gamma_T_save(k_tmp) = mu_gamma_T;
  s2_gamma_T_save(k_tmp) = s2_gamma_T;
  xi2_T_save.col(k_tmp) = xi2_T;
  mu_xi2_T_save(k_tmp) = mu_xi2_T;
  s2_xi2_T_save(k_tmp) = s2_xi2_T;
  gamma_P_save.col(k_tmp) = gamma_P;
  s2_gamma_P_save(k_tmp) = s2_gamma_P;
  mu_gamma_P_save(k_tmp) = mu_gamma_P;
  xi2_P_save.col(k_tmp) = xi2_P;
  mu_xi2_P_save(k_tmp) = mu_xi2_P;
  s2_xi2_P_save(k_tmp) = s2_xi2_P;
  W_save.slice(k_tmp) = W;
  x_save.slice(k_tmp) = x;
      }
    }    
  }    
  return(List::create(
          _["calibration"] = List::create(
                  _["beta_0"] = beta_0_save,
                  _["mu_beta_0"] = mu_beta_0_save,
                  _["s2_beta_0"] = s2_beta_0_save,
                  _["beta_1"] = beta_1_save,
                  _["mu_beta_1"] = mu_beta_1_save,
                  _["s2_beta_1"] = s2_beta_1_save,
                  _["s2"] = s2_save,
                  _["beta_0_tilde"] = beta_0_tilde_save,
                  _["mu_beta_0_tilde"] = mu_beta_0_tilde_save,
                  _["s2_beta_0_tilde"] = s2_beta_0_tilde_save,
                  _["beta_1_tilde"] = beta_1_tilde_save,
                  _["mu_beta_1_tilde"] = mu_beta_1_tilde_save,
                  _["s2_beta_1_tilde"] = s2_beta_1_tilde_save,
                  _["s2_tilde"] = s2_tilde_save
                  ),
          _["growth"] = List::create(
             _["Temp_min"] = Temp_min_save,
             _["Temp_max"] = Temp_max_save,
             _["P_min"] = P_min_save,
             _["P_max"] = P_max_save,
             _["gamma_T"] = gamma_T_save,
             _["mu_gamma_T"] = mu_gamma_T_save,
             _["s2_gamma_T"] = s2_gamma_T_save,
             _["xi2_T"] = xi2_T_save,   
             _["mu_xi2_T"] = mu_xi2_T_save,
             _["s2_xi2_T"] = s2_xi2_T_save,
             _["gamma_P"] = gamma_P_save,
             _["mu_gamma_P"] = mu_gamma_P_save,
             _["s2_gamma_P"] = s2_gamma_P_save,
             _["xi2_P"] = xi2_P_save,
             _["mu_xi2_P"] = mu_xi2_P_save,
             _["s2_xi2_P"] = s2_xi2_P_save),
          _["W"] = W_save,
          _["x"] = x_save,
          _["accept"] = List::create(
             _["Temp_min_accept"] = Temp_min_accept,
             _["Temp_max_accept"] = Temp_max_accept,
             _["P_min_accept"] = P_min_accept,
             _["P_max_accept"] = P_max_accept,
             _["gamma_T_accept"] = gamma_T_accept,
             _["xi2_T_accept"] = xi2_T_accept,
             _["gamma_P_accept"] = gamma_P_accept,
             _["xi2_P_accept"] = xi2_P_accept,
             _["W_T_accept"] = W_T_accept,
             _["W_P_accept"] = W_P_accept
             )
          ));
}  
    
