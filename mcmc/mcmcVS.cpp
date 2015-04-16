#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, myFunctions, bayesTreeRing)]]
#include "mcmc.h"

using namespace Rcpp;
using namespace arma;
 
// [[Rcpp::export]]
List makeMCMC(mat y, mat Temp, mat P, vec species, List params, List process,
              bool sim = false){
        	
  // priors for calibration model
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
  
  // priors for VS-Lite growth model
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
  
  int num_species = as<int>(params["num_species"]);
  vec day_len = as<vec>(params["day_len"]);
  double W_T_tune = as<double>(params["W_T_tune"]);
  double W_P_tune = as<double>(params["W_P_tune"]);
  uvec T_obs_idx = as<uvec>(params["T_obs_idx"]);
  uvec P_obs_idx = as<uvec>(params["P_obs_idx"]);
  int n_mcmc = as<int>(params["n_mcmc"]);
  int n_thin = as<int>(params["n_thin"]);
  
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
  vec N_obs_species(num_species, fill::zeros);
  for(int i = 0; i < p; i++){
    int species_idx = species(i) - 1;
    N_obs_species(species_idx) += H_col_sums(i); 
  }
     
  int t12 = 12 * t;
  IntegerVector idx = rep(seq_len(12), t);
  vec species_m1 = species - 1;
  mat I_12 = eye(12, 12);
  mat I_t12 = eye(t12, t12);
  mat I_num_species = eye(num_species, num_species);
  int N_T = T_obs_idx.n_elem; 
  int N_P = P_obs_idx.n_elem;
  
  vec J_num_species(num_species, fill::ones);
  vec J(24);
  for(int i = 0; i < 24; i++){
    if(i < 12){
      J(i) = 1;
    } else {
      J(i) = 0;
    }
  }
  
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
  double s2_beta_0 = 1.0 / R::rgamma(alpha_s2_beta_0, 1.0 / beta_s2_beta_0);
  double mu_beta_1 = R::rnorm(mu_mu_beta_1, sqrt(s2_mu_beta_1));
  double s2_beta_1 = 1.0 / R::rgamma(alpha_s2_beta_1, 1.0 / beta_s2_beta_1);
  
  vec beta_0(num_species);
  vec beta_1(num_species);
  vec s2(num_species); 
  vec Temp_min(num_species);
  vec Temp_max(num_species);
  vec P_min(num_species);
  vec P_max(num_species);
  
  for(int i = 0; i < num_species; i++){
    beta_0(i) = R::rnorm(mu_beta_0, sqrt(s2_beta_0));
    beta_1(i) = R::rnorm(mu_beta_1, sqrt(s2_beta_1));
    s2(i) = 1.0 / R::rgamma(alpha_s2, 1.0 / beta_s2);
    Temp_min(i) = R::rbeta(alpha_Temp_min, beta_Temp_min) * 
      (Temp_min_upper - Temp_min_lower) + Temp_min_lower;
    Temp_max(i) = R::rbeta(alpha_Temp_max, beta_Temp_max) * 
      (Temp_max_upper - Temp_max_lower) + Temp_max_lower;
    P_min(i) = R::rbeta(alpha_P_min, beta_P_max) *
      (P_min_upper - P_min_lower) + P_min_lower;
    P_max(i) = R::rbeta(alpha_P_max, beta_P_max) *
      (P_max_upper - P_max_lower) + P_max_lower;
  }

  // adaptive tuning
  vec Temp_min_tune = Temp_min_tune_tmp * J_num_species;
  vec Temp_max_tune = Temp_max_tune_tmp * J_num_species;
  vec P_min_tune = P_min_tune_tmp * J_num_species;
  vec P_max_tune = P_max_tune_tmp * J_num_species;

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
  
  mat zeta = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max, 
                      species_m1);
  
  //
  // Initialize save variables
  //
  
  int n_burn = n_mcmc / 2;
  int n_save = (n_mcmc - n_burn) / n_thin;
  mat beta_0_save(num_species, n_save, fill::zeros);
  vec mu_beta_0_save(n_save, fill::zeros);
  vec s2_beta_0_save(n_save, fill::zeros);
  mat beta_1_save(num_species, n_save, fill::zeros);
  vec mu_beta_1_save(n_save, fill::zeros);
  vec s2_beta_1_save(n_save, fill::zeros);
  mat s2_save(num_species, n_save, fill::zeros);
  mat Temp_min_save(num_species, n_save, fill::zeros);
  mat Temp_max_save(num_species, n_save, fill::zeros);
  mat P_min_save(num_species, n_save, fill::zeros);
  mat P_max_save(num_species, n_save, fill::zeros);
  cube W_save(24, t, n_save, fill::zeros);
  vec Temp_min_accept(num_species, fill::zeros);
  vec Temp_max_accept(num_species, fill::zeros);
  vec P_min_accept(num_species, fill::zeros);
  vec P_max_accept(num_species, fill::zeros);
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
  
  Rcout << "\n\n" << "Starting Hierarchical VS-Lite model fit, will run for " << n_mcmc <<
            " iterations\n\n";
  for(int k = 0; k < n_mcmc; k++){
    if((k + 1) % 100 == 0){
      Rcout << " " << k + 1;
    }
    
    Rcpp::checkUserInterrupt();
    
    //
    // sample W_T
    //
    
    vec like_T = makeLikelihood(y, H, beta_0, beta_1, zeta, t, p, 
                                 species_m1, s2);
    // sample for t = t-N_T to 2
    for(int i = t - N_T - 1; i > 0; i--){
      vec W_tilde_star_T = W_tilde.col(i);
      vec W_star_vec_T = W.col(i);
      for(int l = 0; l < 12; l++){
        W_tilde_star_T(l)  =  W_tilde_star_T(l) + R::rnorm(0, W_T_tune);
        W_star_vec_T(l) = W_tilde_star_T(l) * s_T(l) + mu_T(l);
      }
      vec zeta_star_T = makeZetaIndividualVS(p, day_len, W_star_vec_T, Temp_min,
                                           Temp_max, P_min, P_max, species_m1);
      vec tmp1_star_T(24);
      vec tmp2_star_T(24);
      vec tmp1_T(24);
      vec tmp2_T(24);
      if (i == t - N_T - 1){
        tmp1_star_T = (W_tilde_star_T - trend_0 * J - trend_1 * J) - phi_vec %
                       (W_tilde.col(i - 1) - trend_0 * J  - trend_1 * J);
        tmp2_star_T = (W_tilde.col(i + 1) - trend_0 * J - trend_1 * J) - 
                       phi_vec % (W_tilde_star_T - trend_0 * J - trend_1 * J);
        tmp1_T = (W_tilde.col(i) - trend_0 * J - trend_1 * J) - phi_vec % 
                 (W_tilde.col(i - 1) - trend_0 * J - trend_1 * J);
        tmp2_T = (W_tilde.col(i + 1) - trend_0 * J - trend_1 * J) - phi_vec % 
                 (W_tilde.col(i) - trend_0 * J - trend_1 * J);
      } else {
        tmp1_star_T = (W_tilde_star_T - trend_0 * J) - phi_vec % 
                       (W_tilde.col(i - 1) - trend_0 * J);
        tmp2_star_T = (W_tilde.col(i + 1) - trend_0 * J) - phi_vec % 
                       (W_tilde_star_T - trend_0 * J);
        tmp1_T = (W_tilde.col(i) - trend_0 * J) - phi_vec % 
                 (W_tilde.col(i - 1) - trend_0 * J);
        tmp2_T = (W_tilde.col(i + 1) - trend_0 * J) - phi_vec % 
                 (W_tilde.col(i) - trend_0 * J);
      }
      vec mh1_T = - 0.5 * trans(tmp1_star_T) * Sigma_inv * tmp1_star_T - 
                  0.5 * trans(tmp2_star_T) * Sigma_inv * tmp2_star_T + 
                  makeLikelihoodIndividual(y, H, beta_0, beta_1, zeta_star_T,
                                             i, p, species_m1, s2);
      vec mh2_T = - 0.5 * trans(tmp1_T) * Sigma_inv * tmp1_T - 
                  0.5 * trans(tmp2_T) * Sigma_inv * tmp2_T + like_T(i);
      double mh_T = exp(mh1_T(0) - mh2_T(0));
      if(mh_T >   R::runif(0, 1)){
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
    vec zeta_star_T = makeZetaIndividualVS(p, day_len, W_star_vec_T, Temp_min, 
                                         Temp_max, P_min, P_max, species_m1);
    vec tmp_star_T = (W_tilde.col(1) - trend_0 * J) - 
                       phi_vec % (W_tilde_star_T - trend_0 * J);
    int idx_val_T = 0;
    vec mh1_T = - 0.5 * tmp_star_T.t() * Sigma_inv * tmp_star_T + 
                 makeLikelihoodIndividual(y, H, beta_0, beta_1, zeta_star_T,
                                           idx_val_T, p, species_m1, s2);
    vec tmp_T = (W_tilde.col(1) - trend_0 * J) - phi_vec % 
                 (W_tilde.col(0) - trend_0 * J);
    vec mh2_T = - 0.5 * tmp_T.t() * Sigma_inv * tmp_T + like_T(0);
    double mh_T = exp(mh1_T(0) - mh2_T(0));
    if(mh_T >   R::runif(0, 1)){
      W_tilde.col(0) = W_tilde_star_T;
      W.col(0) = W_star_vec_T;
      if(k > n_burn){
          W_T_accept += 1.0  / (n_burn * (t - N_T));
        }
        W_T_accept_tmp += 1.0  / (50 * (t - N_T));
    } 
    zeta = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max,
                    species_m1);
    
    //
    // sample W_P
    //
    
    vec like_P = makeLikelihood(y, H, beta_0, beta_1, zeta, t, p, 
                                species_m1, s2);
    // sample for t = t-N_T to 2
    for(int i = t - N_T - 1; i > 0; i--){
      vec W_tilde_star = W_tilde.col(i);
      vec W_star_vec = W.col(i);
      for(int l = 12; l < 24; l++){
        W_tilde_star(l)  =  W_tilde_star(l) + R::rnorm(0, W_P_tune);
        W_star_vec(l) = W_tilde_star(l) * s_P(l - 12) + mu_P(l - 12);
      }
      vec zeta_star_P = makeZetaIndividualVS(p, day_len, W_star_vec, Temp_min, 
                                           Temp_max, P_min, P_max, species_m1);
      vec tmp1_star(24);
      vec tmp2_star(24);
      tmp1_star = W_tilde_star - phi_vec % W_tilde.col(i - 1);
      tmp2_star = W_tilde.col(i + 1) - phi_vec % W_tilde_star;
      vec mh1 = - 0.5 * trans(tmp1_star) * Sigma_inv * tmp1_star - 
                0.5 * trans(tmp2_star) * Sigma_inv * tmp2_star + 
                makeLikelihoodIndividual(y, H, beta_0, beta_1, 
                                           zeta_star_P, i, p, 
                                           species_m1, s2);
      vec tmp1 = W_tilde.col(i) - phi_vec % W_tilde.col(i - 1);
      vec tmp2 = W_tilde.col(i + 1) - phi_vec % W_tilde.col(i);
      vec mh2 = - 0.5 * trans(tmp1) * Sigma_inv * tmp1 - 
                0.5 * trans(tmp2) * Sigma_inv * tmp2 + like_P(i);
      double mh = exp(mh1(0) - mh2(0));
      if(mh >   R::runif(0, 1)){
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
    vec zeta_star_P = makeZetaIndividualVS(p, day_len, W_star_vec, Temp_min,
                                         Temp_max, P_min, P_max, species_m1);
    vec tmp_star = W_tilde.col(1) - phi_vec % W_tilde_star;
    int idx_val = 0;
    vec mh1 = - 0.5 * tmp_star.t() * Sigma_inv * tmp_star +
              makeLikelihoodIndividual(y, H, beta_0, beta_1, zeta_star_P, 
              idx_val, p, species_m1, s2);
    vec tmp = W_tilde.col(1) - phi_vec % W_tilde.col(0);
    vec mh2 = - 0.5 * tmp.t() * Sigma_inv * tmp + like_P(0);
    double mh = exp(mh1(0) - mh2(0));
    if(mh >   R::runif(0, 1)){
      W_tilde.col(0) = W_tilde_star;
      W.col(0) = W_star_vec;
      if(k > n_burn){
          W_P_accept += 1.0  / (n_burn * (t - N_P));
        }
        W_P_accept_tmp += 1.0  / (50 * (t - N_P));
    } 
    zeta = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max,
                    species_m1);
    
    //
    // sample beta_0
    //
    
    mat A0(num_species, num_species, fill::zeros);
    vec A_vec0(num_species, fill::zeros);
    vec b0 = makeBeta0(y, H, beta_1, zeta, species_m1, t, p, num_species, s2);
    for(int i = 0; i < num_species; i++){
      A_vec0(i) = N_obs_species(i) / s2(i) + 1.0 / s2_beta_0;
      b0(i) += mu_beta_0 / s2_beta_0;
    }
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
    vec A_vec1 = makeBeta1Var(H, zeta, species_m1, t, p, num_species, s2);
    vec b1 = makeBeta1Mean(y, H, beta_0, zeta, species_m1, t, p, num_species, s2);
    for(int i = 0; i < num_species; i++){
      A_vec1(i) += 1.0 / s2_beta_1;
      b1(i) += mu_beta_1 / s2_beta_1;
    }
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
    
    mat zeta_star = makeZetaVS(t, p, day_len, W, Temp_min_star, Temp_max, P_min, 
                               P_max, species_m1);
    vec like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta_star,
                                          species_m1, t, p, num_species, s2);
    vec like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, species_m1, 
                                      t, p, num_species, s2);
    for(int i = 0; i < num_species; i++){
      double x_star = (Temp_min_star(i) - Temp_min_lower) / 
                      (Temp_min_upper - Temp_min_lower);
      double x = (Temp_min(i) - Temp_min_lower) / 
                  (Temp_min_upper - Temp_min_lower);
      double mh1 = like_star(i) + // conditional likelihood under proposed value
                    // density of proposal under prior
                    (alpha_Temp_min - 1.0) * log(x_star) + 
                    (beta_Temp_min - 1.0) * log(1.0 - x_star);
      double mh2 = like(i) + // conditional likelihood under current value
                    // density of current value under prior
                    (alpha_Temp_min - 1.0) * log(x) + 
                    (beta_Temp_min - 1.0) * log(1.0 - x);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
        Temp_min(i) = Temp_min_star(i);
        if(k > n_burn){
        	Temp_min_accept(i) += 1.0 / n_burn;
        }
       	Temp_min_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max, 
        species_m1);
    
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
    
    zeta_star = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max_star, P_min, 
                           P_max, species_m1);
    like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta_star, 
                                      species_m1, t, p, num_species, s2);
    like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, species_m1, t, p,
                                  num_species, s2);
    
    for(int i = 0; i < num_species; i++){
      double x_star = (Temp_max_star(i) - Temp_max_lower) / 
                      (Temp_max_upper - Temp_max_lower);
      double x = (Temp_max(i) - Temp_max_lower) / 
                  (Temp_max_upper - Temp_max_lower);
      double mh1 = like_star(i) + // conditional likelihood under proposed value
                    // density of proposal under prior
                    (alpha_Temp_max - 1.0) * log(x_star) + 
                    (beta_Temp_max - 1.0) * log(1.0 - x_star);                
      double mh2 = like(i) + // conditional likelihood under proposed value
                    // density of proposal under prior
                    (alpha_Temp_max - 1.0) * log(x) + 
                    (beta_Temp_max - 1.0) * log(1.0 - x);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
        Temp_max(i) = Temp_max_star(i);
        if(k > n_burn){
        	Temp_max_accept(i) += 1.0 / n_burn;
        }
       	Temp_max_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max, 
        species_m1);
    
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
    zeta_star = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min_star, 
                           P_max, species_m1);
    like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta_star, 
                                      species_m1, t, p, num_species, s2);
    like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, species_m1, t, p,
                                  num_species, s2);
    
    for(int i = 0; i < num_species; i++){
      double x_star = (P_min_star(i) - P_min_lower) / 
                       (P_min_upper - P_min_lower);
      double x = (P_min(i) - P_min_lower) / 
                  (P_min_upper - P_min_lower);
      double mh1 = like_star(i) + (alpha_P_min - 1.0) * 
  log(x_star) + (beta_P_min - 1.0) * 
  log(1.0 - x_star);             
      double mh2 = like(i) + (alpha_P_min - 1.0) * 
  log(x) + (beta_P_min - 1.0) * 
  log(1.0 - x);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
  P_min(i) = P_min_star(i);
  if(k > n_burn){
        	P_min_accept(i) += 1.0 / n_burn;
        }
       	P_min_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max, 
        species_m1);
    
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
    
    zeta_star = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, 
       P_max_star, species_m1);
    like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1,
            zeta_star, species_m1, t, p,
            num_species, s2);
    like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, 
             species_m1, t, p, num_species, s2);
    for(int i = 0; i < num_species; i++){
      double x_star = (P_max_star(i) - P_max_lower) / 
  (P_max_upper - P_max_lower);
      double x = (P_max(i) - P_max_lower) / 
  (P_max_upper - P_max_lower);
      double mh1 = like_star(i) + (alpha_P_max - 1.0) * 
  log(x_star) + (beta_P_max - 1.0) * 
  log(1.0 - x_star);             
      double mh2 = like(i) + (alpha_P_max - 1.0) * log(x) + 
                  (beta_P_max - 1.0) * log(1.0 - x);   
      double mh = exp(mh1 - mh2);    
      if(mh > R::runif(0, 1)){
        P_max(i) = P_max_star(i);
        if(k > n_burn){
        	P_max_accept(i) += 1.0 / n_burn;
        }
       	P_max_accept_tmp(i) += 1.0 / 50;
      }
    }
    zeta = makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max, 
                     species_m1);
    
    //
    // sample s2
    //
    
		vec s2_beta_vec = makeS2Beta(y, H, beta_0, beta_1, zeta, species_m1, t, p, 
		                             num_species);
		for(int i = 0; i < num_species; i++){
			s2(i) = 1.0 / R::rgamma(alpha_s2 + 0.5 * N_obs_species(i), 
			                        1.0 / (beta_s2 + 0.5 * s2_beta_vec(i)));
		}

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
        // Adaptive tuning for W_T
        updateTuning(k, W_T_accept_tmp, W_T_tune);
        // Adaptive tuning for W_P
        updateTuning(k, W_P_accept_tmp, W_P_tune);
    }

		//
		// save variables
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
				s2_save.col(k_tmp) = s2;
        Temp_min_save.col(k_tmp) = Temp_min;
        Temp_max_save.col(k_tmp) = Temp_max;
        P_min_save.col(k_tmp) = P_min;
        P_max_save.col(k_tmp) = P_max;
        W_save.slice(k_tmp) = W;
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
						      _["s2"] = s2_save),
          _["Temp_min"] = Temp_min_save,
          _["Temp_max"] = Temp_max_save,
          _["P_min"] = P_min_save,
          _["P_max"] = P_max_save,
          _["W"] = W_save,
          _["accept"] = List::create(
            _["Temp_min_accept"] = Temp_min_accept,
            _["Temp_max_accept"] = Temp_max_accept,
            _["P_min_accept"] = P_min_accept,
            _["P_max_accept"] = P_max_accept,
            _["W_T_accept"] = W_T_accept,
            _["W_P_accept"] = W_P_accept
            	)
          ));
}
