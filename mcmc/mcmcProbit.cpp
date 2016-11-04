#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, myFunctions, bayesTreeRing)]]
#include "mcmc.h"

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
  
  // priors for probit growth model
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
  for(int i = 0; i < p; i++){
    for(int k = 0; k < t; k++){
      if(y(k, i) == 0){
        H(k, i) = 0;
      } else {
        H(k, i) = 1;
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
     
  vec species_m1 = species - 1;
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
  double s2_beta_0 = 1.0 ;
//  / R::rgamma(alpha_s2_beta_0, 1.0 / beta_s2_beta_0);
  double mu_beta_1 = R::rnorm(mu_mu_beta_1, sqrt(s2_mu_beta_1));
  double s2_beta_1 = 1.0 ;
//  / R::rgamma(alpha_s2_beta_1, 1.0 / beta_s2_beta_1);
  
  vec beta_0(num_species);
  vec beta_1(num_species);
  vec s2(num_species); 
  for(int i = 0; i < num_species; i++){
    beta_0(i) = R::rnorm(mu_beta_0, sqrt(s2_beta_0));
    beta_1(i) = R::rnorm(mu_beta_1, sqrt(s2_beta_1));
    s2(i) = 1.0 ;
//    / R::rgamma(alpha_s2, 1.0 / beta_s2);
  }

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
  
  vec gamma_T(num_species);
  vec xi2_T(num_species);
  vec gamma_P(num_species);
  vec xi2_P(num_species);
  
  for(int i = 0; i < num_species; i++){
    gamma_T(i) = R::rnorm(mu_gamma_T, sqrt(s2_gamma_T));
    xi2_T(i) = R::rlnorm(mu_xi2_T, sqrt(s2_xi2_T));
    if(xi2_T(i) > 5){
    	xi2_T(i) = 5;
    }
    gamma_P(i) = R::rnorm(mu_gamma_P, sqrt(s2_gamma_P));
    xi2_P(i) = R::rlnorm(mu_xi2_P, sqrt(s2_xi2_P));
    if(xi2_P(i) > 5){
    	xi2_P(i) = 5;
    }
  }
  vec xi_T = sqrt(xi2_T);
  vec xi_P = sqrt(xi2_P);  

    // adaptive tuning
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
  
  mat zeta = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, species_m1);
  
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
  vec gamma_T_accept(num_species, fill::zeros);
  vec xi2_T_accept(num_species, fill::zeros);
  vec gamma_P_accept(num_species, fill::zeros);
  vec xi2_P_accept(num_species, fill::zeros);
  double W_T_accept = 0;
  double W_P_accept = 0;
  vec gamma_T_accept_tmp(num_species, fill::zeros);
  vec xi2_T_accept_tmp(num_species, fill::zeros);
  vec gamma_P_accept_tmp(num_species, fill::zeros);
  vec xi2_P_accept_tmp(num_species, fill::zeros);
  double W_T_accept_tmp = 0;
  double W_P_accept_tmp = 0;
  
  //
  // start MCMC chain
  //
  
  Rcout << "\n\n" << "Starting Hierarchcial Probit model fit with random growth parameters, will run for " << n_mcmc << " iterations\n\n";
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
      vec zeta_star_T = makeZetaIndividualPro(p, day_len, W_star_vec_T, gamma_T, 
					   xi_T, gamma_P, xi_P, species_m1);
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
    vec zeta_star_T = makeZetaIndividualPro(p, day_len, W_star_vec_T, gamma_T, 
                                         xi_T, gamma_P, xi_P, species_m1);
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
    if(mh_T > R::runif(0, 1)){
      W_tilde.col(0) = W_tilde_star_T;
      W.col(0) = W_star_vec_T;
      if(k > n_burn){
      	W_T_accept += 1.0  / (n_burn * (t - N_T));
      }
      W_T_accept_tmp += 1.0  / (50 * (t - N_T));
    } 
    zeta = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, 
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
      vec zeta_star_P = makeZetaIndividualPro(p, day_len, W_star_vec,       
					   gamma_T, xi_T, gamma_P, 
					   xi_P, species_m1);
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
    vec zeta_star_P = makeZetaIndividualPro(p, day_len, W_star_vec, gamma_T,
					 xi_T, gamma_P, xi_P, species_m1);
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
    zeta = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, 
                    species_m1);    
     
    //
    // sample beta_0
    //
    
    mat A0(num_species, num_species, fill::zeros);
    vec A_vec0(num_species);
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
    vec b1 = makeBeta1Mean(y, H, beta_0, zeta, species_m1, t, p, num_species,
			   s2);
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
    // sample probit paramters
    //
    
    // sample gamma_T
    vec gamma_T_star(num_species); 
    for(int i = 0; i < num_species; i++){
      gamma_T_star(i) = R::rnorm(gamma_T(i), gamma_T_tune(i));
    }
    mat zeta_star = makeZetaPro(t, p, day_len, W, gamma_T_star, xi_T, gamma_P,
			     xi_P, species_m1);
    vec like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta_star,
					  species_m1, t, p, num_species, s2);
    vec like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, species_m1,
				     t, p, num_species, s2);
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
    zeta = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, species_m1);
    
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
      zeta_star = makeZetaPro(t, p, day_len, W, gamma_T, xi_T_star, gamma_P, xi_P, 
			 species_m1);
    // growth model likelihood by species for proposed value
    like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta_star,
				      species_m1, t, p, num_species, s2);
    // growth model likelihood by species for current value
    like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, species_m1, t,
				 p, num_species, s2);
    for(int i = 0; i < num_species; i++){
    	if(xi2_T_star(i) > 0){
        double mh1 = like_star(i) + // conditional likelihood under proposal
//  	// likelihood of proposal under IG prior
//	  - (alpha_xi2_T + 1.0) * log(xi2_T_star(i)) - beta_xi2_T / xi2_T_star(i);
		  // likelihood of proposal under Lognormal prior
  		- log(xi2_T_star(i)) - 0.5 * pow(log(xi2_T_star(i)) - mu_xi2_T, 2) / s2_xi2_T;
        double mh2 = like(i) + // conditional likelihood under current value
//	  // likelihood of current value under IG prior
//  	- (alpha_xi2_T + 1.0) * log(xi2_T(i)) - beta_xi2_T / xi2_T(i);
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
    zeta = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, species_m1);

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
    zeta_star = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P_star, xi_P, 
			 species_m1);
    like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta_star,
				      species_m1, t, p, num_species, s2);
    like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, species_m1, t, p,
				 num_species, s2);
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
    zeta = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P,
		    species_m1);
    
    // sample mu_gamma_P
    double a_mu_gamma_P = num_species / s2_gamma_P + 1.0 / s2_mu_gamma_P;
    double b_mu_gamma_P = sum(gamma_P) / s2_gamma_P + mu_mu_gamma_P / s2_mu_gamma_P;
    mu_gamma_P = rMVNArmaScalar(a_mu_gamma_P, b_mu_gamma_P);
        
    // sample s2_gamma_P
    s2_gamma_P = 1.0 / R::rgamma(alpha_s2_gamma_P + 0.5 * num_species, 1.0 / 
				 (beta_s2_gamma_P + 0.5 * sum(pow(gamma_P - 
								  mu_gamma_P, 2))));
    
   // sample xi2_P
    vec xi2_P_star(num_species);
    vec xi_P_star = xi_P;
    for(int i = 0; i < num_species; i++){
      xi2_P_star(i) = R::rnorm(xi2_P(i), xi2_P_tune(i));
      if(xi2_P_star(i) > 0){
      	xi_P_star(i) = sqrt(xi2_P_star(i));
      }
    }
    zeta_star = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P_star, 
			 species_m1);
    // growth model likelihood by species for propose value
    like_star = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta_star,
				      species_m1, t, p, num_species, s2);
    // growth model likelihood by species for current value
    like = makeLikelihoodSpecies(y, H, beta_0, beta_1, zeta, species_m1, t, p,
				 num_species, s2);
    for(int i = 0; i < num_species; i++){
    	if(xi2_P_star(i) > 0){
        double mh1 = like_star(i) + // conditional likelihood under proposal
//  	// likelihood of proposal under IG prior
//	  - (alpha_xi2_P + 1.0) * log(xi2_P_star(i)) - beta_xi2_P / xi2_P_star(i);
  	  // likelihood of proposal under Lognormal prior
	    - log(xi2_P_star(i)) - 0.5 * pow(log(xi2_P_star(i)) - mu_xi2_P, 2) / s2_xi2_P;
        double mh2 = like(i) + // conditional likelihood under current value
//  	// likelihood of current value under IG prior
//	  - (alpha_xi2_P + 1.0) * log(xi2_P(i)) - beta_xi2_P / xi2_P(i);
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
    zeta = makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, species_m1);

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
    // save variables
    //
    
    if(k >= n_burn){
      if(k % n_thin == 0){
        int k_tmp = (k - n_burn) / n_thin;
        beta_0_save.col(k_tmp) = beta_0;
        mu_beta_0_save(k_tmp) = mu_beta_0;
        s2_beta_0_save(k_tmp) = s2_beta_0;
        beta_1_save.col(k_tmp) = beta_1;
        mu_beta_1_save(k_tmp) = mu_beta_1;
        s2_beta_1_save(k_tmp) = s2_beta_1;
        s2_save.col(k_tmp) = s2;
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
		      _["s2_xi2_P"] = s2_xi2_P_save,
		      _["W"] = W_save,
		      _["accept"] = List::create(
						 _["gamma_T_accept"] = gamma_T_accept,
						 _["xi2_T_accept"] = xi2_T_accept,
						 _["gamma_P_accept"] = gamma_P_accept,
						 _["xi2_P_accept"] = xi2_P_accept,
						 _["W_T_accept"] = W_T_accept,
						 _["W_P_accept"] = W_P_accept)
						 ));
}
