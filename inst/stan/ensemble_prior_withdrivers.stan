functions{
  /**
   * Rescaling [-1,1] -> [0,1] to make the Beta distribution applicable for correlations
   */
  real As(real Rho){
    return 1/pi() * atan(Rho/sqrt(1-Rho^2))+0.5;
  }

  /**
   * Beta prior for correlation matrices
   */
  real priors_cor_beta(matrix Rho, int N, matrix beta_1, matrix beta_2) {
    real log_prior;
    log_prior =0;
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        log_prior += beta_lpdf(As(Rho[i,j])|beta_1[i,j],beta_2[i,j]);
      }
    }
    return log_prior;
  }


  /**
   * Hierarchical priors...
   */
  real priors_cor_hierarchical_beta(matrix ind_st_cor, int N, matrix M){
    real log_prior = 0;
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        log_prior += beta_lpdf(0.5*(ind_st_cor[i, j] + 1)| M[i, j], M[j, i] );
      }
    }
    return log_prior;
  }

      /**
   * Extra functions...
   */


  int sq_int(int[] model_num_species, int M){
    int ret = 0;
    for (i in 1:M){
	    ret += model_num_species[i] * model_num_species[i];
	  }
	  return ret;
  }

  int which_fun_subset(int i, int j, int [,] mod_to_dri){
    int ret = 0;
    int num = 0;
    while (num < j){
      ret += 1;
      if (mod_to_dri[i,ret] == ret){
        num += 1;
      }
    }
    return ret;
  }
}



data{
  int <lower=0> N;   // Number of variables
  //int <lower=0> time;// How long the model is run for

  /**
   * Observations
   */
  //array[time] vector [N] observations;
  //cov_matrix [N] obs_covariances;

  /**
  * Simulators and drivers
  */
  int<lower=0> M; // number of models
  int<lower=0> MM; // number of drivers
  //array[M] int<lower=0> n_dri; // number of drivers per simulator
  //array[M, MM] int<lower=0> mod_to_dri; //array assigning driver index to each simulator. Entries equal to zero corresponding to non-assigned drivers.
  //array[sum(n_dri)] int<lower=0> model_num_species;
  //matrix [time,M+1] observation_times; // times that observations happen: first column is the observations the rest is for the simulators and drivers.

  //matrix[sum(model_num_species), N] Ms;
  //matrix [time, sum(model_num_species)] model_outputs;
  //vector [sq_int(model_num_species, sum(n_dri))] model_covariances;

	/**
	 * Prior choice paramters:
	 * These determine the form of the prior parametrisation using an integar representation.
	 * Choices 0,1,2 use a decomposition into a diagonal variance matrix and a correlation matrix,
	 * with inverse-gamma distributions on the variance terms.
	 *      0 - LKJ correlation matrix
	 *      1 - Inverse Wishart correlation matrix
	 *      2 - Beta distributions on correlation matrix entries.
	 *      ONLY IMPLEMENTED FOR SHORT-TERM: 3 - Hierarchical beta priors
	 *
	 */
	 int<lower=0, upper=3> form_prior_ind_st;
	 int<lower=0, upper=2> form_prior_ind_lt;
	 int<lower=0, upper=2> form_prior_sha_st;


  /**
   * Prior parameters
   */
  //Individual short-term
  vector [N] prior_ind_st_var_a; // shape parameter (alpha) of inverse gamma
	vector [N] prior_ind_st_var_b; // scale parameter (beta) of inverse gamma
	array[form_prior_ind_st == 0 ? 1 : 0] real prior_ind_st_cor_lkj; // LKJ shape parameter
  matrix[form_prior_ind_st == 1 ? N : 0, form_prior_ind_st == 1 ? N : 0] prior_ind_st_cor_wish_sigma;//inverse wishart
	array[form_prior_ind_st == 1 ? 1 : 0] real<lower=N-1>	prior_ind_st_cor_wish_nu; //inverse wishart
	matrix [form_prior_ind_st == 2 ? N : 0, form_prior_ind_st == 2 ? N : 0] prior_ind_st_cor_beta_1; // alpha shape parameter for Beta distribution
  matrix [form_prior_ind_st == 2 ? N : 0, form_prior_ind_st == 2 ? N : 0] prior_ind_st_cor_beta_2; // beta shape parameter for Beta distribution

  //JM 06/05: Adding in hierarchical prior options.
  //Correlation: Each correlation matrix element is a Beta(a, b) with a~ Gamma(l,m), b ~ Gamma(n, o) The elements of this vector are: (1) l (2) m (3) n (4) o
  //Variance: This is done via variance ~ Gamma(a, b) with a~ Gamma(l,m), b ~ Gamma(n, o) The elements of this vector are: (1) l (2) m (3) n (4) o
  vector [form_prior_ind_st == 3 ? 4 : 0] prior_ind_st_cor_hierarchical_beta_hyper_params;
  vector [form_prior_ind_st == 3 ? 4 : 0] prior_ind_st_var_hierarchical_hyperparams;

  //JM 22/07 Now have beta priors on the AR parameters
  real<lower=0> prior_ind_st_ar_alpha;
  real<lower=0> prior_ind_st_ar_beta;


  //Individual long-term
  vector [N] prior_ind_lt_var_a ; // shape parameter (alpha) of inverse gamma
  vector [N] prior_ind_lt_var_b ; // scale parameter (beta) of inverse gamma
  array[form_prior_ind_lt == 0 ? 1 : 0] real prior_ind_lt_cor_lkj; // LKJ shape parameter
  matrix[form_prior_ind_lt == 1 ? N : 0, form_prior_ind_lt == 1 ? N : 0] prior_ind_lt_cor_wish_sigma;//inverse wishart
  array[form_prior_ind_lt == 1 ? 1 : 0] real<lower=N-1>	prior_ind_lt_cor_wish_nu; //inverse wishart
  matrix [form_prior_ind_lt == 2 ? N : 0, form_prior_ind_lt == 2 ? N : 0] prior_ind_lt_cor_beta_1; // alpha shape parameter for Beta distribution
  matrix [form_prior_ind_lt == 2 ? N : 0, form_prior_ind_lt == 2 ? N : 0] prior_ind_lt_cor_beta_2; // beta shape parameter for Beta distribution

  //Shared short-term
	//real<lower=0> prior_sha_st_var_exp; // Scale parameter of exponential
	vector [N] prior_sha_st_var_a ; // shape parameter (alpha) of inverse gamma
  vector [N] prior_sha_st_var_b ; // scale parameter (beta) of inverse gamma
	array[form_prior_sha_st == 0 ? 1: 0] real prior_sha_st_cor_lkj; // LKJ shape parameter
	matrix[form_prior_sha_st == 1 ? N: 0,form_prior_sha_st == 1 ? N: 0] prior_sha_st_cor_wish_sigma;//inverse wishart
	array[form_prior_sha_st == 1 ? 1: 0] real<lower=N-1>	prior_sha_st_cor_wish_nu; //inverse wishart
	matrix [form_prior_sha_st == 2 ? N: 0,form_prior_sha_st == 2 ? N: 0] prior_sha_st_cor_beta_1; // alpha shape parameter for Beta distribution
  matrix [form_prior_sha_st == 2 ? N: 0,form_prior_sha_st == 2 ? N: 0] prior_sha_st_cor_beta_2; // beta shape parameter for Beta distribution
  //JM 22/07 Now have beta priors on the AR parameters
  real<lower=0> prior_sha_st_ar_alpha;
  real<lower=0> prior_sha_st_ar_beta;

 //Shared long-term
	vector <lower=0> [N] prior_sha_lt_sd; //sd for prior on error

  //Random walk on y
  vector <lower=0> [N] prior_y_init_mean;
  vector [N] prior_y_init_var;
  real<lower=N-1>	prior_sigma_t_inv_wish_nu; //inverse wishart
  matrix[N, N] prior_sigma_t_inv_wish_sigma;//inverse wishart

}

parameters{
  /**
   * Simulator discrepancies
   */

  // individual driver
  array [MM] vector <lower=0,upper=1>[N] ind_st_ar_param_dri_raw;
  //array [MM] vector [N] log_ind_st_var_dri;
  array [MM] vector<lower=0>[N] ind_st_var_dri;
  array [MM] corr_matrix [N] ind_st_cor_dri;
  array [MM] vector[N] ind_lt_raw_dri;
  vector <lower=0> [N] ind_lt_var_dri;
  corr_matrix [N] ind_lt_cor_dri;

  // Individual simulator
  array[M] vector <lower=0,upper=1>[N] ind_st_ar_param_raw;
  //array[M] vector [N] log_ind_st_var;
  array [M] vector<lower=0>[N] ind_st_var;
  array[M] corr_matrix [N] ind_st_cor;
  array[M] vector[N] ind_lt_raw;
  vector <lower=0> [N] ind_lt_var;
  corr_matrix [N] ind_lt_cor;
  // Shared
  vector <lower=0,upper=1>[N] sha_st_ar_param_raw;
  vector <lower=0> [N] sha_st_var;
  corr_matrix [N] sha_st_cor;
  vector [N] sha_lt_raw;

  /**
   * Random walk on y
   */
  cov_matrix [N] SIGMA_t;

  /**
   * Hierarchical parameters
   *
   * Correlations:
   * The hierarchical beta priors on individual short term discrepancy correlations are encoded by a matrix M and a vector v. In this scheme, the correlations c[i,j] have c[i, j] ~ Beta(a[i, j], b[i, j]). Since correlation matrices are symmetric we need only care about one section of the matrix when making calculations. So the matrix M has:
   * M[i, j] = a[i, j] for i < j
   * M[i, j] = b[i, j] for i > j
   */
   	//matrix<lower = 0>[N , N ] prior_ind_st_cor_hierarchical_beta_params; // a shape parameters for Beta distribution
    //vector [N ] prior_ind_st_var_hierarchical_mu_params;

	//matrix<lower = 0>[N , N ] prior_ind_st_cor_hierarchical_beta_params_dri; // a shape parameters for Beta distribution
    //vector [N ] prior_ind_st_var_hierarchical_mu_params_dri;

}
transformed parameters{
  array[M] matrix [N,N] SIGMA_x;
  array[M] vector [N] ind_st_sd;
  vector [N] sha_lt = prior_sha_lt_sd .* sha_lt_raw;
  array[M] vector [N] ind_lt;
  vector [N] ind_lt_sd = sqrt(ind_lt_var);
  matrix [N,N] ind_lt_covar = diag_post_multiply(diag_pre_multiply(ind_lt_sd,ind_lt_cor),ind_lt_sd);
  matrix [N,N] ind_lt_cov_cholesky = cholesky_decompose(ind_lt_covar);
  // Drivers
  array[MM] matrix [N,N] SIGMA_x_dri;
  array[MM] vector [N] ind_st_sd_dri;
  array[MM] vector [N] ind_lt_dri;
  vector [N] ind_lt_sd_dri = sqrt(ind_lt_var_dri);
  matrix [N,N] ind_lt_covar_dri = diag_post_multiply(diag_pre_multiply(ind_lt_sd_dri,ind_lt_cor_dri),ind_lt_sd_dri);
  matrix [N,N] ind_lt_cov_cholesky_dri = cholesky_decompose(ind_lt_covar_dri);

  array [MM] vector <lower=-1,upper=1>[N] ind_st_ar_param_dri;
  array[M] vector <lower=-1,upper=1>[N] ind_st_ar_param;
  vector <lower=-1,upper=1>[N] sha_st_ar_param = 2.0 * sha_st_ar_param_raw - 1.0;

  //JM 25-08-2022: This is the hierarchical version, we have to do parametrisations after it.
  //array[M] vector <lower=0>[N] ind_st_var_2;
  //array[MM] vector <lower=0>[N] ind_st_var_2_dri;

  vector [(MM + M+2) * N] x_hat = append_row(prior_y_init_mean, rep_vector(0.0, N * (MM + M + 1))); // added columns for drivers
  matrix [(MM + M+2) * N, (MM + M+2) * N] SIGMA_init = rep_matrix(0,(MM + M+2) * N,(MM + M+2) * N ); // The precision matrix will have a lot of zeros

  vector [N] sha_st_sd = sqrt(sha_st_var);
  matrix [N,N] SIGMA_mu = diag_post_multiply(diag_pre_multiply(sha_st_sd, sha_st_cor), sha_st_sd);

  /**
  *  Kalman filter parameters:
  *  In each case, the ordering to be passed through to the Kalman Filter is:
  *  (1) Observations (2) Consensus (3) Models (4) Drivers
  *     [if (1) is available, otherwise this is ignored]
  *
  */
  matrix[(MM + M + 2) * N , (MM + M + 2) * N] SIGMA = rep_matrix(0,(MM + M + 2) * N,(MM + M + 2) * N); // a lot of zeros
  vector[(MM + M + 2) * N ] lt_discrepancies;
  vector[(MM + M + 2) * N] AR_params;

  //SIGMA

  //SIGMA[1:N, 1:N ] = SIGMA_t;
  //SIGMA[(N + 1):(2*N), (N + 1):(2*N) ] = SIGMA_mu;
  //for (i in 1:M){
    //ind_st_var_2[i] = exp(prior_ind_st_var_hierarchical_mu_params + sqrt(diagonal(prior_ind_st_cor_hierarchical_beta_params)) .* log_ind_st_var[i]);
    //ind_st_sd[i] = sqrt(ind_st_var_2[i]);
    //SIGMA_x[i] = diag_post_multiply(diag_pre_multiply(ind_st_sd[i],ind_st_cor[i]),ind_st_sd[i]);
	  //SIGMA[((i+1) * N + 1):((i+2)*N ),((i + 1) * N + 1):((i+2)*N )] = SIGMA_x[i];
  //}
  //for (i in 1:MM){
	  //ind_st_var_2_dri[i] = exp(prior_ind_st_var_hierarchical_mu_params_dri + sqrt(diagonal(prior_ind_st_cor_hierarchical_beta_params_dri)) .* log_ind_st_var_dri[i]);
	  //ind_st_sd_dri[i] = sqrt(ind_st_var_2_dri[i]);
	  //SIGMA_x_dri[i] = diag_post_multiply(diag_pre_multiply(ind_st_sd_dri[i],ind_st_cor_dri[i]),ind_st_sd_dri[i]);
	  //SIGMA[((i+ M + 1) * N + 1):((i + M +2)*N ),((i + M + 1) * N + 1):((i + M + 2)*N )] = SIGMA_x_dri[i];
  //}

  SIGMA[1:N, 1:N ] = SIGMA_t;
  SIGMA[(N + 1):(2*N), (N + 1):(2*N) ] = SIGMA_mu;
  for (i in 1:M){
    ind_st_sd[i] = sqrt(ind_st_var[i]);
    SIGMA_x[i] = diag_post_multiply(diag_pre_multiply(ind_st_sd[i],ind_st_cor[i]),ind_st_sd[i]);
	  SIGMA[((i+1) * N + 1):((i+2)*N ),((i + 1) * N + 1):((i+2)*N )] = SIGMA_x[i];
  }
  for (i in 1:MM){
    ind_st_sd_dri[i] = sqrt(ind_st_var_dri[i]);
    SIGMA_x_dri[i] = diag_post_multiply(diag_pre_multiply(ind_st_sd_dri[i],ind_st_cor_dri[i]),ind_st_sd_dri[i]);
	  SIGMA[((i+ M + 1) * N + 1):((i + M +2)*N ),((i + M + 1) * N + 1):((i + M + 2)*N )] = SIGMA_x_dri[i];
  }

    //

  //SIGMA_init
  SIGMA_init[1:N,1:N] = diag_matrix(prior_y_init_var);
  SIGMA_init[(N + 1):(2*N), (N + 1):(2*N) ] = SIGMA_mu ./ (1 - sha_st_ar_param * sha_st_ar_param');
  for (i in 1:M){
    ind_st_ar_param[i] = 2.0 * ind_st_ar_param_raw[i] - 1.0;
    SIGMA_init[((i+1) * N + 1):((i+2)*N ),((i+1) * N + 1):((i+2)*N )] = SIGMA_x[i] ./ (1 - ind_st_ar_param[i] * ind_st_ar_param[i]');
  }
  for (i in 1:MM){
    ind_st_ar_param_dri[i] = 2.0 * ind_st_ar_param_dri_raw[i] - 1.0;
    SIGMA_init[((i + M + 1) * N + 1):((i + M + 2)*N ),((i + M + 1) * N + 1):((i + M + 2)*N )] = SIGMA_x_dri[i] ./ (1 - ind_st_ar_param_dri[i] * ind_st_ar_param_dri[i]');
  }


  lt_discrepancies[1:(2 * N)] = append_row(rep_vector(0.0,N), sha_lt);
  AR_params[1:(2 * N)] = append_row(rep_vector(1.0,N), sha_st_ar_param);
  for (i in 1:M){
    ind_lt[i] = ind_lt_cov_cholesky*ind_lt_raw[i];
    lt_discrepancies[((i+1) * N + 1):((i+2)*N )] = ind_lt[i];
	  AR_params[((i+1) * N + 1):((i+2)*N )] = ind_st_ar_param[i];
  }
  for (i in 1:MM){
    ind_lt_dri[i] = ind_lt_cov_cholesky_dri*ind_lt_raw_dri[i];
    lt_discrepancies[((i+ M + 1) * N + 1):((i + M + 2)*N )] = ind_lt_dri[i];
	  AR_params[((i + M + 1) * N + 1):((i + M + 2)*N )] = ind_st_ar_param_dri[i];
  }
}
model{
  /**
  * Priors
  */
  //Random walk on y
  SIGMA_t ~ inv_wishart(prior_sigma_t_inv_wish_nu, prior_sigma_t_inv_wish_sigma); // the random walk of y //Same MS


  // Shared discrepancies
  sha_lt_raw ~ std_normal();
  sha_st_var ~ gamma(prior_sha_st_var_a,prior_sha_st_var_b); // Variance
  //JM 22/07: Beta priors on the AR parameters
  target += beta_lpdf((sha_st_ar_param + 1)/2 | prior_sha_st_ar_alpha, prior_sha_st_ar_beta);  // same

  // Correlation matrix
  if(form_prior_sha_st == 0){
    sha_st_cor ~ lkj_corr(prior_sha_st_cor_lkj[1]);
  } else if(form_prior_sha_st == 1){
    sha_st_cor ~ inv_wishart(prior_sha_st_cor_wish_nu[1], prior_sha_st_cor_wish_sigma);
  } else {
    target += priors_cor_beta(sha_st_cor, N, prior_sha_st_cor_beta_1, prior_sha_st_cor_beta_2);
  }


  // Individual discrepancies
  // Note that we're assuming long-term discrepancies are drawn from a N(0,C) distribution
  // where C is independent of the simulators. This means we treat C outside the for loop.
  ind_lt_var ~ gamma(prior_ind_lt_var_a,prior_ind_lt_var_b); // Variance
  ind_lt_var_dri ~ gamma(prior_ind_lt_var_a,prior_ind_lt_var_b);
  //Long term correlations
  if(form_prior_ind_lt == 0){  // recommend lkj_corr or beta
    ind_lt_cor ~ lkj_corr(prior_ind_lt_cor_lkj[1]);
  } else if(form_prior_ind_lt == 1){
    ind_lt_cor ~ inv_wishart(prior_ind_lt_cor_wish_nu[1], prior_ind_lt_cor_wish_sigma);
  } else{
    target += priors_cor_beta(ind_lt_cor, N, prior_ind_lt_cor_beta_1, prior_ind_lt_cor_beta_2);
  }

  //Long term correlations -- Drivers
  if(form_prior_ind_lt == 0){  // recommend lkj_corr or beta
    ind_lt_cor_dri ~ lkj_corr(prior_ind_lt_cor_lkj[1]);
  } else if(form_prior_ind_lt == 1){
    ind_lt_cor_dri ~ inv_wishart(prior_ind_lt_cor_wish_nu[1], prior_ind_lt_cor_wish_sigma);
  } else{
    target += priors_cor_beta(ind_lt_cor_dri, N, prior_ind_lt_cor_beta_1, prior_ind_lt_cor_beta_2);
  }


  for(i in 1:M){
    //AR Parameters
    target += beta_lpdf((ind_st_ar_param[i] + 1)/2 | prior_ind_st_ar_alpha, prior_ind_st_ar_beta);

    ind_st_var[i] ~ gamma(prior_ind_st_var_a, prior_ind_st_var_b);
    ind_lt_raw[i] ~ std_normal();

    //ind_lt_raw[i] ~ std_normal();

	  //log_ind_st_var[i] ~ std_normal();

    // Correlation matrix
  	//target += priors_cor_hierarchical_beta(ind_st_cor[i], N, prior_ind_st_cor_hierarchical_beta_params);

  	if (form_prior_ind_st == 0){
  	  ind_st_cor[i] ~ lkj_corr(prior_ind_st_cor_lkj[1]);
    } else if(form_prior_ind_st == 1){
  	  ind_st_cor[i] ~ inv_wishart(prior_ind_st_cor_wish_nu[1], prior_ind_st_cor_wish_sigma);
  	} else{
  	  target += priors_cor_beta(ind_st_cor[i], N, prior_ind_st_cor_beta_1, prior_ind_st_cor_beta_2);
  	}


  }


  for(i in 1:MM){
    //AR Parameters
    target += beta_lpdf((ind_st_ar_param_dri[i] + 1)/2 | prior_ind_st_ar_alpha, prior_ind_st_ar_beta);

    ind_st_var_dri[i] ~ gamma(prior_ind_st_var_a, prior_ind_st_var_b);
    ind_lt_raw_dri[i] ~ std_normal();

    //ind_lt_raw_dri[i] ~ std_normal();

	  //log_ind_st_var_dri[i] ~ std_normal();




    // Correlation matrix
	  //target += priors_cor_hierarchical_beta(ind_st_cor_dri[i], N, prior_ind_st_cor_hierarchical_beta_params_dri);

  	if (form_prior_ind_st == 0){
  	  ind_st_cor_dri[i] ~ lkj_corr(prior_ind_st_cor_lkj[1]);
    } else if(form_prior_ind_st == 1){
  	  ind_st_cor_dri[i] ~ inv_wishart(prior_ind_st_cor_wish_nu[1], prior_ind_st_cor_wish_sigma);
  	} else{
  	  target += priors_cor_beta(ind_st_cor_dri[i], N, prior_ind_st_cor_beta_1, prior_ind_st_cor_beta_2);
  	}


  }

}











