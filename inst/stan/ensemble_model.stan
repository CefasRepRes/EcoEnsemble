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
    for (i in 1:(N-1))  {
      for (j in (i+1):N){
        log_prior += beta_lpdf(As(Rho[i,j])|beta_1[i,j],beta_2[i,j]);
      }
    }
    return log_prior;
  }

  /**
   * The Kalman filter - likelihood
   */
  real KalmanFilter_seq_em(vector rhos, vector lt_discrepancies,
                           vector R, matrix Q, matrix C,// parameters
                           matrix P, vector xhat, int Time, matrix y,matrix obs){
  real log_like = 0;
  matrix [num_elements(xhat),num_elements(xhat)] P_ = P;
  matrix [num_elements(xhat),num_elements(xhat)] A = rhos * rhos';
  vector [num_elements(xhat)] xhat_=xhat;
  matrix [num_elements(xhat),cols(y)] G;
  matrix [num_elements(xhat),cols(y)] G_;
  real er;
  real est;
  real Qs_inv;
  for (i in 1:Time)
  {
    //Prediction step
    xhat_ = rhos .* xhat_; // State estimate
    P_ = P_ .* A + Q;      // Covariance estimate

    for (j in 1:cols(y))
    {
      if (obs[i,j]==1.0)
      {
        G_[,j] = P_ * C[j,]';
        Qs_inv = 1 / (C[j,] * G_[,j] + R[j]);

        est = dot_product((xhat_ + lt_discrepancies)' ,  C[j,]) ; // careful we are required to make a transformation in dee (the discrepency) space!! bigM * dee
        er = y[i,j] - est; // pre-fit residuals

        G[,j] = Qs_inv * G_[,j];
        xhat_ += G[,j] * (y[i,j] - est);
        P_ -=  G[,j] * G_[,j]';

        log_like -= 0.5 * (- log(Qs_inv) + er^2 * Qs_inv);
      }
    }
  }
  return log_like;
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

  int sq_int(int[] model_num_species, int M){
    int ret = 0;
    for (i in 1:M){
	    ret += model_num_species[i] * model_num_species[i];
	  }
	  return ret;
  }
}

data{
  int <lower=0> N;   // Number of variables
  int <lower=0> time;// How long the model is run for

  /**
   * Observations
   */
  vector [N] observations[time];
  cov_matrix [N] obs_covariances;

  /**
  * Simulators
  */
  int<lower=0> M; // number of models
  matrix [time,M+1] observation_times; // times that observations happen
  int<lower=0> model_num_species[M]; // might have to be one
  matrix[sum(model_num_species),N] Ms;  // the Ms -- this assumes that the observations are always the same
  matrix [time ,sum(model_num_species)] model_outputs; // vector of observations from the models
  vector [sq_int(model_num_species, M)] model_covariances; // vector of covariance matrices -- this needs to be checked that the individual matrices are positive definite!

	/**
	 * Prior choice paramters:
	 * These determine the form of the prior parametrisation using an integar representation.
	 * Choices 0,1,2 use a decomposition into a diagonal variance matrix and a correlation matrix,
	 * with inverse-gamma distributions on the variance terms.
	 *      0 - LKJ correlation matrix
	 *      1 - Inverse Wishart correlation matrix
	 *      2 - Beta distributions on correlation matrix entries.
	 *      ONLY IMPLEMENTED FOR SHORT-TERM: 3 - Hierarchical beta priors
	 *      NOT IMPLEMENTED: 4 - Inverse Wishart covariance matrix
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
	real prior_ind_st_cor_lkj[form_prior_ind_st == 0 ? 1 : 0]; // LKJ shape parameter
  matrix[form_prior_ind_st == 1 ? N : 0, form_prior_ind_st == 1 ? N : 0] prior_ind_st_cor_wish_sigma;//inverse wishart
	real<lower=N-1>	prior_ind_st_cor_wish_nu[form_prior_ind_st == 1 ? 1 : 0]; //inverse wishart
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
  real prior_ind_lt_cor_lkj[form_prior_ind_lt == 0 ? 1 : 0]; // LKJ shape parameter
  matrix[form_prior_ind_lt == 1 ? N : 0, form_prior_ind_lt == 1 ? N : 0] prior_ind_lt_cor_wish_sigma;//inverse wishart
	real<lower=N-1>	prior_ind_lt_cor_wish_nu[form_prior_ind_lt == 1 ? 1 : 0]; //inverse wishart
	matrix [form_prior_ind_lt == 2 ? N : 0, form_prior_ind_lt == 2 ? N : 0] prior_ind_lt_cor_beta_1; // alpha shape parameter for Beta distribution
  matrix [form_prior_ind_lt == 2 ? N : 0, form_prior_ind_lt == 2 ? N : 0] prior_ind_lt_cor_beta_2; // beta shape parameter for Beta distribution


	//Shared short-term
	//real<lower=0> prior_sha_st_var_exp; // Scale parameter of exponential
	vector [N] prior_sha_st_var_a ; // shape parameter (alpha) of inverse gamma
  vector [N] prior_sha_st_var_b ; // scale parameter (beta) of inverse gamma
	real prior_sha_st_cor_lkj[form_prior_sha_st == 0 ? 1: 0]; // LKJ shape parameter
	matrix[form_prior_sha_st == 1 ? N: 0,form_prior_sha_st == 1 ? N: 0] prior_sha_st_cor_wish_sigma;//inverse wishart
	real<lower=N-1>	prior_sha_st_cor_wish_nu[form_prior_sha_st == 1 ? 1: 0]; //inverse wishart
	matrix [form_prior_sha_st == 2 ? N: 0,form_prior_sha_st == 2 ? N: 0] prior_sha_st_cor_beta_1; // alpha shape parameter for Beta distribution
  matrix [form_prior_sha_st == 2 ? N: 0,form_prior_sha_st == 2 ? N: 0] prior_sha_st_cor_beta_2; // beta shape parameter for Beta distribution
  //JM 22/07 Now have beta priors on the AR parameters
  real<lower=0> prior_sha_st_ar_alpha;
  real<lower=0> prior_sha_st_ar_beta;

 //Shared long-term
	vector <lower=0> [N] prior_sha_lt_sd; //sd for prior on error

	//Random walk on y
	vector [N] prior_y_init_mean;
	vector <lower=0> [N] prior_y_init_var;
	real<lower=N-1>	prior_sigma_t_inv_wish_nu; //inverse wishart
	matrix[N, N] prior_sigma_t_inv_wish_sigma;//inverse wishart

}
transformed data{

  matrix[N + sum(model_num_species) , (M+2) * N ] bigM = rep_matrix(0,N + sum(model_num_species) , (M+2) * N );
  // transfomations used for sequential modelling
  vector[N + sum(model_num_species) ] all_eigenvalues_cov;
  matrix[N + sum(model_num_species)  , N + sum(model_num_species) ] all_eigenvectors_cov = rep_matrix(0,N + sum(model_num_species) , N + sum(model_num_species) );
  // vector
  matrix[time , N + sum(model_num_species) ] observation_available;
  matrix[time , N + sum(model_num_species) ] new_data;

  bigM[1:N,1:N] = diag_matrix(rep_vector(1.0,N));


  all_eigenvalues_cov[1:N] = eigenvalues_sym(obs_covariances);
  all_eigenvectors_cov[1:N,1:N] = eigenvectors_sym(obs_covariances);
  if (M > 0){
    bigM[(N + 1):(N + model_num_species[1] ),1:N] = Ms[1:model_num_species[1],];
    bigM[(N + 1):(N + model_num_species[1] ),(N+1):(2*N)] = Ms[1:model_num_species[1],];
    bigM[(N + 1):(N + model_num_species[1] ),(2*N + 1):(3 * N)] = Ms[1:model_num_species[1],];
	  all_eigenvalues_cov[(N + 1):(N + model_num_species[1] )] = eigenvalues_sym(to_matrix(model_covariances[1:sq_int(model_num_species,1)],model_num_species[1],model_num_species[1]));
	  all_eigenvectors_cov[(N + 1):(N + model_num_species[1]),(N + 1):(N + model_num_species[1])] = eigenvectors_sym(to_matrix(model_covariances[1:sq_int(model_num_species,1)],model_num_species[1],model_num_species[1]));
	  if(M > 1){
		  for(i in 2:M){
		    int start_index =  N + sum(model_num_species[1:(i-1)]) + 1;
		    int end_index = N + sum(model_num_species[1:i]);
		    matrix[model_num_species[i], N ] Ms_transformation = Ms[(sum(model_num_species[1:(i-1)]) + 1):sum(model_num_species[1:i]),];
		    bigM[start_index:end_index,1:N] = Ms_transformation;
  		  bigM[start_index:end_index,(N+1):(2*N)] = Ms_transformation;
	   	  bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] = Ms_transformation;
		    all_eigenvalues_cov[start_index:end_index] = eigenvalues_sym(to_matrix(model_covariances[(sq_int(model_num_species,i-1) + 1):(sq_int(model_num_species,i))],model_num_species[i],model_num_species[i]));
		    all_eigenvectors_cov[start_index:end_index, start_index:end_index] = eigenvectors_sym(to_matrix(model_covariances[(sq_int(model_num_species,i-1) + 1):(sq_int(model_num_species,i))],model_num_species[i],model_num_species[i]));
		  }
	  }
  }

  //Expand the observation times to be observations for each model for each species for each time step
  for (i in 1:time){
    //Data observations
    observation_available[i,1:N] = rep_row_vector(observation_times[i,1],N);
	  if (M > 0){
	    //Model observations
		  observation_available[i,(N+1):(N + model_num_species[1])] = rep_row_vector(observation_times[i,2],model_num_species[1]);
		  for (j in 2:M){
			  observation_available[i,(N + sum(model_num_species[1:(j-1)]) + 1):(N + sum(model_num_species[1:j]))] = rep_row_vector(observation_times[i,j + 1],model_num_species[j]);
		  }
	  }

	  new_data[i,] = (all_eigenvectors_cov' * append_row(observations[i],model_outputs[i,]'))';
  }

  bigM = all_eigenvectors_cov' * bigM;
}

parameters{
  /**
   * Simulator discrepancies
   */
  // Individual
  vector <lower=-1,upper=1>[N] ind_st_ar_param[M];
  vector<lower=0>[N] ind_st_var[M];
  corr_matrix [N] ind_st_cor[M];
  vector[N] ind_lt_raw[M];
  vector <lower=0> [N] ind_lt_var;
  corr_matrix [N] ind_lt_cor;
  // Shared
  vector <lower=-1,upper=1>[N] sha_st_ar_param;
  vector <lower=0> [N] sha_st_var;
  corr_matrix [N] sha_st_cor;
  vector [N] sha_lt_raw;

  /**
   * Random walk on y
   */
  cov_matrix [N] SIGMA_t;

}
transformed parameters{
  matrix [N,N] SIGMA_x[M];
  vector [N] ind_st_sd[M];
  vector [N] sha_lt = prior_sha_lt_sd .* sha_lt_raw;
  vector [N] ind_lt[M];
  vector [N] ind_lt_sd = sqrt(ind_lt_var);
  matrix [N,N] ind_lt_covar = diag_post_multiply(diag_pre_multiply(ind_lt_sd,ind_lt_cor),ind_lt_sd);
  matrix [N,N] ind_lt_cov_cholesky = cholesky_decompose(ind_lt_covar);


  vector [(M+2) * N] x_hat = append_row(prior_y_init_mean,rep_vector(0.0,N * (M + 1)));
  matrix [(M+2) * N,(M+2) * N] SIGMA_init = rep_matrix(0,(M+2) * N,(M+2) * N );

  vector [N] sha_st_sd = sqrt(sha_st_var);
  matrix [N,N] SIGMA_mu = diag_post_multiply(diag_pre_multiply(sha_st_sd, sha_st_cor), sha_st_sd);

  /**
  *  Kalman filter parameters:
  *  In each case, the ordering to be passed through to the Kalman Filter is:
  *  (1) Observations (2) Consensus (3) Models
  *     [if (1) is available, otherwise this is ignored]
  *
  */
  matrix[(M+2) * N , (M+2) * N] SIGMA = rep_matrix(0,(M+2) * N,(M+2) * N );
  vector[(M+2) * N ] lt_discrepancies;
  vector[(M+2) * N] AR_params;

  //SIGMA
  SIGMA[1:N, 1:N ] = SIGMA_t;
  SIGMA[(N + 1):(2*N), (N + 1):(2*N) ] = SIGMA_mu;
  for (i in 1:M){
    ind_st_sd[i] = sqrt(ind_st_var[i]);

    SIGMA_x[i] = diag_post_multiply(diag_pre_multiply(ind_st_sd[i],ind_st_cor[i]),ind_st_sd[i]);
	  SIGMA[((i+1) * N + 1):((i+2)*N ),((i + 1) * N + 1):((i+2)*N )] = SIGMA_x[i];
  }

  //SIGMA_init
  SIGMA_init[1:N,1:N] = diag_matrix(prior_y_init_var);
  SIGMA_init[(N + 1):(2*N), (N + 1):(2*N) ] = SIGMA_mu ./ (1 - sha_st_ar_param * sha_st_ar_param');;
  for (i in 1:M){
    SIGMA_init[((i+1) * N + 1):((i+2)*N ),((i+1) * N + 1):((i+2)*N )] = SIGMA_x[i] ./ (1 - ind_st_ar_param[i] * ind_st_ar_param[i]');

  }


  lt_discrepancies[1:(2 * N)] = append_row(rep_vector(0.0,N), sha_lt);
  AR_params[1:(2 * N)] = append_row(rep_vector(1.0,N), sha_st_ar_param);
  for (i in 1:M){
    ind_lt[i] = ind_lt_cov_cholesky*ind_lt_raw[i];
    lt_discrepancies[((i+1) * N + 1):((i+2)*N )] = ind_lt[i];
	  AR_params[((i+1) * N + 1):((i+2)*N )] = ind_st_ar_param[i];
  }
}
model{
  /**
  * Priors
  */
  //Random walk on y
  SIGMA_t ~ inv_wishart(prior_sigma_t_inv_wish_nu, prior_sigma_t_inv_wish_sigma); // the random walk of y


  // Shared discrepancies
  sha_lt_raw ~ std_normal();
  sha_st_var ~ gamma(prior_sha_st_var_a,prior_sha_st_var_b); // Variance
  //JM 22/07: Beta priors on the AR parameters
  target += beta_lpdf((sha_st_ar_param + 1)/2 | prior_sha_st_ar_alpha, prior_sha_st_ar_beta);

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
  //Long term correlations
  if(form_prior_ind_lt == 0){
    ind_lt_cor ~ lkj_corr(prior_ind_lt_cor_lkj[1]);
  } else if(form_prior_ind_lt == 1){
    ind_lt_cor ~ inv_wishart(prior_ind_lt_cor_wish_nu[1], prior_ind_lt_cor_wish_sigma);
  } else{
	  target += priors_cor_beta(ind_lt_cor, N, prior_ind_lt_cor_beta_1, prior_ind_lt_cor_beta_2);
  }


  for(i in 1:M){
    //AR Parameters
    target += beta_lpdf((ind_st_ar_param[i] + 1)/2 | prior_ind_st_ar_alpha, prior_ind_st_ar_beta);

    ind_st_var[i] ~ gamma(prior_ind_st_var_a, prior_ind_st_var_b);
    ind_lt_raw[i] ~ std_normal();
	if (form_prior_ind_st == 0){
	  ind_st_cor[i] ~ lkj_corr(prior_ind_st_cor_lkj[1]);
    } else if(form_prior_ind_st == 1){
	  ind_st_cor[i] ~ inv_wishart(prior_ind_st_cor_wish_nu[1], prior_ind_st_cor_wish_sigma);
	} else{
	  target += priors_cor_beta(ind_st_cor[i], N, prior_ind_st_cor_beta_1, prior_ind_st_cor_beta_2);
	}
  }


  /**
   * Likelihood
   */
  // Kalman filter
  target += KalmanFilter_seq_em(AR_params, lt_discrepancies, all_eigenvalues_cov, SIGMA, bigM,
                                SIGMA_init, x_hat, time, new_data, observation_available);
}

