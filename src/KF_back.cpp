#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//' @title Backwards Kalman filter
//' @description Finds the most likely path through a dynamical linear model.
//' @param rhos A \code{numeric} containing the diagonal elements of the transition matrix of the evolution equation.
//' @param dee A \code{numeric} of the same length as \code{rho} containing the discrepancies or biases in the observation process.
//' @param R A \code{numeric} representing the variances of the observation process.
//' @param Q A \code{matrix} of dimensions \code{length(rho)} and \code{length(rho)} representing the covariance of the evolution process.
//' @param C A a \code{matrix} of dimensions \code{length(rho)} and \code{length(R)} representing the observation operator of the observation process.
//' @param P A \code{matrix} of dimensions \code{length(rho)} and \code{length(rho)} representing the covariance matrix of the first state of the evolution process.
//' @param xhat A \code{numeric} of the same length as \code{rho} representing the expectation of the first state of the evolution process.
//' @param Time A numeric The length of time of the dynamical linear model.
//' @param y A \code{matrix} of dimensions \code{Time} and \code{length(R)} of observations of the observation process.
//' @param obs A \code{matrix} of dimensions \code{Time} and \code{length(R)}. \code{1} means in the \code{i},\code{j}th element means that the \code{j}th output is observed at tim \code{i}.
//' @details For the model with the evolution process \deqn{x_{t+1}\sim{}N(\rho{}\cdot{}x_{t},Q)} and observation process \deqn{y_{t}\sim{}N(\rho{}(x_{t} + \delta),diag(R))}.
//' @details Using the sequential Kalman filter, the function gives the mostly path of \eqn{x_{t}} for all \eqn{t}.
//' @return A \code{matrix} with dimensions \code{nrow(time)} and \code{length(xhat)} representing the most likely values of the latent variables.
//' @references Chui, C.K. & Chen, G. (2009) Kalman Filtering with Real-Time Applications. Springer, Berlin, Heidelberg, Fourth Edtion.
//' @references Kalman, R. E. (1960) A new approach to linear filtering and prediction problems. Trans. ASME, J. Basic Eng., 82, pp. 35-45.
//' @examples
//'\donttest{
//' fit <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
//'                simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
//'                                  list(SSB_fs,  Sigma_fs, "FishSUMS"),
//'                                  list(SSB_lm,  Sigma_lm, "LeMans"),
//'                                  list(SSB_miz, Sigma_miz, "Mizer")),
//'                priors = EnsemblePrior(4,
//'                ind_st_params = IndSTPrior(parametrisation_form = "lkj",
//'                var_params= list(1,1), cor_params = 10, AR_params = c(2, 2))),
//'                full_sample = FALSE) #Only optimise in this case
//' transformed_data <- get_transformed_data(fit)
//' ex.fit <- fit@point_estimate$par
//' params <- get_parameters(ex.fit)
//' ret <- KalmanFilter_back(params$AR_params, params$lt_discrepancies,
//'                           transformed_data$all_eigenvalues_cov,params$SIGMA,
//'                           transformed_data$bigM, params$SIGMA_init, params$x_hat,
//'                           fit@ensemble_data@stan_input$time,transformed_data$new_data,
//'                           transformed_data$observation_available)
//'}
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd KalmanFilter_back(Eigen::VectorXd rhos, Eigen::VectorXd dee, Eigen::VectorXd R, Eigen::MatrixXd Q, Eigen::MatrixXd C,// parameters
                            Eigen::MatrixXd P, Eigen::VectorXd xhat, int Time, Eigen::MatrixXd y, Eigen::MatrixXd obs){

  Eigen::MatrixXd xhat_b = Eigen::MatrixXd::Zero(Time, xhat.size());
  Eigen::MatrixXd P_ = P;
  // DO NOT CHANGE THE BELOW. For linux distributions, the below is not a valid definition
  //  because Time is not know at compile time. To have dynamically sized arrays we must use std::vector.
  //  This is not a problem for windows distributions, but will cause problems when you submit to CRAN.
  //Eigen::MatrixXd P_s[Time];
  //Eigen::MatrixXd G[Time];
  std::vector<Eigen::MatrixXd> P_s(Time, Eigen::MatrixXd(xhat.size(), y.cols()));
  std::vector<Eigen::MatrixXd> G(Time, Eigen::MatrixXd(xhat.size(), y.cols()));
  Eigen::MatrixXd A = rhos * rhos.transpose();
  Eigen::VectorXd xhat_ = xhat;
  Eigen::MatrixXd xhat_s = Eigen::MatrixXd::Zero(Time, xhat.size());
  Eigen::MatrixXd G_ = Eigen::MatrixXd::Zero(xhat.size(), y.cols());
  Eigen::MatrixXd er = Eigen::MatrixXd::Zero(Time, y.cols());
  double est;
  Eigen::MatrixXd Qs_inv = Eigen::MatrixXd::Zero(Time, y.cols());
  /// smoothing bit
  Eigen::MatrixXd Identity = Eigen::MatrixXd::Identity(xhat.size(), xhat.size());
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(xhat.size(), xhat.size());
  Eigen::VectorXd r = Eigen::VectorXd::Zero(xhat.size());
  int k_i;
  int l_j;

  for (int i = 0; i < Time; ++i)
  {
    //JM: Need to initialise some of the other variables
    //G[i] = Eigen::MatrixXd(xhat.size(), y.cols());
    //G.push_back(Eigen::MatrixXd(xhat.size(), y.cols()));

    P_ = P_.cwiseProduct(A) + Q;
    xhat_ = rhos.cwiseProduct(xhat_);

    //P_s[i] = P_ ;
    P_s.at(i) = P_;

    xhat_s.row(i) = xhat_.transpose() ;
    for (int j = 0; j < y.cols(); ++j) // change this to some input that lists the numbers of the elements of y that we are interested in, maybe another value that says how many we are interested in.
    {
      if (obs(i,j)==1.0)
      {
        G_.col(j) = P_ * C.row(j).transpose();
        double Qs = (C.row(j) * G_.col(j) + R(j));
        Qs_inv(i, j) =  1 /Qs;
        est = (xhat_ + dee).dot(C.row(j)) ;
        er(i,j) = y(i,j) - est;
        G[i].col(j) = Qs_inv(i,j) * G_.col(j);
        xhat_ = xhat_ + G[i].col(j) * er(i,j);
        P_ -=  G[i].col(j) * G_.col(j).transpose();
      }
    }
  }
  for (int i = 0; i < Time; ++i)
  {
    //k_i = Time + 1 - i;
    k_i = Time - i - 1;
    int num_cols = y.cols();
    for (int j = 0; j < num_cols; ++j)
    {
      //l_j = num_cols + 1 - j;
      l_j = num_cols - j - 1;
      if (obs(k_i,l_j) == 1.0)
      {
        L = Identity - G[k_i].col(l_j) * C.row(l_j);
        r = C.row(l_j).transpose()  * Qs_inv(k_i,l_j) * er(k_i,l_j) + L.transpose() * r;

      }
    }

    xhat_b.row(k_i) = xhat_s.row(k_i) + (P_s[k_i] * r).transpose();
    r = rhos.cwiseProduct(r);
  }
  return xhat_b;
}
