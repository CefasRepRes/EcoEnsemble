functions{
  matrix KalmanFilter_back(vector rhos, vector dee, vector R, matrix Q, matrix C,// parameters
  matrix P, vector xhat, int Time, matrix y,matrix obs){
  matrix [Time,num_elements(xhat)] xhat_b;
  matrix [num_elements(xhat),num_elements(xhat)] P_ = P;
  matrix [num_elements(xhat),num_elements(xhat)] P_s [Time];
  matrix [num_elements(xhat),num_elements(xhat)] A = rhos * rhos';
  vector [num_elements(xhat)] xhat_=xhat;
  matrix [Time,num_elements(xhat)] xhat_s;
  matrix [num_elements(xhat),cols(y)] G[Time];
  matrix [num_elements(xhat),cols(y)] G_;
  matrix [Time,cols(y)] er;
  real est;
  matrix[Time,cols(y)] Qs_inv;
  /// smoothing bit
  matrix [num_elements(xhat),num_elements(xhat)] Identity = diag_matrix(rep_vector(1.0,num_elements(xhat)));
  matrix [num_elements(xhat),num_elements(xhat)] L;
  vector [num_elements(xhat)] r = rep_vector(0,num_elements(xhat));
  int k_i;
  int l_j;
  for (i in 1:Time)
  {
    P_ = P_ .* A + Q; // l 455
    xhat_ = rhos .* xhat_; // l 454
    P_s[i,,] = P_ ;
    xhat_s[i,] = xhat_' ;
    for (j in 1:cols(y)) // change this to some input that lists the numbers of the elements of y that we are interested in, maybe another value that says how many we are interested in.
    {
      if (obs[i,j]==1.0)
      {
        G_[,j] = P_ * C[j,]';
        Qs_inv[i,j] = 1 / (C[j,] * G_[,j] + R[j]);
        est = dot_product((xhat_ + dee)' ,  C[j,]) ;
        er[i,j] = y[i,j] - est;
        G[i,,j] = Qs_inv[i,j] * G_[,j];
        xhat_ += G[i,,j] * er[i,j];
        P_ -=  G[i,,j] * G_[,j]';
      }
    }
  }
  for (i in 1:Time)
  {
    k_i = Time + 1 - i;
    for (j in 1:cols(y))
    {
      l_j = cols(y) + 1 - j;
      if (obs[k_i,l_j] == 1.0)
      {
        L = Identity - G[k_i,,l_j] * C[l_j,];
        r = C[l_j,]'  * Qs_inv[k_i,l_j] * er[k_i,l_j] + L' * r;
      }
    }
    xhat_b[k_i,] = xhat_s[k_i,] + (P_s[k_i,,] * r)';
    r = rhos .* r;
  }
  return xhat_b;
  }
}
