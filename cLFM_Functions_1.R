###############################################################################
## Description: Functions for conducting the proposed EM estimation algorithm   
##              described in Algorithm 1 of 'Contrastive Latent Functional Model'.
###############################################################################
## Functions included:
## Main function:
##    run_em: Function for fitting a cLFM model and executing the 
##            proposed EM estimation algorithm described in Algorithm 1.
## Supporting functions used by the main function:
##    1. FPCA_vf: Function for applying FPCA to the input data. Number of retained eigenfunctions 
##                is based on the preset variation threshold.
##    2. Cov_mat: Function for calculating the empirical or 2D smoothed covariance 
##                surface of the given data.
##    3. adap: Function for adapting the input functional matrix or vector by retaining their 
##             information only at the observed time points of the specified data subject with
##             subject-specific adaptation matrix A_i^x/A_j^y.
##    4. smooth_est: Function used in Step 7 for deriving the smoothed estimates of latent components 
##                    and score variances from the estimated shared/unique covariance spaces.
###############################################################################

run_em = function(X, Y, L, K, R, tol, verbose){
  
  #############################################################################
  ## Description: (main function) This function executes the proposed EM estimation algorithm described in Algorithm 1.
  ##              The two input data grids (X, Y) are constructed with N columns, representing measurements at the 
  ##              collection of total N distinct time points (t_1, ..., t_N) across all subjects from both 
  ##              group, where NA is used for missing data.
  ## Args:        X: data grid of the first group of data (matrix, n_x*N)
  ##              Y: data grid of the second group of data (matrix, n_y*N)
  ##              L: number of shared latent components, L (integer).
  ##              K: number of unique latent components in X, K (integer).
  ##              R: number of unique latent components in Y, R (integer).
  ##              tol: predefined tolerance level to stop the EM iteration, (scalar, varepsilon in Step 6)
  ##              verbose: indicator whether to display the EM estimation progress. (logical)
  ## Returns:     a list containing:  
  ##                mu_x: estimate of the over-all mean function mu^x(t) (vector, N*1)
  ##                mu_y: estimate of the over-all mean function mu^y(t) (vector, N*1)
  ##                s2_x: estimate of the error variance sigma^2_X (scalar)
  ##                s2_y: estimate of the error variance sigma^2_Y (scalar)
  ##                Psi: smoothed estimate of the shared latent components {psi_l(t), l = 1,...,L} (matrix, N*L).
  ##                Phi: smoothed estimate of the unique latent components {phi_k(t), k = 1,...,K} (matrix, N*K).
  ##                Gamma: smoothed estimate of the unique latent components {gamma_r(t), r = 1,...,R} (matrix, N*R).
  ##                Omega_x: diagonal matrix of estimates of the score variances {omega^x_l, l = 1,...,L} (matrix, L*L).
  ##                Omega_y: diagonal matrix of estimates of the score variances {omega^y_l, l = 1,...,L} (matrix, L*L).
  ##                Theta: diagonal matrix of estimates of the score variances {theta_k, k = 1,...,K} (matrix, K*K).
  ##                Lambda: diagonal matrix of estimates of the score variance {lambda_r, r = 1,...,R} (matrix, R*R).
  #############################################################################
  

  if(ncol(X) != ncol(Y)){
    stop(paste("The total distinct time points are inconsistent between two",
               "groups of data.") )
  }else{
    N = ncol(X) # total number of time points
    tobs = seq(0, 1, length.out = N)  # total time grid T (vector, N*1).
  }
  
  if(L < 0){
    stop("The dimension of shared variation space is negative.")
  }else if(K < 0){
    stop("The dimension of unique variation space of X is negative.")
  }else if(R < 0){
    stop("The dimension of unique variation space of Y is negative.")
  }else if(L + K == 0){
    stop("The dimension of total variation space of X is non-positive.")
  }else if(L + R == 0){
    stop("The dimension of total variation space of Y is non-positive.")
  }
  
  
  
  #################################
  # Step 1. Estimate the overall mean functions mu^x(t) and mu^y(t) (by penalized smoothing spline) 
  #         and error variances sigma^2_X and sigma^2_Y, and centralize the input data.
  #################################
  
  ## mu^x(t)
  n_x = nrow(X) # number of subjects within the first group
  data_tab = data.frame(na.omit(cbind(rep(1 : n_x, each = N), rep(tobs, n_x), 
                                       as.vector(t(X))))) 
  colnames(data_tab) = c("id", "timePoints", "X")
  est_mean = smooth.spline(data_tab$timePoints, data_tab$X, cv = F, spar = NULL) 
    # `spar` is a scalar between 0 to 1 working as a smoothing parameter controlling this spline 
    # estimation. When it is 'NULL', it will be estimated by generalized cross-validation. It
    # can be specified by users when the it is over- or under-smoothed. 
  mu_x = matrix(data = est_mean$y, nrow = length(tobs), ncol = 1)
  X1.c = sweep(X, 2, mu_x) # centralize the first group of data, (X(t) - mu^x(t))
  
  
  ## mu^y
  n_y = nrow(Y) # number of subjects within the second group
  data_tab = data.frame(na.omit(cbind(rep(1:n_y, each = N), rep(tobs, n_y), 
                                       as.vector(t(Y))))) 
  colnames(data_tab) = c("id", "timePoints", "Y")
  est_mean = smooth.spline(data_tab$timePoints, data_tab$Y, cv = F, spar = NULL)
  mu_y = matrix(data = est_mean$y, nrow = length(tobs), ncol = 1)
  Y1.c = sweep(Y, 2, mu_y) # centralize the second group of data (Y(t) - mu^y(t))
  
  
  ## sigma^2_X
  empCov = Cov_mat(X1.c, tobs = tobs, smooth_2D = FALSE) # empirical covariance matrix
  smoothCov = Cov_mat(X1.c, tobs = tobs, smooth_2D = TRUE, 
                      user_k = NULL) 
    # 2D smoothed covariance matrix without contamination from errors within the diagonal elements.
    # The default number of basis functions used in smoothing process is 10, and it can be 
    # specified by users through `user_k`.
  s2_x = mean(diag(empCov) - diag(smoothCov))
  if(s2_x <= 0){
    stop("The estimated error variance from the first sample is non-positive.")
  }
  
  
  ## sigma^2_Y
  empCov = Cov_mat(Y1.c, tobs = tobs, smooth_2D = FALSE)
  smoothCov = Cov_mat(Y1.c, tobs = tobs, smooth_2D = TRUE, 
                      user_k = NULL)
  s2_y = mean(diag(empCov) - diag(smoothCov))
  if(s2_y <= 0){
    stop("The estimated error variance from the second sample is non-positive.")
  }

  
  Psi = NULL # Shared latent component matrix
  Phi = NULL # Unique latent component matrix of X
  Gamma = NULL # Unique latent component matrix of Y
  Omega_x = NULL # Score variance matrix of shared components in X
  Omega_y = NULL # Score variance matrix of shared components in Y
  Theta = NULL # Score variance matrix of unique components in X
  Lambda = NULL # Score variance matrix of unique components in Y
  G_Psi_x = NULL # Shared covariance surface of X
  G_Psi_y = NULL # Shared covariance surface of Y
  G_Phi = NULL # Unique covariance surface of X
  G_Gamma = NULL # Unique covariance surface of Y
  stop_flag = 0 # Flag if it fails to converge
  repeat{ # iterations between Step 3 and Step 8 until the percent change in covariance surface estimates < 10%.
    
    #################################
    # Step 2. Initialize the model parameters: latent components (Psi, Phi, Gamma) and
    #         their respective score variances (Omega^x, Omega^y, Theta, Lambda).
    #################################
    
    ## Generate a N*N-dim identity matrix. 
    id_mat = matrix(data = 0, nrow = N, ncol = N) 
    diag(id_mat) = rep(1, N)
    
    
    ## Psi, Omega^x, Omega^y
    if(L > 0){
      if(is.null(Psi)){
        # Init Psi with the first N*L sub-matrix 
        Psi = matrix(data = id_mat[1:N, 1:L], ncol = L) 
      }
      if(is.null(Omega_x)){Omega_x = diag(L)}
      if(is.null(Omega_y)){Omega_y = diag(L)}
    }else{ # when L = 0
      Psi = NULL
      Omega_x = NULL
      Omega_y = NULL
    } 
    
    
    ## Phi, Theta
    if(K > 0){
      if(is.null(Phi)){
        # Init Phi with the following N*K sub-matrix
        Phi = matrix(data = id_mat[1:N, (L+1):(L+K)], ncol = K)
      }
      if(is.null(Theta)){Theta = diag(K)}
    }else{ # when K = 0
      Phi = NULL
      Theta = NULL
    } 
    
    
    ## Gamma, Lambda
    if(R > 0){
      if(is.null(Gamma)){
        # Init Gamma with the next N*R sub-matrix
        Gamma = matrix(data = id_mat[1:N, (L+K+1):(L+K+R)], ncol = R)
      }
      if(is.null(Lambda)){Lambda = diag(R)}
    }else{ # when R = 0
      Gamma = NULL
      Lambda = NULL
    } 

    
    # Start EM iterations (Step 3 to Step 6) to estimate the latent component matrices and score
    # variance matrices. The estimation process is stopped if the percent change in log-likelihood 
    # is below the predefined tolerance level. The default maximum of this iteration process is set to 200.
    
    ## Four indicators denoting the origin of the over-fitting issue if the estimation fails to converge.
    L_x_of = FALSE
    L_y_of = FALSE
    K_of = FALSE
    R_of = FALSE
    
    llh_old = NULL 
      # The log-likelihood calculated in the last iteration and compared with the latest one 
      # to determine the convergence of EM.
    for (ite in c(1:200)) {
      
      #################################
      # Step 3. (E-step) Calculate the conditional quantities and the current log-likelihood value.
      #################################
      
      ## Pre-compute the quantities that will be used multiple times in the following computation
      if(L > 0){
        PsiO_x = Psi %*% Omega_x
        PsiO_y = Psi %*% Omega_y
      }
      
      if(K > 0){
        PhiT = Phi %*% Theta
      }
      
      if(R > 0){
        GammaL = Gamma %*% Lambda
      }
      
      
      ## Conditional expectations
      if(L > 0){
        E_eta_x = matrix(data = NA, nrow = n_x, ncol = L) 
          # Each row saves the conditional expectation E(eta^x_i|X_i), i = 1,...,n_x. (vector, 1*L)
        for (i in c(1 : n_x)) {
          x = matrix(X1.c[i,], nrow = N, ncol = 1)
          E_eta_x[i,] = t(adap(x, x)) %*% adap(PsiO_x, x) %*% 
            t(solve(s2_x * diag(L) + Omega_x %*% crossprod(adap(Psi, x))))
        }
        
        E_eta_y = matrix(data = NA, nrow = n_y, ncol = L) 
          # Each row saves the conditional expectation E(eta^y_j|Y_j), j = 1,...,n_y. (vector, 1*L)
        for (i in c(1 : n_y)) {
          y = matrix(Y1.c[i,], nrow = N, ncol = 1)
          E_eta_y[i,] = t(adap(y, y)) %*% adap(PsiO_y, y) %*% 
            t(solve(s2_y * diag(L) + Omega_y %*% crossprod(adap(Psi, y))))
        }
      }
      
      if(K > 0){
        E_xi = matrix(data = NA, nrow = n_x, ncol = K) 
          # Each row saves the conditional expectation E(xi_i|X_i), i = 1,...,n_x. (vector, 1*K)
        for (i in c(1 : n_x)) {
          x = matrix(X1.c[i,], nrow = N, ncol = 1)
          E_xi[i,] = t(adap(x, x)) %*% adap(PhiT, x) %*% 
            t(solve(s2_x * diag(K) + Theta %*% crossprod(adap(Phi, x))))
        }
      }
      
      if(R > 0){
        E_zeta = matrix(data = NA, nrow = n_y, ncol = R) 
          # Each row saves the conditional expectation E(zeta_j|Y_j), j = 1,...,n_y. (vector, 1*R)
        for (i in c(1 : n_y)) {
          y = matrix(Y1.c[i,], nrow = N, ncol = 1)
          E_zeta[i,] = t(adap(y, y)) %*% adap(GammaL, y) %*% 
            t(solve(s2_y * diag(R) + Lambda %*% crossprod(adap(Gamma, y))))
        }
      }
      
      
      ## Conditional variances and conditional 2nd moments
      if(L > 0){
        E_ee_x = matrix(data = NA, nrow = n_x, ncol = (L^2)) 
          # Each row captures the conditional second moments E(eta^x_i (eta^x_i)^T|X_i) 
          # (matrix, L*L) flattened to a 1*L^2 vector, i = 1,...,n_x.
        for (i in c(1 : n_x)) {
          x = matrix(X1.c[i,], nrow = N, ncol = 1) # column vector of the centralized subject data (X_i-mu^x)
          E_sh = matrix(E_eta_x[i,], nrow = 1) # transpose of the conditional expectation vector E(eta^x_i|X_i)
          Var_sh = s2_x * 
            solve(s2_x * diag(L) + Omega_x %*% crossprod(adap(Psi, x))) %*% 
            Omega_x # conditional variance matrix Var(eta^x_i|X_i) 
          E_ee_x[i,] = t(as.vector(crossprod(E_sh) + Var_sh))
        }
        
        
        E_ee_y = matrix(data = NA, nrow = n_y, ncol = (L^2)) 
          # Each row captures the conditional second moments E(eta^y_j (eta^y_j)^T|Y_j) 
          # (L*L matrix) flattened to a 1*L^2 vector, j = 1,...,n_y.
        for (j in c(1 : n_y)) {
          y = matrix(Y1.c[j,], nrow = N, ncol = 1) # column vector of the centralized subject data (Y_j-mu^y)
          E_sh = matrix(E_eta_y[j,], nrow = 1) # transpose of the conditional expectation vector E(eta^y_j|Y_j) 
          Var_sh = s2_y * 
            solve(s2_y * diag(L) + Omega_y %*% crossprod(adap(Psi, y))) %*% 
            Omega_y # conditional variance matrix Var(eta^y_j|Y_j) 
          E_ee_y[j,] = t(as.vector(crossprod(E_sh) + Var_sh))
        }
      }
      
      
      if(K > 0){
        E_xx = matrix(data = NA, nrow = n_x, ncol = (K^2)) 
          # Each row captures the conditional second moments E(xi_i (xi_i)^T|X_i) 
          # (K*K matrix) flattened to a 1*K^2 vector, i = 1,...,n_x.
        for (i in c(1 : n_x)) {
          x = matrix(X1.c[i,], nrow = N, ncol = 1) # column vector of the centralized subject data (X_i-mu^x)
          E_uni = matrix(E_xi[i,], nrow = 1) # transpose of the conditional expectation vector E(xi_i|X_i) 
          Var_uni = s2_x * 
            solve(s2_x * diag(K) + Theta %*% crossprod(adap(Phi, x))) %*% 
            Theta # conditional variance matrix Var(xi_i|X_i) 
          E_xx[i,] = t(as.vector(crossprod(E_uni) + Var_uni))
        }
      }
      
      
      if(R > 0){
        E_zz = matrix(data = NA, nrow = n_y, ncol = (R^2)) 
          # Each row captures the conditional second moments E(zeta_j (zeta_j)^T|Y_j) 
          # (R*R matrix) flattened to a 1*R^2 vector, j = 1,...,n_y.
        for (j in c(1 : n_y)) {
          y = matrix(Y1.c[j,], nrow = N, ncol = 1) # column vector of the centralized subject data (Y_j-mu^y)
          E_uni = matrix(E_zeta[j,], nrow = 1) # transpose of the conditional expectation vector E(zeta_j|Y_j) 
          Var_uni = s2_y * 
            solve(s2_y * diag(R) + Lambda %*% crossprod(adap(Gamma, y))) %*% 
            Lambda # conditional variance matrix Var(zeta_j|Y_j)
          E_zz[j,] = t(as.vector(crossprod(E_uni) + Var_uni))
        }
      }
      
      
      
      #################################
      # Step 4. (M-step) Derive the raw updates of shared/unique latent components and score variances.
      #################################
      
      ## Update latent component matrices Psi, Phi and Gamma at each time t (row-by-row).
      
      ### Psi
      if(L > 0){
        Psi_2 = Psi # create a new matrix to save the updated Psi
        for (t in 1 : N) {
          part1 = 0 # the sum of X part in the right square brackets of the updating formula of Psi
          part2 = 0 # the sum of X part in the left inverse curly brackets of the updating formula of Psi
          for (i in 1 : n_x) {
            if(is.na(X1.c[i,t])){ # skip to involve (X_i(t) - mu^x(t)) if missing
              next
            }else{
              if(K > 0){
                part1 = part1 + 
                  (X1.c[i,t] - Phi[t,] %*% matrix(E_xi[i,], ncol = 1)) %*% 
                  matrix(E_eta_x[i,], nrow = 1)
                part2 = part2 + 
                  matrix(E_ee_x[i,], nrow = L, ncol = L, byrow = FALSE)
                  # reform the conditional second moments from its flattened vector.
              }else{
                part1 = part1 + (X1.c[i,t]) %*% matrix(E_eta_x[i,], nrow = 1)
                part2 = part2 + 
                  matrix(E_ee_x[i,], nrow = L, ncol = L, byrow = FALSE)
              }
            }
          }
          
          part3 = 0 # the sum of Y part in the right square brackets of the updating formula of Psi
          part4 = 0 # the sum of Y part in the left inverse curly brackets of the updating formula of Psi
          for (j in 1 : n_y) {
            if(is.na(Y1.c[j,t])){ # skip to involve Y_j(t) if missing
              next
            }else{
              if(R > 0){
                part3 = part3 + 
                  (Y1.c[j,t] - Gamma[t,] %*% matrix(E_zeta[j,], ncol = 1)) %*% 
                  matrix(E_eta_y[j,], nrow = 1)
                part4 = part4 + 
                  matrix(E_ee_y[j,], nrow = L, ncol = L, byrow = FALSE)
              }else{
                part3 = part3 + (Y1.c[j,t]) %*% matrix(E_eta_y[j,], nrow = 1)
                part4 = part4 + 
                  matrix(E_ee_y[j,], nrow = L, ncol = L, byrow = FALSE)
              }
            }
          }
        
          eigen_res = eigen(part2 / s2_x + part4 / s2_y) # eigen-decomposition of the input square matrix 
          eigenval = eigen_res$values
          eigenfun = eigen_res$vector
          eigenval = eigenval[which(eigenval > 0)] # retain positive eigenvalues only. 
          if(length(eigenval) > 0){
            component = length(eigenval)
            inv_mat = 0
            for (ind in 1 : component) {
              inv_mat = inv_mat + (1 / eigenval[ind]) * 
                tcrossprod(eigenfun[, ind])
            }
          }else{
            stop("Error, no positive eigenvalues are available for updating 
                 Psi.")
          }
          
          Psi_2[t,] = (part1 / s2_x + part3 / s2_y) %*% inv_mat
        }
      }else{
        Psi_2 = NULL
      }
      
      
      ### Phi
      if(K > 0){
        Phi_2 = Phi # create a new matrix to save the updated Phi
        for (t in 1 : N) {
          part1 = 0 # the sum in the right square brackets of the updating formula of Phi
          part2 = 0 # the sum in the left inverse curly brackets of the updating formula of Phi
          for (i in 1 : n_x) {
            if(is.na(X1.c[i,t])){
              next
            }else{
              if(L > 0){
                part1 = part1 + 
                  (X1.c[i,t] - Psi[t,] %*% matrix(E_eta_x[i,], ncol = 1)) %*% 
                  matrix(E_xi[i,], nrow = 1)
                part2 = part2 + 
                  matrix(E_xx[i,], nrow = K, ncol = K, byrow = FALSE)
              }else{
                part1 = part1 + (X1.c[i,t]) %*% matrix(E_xi[i,], nrow = 1)
                part2 = part2 + 
                  matrix(E_xx[i,], nrow = K, ncol = K, byrow = FALSE)
              }
            }
          }
          
          eigen_res = eigen(part2) 
          eigenval = eigen_res$values
          eigenfun = eigen_res$vector
          eigenval = eigenval[which(eigenval > 0)]
          if(length(eigenval) > 0){
            inv_mat = 0
            component = length(eigenval)
            for (ind in 1 : component) {
              inv_mat = inv_mat + (1 / eigenval[ind]) * 
                tcrossprod(eigenfun[, ind])
            }
          }else{
            stop("Error, no positive eigenvalues are available for updating 
                 Phi.")
          }
          
          Phi_2[t,] = part1 %*% inv_mat
        }
      }else{
        Phi_2 = NULL
      }
      
      
      ### Gamma
      if(R > 0){
        Gamma_2 = Gamma # create a new matrix to save the updated Gamma
        for (t in 1 : N) {
          part1 = 0 # the sum in the right square brackets of the updating formula of Gamma
          part2 = 0 # the sum in the left inverse curly brackets of the updating formula of Gamma
          for (j in 1 : n_y) {
            if(is.na(Y1.c[j,t])){
              next
            }else{
              if(L > 0){
                part1 = part1 + 
                  (Y1.c[j,t] - Psi[t,] %*% matrix(E_eta_y[j,], ncol = 1)) %*% 
                  matrix(E_zeta[j,], nrow = 1)
                part2 = part2 + 
                  matrix(E_zz[j,], nrow = R, ncol = R, byrow = FALSE)
              }else{
                part1 = part1 + (Y1.c[j,t]) %*% matrix(E_zeta[j,], nrow = 1)
                part2 = part2 + 
                  matrix(E_zz[j,], nrow = R, ncol = R, byrow = FALSE)
              }
            }
          }
          
          eigen_res = eigen(part2) 
          eigenval = eigen_res$values
          eigenfun = eigen_res$vector
          eigenval = eigenval[which(eigenval > 0)]
          if(length(eigenval) > 0){
            inv_mat = 0
            component = length(eigenval)
            for (ind in 1 : component) {
              inv_mat = inv_mat + (1 / eigenval[ind]) * 
                tcrossprod(eigenfun[, ind])
            }
          }else{
            stop("Error, no positive eigenvalues are available for updating
                  Gamma.")
          }
          
          Gamma_2[t,] = part1 %*% inv_mat
        }
      }else{
        Gamma_2 = NULL
      }
      
      
      ## Update the score variance matrices Omega_x, Omega_y, Theta and Lambda
      if(L > 1){
        Omega_x_diag = diag(matrix(colSums(E_ee_x), nrow = L, 
                                   ncol = L, byrow = FALSE)) / n_x 
          # obtain the update of diagonal elements of Omega^x
        Omega_x_2 = diag(x = Omega_x_diag) # reconstruct the updated matrix Omega^x (L*L matrix)
        
        Omega_y_diag = diag(matrix(colSums(E_ee_y), nrow = L, 
                                   ncol = L, byrow = FALSE)) / n_y 
          # obtain the update of diagonal elements of Omega^y
        Omega_y_2 = diag(x = Omega_y_diag) # reconstruct the updated matrix Omega^y (L*L matrix)
      }else if(L == 1){
        Omega_x_2 = matrix(sum(E_ee_x) / n_x) # update Omega^x (1*1 matrix)
        Omega_y_2 = matrix(sum(E_ee_y) / n_y) # update Omega^y (1*1 matrix)
      }else{
        Omega_x_2 = NULL
        Omega_y_2 = NULL
      }
      
      
      if(K > 1){
        Theta_diag = diag(matrix(colSums(E_xx), nrow = K, 
                                 ncol = K, byrow = FALSE)) / n_x 
          # obtain the update of diagonal elements of Theta
        Theta_2 = diag(x = Theta_diag) # reconstruct the updated matrix Theta (K*K matrix)
      }else if(K == 1){
        Theta_2 = matrix(sum(E_xx) / n_x) # updated Theta (1*1 matrix)
      }else{
        Theta_2 = NULL
      }
      
      
      if(R > 1){
        Lambda_diag = diag(matrix(colSums(E_zz), nrow = R, 
                                  ncol = R, byrow = FALSE)) / n_y 
        # obtain the update of diagonal elements of Lambda
        Lambda_2 = diag(x = Lambda_diag) # reconstruct the updated matrix Lambda (R*R matrix)
      }else if(R == 1){
        Lambda_2 = matrix(sum(E_zz) / n_y) # updated Lambda (1*1 matrix)
      }else{
        Lambda_2 = NULL
      }
      
      
      #################################
      # Step 5. Gram-Schmidt process to the derived raw estimates
      #################################
      
      if(L > 1){
        for(e in 1 : L){
          normal.factor <- trapz(tobs, Psi_2[,e]^2) 
            # captures the residual variation contained in the latent component estimates.
          Psi_2[,e] <- (Psi_2[,e] / sqrt(normal.factor)) 
            # ensures unit integration of the squared latent component over the total observed times T.
          Omega_x_2[e,e] <- Omega_x_2[e,e] * normal.factor 
          Omega_y_2[e,e] <- Omega_y_2[e,e] * normal.factor
            # adjusts the corresponding score variances to ensure consistency of the original covariance surface.
        }
        # Gram-Schmidt process
        if(!inherits(try(suppressMessages(gramSchmidt(Psi_2)), silent = TRUE), 
                     "try-error")){
          gs <- suppressMessages(gramSchmidt(Psi_2))
          Psi_2 <- gs$Q
          Omega_x_2 <- diag(x = diag(gs$R %*% Omega_x_2 %*% t(gs$R)), nrow = L)
          Omega_y_2 <- diag(x = diag(gs$R %*% Omega_y_2 %*% t(gs$R)), nrow = L)
        }
      }else if(L == 1){
        normal.factor <- trapz(tobs, Psi_2^2)
        Psi_2 <- matrix(data = (Psi_2 / sqrt(normal.factor)), ncol = 1)
        Omega_x_2 <- Omega_x_2 * normal.factor
        Omega_y_2 <- Omega_y_2 * normal.factor
      }else{}
      
      
      if(K > 1){
        for(e in 1 : K){
          normal.factor <- trapz(tobs, Phi_2[,e]^2)
          Phi_2[,e] <- (Phi_2[,e] / sqrt(normal.factor))
          Theta_2[e,e] <- Theta_2[e,e] * normal.factor
        }
        # Gram-Schmidt process
        if(!inherits(try(suppressMessages(gramSchmidt(Phi_2)), silent = TRUE), 
                     "try-error")){
          gs <- suppressMessages(gramSchmidt(Phi_2))
          Phi_2 <- gs$Q
          Theta_2 <- diag(x = diag(gs$R %*% Theta_2 %*% t(gs$R)), nrow = K)
        }
      }else if(K == 1){
        normal.factor <- trapz(tobs, Phi_2^2)
        Phi_2 <- matrix(data = (Phi_2 / sqrt(normal.factor)), ncol = 1)
        Theta_2 <- Theta_2 * normal.factor
      }else{}
      
      
      if(R > 1){
        for(e in 1 : R){
          normal.factor <- trapz(tobs, Gamma_2[,e]^2)
          Gamma_2[,e] <- (Gamma_2[,e] / sqrt(normal.factor))
          Lambda_2[e,e] <- Lambda_2[e,e] * normal.factor
        }
        # Gram-Schmidt process
        if(!inherits(try(suppressMessages(gramSchmidt(Gamma_2)), silent = TRUE), 
                     "try-error")){
          gs <- suppressMessages(gramSchmidt(Gamma_2))
          Gamma_2 <- gs$Q
          Lambda_2 <- diag(x = diag(gs$R %*% Lambda_2 %*% t(gs$R)), nrow = R)
        }
      }else if(R == 1){
        normal.factor <- trapz(tobs, Gamma_2^2)
        Gamma_2 <- matrix(data = (Gamma_2 / sqrt(normal.factor)), ncol = 1)
        Lambda_2 <- Lambda_2 * normal.factor
      }else{}
      
      ## Pass the orthogonalized updates to the original ones and continue the EM iteration process.
      Psi <- Psi_2 # Psi
      Phi <- Phi_2 # Phi
      Gamma <- Gamma_2 # Gamma
      Omega_x <- Omega_x_2 # Omega_x
      Omega_y <- Omega_y_2 # Omega_y
      Theta <- Theta_2 # Theta
      Lambda <- Lambda_2 # Lambda
      
      
      #################################
      # Stopping criteria (Step 6)
      ## The estimation algorithm may fail to converge due to over-fitted shared or unique components. 
      ## The estimated score variances of over-fitted shared/unique components are much
      ## close to 0, indicating their ignorable explained variation within the given data. Hence, the 
      ## three dimensions L, K and R will be adjusted based on the current estimate of score variances 
      ## and cLFM model will be reset and refitted with the new L, K and R immediately.
      #################################
      
      ## Likelihood calculation
      if(L > 0){
        PsiO_x = Psi %*% Omega_x
        PsiO_y = Psi %*% Omega_y
      }
      
      if(K > 0){
        PhiT = Phi %*% Theta
      }
      
      if(R > 0){
        GammaL = Gamma %*% Lambda
      }
      
      part1 = 0 # part 1 involves the first group of data X only
      for (i in c(1 : n_x)) {
        x = matrix(X1.c[i,], nrow = N, ncol = 1) 
        non_ms = sum(!is.na(x)) # record the number of non-missing measurements of the subject X_i
        if(L > 0 & K > 0){
          E_sh = matrix((t(adap(x,x)) %*% adap(PsiO_x, x) %*% 
                           t(solve(s2_x * diag(L) + Omega_x %*% 
                                     crossprod(adap(Psi, x))))), nrow = 1)
          Var_sh = s2_x * 
            solve(s2_x * diag(L) + Omega_x %*% crossprod(adap(Psi, x))) %*% 
            Omega_x
          
          E_uni = matrix((t(adap(x,x)) %*% adap(PhiT, x) %*% 
                            t(solve(s2_x * diag(K) + Theta %*% 
                                      crossprod(adap(Phi, x))))), nrow = 1)
          Var_uni = s2_x * 
            solve(s2_x * diag(K) + Theta %*% crossprod(adap(Phi, x))) %*% 
            Theta
          
          
          part1 = part1 - (non_ms + L + K) / 2 * log(2 * pi) - 
            (non_ms / 2) * log(s2_x) - 
            (crossprod(adap(x, x)) + (tcrossprod(E_sh %*% t(adap(Psi, x))) + 
                                        tr(adap(Psi, x) %*% Var_sh %*% 
                                             t(adap(Psi, x)))) + 
               (tcrossprod(E_uni %*% t(adap(Phi, x))) + 
                  tr(adap(Phi, x) %*% Var_uni %*% t(adap(Phi, x)))) - 
               2 * (t(adap(x, x)) %*% adap(Psi, x) %*% t(E_sh)) - 
               2 * (t(adap(x,x)) %*% adap(Phi, x) %*% t(E_uni)) + 
               2 * (E_sh %*% t(adap(Psi,x)) %*% adap(Phi,x) %*% t(E_uni))) /
            (2 * s2_x) - 
            0.5 * (sum(log(diag(Omega_x))) + 
                     (1 / diag(Omega_x)) %*% (t(E_sh^2) + diag(Var_sh))) - 
            0.5 * (sum(log(diag(Theta))) + 
                     (1 / diag(Theta)) %*% (t(E_uni^2) + diag(Var_uni)))
        }else if(L > 0 & K == 0){
          E_sh = matrix((t(adap(x,x)) %*% adap(PsiO_x, x) %*% 
                           t(solve(s2_x * diag(L) + Omega_x %*% 
                                     crossprod(adap(Psi, x))))), nrow = 1)
          Var_sh = s2_x * 
            solve(s2_x * diag(L) + Omega_x %*% crossprod(adap(Psi, x))) %*% 
            Omega_x
          
          part1 = part1 - (non_ms + L) / 2 * log(2 * pi) - 
            (non_ms / 2) * log(s2_x) - 
            (crossprod(adap(x, x)) + (tcrossprod(E_sh %*% t(adap(Psi, x))) + 
                                        tr(adap(Psi, x) %*% Var_sh %*% 
                                             t(adap(Psi, x)))) - 
               2 * (t(adap(x, x)) %*% adap(Psi, x) %*% t(E_sh))) / 
            (2 * s2_x) - 0.5 * (sum(log(diag(Omega_x))) + 
                                  (1 / diag(Omega_x)) %*% (t(E_sh^2) + 
                                                             diag(Var_sh)))
        }else if(L == 0 & K > 0){
          E_uni = matrix((t(adap(x,x)) %*% adap(PhiT, x) %*% 
                            t(solve(s2_x * diag(K) + Theta %*% 
                                      crossprod(adap(Phi, x))))), nrow = 1)
          Var_uni = s2_x * 
            solve(s2_x * diag(K) + Theta %*% crossprod(adap(Phi, x))) %*% Theta
          
          part1 = part1 - (non_ms + K) / 2 * log(2 * pi) - 
            (non_ms / 2) * log(s2_x) - 
            (crossprod(adap(x, x)) + (tcrossprod(E_uni %*% t(adap(Phi, x))) + 
                                        tr(adap(Phi, x) %*% Var_uni %*% 
                                             t(adap(Phi, x)))) - 
               2 * (t(adap(x,x)) %*% adap(Phi, x) %*% t(E_uni))) / (2 * s2_x) - 
            0.5 * (sum(log(diag(Theta))) + (1 / diag(Theta)) %*% 
                     (t(E_uni^2) + diag(Var_uni)))
        }else{
          stop("L, K are both 0.")
        }
      }
      
      
      part2 = 0 # part 2 involves the second group of data Y only
      for (j in c(1 : n_y)) {
        y = matrix(Y1.c[j,], nrow = N, ncol = 1) 
        non_ms = sum(!is.na(y))
        if(L > 0 & R > 0){
          E_sh = matrix((t(adap(y,y)) %*% adap(PsiO_y, y) %*% 
                           t(solve(s2_y * diag(L) + Omega_y %*% 
                                     crossprod(adap(Psi, y))))), nrow = 1) 
          Var_sh = s2_y * 
            solve(s2_y * diag(L) + Omega_y %*% crossprod(adap(Psi, y))) %*% 
            Omega_y
          
          E_uni = matrix((t(adap(y,y)) %*% adap(GammaL, y) %*% 
                            t(solve(s2_y * diag(R) + Lambda %*% 
                                      crossprod(adap(Gamma, y))))), nrow = 1) 
          Var_uni = s2_y * 
            solve(s2_y * diag(R) + Lambda %*% crossprod(adap(Gamma, y))) %*% 
            Lambda
          
          part2 = part2 - (non_ms + L + R) / 2 * log(2 * pi) - 
            (non_ms / 2) * log(s2_y) - 
            (crossprod(adap(y, y)) + (tcrossprod(E_sh %*% t(adap(Psi, y))) + 
                                        tr(adap(Psi, y) %*% Var_sh %*% 
                                             t(adap(Psi, y)))) + 
               (tcrossprod(E_uni %*% t(adap(Gamma, y))) + 
                  tr(adap(Gamma, y) %*% Var_uni %*% t(adap(Gamma, y)))) - 
               2 * (t(adap(y, y)) %*% adap(Psi, y) %*% t(E_sh)) - 
               2 * (t(adap(y,y)) %*% adap(Gamma, y) %*% t(E_uni)) + 
               2 * (E_sh %*% t(adap(Psi,y)) %*% adap(Gamma,y) %*% t(E_uni))) / 
            (2 * s2_y) - 
            0.5 * (sum(log(diag(Omega_y))) + 
                     (1 / diag(Omega_y)) %*% (t(E_sh^2) + diag(Var_sh))) - 
            0.5 * (sum(log(diag(Lambda))) + 
                     (1 / diag(Lambda)) %*% (t(E_uni^2) + diag(Var_uni)))
        }else if(L > 0 & R == 0){
          E_sh = matrix((t(adap(y,y)) %*% adap(PsiO_y, y) %*% 
                           t(solve(s2_y * diag(L) + Omega_y %*% 
                                     crossprod(adap(Psi, y))))), nrow = 1) 
          Var_sh = s2_y * 
            solve(s2_y * diag(L) + Omega_y %*% crossprod(adap(Psi, y))) %*% 
            Omega_y
          
          part2 = part2 - (non_ms + R) / 2 * log(2 * pi) - 
            (non_ms / 2) * log(s2_y) - 
            (crossprod(adap(y, y)) + (tcrossprod(E_sh %*% t(adap(Psi, y))) + 
                                        tr(adap(Psi, y) %*% Var_sh %*% 
                                             t(adap(Psi, y)))) - 
               2 * (t(adap(y, y)) %*% adap(Psi, y) %*% t(E_sh))) / 
            (2 * s2_y) - 
            0.5 * (sum(log(diag(Omega_y))) + 
                     (1 / diag(Omega_y)) %*% (t(E_sh^2) + diag(Var_sh)))
        }else if (L == 0 & R > 0){
          E_uni = matrix((t(adap(y,y)) %*% adap(GammaL, y) %*% 
                            t(solve(s2_y * diag(R) + Lambda %*% 
                                      crossprod(adap(Gamma, y))))), nrow = 1) 
          Var_uni = s2_y * 
            solve(s2_y * diag(R) + Lambda %*% crossprod(adap(Gamma, y))) %*% 
            Lambda
          
          part2 = part2 - (non_ms + L + R) / 2 * log(2 * pi) - 
            (non_ms / 2) * log(s2_y) - 
            (crossprod(adap(y, y)) + (tcrossprod(E_uni %*% t(adap(Gamma, y))) + 
                                        tr(adap(Gamma, y) %*% Var_uni %*% 
                                             t(adap(Gamma, y)))) - 
               2 * (t(adap(y,y)) %*% adap(Gamma, y) %*% t(E_uni))) / 
            (2 * s2_y) - 
            0.5 * (sum(log(diag(Lambda))) + 
                     (1 / diag(Lambda)) %*% (t(E_uni^2) + diag(Var_uni)))
        }else{
          stop("L, R are both 0.")
        }
      }
      
      
      llh = part1 + part2 # the complete log-likelihood is derived as the sum of two parts
      if(verbose){print(paste0("Log-likelihood: ", round(llh, digits = 3)))}
      
      
      
      if(is.nan(llh) | ite == 200 | stop_flag == 20){ # detection of convergence failure
        # Detection in X
        ## If any shared/unique component estimates in X explain less than 1%, then L/K will be reduced.
        if(L > 0 & K > 0){
          scoreVar = c(diag(Omega_x), diag(Theta))
          prop = scoreVar / sum(scoreVar)
          for(i in 1 : length(prop)){
            if(prop[i] < .01){  
              if(i <= L){
                L_x_of = TRUE
              }else{
                K_of = TRUE
              }
            }
          }
        }else if(L > 0){
          scoreVar = c(diag(Omega_x))
          prop = scoreVar / sum(scoreVar)
          for(i in 1 : length(prop)){
            if(prop[i] < .01){
              L_x_of = TRUE
            }
          }
        }else{ # K > 0
          scoreVar = c(diag(Theta))
          prop = scoreVar / sum(scoreVar)
          for(i in 1 : length(prop)){
            if(prop[i] < .01){
              K_of = TRUE
            }
          }
        }
        
        
        # Detection in Y
        ## If any shared/unique component estimates in Y explain less than 1%, then L/R will be reduced.
        if(L > 0 & R > 0){
          scoreVar = c(diag(Omega_y), diag(Lambda))
          prop = scoreVar / sum(scoreVar)
          for(i in 1 : length(prop)){
            if(prop[i] < .01){
              if(i <= L){
                L_y_of = TRUE
              }else{
                R_of = TRUE
              }
            }
          }
        }else if(L > 0){
          scoreVar = c(diag(Omega_y))
          prop = scoreVar / sum(scoreVar)
          for(i in 1 : length(prop)){
            if(prop[i] < .01){
              L_y_of = TRUE
            }
          }
        }else{ # R > 0
          scoreVar = c(diag(Lambda))
          prop = scoreVar / sum(scoreVar)
          for(i in 1 : length(prop)){
            if(prop[i] < .01){
              R_of = TRUE
            }
          }
        }
        
        
        ## If the convergence failure cannot be solved by adjustment on shared/unique 
        ## space dimensions (L, K, R), an inspection is needed to address this issue.
        if(!L_x_of & !L_y_of & !K_of & !R_of){
          stop(paste("The estimation cannot converge with current setting of", 
                     "(L, K, R). Please check the input quantities and",
                     "smoothing parameters."))
        }else{
          break
        }
      }else{ # compare the percent change in log-likelihood to a predetermined tolerance level.
        if(!is.null(llh_old) && (abs(llh - llh_old) < tol * abs(llh_old))){
          break
        }else{
          llh_old = llh
        }
      }
      
      if(verbose){print(paste0(ite, "th round is finished."))}
    }
    
    
    #################################
    # Detection of over-fitting issue
    #################################
    
    # If the initial values of (L, K, R) lead to overfitting, variance components within Omega^x$,
    # Omega^y, Theta and Lambda explaining less than 1% of the variation will be discarded 
    # and (L, K, R) are adjusted accordingly. 
    
    L_old = L
    K_old = K
    R_old = R
    
    if(L_x_of & L_y_of){ 
      L = L - 1
    }else if(L_x_of & !L_y_of){
      if(!R_of){ # R is not large enough
        L = L - 1
        R = R + 1
      }else{
        L = L - 1
      }
    }else if(!L_x_of & L_y_of){ 
      if(!K_of){
        L = L - 1
        K = K + 1
      }else{
        L = L - 1
      }
    }else if(K_of & R_of){
      K = K - 1
      R = R - 1
    }else if(!K_of & R_of){
      R = R - 1
    }else if(K_of & !R_of){
      K = K - 1
    }else{}
    
    
    if(L != L_old | K != K_old | R != R_old){ 
      if(verbose){
        print(paste0("The proposed EM estimation algorithm fails to converge ", 
                     "with the initial setting of (L, K, R) = (", L_old, ", ", 
                     K_old, ", ", R_old, ")", ". cLFM model is refitted ", 
                     "with the adjusted (L, K, R) = (", L, ", ", K, ", ", R, 
                     ")."))
      }
      return(run_em(X, Y, L, K, R, tol, verbose = verbose)) # refit the model with new (L, K, R)
    }
    
    
    
    #################################
    # Step 7. Smoothing out the estimates of latent components and score variances.
    #################################
    
    # Create a convergence indicator if the percent change in  
    # all of the shared and unique covariance spaces is below 10%.
    conv_flag = TRUE 
    
    if(K > 0){
      smooth_est = smooth_est(Phi, Theta, K, tobs = tobs) # G_phi
      Phi = smooth_est[[1]]
      Theta = smooth_est[[2]]
      if(conv_flag & !is.null(G_Phi)){
        conv_flag = sum((G_Phi - smooth_est[[3]])^2) / sum((G_Phi)^2) < .10
        if(verbose){print(paste("G_Phi", round(
          sum((G_Phi - smooth_est[[3]])^2) / sum((G_Phi)^2), digits = 3)))}
      }else{ # if no prior estimate of G_phi
        conv_flag = FALSE
      }
      G_Phi = smooth_est[[3]]
    }else{
      Phi = NULL
      Theta = NULL
    }
    
    if(R > 0){
      smooth_est = smooth_est(Gamma, Lambda, R, tobs = tobs, user_k = 6) # G_gamma
      Gamma = smooth_est[[1]]
      Lambda = smooth_est[[2]]
      if(conv_flag & !is.null(G_Gamma)){
        conv_flag = sum((G_Gamma - smooth_est[[3]])^2) / sum((G_Gamma)^2) < .10
        if(verbose){print(paste("G_Gamma", round(
          sum((G_Gamma - smooth_est[[3]])^2) / sum((G_Gamma)^2), digits = 3)))}
      }else{ # if no prior estimate of G_gamma
        conv_flag = FALSE
      }
      G_Gamma = smooth_est[[3]]
    }else{
      Gamma = NULL
      Lambda = NULL
    }
    
    if(L > 0){
      smooth_est = smooth_est(Psi, Omega_x, L, tobs = tobs) # G^x_psi
      Psi_sm = smooth_est[[1]] # smoothed estimate of Psi
      Omega_x = smooth_est[[2]] # smoothed estimate of Omega^x
      
      if(conv_flag & !is.null(G_Psi_x)){
        conv_flag = sum((G_Psi_x - smooth_est[[3]])^2) / sum((G_Psi_x)^2) < .10
        if(verbose){print(paste("G_Psi_x", round(
          sum((G_Psi_x - smooth_est[[3]])^2) / sum((G_Psi_x)^2), digits = 3)))}
      }else{ # if no prior estimate of G^x_psi
        conv_flag = FALSE
      }
      G_Psi_x = smooth_est[[3]]
      
      smooth_est = smooth_est(Psi, Omega_y, L, tobs = tobs) # G^y_psi
      if(conv_flag & !is.null(G_Psi_y)){
        conv_flag = sum((G_Psi_y - smooth_est[[3]])^2) / sum((G_Psi_y)^2) < .10
        if(verbose){print(paste("G_Psi_y", round(
          sum((G_Psi_y - smooth_est[[3]])^2) / sum((G_Psi_y)^2), digits = 3)))}
      }else{ # if no prior estimate of G^y_psi
        conv_flag = FALSE
      }
      G_Psi_y = smooth_est[[3]]
      
      
      # When L > 1, the correspondence of the smoothed estimate of Psi derived from X 
      # and Y is determined by their MSDD.
      ind_Psi = c()
      for(l_1 in 1:L){
        MSDD_ls = c()
        for(l_2 in 1:L){
          MSDD_ls = append(MSDD_ls, min(trapz(tobs, (smooth_est[[1]][,l_1] 
                                                     - Psi_sm[,l_2])^2),
                                        trapz(tobs, (smooth_est[[1]][,l_1] 
                                                     + Psi_sm[,l_2])^2)))
        }
        if(which.min(MSDD_ls) %in% ind_Psi){
          ind_Psi = append(ind_Psi, NA)
        }else{
          ind_Psi = append(ind_Psi, which.min(MSDD_ls))
        }
      }
      
      if(sum(is.na(ind_Psi)) == 0){
        Omega_y = diag(x = diag(smooth_est[[2]])[ind_Psi], ncol = L)
        Psi = Psi_sm
      }else{ # Psi derived from X is different with that derived from Y, suggesting non-converged result of Psi.
        Omega_y = smooth_est[[2]]
        Psi = Psi_sm 
        conv_flag = FALSE # the estimation continues
      }
    }else{
      Psi = NULL
      Omega_x = NULL
      Omega_y = NULL
    }
    
    # Convergence flag with percent change in covariance spaces is all below 10%.
    if(conv_flag){
      break
    }else{
      stop_flag = stop_flag + 1
    }
  }

  
  ls_res = list(mu_x, mu_y, s2_x, s2_y, Psi, Phi, Gamma, Omega_x, Omega_y, 
                Theta, Lambda)
  return(ls_res)
}



FPCA_vf = function(data_grid, spar = NULL, elbow_plot = FALSE, vf_th = .90, 
                   user_k = NULL, verbose = TRUE){
  
  #############################################################################
  ## Description: (supporting function) This function applies FPCA to the input data and
  ##              the number of retained eigenfunctions is based on the preset variation threshold.
  ## Args:        data_grid: matrix of subjects' data within X or Y (matrix, n_x*N or n_y*N).
  ##              spar: smoothing parameter specified for spline estimation. (scalar, in (0,1], 
  ##                    'NULL' means it is generated by 'generalized' cross-validation).
  ##              elbow_plot: indicator whether to display elbow plots of variation proportions 
  ##                          within FPCA result. (logical)
  ##              vf_th: a threshold to determine the number of eigenfunctions retained in FPCA. (scalar)
  ##              user_k: the dimension of bases used in smoothing covariance surfaces 
  ##                               if specified by users. (integer)
  ##              verbose: indicator whether to display the diagnostic messages. (logical)
  ## Returns:     a list containing:
  ##                mu_e: estimate of the overall mean function (vector, N*1)
  ##                eigenval: an array saves retained eigenvalues. (array)
  ##                eigenfun: a matrix saves retained eigenfunctions. (matrix)
  #############################################################################
  
  # Calculation of smoothed covariance surface
  N <- ncol(data_grid) # total number of time points
  tobs = seq(0, 1, length.out = N)  # total time grid T (vector, N*1).
  n_data = nrow(data_grid) # number of subjects within the input data set
  data_tab <- data.frame(na.omit(cbind(
    rep(1:n_data, each = N), rep(tobs, n_data), as.vector(t(data_grid))))) 
  colnames(data_tab) <- c("id", "timePoints", "data_grid")
  # Calculate the mean
  est_ss <- smooth.spline(data_tab$timePoints, data_tab$data_grid, cv = F, 
                          spar = spar)
  if(verbose){
    print(c("The smoothing parameter of the estimated mean function is ", 
            round(est_ss$spar, digits = 2)))
  }
  
  mu_e = matrix(data = est_ss$y, nrow = length(tobs), ncol = 1)
  data.c <- sweep(data_grid, 2, mu_e)
  if(!is.null(user_k)){
    cov_smooth <- Cov_mat(data.c, tobs = tobs, smooth_2D = TRUE, 
                          user_k = user_k)
  }else{
    cov_smooth <- Cov_mat(data.c, tobs = tobs, smooth_2D = TRUE)
  }
  
  
  # Eigen-decomposition of the smoothed covariance surface.
  eigen_res <- eigen(cov_smooth, symmetric = TRUE)
  eigen_res$values <- eigen_res$values[which(eigen_res$values > 0)]  # retain positive eigenvalues only.
  if(length(eigen_res$values) > 0){
    eigen_res$vectors <- eigen_res$vectors[, 1:length(eigen_res$values)]
  }else{
    stop("No eigenfunctions with positive eigenvalues are returned.")
  }
  
  
  # Normalization
  eigenfun <- NULL # save the final normalized eigen-functions. (matrix)
  eigenval <- NULL # save the corresponding adjusted eigen-values. (array)
  if(length(eigen_res$values) == 1){
    normal.factor <- trapz(tobs, eigen_res$vectors^2)
    eigen_res$vectors <- eigen_res$vectors / sqrt(normal.factor) 
    eigen_res$values <- eigen_res$values * normal.factor 
    eigenfun = eigen_res$vectors 
    eigenval = eigen_res$values
    print(paste("Only one eigenfunction with positive eigenvalue is retained,", 
                "and no elbow plot is displayed."))
  }else{
    for(e in 1:length(eigen_res$values)){
      normal.factor <- trapz(tobs, eigen_res$vectors[,e]^2)
      eigen_res$vectors[,e] <- eigen_res$vectors[,e] / sqrt(normal.factor) 
      eigen_res$values[e] <- eigen_res$values[e] * normal.factor 
    }
    eigenfun = eigen_res$vectors 
    eigenval = eigen_res$values
    
    ## Elbow plot
    if(elbow_plot){
      plot(eigenval[1:min(10, length(eigenval))], type = "b", 
           xlab = "Component number", ylab = "Eigenvalues", lwd = 3,  
           cex.lab = 1.5, cex.axis = 1.5) 
      title(main = list("Elbow plot", cex = 1.5))
    }
    vf_prop = eigenval*100/sum(eigenval) # variation proportion explained by each eigenfunction.
    n_components = match(TRUE, (cumsum(eigenval)/sum(eigenval)) > vf_th) 
    # the number of eigenfunctions being retained based on the preset threshold.
    if(n_components == 1){
      print(paste0("The first eigenfunction has already explained over ", 
                   (vf_th*100), 
                   "% of the total variation within this data set."))
    }
    eigenval = eigenval[1:n_components] 
    eigenfun = eigenfun[,1:n_components]
  }
  
  if(verbose){
    print("------------")
    print(c(paste("The variation fraction explained by the eigenfunctions", 
                  "retained from FPCA is"), 
            paste0(round(vf_prop[1:n_components], digits = 2), "%")))
    print(c("------------"))
  }
  
  return(list(mu_e, eigenval, eigenfun))
}



Cov_mat = function(data.c, tobs, smooth_2D, user_k = NULL){
  
  #############################################################################
  ## Description: (supporting function) This function calculates the empirical or 2D smoothed covariance 
  ##              surface of a given datagrid.
  ## Args:        data.c: matrix of the subjects' data from one group after centralization (matrix, n_x*N or n_y*N).
  ##              tobs: total time grid T (vector, N*1).
  ##              smooth_2D: indicator whether to conduct 2D smoothing process (TRUE: apply 
  ##                         2D smoothing; FALSE: do not apply)
  ##              user_k: the dimension(s) of bases used in smoothing process specified by users. (integer)
  ## Returns:     An empirical or 2D smoothed covariance surface matrix (matrix, N*N).
  #############################################################################
  
  N <- length(tobs) # total number of time points
  n_data = nrow(data.c) # number of subjects
  
  # Calculation of empirical covariance surface
  empCov <- matrix(nrow = N, ncol = N)
  for (t1 in 1:N) {
    for(t2 in 1:N){
      pp = na.omit(data.c[,t1] * data.c[,t2]) 
        # point-wise product of the input data at time t1 and t2 omitting missing values.
      empCov[t1,t2] = sum(pp) / (length(pp) - 1) 
        # calculates the estimated covariance using available data at time t1, t2.
    }
  }
  
  
  # 2D smoothing process (skipped if the empirical one is required)
  if(smooth_2D){
    x0 <- rep(tobs, each = N) 
    x1 <- rep(tobs, N)
     
    # Remove the diagonal elements to avoid effects from errors.
    for(i in 1:nrow(empCov)){
      empCov[i,i] <- NA
    }
    
    # Smooth the covariance surface by 2D penalized smoothing splines, with smoothing 
    # parameters chosen by REML and default dimension of basis function is k = 10.
    if(is.null(user_k)){
      smooth_cov <- gam(as.vector(empCov) ~ te(x0, x1, k = 10, bs = "bs"), 
                        method = "REML")
    }else{
      smooth_cov <- gam(as.vector(empCov) ~ te(x0, x1, k = user_k, 
                                               bs = "bs"), method = "REML")
    }
    
    grids2d <- data.frame(x0 = x0, x1 = x1) 
    covSm <- matrix(predict(smooth_cov, newdata = grids2d), nrow = N) # derive the 2D smoothed estimates.
    cov_return <- (covSm + t(covSm)) / 2 # ensure the returned matrix is symmetric.
  }else{
    cov_return <- empCov
  }
  
  return(cov_return)
}



adap = function(mat, data){
  
  #############################################################################
  ## Description: (supporting function) This function adapt the input functional matrix or vector `mat` 
  ##              as described in Step 3, by retaining their information only at the observed time 
  ##              points of the specified data subject `data` with subject-specific adaptation matrix A_i^x/A_j^y.
  ## Args:        mat: a functional matrix or vector being adapted to non-missing time points of the 
  ##                   input subject data. (matrix, N*(L/K/R))
  ##              data: a generic subject vector X_i or Y_j (N*1 vector with N^x_i/N^y_j non-missing observations)
  ## Returns:     An adapted latent component matrix (matrix, (N^x_i/N^y_j)*(L/K/R)).
  #############################################################################
  
  if(ncol(mat) > 1){
    if(sum(!is.na(data)) == 1){ 
      mat_adp <- t(mat[!is.na(data),]) # return a 1*(L/K/R) matrix
    }else{
      mat_adp <- mat[!is.na(data),] 
    }
  }else{
    mat_adp <- mat[!is.na(data)]
  }
  return(mat_adp)
}



smooth_est = function(latent_mat, scoreVar, component, tobs, 
                       user_k = NULL){
  
  #############################################################################
  ## Description: (supporting function) This function derives the smoothed estimates of latent components 
  ##              and score variances from the estimated shared/unique covariance spaces in Step 7 of Algorithm 1.
  ## Args:        latent_mat: a non-smoothed shared/unique latent component matrix (matrix, N*(L/K/R)).
  ##              scoreVar: a non-smoothed score variance matrix (matrix, L*L/K*K/R*R)
  ##              component: number of shared/unique latent components (integer)
  ##              tobs: total time grid T (vector, N*1).
  ##              user_k: the dimension(s) of bases used in smoothing process. (integer)
  ## Returns:     a list containing:
  ##                latent_mat_sm: smoothed estimate of the input latent component matrix. (matrix, same 
  ##                               dimension as `latent_mat`)
  ##                scoreVar_sm: smoothed estimate of the input score variance matrix. (matrix, same  
  ##                             dimension as `scoreVar`)
  ##                Ghat_sm: estimated shared/unique covariance space. (matrix, N*N)
  #############################################################################
  
  # Reconstruct the shared/unique covariance space (hat(G)^x_psi, hat(G)^y_psi, hat(G)_phi or hat(G)_gamma)
  # using the input non-smoothed latent component and score variance matrices.
  Ghat = latent_mat %*% scoreVar %*% t(latent_mat) 
  
  
  # Conduct FPCA on the derived covariance space and retain eigenfunctions and eigenvalues based on the
  # preset component number (L, K or R).
  N <- length(tobs)
  x0 <- rep(tobs, each = N)
  x1 <- rep(tobs, N)
  
  ## 2D smoothing of the input covariance surface
  if(is.null(user_k)){
    smooth_cov <- gam(as.vector(Ghat) ~ te(x0, x1, bs = "bs", k = 5, 
                                           sp = c(-1, -1)), method = "REML")
  }else{
    smooth_cov <- gam(as.vector(Ghat) ~ te(x0, x1, bs = "bs", 
                                           k = user_k, sp = c(-1, -1)), 
                      method = "REML")
  }
  
  grids2d <- data.frame(x0 = x0, x1 = x1)
  covmat <- matrix(predict(smooth_cov, newdata = grids2d), nrow = N)
  cov_smooth <- (covmat + t(covmat)) / 2 #  Symmetrize covariance surface
  
  ## Eigenfunction decomposition
  eigen_res <- eigen(cov_smooth, symmetric = TRUE)
  eigen_res$values <- eigen_res$values[which(eigen_res$values > 0)]  # keep positive eigenvalues only
  eigen_res$vectors <- eigen_res$vectors[, 1:length(eigen_res$values)] # keep corresponding eigenfunctions.
  if(length(eigen_res$values) < component){
    stop("The number of required components is greater than the number of  
         returned eigenfunctions.")
  }else{
    ## Normalize the result on the given times T.
    for(e in 1:component){
      normal.factor <- trapz(tobs, eigen_res$vectors[, e]^2) 
      eigen_res$vectors[, e] <- (eigen_res$vectors[, e] / sqrt(normal.factor)) 
      eigen_res$values[e] <- eigen_res$values[e] * normal.factor 
    }
  }
  
  if(component > 1){
    latent_mat_sm = eigen_res$vectors[, 1:component]
    scoreVar_sm = diag(x = eigen_res$values[1:component])
  }else{
    latent_mat_sm = matrix(data = eigen_res$vectors[, 1:component], ncol = 1)
    scoreVar_sm = matrix(data = eigen_res$values[1:component], ncol = 1)
  }
  
  # Reconstruct the shared/unique covariance space using derived smoothed estimates.
  Ghat_sm = latent_mat_sm %*% scoreVar_sm %*% t(latent_mat_sm) 
  
  return(list(latent_mat_sm, scoreVar_sm, Ghat_sm))  
}



