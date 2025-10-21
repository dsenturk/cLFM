###############################################################################
## Description: Functions for generating contrastive functional data sets with 
##              latent components and score variances of Scenario 1 described in   
##              the simulation section in 'Contrastive Latent Functional Model'.
###############################################################################

data_gen = function(n_x, n_y, s2_x, s2_y, case = 1){
  
  #############################################################################
  ## Description: This function generates the data sets described in the first 
  ##              simulation scenario of the simulation section of the paper.
  ## Args:        n_x: number of subjects within the first data set (integer).
  ##              n_y: number of subjects within the second data set (integer).
  ##              s2_x: error variance for the first data set, sigma^2_X. (scalar)
  ##              s2_y: error variance for the second data set, sigma^2_Y. (scalar)
  ##              case: index of simulation scenarios. (integer)
  ## Returns:     A list containing generated data pairs (X, Y) (matrix, n_x*30, n_y*30).
  #############################################################################
  
  N = 30 # total number of time points
  tobs = seq(0, 1, length.out = N) # vector of total time points
  
  # Construct the shared/unique latent component matrices
  if(case <= 2){
    psi = suppressWarnings(create.fourier.basis(nbasis = 3, dropind = c(1,2))) # psi_1(t) = sqrt(2)*cos(2pi*t)
    Psi = eval.basis(tobs, psi)
    
    phi = suppressWarnings(create.fourier.basis(nbasis = 3, dropind = c(1,3))) # phi_1(t) = sqrt(2)*sin(2pi*t)
    Phi = eval.basis(tobs, phi)
    
    gamma = suppressWarnings(create.fourier.basis(nbasis = 5, 
                                                  dropind = c(1,2,3,5))) # gamma_1(t) = sqrt(2)*sin(4pi*t)
    Gamma = eval.basis(tobs, gamma)
  }else if(case == 3){
    psi = suppressWarnings(create.fourier.basis(nbasis = 3, dropind = c(1))) # psi_1(t) = sqrt(2)*cos(2pi*t), psi_2(t) = sqrt(2)*sin(2pi*t)
    Psi = eval.basis(tobs, psi)
  }else if(case == 4){
    phi = suppressWarnings(create.fourier.basis(nbasis = 5, dropind = c(1,3,4))) # phi_1(t) = sqrt(2)*sin(2pi*t), phi_2(t) = sqrt(2)*cos(4pi*t)
    Phi = eval.basis(tobs, phi)
    
    gamma = suppressWarnings(create.fourier.basis(nbasis = 5, 
                                                  dropind = c(1,2,5))) # gamma_1(t) = sqrt(2)*cos(2pi*t), gamma_2(t) = sqrt(2)*sin(4pi*t)
    Gamma = eval.basis(tobs, gamma)
  }else{
    stop("Wrong scenario index.")
  }
  
  
  # Construct the score variance matrices
  if(case == 1){
    Omega_x = matrix(14)
    Omega_y = matrix(11.5)
    Theta = matrix(6)
    Lambda = matrix(8.5)
    L = 1 # dimension of the shared space
    K = 1 # dimension of the unique space of the first group of data X
    R = 1 # dimension of the unique space of the second group of data Y
  }else if(case == 2){
    Omega_x = matrix(9)
    Omega_y = matrix(9)
    Theta = matrix(14)
    Lambda = matrix(14)
    L = 1 
    K = 1 
    R = 1
  }else if(case == 3){
    Omega_x = diag(c(14, 6))
    Omega_y = diag(c(8, 12))
    L = 2 
    K = 0 
    R = 0
  }else if(case == 4){
    Theta = diag(c(14, 6))
    Lambda = diag(c(12, 8))
    L = 0
    K = 2 
    R = 2
  }else{
    stop("Wrong case index.")
  }
  
  
  # Generate X
  mu_x = matrix(data = sapply(tobs, 
                              function(x){return(4.5 * x + 1.5 * exp(x ^ 2))}), 
                nrow = N, ncol = 1) # generate the overall mean vector mu^x 
  X = matrix(data = NA, nrow = n_x, ncol = N) # create an empty matrix to store generated data
  for (i in 1 : n_x) {
    if(L > 0){
      eta = matrix(data = rnorm(L, mean = rep(0, L), 
                                sd = sqrt(diag(Omega_x))), ncol = 1) # subject-specific score of shared components, eta^x_i
    }
    if(K > 0){
      xi = matrix(data = rnorm(K, mean = rep(0, K), 
                               sd = sqrt(diag(Theta))), ncol = 1) # subject-specific score of unique components, xi_i
    }
    epsilon = matrix(data = rnorm(N, mean = 0, sd = sqrt(s2_x)), ncol = 1) # noise term, epsilon^x_i
    if(L > 0 & K > 0){
      X[i,] = t(mu_x + Psi %*% eta + Phi %*% xi + epsilon)
    }else if(L > 0){
      X[i,] = t(mu_x + Psi %*% eta + epsilon)
    }else if(K > 0){
      X[i,] = t(mu_x + Phi %*% xi + epsilon)
    }else{
      stop("No components are used in generating data X.")
    }
  }
  
  
  # Generate Y
  mu_y = matrix(data = sapply(tobs, function(x){
    return(-10 * x + 5 * exp(-50 * (x - 0.7) ^ 2) + 5)}), nrow = N, ncol = 1) # generate the overall mean vector mu^y 
  Y = matrix(data = NA, nrow = n_y, ncol = N)
  for (j in 1 : n_y) {
    if(L > 0){
      eta = matrix(data = rnorm(L, mean = rep(0, L), 
                                sd = sqrt(diag(Omega_y))), ncol = 1) # subject-specific score of shared components, eta^y_j
    }
    if(R > 0){
      zeta = matrix(data = rnorm(R, mean = rep(0, R), 
                                 sd = sqrt(diag(Lambda))), ncol = 1) # subject-specific score of unique components, zeta_j
    }
    epsilon = matrix(data = rnorm(N, mean = 0, sd = sqrt(s2_y)), ncol = 1) # noise term, epsilon^y_j
    if(L > 0 & R > 0){
      Y[j,] = t(mu_y + Psi %*% eta + Gamma %*% zeta + epsilon)
    }else if(L > 0){
      Y[j,] = t(mu_y + Psi %*% eta + epsilon)
    }else if(R > 0){
      Y[j,] = t(mu_y + Gamma %*% zeta + epsilon)
    }
  }
  
  
  # Add missing data into the generated data pairs
  missing_prop = 0.2 # predefined missing proportion within the generated datasets. (proportion, from 0 to 1)
  if(missing_prop > 0){
    for (ind in 1:n_x) {
      time_missing = sample(1 : N, size = (round(N * missing_prop)), 
                            replace = FALSE) # sample time points with a preset proportion
      X[ind, time_missing] = NA # replace the generated data with missing values at the selected time points.
    }
    
    for (ind in 1:n_y) {
      time_missing = sample(1 : N, size = (round(N * missing_prop)), 
                            replace = FALSE)
      Y[ind, time_missing] = NA
    }
  }
  
  return(list(X, Y))
}
