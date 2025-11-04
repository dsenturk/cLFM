###############################################################################
## Description: Functions for implementing cLVM modeling with the trimming-refinement process   
##              described in Algorithm 2 of 'Contrastive Latent Functional Model'.
###############################################################################

cLFM = function(X, Y, alpha_L, alpha_KR, alpha_refine, verbose = TRUE){
  
  #############################################################################
  ## Description: This function constructs the contrastive latent functional model by running the EM estimation (Algorithm 1)
  ##              and applying trimming-refinement process (Algorithm 2) to the three quantities (L, K, R) based on shared and
  ##              unique variation proportions. 
  ## Args:        X: data grid of the first group of data (matrix, n_x*N)
  ##              Y: data grid of the second group of data (matrix, n_y*N)
  ##              alpha_L: threshold of variation used in the trimming stage to retain the shared component estimates 
  ##                       within the two shared covariance spaces (scalar, (0, 1))
  ##              alpha_KR: threshold of variation used in the trimming stage to retain the unique component  
  ##                        estimates within the two total variation spaces (scalar, (0, 1))
  ##              alpha_refine: threshold of variation used in the refinement stage to retain the shared component  
  ##                            estimates within the two total covariance spaces (scalar, (0, 1))
  ##              verbose: indicator whether to display the EM estimation progress, and the
  ##                       trimming-refinement details of (L, K, R). (logical)
  ## Returns:     A result list by fitting the proposed cLFM containing:
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
  
  # Initialize the dimensions of shared and unique variation spaces based on group FPCA results.
  K = length(FPCA_vf(X, verbose = verbose)[[2]])
  R = length(FPCA_vf(Y, verbose = verbose)[[2]])
  L = min(K, R)
  user_tol = 1e-3
  
  
  results = run_em(X, Y, L, K, R, tol = user_tol, verbose = verbose) 
  
  Phi = results[[6]] 
  Gamma = results[[7]]
  Omega_x = results[[8]] 
  Omega_y = results[[9]] 
  Theta = results[[10]]
  Lambda = results[[11]]
  
  
  #################################
  # Trimming stage
  #################################
  
  stop_flag = 0
  
  repeat{
    # Read the untrimmed values of (L, K, R)
    L = ifelse(is.null(Omega_x), 0, nrow(Omega_x)) 
    K = ifelse(is.null(Theta), 0, nrow(Theta))
    R = ifelse(is.null(Lambda), 0, nrow(Lambda))
    
    
    # Firstly, we trim the quantity L with respect to the two shared covariance spaces.
    if(L <= 1){ # when L = 0 or 1, we skip trimming L and further adjustments are delayed to the refinement stage.
      L_new = L
      K1 = K # intermediate values of K
      R1 = R # intermediate values of R
    }else{
      ## Trim L within the shared variation space of the first group of data X.
      ### 1. calculate the variation constitution
      scoreVar = diag(Omega_x) # extract the score variances from the diagonal of the matrix Omega^x. 
      vf = scoreVar * 100 / sum(scoreVar)
      vf_print = round(vf[1], digits = 2)
      if(length(vf) > 1){
        for(ind in 2:length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the estimated",
                    "shared components within the shared space of X are:", 
                    vf_print))
      }
      
      
      ### 2. find the least number of shared components that explain over `alpha_L` of the shared variation within X
      scoreVar_s <- sort(scoreVar, decreasing = TRUE) # sort the score variances in a decreasing order
      L_x = match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > alpha_L) 
      
      ## Trim L within the shared variation space of the second group of data Y.
      ### 1. Calculate the variation constitution
      scoreVar = diag(Omega_y) # extract the score variances from the diagonal of the matrix Omega^y. 
      vf = scoreVar * 100 / sum(scoreVar)
      vf_print = round(vf[1], digits = 2)
      if(length(vf) > 1){
        for(ind in 2:length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the estimated", 
                    "shared components within the shared space of Y are:", 
                    vf_print))
      }
      
      
      ### 2. Find the least number of shared components that explain over `alpha_L` of the shared variation within Y
      scoreVar_s <- sort(scoreVar, decreasing = TRUE) # sort the score variances in a decreasing order
      L_y = match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > alpha_L) 
      
      ## Update L with the minimal value among two shared variation spaces.
      L_new = min(L_x, L_y)
      if(L_new < L){# if L is trimmed, we refit the cLFM with new (L, K, R).
        if(verbose){
          print(paste0("L is trimmed based on the constitution of the two ",
                       "shared variation spaces. The cLFM is refitted with (", 
                       L_new, ", ", K, ", ", R, ")."))
        }
        results = run_em(X, Y, L_new, K, R, tol = user_tol, verbose = verbose)
        Phi = results[[6]]
        Gamma = results[[7]]
        Omega_x = results[[8]]
        Omega_y = results[[9]]
        Theta = results[[10]] 
        Lambda = results[[11]]
        L_new = ifelse(is.null(Omega_x), 0, nrow(Omega_x)) 
        # K1 and R1 are used to capture the adjusted values of K, R in case of convergence 
        # failure of EM estimation with (L_new, K, R)
        K1 = ifelse(is.null(Theta), 0, nrow(Theta))  
        R1 = ifelse(is.null(Lambda), 0, nrow(Lambda))
      }else{# no updates on L, i.e. L_new = L
        K1 = K
        R1 = R
      } 
    }

    # Secondly, we trim the quantities K and R with respect to the two total covariance spaces.
    ## Adjustment on K
    if(K1 > 0){ # skipped if current K = 0
      if(L_new > 0){
        scoreVar = diag(Omega_x) # extract the score variances for the shared components from Omega^x. 
      }else{
        scoreVar = c()
      }
      scoreVar_total = append(scoreVar, diag(Theta))
      vf = scoreVar_total * 100 / sum(scoreVar_total)
      vf_print = round(vf[L_new + 1], digits = 2)
      if(K1 > 1){
        for(ind in (L_new + 2):length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the", 
                    "estimated unique components within the total variation", 
                    "space of X are:", vf_print))
      }
      
      ## Find the least number of unique components that explain over `alpha_KR` of the total 
      ## variation within the first data set X
      scoreVar_s = append(scoreVar, sort(diag(Theta), decreasing = TRUE))
      K_new = match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > alpha_KR) - 
        L_new
      if(K_new < 0){
        K_new <- 0
      }
    }else{
      K_new = 0
    }
    
    
    # Adjustment on R
    if(R1 > 0){ # skipped if current R = 0
      if(L_new > 0){
        scoreVar = diag(Omega_y) # extract the score variances for the shared components from Omega^x. 
      }else{
        scoreVar = c()
      }
      scoreVar_total = append(scoreVar, diag(Lambda))
      vf = scoreVar_total * 100 / sum(scoreVar_total)
      vf_print = round(vf[L_new + 1], digits = 2)
      if(R1 > 1){
        for(ind in (L_new + 2):length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the", 
                    "estimated unique components within the total variation", 
                    "space of Y are:", vf_print))
      }
      
      
      ## Find the least number of unique components that explain over `alpha_KR` of the total variation within 
      ## the second data set Y
      scoreVar_s = append(scoreVar, sort(diag(Lambda), decreasing = TRUE))
      R_new = match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > alpha_KR) - 
        L_new
      if(R_new < 0){
        R_new <- 0
      }
    }else{
      R_new <- 0
    }
    
    
    ## Detect whether K or R are adjusted or not.
    if(K_new < K1 | R_new < R1){
      if(verbose){
        print(paste0("K and R are trimmed based on the constitution of total ", 
                     "variation spaces of the two groups. The cLFM is refitted",
                     "with (", L_new, ", ", K_new, ", ", R_new, ")."))
      }
      
      results = run_em(X, Y, L_new, K_new, R_new, tol = user_tol, 
                       verbose = verbose)
      Phi = results[[6]]
      Gamma = results[[7]]
      Omega_x = results[[8]]
      Omega_y = results[[9]]
      Theta = results[[10]] 
      Lambda = results[[11]]
      L_new = ifelse(is.null(Omega_x), 0, nrow(Omega_x)) 
      K_new = ifelse(is.null(Theta), 0, nrow(Theta))
      R_new = ifelse(is.null(Lambda), 0, nrow(Lambda))
    }
    
    
    # Stop the iterative trimming process until no more changes on (L, K, R).
    if(L_new == L & K_new == K & R_new == R){
      if(verbose){print("Trimming process is done.")}
      break
    }else if(stop_flag > 20){
      stop("Failure in the trimming process.")
    }else{
      stop_flag = stop_flag + 1
    }
  }

  
  #################################
  # Refinement stage
  #################################

  # In the refinement stage, the three dimensions L, K and R get further adjusted by moving
  # comparable unique component pairs to the shared space and trimming all components in 
  # the total variation space.
  stop_flag = 0
  repeat{
    # Detection of comparable unique component pairs (phi_k(t), gamma_r(t))
    ## Once the MSDD of a pair of unique components from the two unique spaces 
    ## is below 10%, the current detection will be terminated and the model will be 
    ## refitted with the adjusted values of new L, K and R. Then the refinement will 
    ## be progressed with the new fitted results.
    L1 = ifelse(is.null(Omega_x), 0, nrow(Omega_x)) # intermediate values of L
    K1 = ifelse(is.null(Theta), 0, nrow(Theta)) # intermediate values of K
    R1 = ifelse(is.null(Lambda), 0, nrow(Lambda)) # intermediate values of R
    if(K1 > 0 & R1 > 0){ # skip it if any unique space is empty
      N <- nrow(Phi) # total number of time points
      tobs = seq(0, 1, length.out = N)  # total time grid T (vector, N*1).
      
      refit = FALSE ## indicator signaling if refitting the cLFM is needed before next refinement
      for (k in 1:K1) {
        for (r in 1:R1) {
          MSDD_unique = min(trapz(tobs, (Phi[,k] - Gamma[,r])^2), 
                            trapz(tobs, (Phi[,k] + Gamma[,r])^2))
          if(MSDD_unique < .15){
            L1 = L1 + 1
            K1 = K1 - 1
            R1 = R1 - 1
            if(verbose){
              print(paste0("MSDD(phi", k, "(t), gamma", r, "(t)) is ", 
                           round(MSDD_unique, digits = 2), ", lower than 15%."))
            }
            refit = TRUE
            break
          }
        }
        if(refit){break} 
      }
      if(refit){
        if(verbose){
          print(paste0("L, K and R are refined due to comparable unique ", 
                       "component pairs. The cLFM is refitted with (", L1, 
                       ", ", K1, ", ", R1, ")."))
        }
        results = run_em(X, Y, L1, K1, R1, tol = user_tol, verbose = verbose)
        Phi = results[[6]]
        Gamma = results[[7]]
        Omega_x = results[[8]]
        Omega_y = results[[9]]
        Theta = results[[10]] 
        Lambda = results[[11]]
        next
      }else{}
    }else{
      if(verbose){print("At least one unique latent component space is empty.")}
    }
    
    
    # Refinement on the shared latent component estimates within the total variation space.
    
    # Total variation space of X
    L_x_adj = FALSE
    if(L1 > 0 & K1 > 0){
      scoreVar = c(diag(Omega_x), diag(Theta))
      vf = scoreVar * 100 / sum(scoreVar)
      vf_print = round(vf[1], digits = 2)
      if((L1 + K1) > 1){
        for(ind in 2:length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), 
                           sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the", 
                    "estimated shared and unique components within X are:", 
                    vf_print))
      }
      
      ## Find the least number of latent components that explain over `alpha_refine` of 
      ## the total variation within X
      scoreVar_s = c(diag(Theta), sort(diag(Omega_x), decreasing = TRUE))
      num_component <- match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > 
                               alpha_refine)
      if(num_component < (L1 + K1)){L_x_adj = TRUE}
    }else if(L1 > 0 & K1 == 0){
      scoreVar = diag(Omega_x)
      vf = scoreVar * 100 / sum(scoreVar)
      vf_print = round(vf[1], digits = 2)
      if(L1 > 1){
        for(ind in 2:length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), 
                           sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the", 
                    "estimated shared and unique components within X are:", 
                    vf_print))
      }
      
      ## Find the least number of latent components that explain over `alpha_refine` of 
      ## the total variation within X
      scoreVar_s = sort(scoreVar, decreasing = TRUE)
      num_component <- match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > 
                               alpha_refine)
      if(num_component < L1){L_x_adj = TRUE}
    }else{
      if(verbose){print("No shared components are estimated for the data X.")}
    }
    
    
    
    # Total variation space of Y
    L_y_adj = FALSE
    if(L1 > 0 & R1 > 0){
      scoreVar = c(diag(Omega_y), diag(Lambda))
      vf = scoreVar * 100 / sum(scoreVar)
      vf_print = round(vf[1], digits = 2)
      if((L1 + R1) > 1){
        for(ind in 2:length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), 
                           sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the", 
                    "estimated shared and unique components within Y are:", 
                    vf_print))
      }
      
      ## Find the least number of latent components that explain over `alpha_refine` of 
      ## the total variation within Y
      scoreVar_s = c(diag(Lambda), sort(diag(Omega_y), decreasing = TRUE))
      num_component <- match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > 
                               alpha_refine)
      if(num_component < (L1 + R1)){L_y_adj = TRUE}
    }else if(L1 > 0 & R1 == 0){
      scoreVar = diag(Omega_y)
      vf = scoreVar * 100 / sum(scoreVar)
      vf_print = round(vf[1], digits = 2)
      if(L1 > 1){
        for(ind in 2:length(vf)){
          vf_print = paste(vf_print, round(vf[ind], digits = 2), 
                           sep = ", ")
        }
      }
      if(verbose){
        print(paste("The variation proportions explained by the", 
                    "estimated shared and unique components within Y are:", 
                    vf_print))
      }
      
      ## Find the least number of latent components that explain over `alpha_refine` of 
      ## the total variation within Y
      scoreVar_s = sort(scoreVar, decreasing = TRUE)
      num_component <- match(TRUE, (cumsum(scoreVar_s) / sum(scoreVar_s)) > 
                               alpha_refine)
      if(num_component < L1){L_y_adj = TRUE}
    }else{
      if(verbose){
        print("No shared components are estimated for the data Y.")
      }
    }
    
    
    # Adjust L, K, R accordingly
    if(L_x_adj){
      L1 <- L1 - 1
      if(!L_y_adj){ # special case, when the disregarded shared component does explain enough variation in Y 
        R1 <- R1 + 1
      }
    }else if(L_y_adj){ # L_x_adj is FALSE
      L1 <- L1 - 1
      K1 <- K1 + 1 # special case for X
    }else{# L_x_adj and L_y_adj are both FALSE
      if(verbose){
        print("Refinement stage is finished.")
        print(paste0(L1, " shared latent components and ", K1, " and ", R1, 
                     " unique latent components for respective X and Y are ", 
                     "estimated by cLFM."))
      }
      return(results)
    }
    
    
    if(verbose){
      print(paste0("L, K and R are further refined within the total variation", 
                   "space, and the cLFM is refitted with (", L1, 
                   ", ", K1, ", ", R1, ")."))
    }
    results = run_em(X, Y, L1, K1, R1, tol = user_tol, verbose = verbose)
    Phi = results[[6]]
    Gamma = results[[7]]
    Omega_x = results[[8]]
    Omega_y = results[[9]]
    Theta = results[[10]] 
    Lambda = results[[11]]
    L1 = ifelse(is.null(Omega_x), 0, nrow(Omega_x)) 
    K1 = ifelse(is.null(Theta), 0, nrow(Theta))
    R1 = ifelse(is.null(Lambda), 0, nrow(Lambda))
    
    
    if(stop_flag > 20){ # if L, K, R do not stop changing, the loop breaks after 20 iterations.
      stop("Failure in the refinement process.")
    }else{
      stop_flag = stop_flag + 1
      next
    }
  }
}
