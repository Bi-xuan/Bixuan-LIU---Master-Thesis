### Functions

# IHT

# grad <- -4 * t((Sigma %*% Lambda_0) %*% (2*(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n)) -  diag(diag(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega*diag(n)))))

# ISTA update with stepsize eta
ft_IHT_update <- function(Lambda_0,Sigma,omega,gma,eta) {
  grad <- -2 * ((Sigma %*% Lambda_0) %*% (Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n))) / frobenius.norm(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n))
  grad_step <- Lambda_0 - eta * grad
  Lambda_updt <- grad_step
  
  mask_diag <- diag(n) == 1        # diagonal entries
  mask_off <- !mask_diag
  
  
  mask_keep <- abs(grad_step) > sqrt(2 * eta * gma)
  # mask_keep[1,3] <- TRUE
  Lambda_updt[mask_off & !mask_keep] <- 0
  
  return (list(Lambda_updt,grad))
}

ft_IHT_bktr <- function(Lambda_0,Sigma,omega,gma,eta_0, eta_inc) {
  # eta_0 is the initial stepsize (before backtracking)
  eta <- list(eta_0)
  
  Lambda <- list(Lambda_0)
  obj <- list(frobenius.norm(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n)) + gma * (sum(abs(c(Lambda_0)) >0)-n))
  obj_gap <- list(obj[[1]])
  
  for (i in 1:1000) {
    eta_t <- eta[[i]]
    l <- 1
    repeat {
      IHT_updt <- ft_IHT_update(Lambda[[i]], Sigma, omega, gma, eta_t)
      Lambda_updt <- IHT_updt[[1]]
      grad_updt <- IHT_updt[[2]]
      
      lhs <- frobenius.norm(Sigma - t(Lambda_updt) %*% Sigma %*% Lambda_updt - omega * diag(n))
      rhs <- frobenius.norm(Sigma - t(Lambda[[i]]) %*% Sigma %*% Lambda[[i]] - omega * diag(n)) + sum((Lambda_updt - Lambda[[i]]) * grad_updt) + (frobenius.norm(Lambda_updt - Lambda[[i]])^2) / (2*eta_t)
      
      if (lhs <= rhs) break
      eta_t <- eta_t / eta_inc
      
      l <- l+1
      if(l > 1000) break
    }
    eta[[i+1]] <- eta_t
    
    IHT_updt <- ft_IHT_update(Lambda[[i]],Sigma,omega,gma,eta[[i+1]])
    Lambda[[i+1]] <- IHT_updt[[1]]
    obj[[i+1]] <- frobenius.norm(Sigma - t(Lambda[[i+1]])%*%Sigma%*%Lambda[[i+1]] - omega * diag(n)) + gma * (sum(abs(c(Lambda[[i+1]]) )>0)-n)
    obj_gap[[i+1]] <- abs(obj[[i+1]] - obj[[i]])
    
    if(obj_gap[[i+1]] < 1e-6) {
      return(list(Lambda[[i+1]],obj[[i+1]],TRUE))
    }
  }
  # print("IHT did not converge within 1000 iterations. Objective gap:")
  # print(obj_gap[[i+1]])
  return(list(Lambda[[i+1]],obj[[i+1]],FALSE))
}


# get_prop_IHT <- function (n,Sigma, omega, Lambda_0, num_exp, eta, eta_inc) {
#   res_pos <- 0
# 
#   for (j in 1:num_exp) {
#     Lambda <- list(matrix(runif(n**2,-1,1),n,n))
#     gma <- list(0)
#     err_edges <- integer(0)
#     res_IHT_bktr <- ft_IHT_bktr(Lambda[[1]],Sigma,omega,gma[[1]],eta,eta_inc)
#     #if(!res_IHT_bktr[[3]]){
#     #  return(list(res_pos / num_exp,FALSE))
#     #}
#     Lambda[[1]] <- res_IHT_bktr[[1]]
#     err_edges[1] <- sum(xor(Lambda_0 != 0, Lambda[[1]] != 0))
# 
#     i <- 1
#     repeat {
#       Lambda[[i+1]] <- Lambda[[i]]
#       t_gma <- gma[[i]]
#       while (all((Lambda[[i]]!=0) == (Lambda[[i+1]]!=0))) {
#         t_gma <- t_gma + 0.01
#         res_IHT_bktr <- ft_IHT_bktr(Lambda[[i]],Sigma,omega,t_gma,eta,eta_inc)
#         # if(!res_IHT_bktr[[3]]){
#         #  return(list(res_pos / num_exp,FALSE))
#         #}
#         Lambda[[i+1]] <- res_IHT_bktr[[1]]
#       }
#       gma[[i+1]] <- t_gma
#       err_edges[i+1] <- sum(xor(Lambda_0 != 0, Lambda[[i+1]] != 0))
#       if(nnzero(Lambda[[i+1]]) == n)  break
#       i <- i+1
#     }
# 
#     if (0 %in% err_edges) {
#       res_pos <- res_pos+1
#     }
#   }
# 
#   return(list(res_pos / num_exp,TRUE))
# }

get_roc_IHT <- function (n,Sigma, omega, Lambda_0, num_exp, eta, eta_inc) {
  TPR <- numeric()
  FPR <- numeric()
  k <- 1 # index for TPR and FPR
  
  off_diag <- row(Lambda_0) != col(Lambda_0)
  num_pos <- sum((Lambda_0 != 0)[off_diag])
  num_neg <- sum((Lambda_0 == 0)[off_diag])
  
  for (j in 1:num_exp) {
    Lambda <- list(matrix(runif(n**2,-1,1),n,n))
    gma <- list(0)
    res_IHT_bktr <- ft_IHT_bktr(Lambda[[1]],Sigma,omega,gma[[1]],eta,eta_inc)
    Lambda[[1]] <- res_IHT_bktr[[1]]
    
    TP <- sum((Lambda_0 != 0 & Lambda[[1]] != 0)[off_diag])
    FP <- sum((Lambda_0 == 0 & Lambda[[1]] != 0)[off_diag])
    TPR[k] <- TP / num_pos
    FPR[k] <- FP / num_neg
    k <- k + 1
    
    if(nnzero(Lambda[[1]]) == n) break
    
    i <- 1
    repeat {
      Lambda[[i+1]] <- Lambda[[i]]
      t_gma <- gma[[i]]
      while (all((Lambda[[i]]!=0) == (Lambda[[i+1]]!=0))) {
        t_gma <- t_gma + 0.01
        res_IHT_bktr <- ft_IHT_bktr(Lambda[[i]],Sigma,omega,t_gma,eta,eta_inc)
        Lambda[[i+1]] <- res_IHT_bktr[[1]]
      }
      gma[[i+1]] <- t_gma
      
      TP <- sum((Lambda_0 != 0 & Lambda[[i+1]] != 0)[off_diag])
      FP <- sum((Lambda_0 == 0 & Lambda[[i+1]] != 0)[off_diag])
      TPR[k] <- TP / num_pos
      FPR[k] <- FP / num_neg
      k <- k+1
      
      if(nnzero(Lambda[[i+1]]) == n) break
      
      i <- i+1
    }
  }
  
  return(list(TPR,FPR))
}


# ISTA
# grad <- -4 * t((Sigma %*% Lambda_0) %*% (2*(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n)) -  diag(diag(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega*diag(n)))))

# ISTA update with stepsize eta
ft_ISTA_update <- function(Lambda_0,Sigma,omega,gma,eta) {
  grad <- -2 * ((Sigma %*% Lambda_0) %*% (Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n))) / frobenius.norm(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n))
  grad_step <- Lambda_0 - eta * grad
  Lambda_updt <- grad_step
  
  mask_diag <- diag(n) == 1
  Lambda_updt[!mask_diag & grad_step < -eta*gma] <- grad_step[!mask_diag & grad_step < -eta*gma] + eta*gma
  Lambda_updt[!mask_diag & grad_step >  eta*gma] <- grad_step[!mask_diag & grad_step >  eta*gma] - eta*gma
  Lambda_updt[!mask_diag & abs(grad_step) <= eta*gma] <- 0
  
  return (list(Lambda_updt,grad))
}

ft_ISTA_bktr <- function(Lambda_0,Sigma,omega,gma,eta_0, eta_inc) {
  # eta_0 is the initial stepsize (before backtracking)
  eta <- list(eta_0)
  
  Lambda <- list(Lambda_0)
  obj <- list(frobenius.norm(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n)) + gma * sum(abs(c(Lambda_0) - c(diag(diag(Lambda_0))))))
  obj_gap <- list(obj[[1]])
  
  for (i in 1:1000) {
    eta_t <- eta[[i]]
    l <- 1
    repeat {
      ISTA_updt <- ft_ISTA_update(Lambda[[i]], Sigma, omega, gma, eta_t)
      Lambda_updt <- ISTA_updt[[1]]
      grad_updt <- ISTA_updt[[2]]
      
      lhs <- frobenius.norm(Sigma - t(Lambda_updt) %*% Sigma %*% Lambda_updt - omega * diag(n))
      rhs <- frobenius.norm(Sigma - t(Lambda[[i]]) %*% Sigma %*% Lambda[[i]] - omega * diag(n)) + sum((Lambda_updt - Lambda[[i]]) * grad_updt) + (frobenius.norm(Lambda_updt - Lambda[[i]])^2) / (2*eta_t)
      
      if (lhs <= rhs) break
      eta_t <- eta_t / eta_inc
      
      l <- l+1
      if(l>1000) break
    }
    eta[[i+1]] <- eta_t
    
    ISTA_updt <- ft_ISTA_update(Lambda[[i]],Sigma,omega,gma,eta[[i+1]])
    Lambda[[i+1]] <- ISTA_updt[[1]]
    obj[[i+1]] <- frobenius.norm(Sigma - t(Lambda[[i+1]])%*%Sigma%*%Lambda[[i+1]] - omega * diag(n)) + gma * sum(abs(c(Lambda[[i+1]]) - c(diag(diag(Lambda[[i+1]])))))
    obj_gap[[i+1]] <- abs(obj[[i+1]] - obj[[i]])
    
    if(obj_gap[[i+1]] < 1e-4) {
      return(list(Lambda[[i+1]],obj[[i+1]],TRUE))
    }
  }
  #print("ISTA did not converge within 10000 iterations. Objective gap:")
  #print(obj_gap[[i+1]])
  return(list(Lambda[[i+1]],obj[[i+1]],FALSE))
}

# get_prop_ISTA <- function (n,Sigma, omega, Lambda_0, num_exp, eta, eta_inc) {
#   res_pos <- 0
#   res_pos_t <- 0
#   
#   for (j in 1:num_exp) {
#     Lambda <- list(matrix(runif(n**2,-1,1),n,n))
#     gma <- list(0)
#     err_edges <- integer(0)
#     err_edges_t <- integer(0)
#     res_ISTA_bktr <- ft_ISTA_bktr(Lambda[[1]],Sigma,omega,gma[[1]],eta,eta_inc)
#     # if(!res_ISTA_bktr[[3]]){
#     #  return(list(res_pos / num_exp,FALSE))
#     # }
#     Lambda[[1]] <- res_ISTA_bktr[[1]]
#     err_edges[1] <- sum(xor(Lambda_0 != 0, Lambda[[1]] != 0))
#     err_edges_t[1] <- sum(xor(t(Lambda_0) != 0, Lambda[[1]] != 0))
#     
#     for (i in 1:10000) {
#       Lambda[[i+1]] <- Lambda[[i]]
#       t_gma <- gma[[i]]
#       while (all((Lambda[[i]]!=0) == (Lambda[[i+1]]!=0))) {
#         t_gma <- t_gma + 0.01
#         res_ISTA_bktr <- ft_ISTA_bktr(Lambda[[i]],Sigma,omega,t_gma,eta,eta_inc)
#         # if(!res_ISTA_bktr[[3]]){
#         #  return(list(res_pos / num_exp,FALSE))
#         # }
#         Lambda[[i+1]] <- res_ISTA_bktr[[1]]
#       }
#       gma[[i+1]] <- t_gma
#       err_edges[i+1] <- sum(xor(Lambda_0 != 0, Lambda[[i+1]] != 0))
#       err_edges_t[i+1] <- sum(xor(t(Lambda_0) != 0, Lambda[[i+1]] != 0))
#       if(nnzero(Lambda[[i+1]]) == n) {
#         break
#       }
#     }
#     
#     if (0 %in% err_edges) {
#       res_pos <- res_pos+1
#     }else if (0 %in% err_edges_t) {
#       res_pos_t <- res_pos_t+1
#     }
#   }
#   
#   return(list(res_pos / num_exp,res_pos_t / num_exp,TRUE))
# }

# get_roc_ISTA <- function (n,Sigma, omega, Lambda_0, num_exp, eta, eta_inc) {
#   TPR <- numeric()
#   FPR <- numeric()
#   TPR_t <- numeric()
#   FPR_t <- numeric()
#   k <- 1
#   
#   off_diag <- row(Lambda_0) != col(Lambda_0)
#   num_pos <- sum((Lambda_0 != 0)[off_diag])
#   num_neg <- sum((Lambda_0 == 0)[off_diag])
#   
#   for (j in 1:num_exp) {
#     Lambda <- list(matrix(runif(n**2,-1,1),n,n))
#     gma <- list(0)
#     res_ISTA_bktr <- ft_ISTA_bktr(Lambda[[1]],Sigma,omega,gma[[1]],eta,eta_inc)
#     Lambda[[1]] <- res_ISTA_bktr[[1]]
#     
#     TP <- sum((Lambda_0 != 0 & Lambda[[1]] != 0)[off_diag])
#     FP <- sum((Lambda_0 == 0 & Lambda[[1]] != 0)[off_diag])
#     TPR[k] <- TP / num_pos
#     FPR[k] <- FP / num_neg
#     
#     TP_t <- sum((Lambda_0 != 0 & t(Lambda[[1]]) != 0)[off_diag])
#     FP_t <- sum((Lambda_0 == 0 & t(Lambda[[1]]) != 0)[off_diag])
#     TPR_t[k] <- TP_t / num_pos
#     FPR_t[k] <- FP_t / num_neg
#     
#     k <- k + 1
#     
#     if(nnzero(Lambda[[1]]) == n) break
#     
#     for (i in 1:10000) {
#       Lambda[[i+1]] <- Lambda[[i]]
#       t_gma <- gma[[i]]
#       while (all((Lambda[[i]]!=0) == (Lambda[[i+1]]!=0))) {
#         t_gma <- t_gma + 0.01
#         res_ISTA_bktr <- ft_ISTA_bktr(Lambda[[i]],Sigma,omega,t_gma,eta,eta_inc)
#         Lambda[[i+1]] <- res_ISTA_bktr[[1]]
#       }
#       gma[[i+1]] <- t_gma
#       
#       TP <- sum((Lambda_0 != 0 & Lambda[[i+1]] != 0)[off_diag])
#       FP <- sum((Lambda_0 == 0 & Lambda[[i+1]] != 0)[off_diag])
#       TPR[k] <- TP / num_pos
#       FPR[k] <- FP / num_neg
#       
#       TP_t <- sum((Lambda_0 != 0 & t(Lambda[[i+1]]) != 0)[off_diag])
#       FP_t <- sum((Lambda_0 == 0 & t(Lambda[[i+1]]) != 0)[off_diag])
#       TPR_t[k] <- TP_t / num_pos
#       FPR_t[k] <- FP_t / num_neg
#       
#       k <- k+1
#       
#       if(nnzero(Lambda[[i+1]]) == n) break
#     }
#   }
#   
#   return(list(TPR,FPR,TPR_t,FPR_t))
# }

get_roc_ISTA <- function (n,Sigma, omega, Lambda_0, Lambda_init, gma_cand, eta, eta_inc) {
  TPR <- numeric()
  FPR <- numeric()
  TPR_t <- numeric()
  FPR_t <- numeric()
  k <- 1
  
  off_diag <- row(Lambda_0) != col(Lambda_0)
  num_pos <- sum((Lambda_0 != 0)[off_diag])
  num_neg <- sum((Lambda_0 == 0)[off_diag])
  
  Lambda <- list(Lambda_init)
  p <- 1
  gma <- c(gma_cand[1])
  res_ISTA_bktr <- ft_ISTA_bktr(Lambda[[1]],Sigma,omega,gma[[1]],eta,eta_inc)
  Lambda[[1]] <- res_ISTA_bktr[[1]]
  
  TP <- sum((Lambda_0 != 0 & Lambda[[1]] != 0)[off_diag])
  FP <- sum((Lambda_0 == 0 & Lambda[[1]] != 0)[off_diag])
  TPR[k] <- TP / num_pos
  FPR[k] <- FP / num_neg
  
  TP_t <- sum((Lambda_0 != 0 & t(Lambda[[1]]) != 0)[off_diag])
  FP_t <- sum((Lambda_0 == 0 & t(Lambda[[1]]) != 0)[off_diag])
  TPR_t[k] <- TP_t / num_pos
  FPR_t[k] <- FP_t / num_neg
  
  k <- k + 1
  
  # if(nnzero(Lambda[[1]]) == n ** 2) break
  if(nnzero(Lambda[[1]]) == n) break
  
  for (i in 1:100) {
    Lambda[[i+1]] <- Lambda[[i]]
    while (all((Lambda[[i]]!=0) == (Lambda[[i+1]]!=0))) {
      p <- p + 1
      t_gma <- gma_cand[p]
      res_ISTA_bktr <- ft_ISTA_bktr(Lambda[[i]],Sigma,omega,t_gma,eta,eta_inc)
      Lambda[[i+1]] <- res_ISTA_bktr[[1]]
    }
    gma[[i+1]] <- t_gma
    
    TP <- sum((Lambda_0 != 0 & Lambda[[i+1]] != 0)[off_diag])
    FP <- sum((Lambda_0 == 0 & Lambda[[i+1]] != 0)[off_diag])
    TPR[k] <- TP / num_pos
    FPR[k] <- FP / num_neg
    
    TP_t <- sum((Lambda_0 != 0 & t(Lambda[[i+1]]) != 0)[off_diag])
    FP_t <- sum((Lambda_0 == 0 & t(Lambda[[i+1]]) != 0)[off_diag])
    TPR_t[k] <- TP_t / num_pos
    FPR_t[k] <- FP_t / num_neg
    
    k <- k+1
    
    if(nnzero(Lambda[[i+1]]) == n) break
  }
  
  return(list(TPR,FPR,TPR_t,FPR_t))
}
