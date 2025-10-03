# experiment_IHT.R
library('fBasics')
library('Matrix')
library("matrixcalc")

library(dplyr)
library(tidyr)

library(reshape2)


library(future)
library(future.apply)
library(parallel)

# --- Parallel Processing Setup ---
plan(multisession, workers = parallel::detectCores() - 1)

source("functions.R")
load("result/Accuracy/res_3.RData")
load("Initializations/inits_lmd_12.RData")
n <- 3

# Lambda_supp <- rbind(c(1,0,1),
#                      c(0,1,1),
#                      c(1,1,1))
# Lambda_0 <- Lambda_supp
# Lambda_0[Lambda_supp == 1] <- runif(sum(Lambda_supp))

Lambda_0 <- Lambda_star[[17]]

omega <- 1

vec_Sigma <- solve(diag(n**2)-(t(Lambda_0) %x% t(Lambda_0))) %*% vec(omega * diag(n))
Sigma <- matrix(vec_Sigma, nrow = n)

# Lambda_init <- lapply(1:10, function(x) matrix(runif(n**2,-1,1), n,n))


# --- Parallelized Main Loop ---
# Define the range to iterate over. The commented-out line shows the original intent.
# We will use future_lapply to run the loop in parallel.
range_k <- 1:(length(Lambda_init))

# The loop is converted into a function that is applied to each element of the range.
# The 'future_lapply' function will execute this function in parallel on the available workers.
df_list <- future_lapply(range_k, function(k) {
  # The code inside this function is the body of the original for-loop.
  # This code will be executed for each value of 'k' on a separate worker.
  
  # Initialize the lists for this iteration
  Lambda <- list(Lambda_init[[k]])
  
  # Define gma and run the first IHT step
  eta <- 1
  eta_inc <- 1.5
  # gma <- log(seq(100,1,-0.1),10)
  # gma <- log(seq(1,100,0.1),10)
  gma <- log(seq(1,100,0.1),10)
  res_IHT_bktr <- ft_IHT_bktr(Lambda[[1]],n,Sigma,omega,gma[[1]],eta,eta_inc)
  Lambda[[1]] <- res_IHT_bktr[[1]]
  obj <- list(res_IHT_bktr[[2]])
  
  for (i in 1:(length(gma) - 1)) {
    res_IHT_bktr <- ft_IHT_bktr(Lambda[[i]],n,Sigma,omega,gma[i+1],eta,eta_inc)
    Lambda[[i+1]] <- res_IHT_bktr[[1]]
    obj[[i+1]] <- res_IHT_bktr[[2]]
    
    if(nnzero(Lambda[[i+1]]) == n) break
  }
  
  # Create the temporary data frame for this iteration
  df_t <- data.frame(
    lmd = array(Lambda),
    gma = gma[1:length(obj)],
    obj = array(unlist(obj)),
    init_id = rep(as.character(k), length(obj))
  )
  
  # Return the data frame, which will be collected by future_lapply
  return(df_t)
})

# --- Combine Results ---
# The result of future_lapply is a list of data frames.
# We use do.call(rbind, ...) to combine them into a single data frame.
df <- do.call(rbind, df_list)

# You can now proceed with further analysis on the combined 'df' data frame.
# For example, print the first few rows to verify the result.
head(df)

# You can stop the parallel workers when you are done.
# plan(sequential)

save(df, file = "results_IHT_Mat17.RData")