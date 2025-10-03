# experiment_ADMM_ISTA.R
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

Lambda_0 <- Lambda_star[[12]]

omega <- 1

vec_Sigma <- solve(diag(n**2)-(t(Lambda_0) %x% t(Lambda_0))) %*% vec(omega * diag(n))
Sigma <- matrix(vec_Sigma, nrow = n)

# Lambda_init <- lapply(1:10, function(x) matrix(runif(n**2,-1,1), n,n))


# --- Parallelized Main Loop ---
# Define the range to iterate over. The commented-out line shows the original intent.
# We will use future_lapply to run the loop in parallel.
range_k <- 1:(length(Lambda_init) - 1)

# The loop is converted into a function that is applied to each element of the range.
# The 'future_lapply' function will execute this function in parallel on the available workers.
df_list <- future_lapply(range_k, function(k) {
  # The code inside this function is the body of the original for-loop.
  # This code will be executed for each value of 'k' on a separate worker.
  
  # Initialize the lists for this iteration
  Lambda_1 <- list(Lambda_init[[k]])
  Lambda_2 <- list(Lambda_init[[k+1]])
  
  # Define gma and run the first ADMM-ISTA step
  gma <- log(seq(100, 1, -0.1), 2)
  res_ADMM_ISTA <- ADMM_ISTA_constant(Lambda_1[[1]], Lambda_2[[1]], Sigma, omega, gma[[1]])
  Lambda_1[[1]] <- res_ADMM_ISTA[[1]]
  Lambda_2[[1]] <- res_ADMM_ISTA[[2]]
  obj <- list(res_ADMM_ISTA[[3]])
  
  # Run the inner loop
  for (i in 1:length(gma)) {
    res_ADMM_ISTA <- ADMM_ISTA_constant(Lambda_1[[i]], Lambda_2[[i]], Sigma, omega, gma[[i]])
    Lambda_1[[i+1]] <- res_ADMM_ISTA[[1]]
    Lambda_2[[i+1]] <- res_ADMM_ISTA[[2]]
    obj[[i+1]] <- res_ADMM_ISTA[[3]]
    
    if (nnzero(Lambda_2[[i+1]]) == n) break
  }
  
  # Create the temporary data frame for this iteration
  df_t <- data.frame(
    lmd = array(Lambda_2),
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

save(df, file = "results_ADMM_ISTA.RData")