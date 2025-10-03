# experiment_run.R
library('fBasics')
library('Matrix')
library("matrixcalc")

library(future.apply)
library(progressr)

source("functions.R")
load("result/Accuracy/res_3.RData")

# enable progress bar globally
handlers(global = TRUE)
handlers("txtprogressbar")   # simple text bar works in background jobs

# -------- Parallel setup --------
plan(multisession, workers = parallel::detectCores() - 1)

n <- 3
# Lambda_supp <- rbind(c(1,0,1),
#                      c(0,1,1),
#                      c(0,0,1))

# n <- 4
# Lambda_supp <- rbind(c(1,0,0,1),
#                      c(0,1,0,1),
#                      c(0,0,1,1),
#                      c(0,0,0,1))

# n <- 5
# Lambda_supp <- rbind(c(1,0,0,0,1),
#                      c(0,1,0,0,1),
#                      c(0,0,1,0,1),
#                      c(0,0,0,1,1),
#                      c(0,0,0,0,1))

# n <- 10
# Lambda_supp <- rbind(c(1,0,0,0,0,0,0,0,0,1),
#                      c(0,1,0,0,0,0,0,0,0,1),
#                      c(0,0,1,0,0,0,0,0,0,1),
#                      c(0,0,0,1,0,0,0,0,0,1),
#                      c(0,0,0,0,1,0,0,0,0,1),
#                      c(0,0,0,0,0,1,0,0,0,1),
#                      c(0,0,0,0,0,0,1,0,0,1),
#                      c(0,0,0,0,0,0,0,1,0,1),
#                      c(0,0,0,0,0,0,0,0,1,1),
#                      c(0,0,0,0,0,0,0,0,0,1))

num_exp  <- 10
num_init <- 1

pairs <- expand.grid(init_id = 1:num_init, exp_id = 1:num_exp)

with_progress({
  p <- progressor(along = 1:nrow(pairs))   
  
  results <- future_lapply(1:num_init, function(i) {
    Lambda_0 <- Lambda_star[[6]]
    omega <- 1
    
    vec_Sigma <- solve(diag(n^2) - (t(Lambda_0) %x% t(Lambda_0))) %*% vec(omega * diag(n))
    Sigma <- matrix(vec_Sigma, nrow = n)
    
    # res_IHT_prop  <- get_prop_IHT(n, Sigma, omega, Lambda_0, num_exp, 1, 1.05)
    # res_ISTA_prop <- get_prop_ISTA(n, Sigma, omega, Lambda_0, num_exp, 1, 1.05)
    
    res_IHT_roc  <- get_roc_IHT(n, Sigma, omega, Lambda_0, num_exp, 1, 1.05)
    res_ISTA_roc <- get_roc_ISTA(n, Sigma, omega, Lambda_0, num_exp, 1, 1.05)
    
    # package results
    out <- list(
      Lambda_0 = Lambda_0,
      res = list(
        data.frame(n=n, algo="IHT",    TPR=res_IHT_roc[[1]], FPR=res_IHT_roc[[2]], init_id = i),
        data.frame(n=n, algo="ISTA",   TPR=res_ISTA_roc[[1]], FPR=res_ISTA_roc[[2]], init_id = i),
        data.frame(n=n, algo="ISTA_t", TPR=res_ISTA_roc[[3]], FPR=res_ISTA_roc[[4]], init_id = i)
      )
    )
    
    # save to disk immediately (append mode)
    saveRDS(out, file = paste0("partial_result/partial_result_", i, ".rds"))
    
    # update progress
    p()
    
    out
  }, future.seed = TRUE)
})

# -------- Collect results --------
df_t <- do.call(rbind, unlist(lapply(results, `[[`, "res"), recursive = FALSE))

save(df_t, file = "res_3_roc_4_all_10.RData")
