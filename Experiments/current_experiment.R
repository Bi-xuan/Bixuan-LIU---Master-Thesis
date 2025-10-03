library('fBasics')
library('Matrix')
library("matrixcalc")

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

Lambda <- list(Lambda_init[[8]])

eta <- 1
eta_inc <- 1.05
# gma_cand <- c(0,log(seq(100,1,-0.1),10))
# gma_cand <- c(0,10 ** (seq(log10(0.36),log10(0.37),0.00001)))
gma_cand <- log(seq(1,100,0.001),10)
gma <- c(gma_cand[1])
# gma_cand <- 10 ** (seq(-5,log10(2.5),0.003))
# gma_cand <- seq(0,2,0.001)
res_IHT_bktr <- ft_IHT_bktr(Lambda[[1]],n,Sigma,omega,gma[1],eta,eta_inc)
# res_IHT_bktr <- ft_IHT_bktr(Lambda_init[[8]],n,Sigma,omega,gma[1],eta,eta_inc)
Lambda[[1]] <- res_IHT_bktr[[1]]
obj <- list(res_IHT_bktr[[2]])

k <- 1

for (i in 1:100) {
  Lambda[[i+1]] <- Lambda[[i]]
  while (all((Lambda[[i]]!=0) == (Lambda[[i+1]]!=0))) {
    k <- k + 1
    t_gma <- gma_cand[k]
    res_IHT_bktr <- ft_IHT_bktr(Lambda[[i]],n,Sigma,omega,t_gma,eta,eta_inc)
    # res_IHT_bktr <- ft_IHT_bktr(Lambda_init[[8]],n,Sigma,omega,t_gma,eta,eta_inc)
    Lambda[[i+1]] <- res_IHT_bktr[[1]]
    
    if(k == length(gma_cand)) break
  }
  obj[[i+1]] <- res_IHT_bktr[[2]]
  gma[i+1] <- t_gma
  if(nnzero(Lambda[[i+1]]) == n + 1 || k == length(gma_cand)) break
}

Lambda

gma