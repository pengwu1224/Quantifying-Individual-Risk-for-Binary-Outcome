# Encoding: UTF-8

rm(list = ls( ))

library(glmnet)

# generate data
source('GenData.R')   

# function of estimate the bound of FNA
source('EstFun.R')

# sample size 
size.lst <- c(500, 1000, 2000)     

# rho
rho.lst <- seq(0, 0.3, 0.1)


set.seed(67)

# true FNA
m <- 1000000
dat.full <- gen.dat5(m)
true.FNA <- sum((dat.full$Y0 == 1) & (dat.full$Y1 == 0))/m
# cor(dat.full$Y0, dat.full$Y1)

# true mu0 and mu1, based on Monte Carlo 
p <- 50
X <- dat.full[, 1:p]
alpha <- sapply(1:p, function(x) 1/(2^(x-1)) )
tmp <- as.matrix(X) %*% alpha

true.mu0 <- rep(NA, m)
true.mu1 <- rep(NA, m)

for(j in 1:m){
  U <- rnorm(1000)
  temp0 <- plogis( tmp[j] + U) 
  temp1 <- plogis( tmp[j] + 1 + U)
  
  true.mu0[j] <- mean(temp0)
  true.mu1[j] <- mean(temp1) 
  cat(j,'\r')
}



RES <- data.frame(matrix(nrow = length(rho.lst), ncol = 14))
names(RES) <- c('rho', 'true.beta', 'bias.500', 'std.500', 'ese.500', 'CP95.500', 
                'bias.1000', 'std.1000', 'ese.1000', 'CP95.1000',
                'bias.2000', 'std.2000', 'ese.2000', 'CP95.2000')

RES$rho <- rho.lst


for(i in 1:length(rho.lst)){
  
  rho <- rho.lst[i]
  
  # true bound for a given rho
  g.eta <- true.mu0*(1-true.mu1) - 
    rho * sqrt(true.mu0*(1-true.mu0)*true.mu1*(1-true.mu1))
  psi <- 1*(g.eta>=0) * g.eta    
  true.beta <- mean(psi)
  
  RES$true.beta[i] <- true.beta
  
  
  for(j in 1:length(size.lst)){
    
    n <- size.lst[j] 
    cat('========== rho =', rho, '& n =', n, ' =========== \n')
    
    # replicate 1000 simulations for each "rho" and "n"
    B <- 500
    est.lst <- rep(NA, B)
    std.lst <- rep(NA, B)
    count.lst <- rep(NA, B)
    
    for(b in 1:B){
      cat(b,'\r')
      
      ## generate data 
      dat.full <- gen.dat5(n) 
      dat <- dat.full[, 1:(p+2)]
      
      ## estimate the FNA for a given $rho$
      res <- EstFun2(dat, rho = rho, K = 2, p = 50)
      
      est.lst[b] <- res$beta.hat 
      std.lst[b] <- res$std
      count.lst[b] <-  (res$beta.hat <= true.beta + 1.96*res$std) &  
        (res$beta.hat >= true.beta - 1.96*res$std)
      
    }
    RES[i, 3+4*(j-1)] <- mean(est.lst - true.beta, na.rm = T)
    RES[i, 4+4*(j-1)] <- sd(est.lst, na.rm = T)
    RES[i, 5+4*(j-1)] <- mean(std.lst, na.rm = T)
    RES[i, 6+4*(j-1)] <- mean(count.lst, na.rm = T)
    
    print(RES)
  }
}  

write.csv(RES, file = 'case5.csv', row.names = F)


