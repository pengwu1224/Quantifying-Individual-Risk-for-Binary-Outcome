# Encoding: UTF-8

rm(list = ls( ))

# generate data
source('GenData.R')   

# function of estimate the bound of FNA
source('EstFun.R')

# rho
rho.lst <- seq(0, 0.4, 0.01)

set.seed(844)

# true FNA
m <- 1000000
dat.full <- gen.dat3(m)
true.FNA <- sum((dat.full$Y0 == 1) & (dat.full$Y1 == 0))/m

# true mu0 and mu1, based on Monte Carlo 
X1 <- dat.full$X1
X2 <- dat.full$X2
true.mu0 <- rep(NA, m)
true.mu1 <- rep(NA, m)


for(j in 1:m){
  U <- rnorm(1000)
  temp0 <- plogis(0.5*X1[j] + 0.5*X2[j] + U) 
  temp1 <- plogis(0.5*X1[j] + 0.5*X2[j] + 1 + 0.5*U)
  
  true.mu0[j] <- mean(temp0)
  true.mu1[j] <- mean(temp1) 
  cat(j,'\r')
}



RES <- data.frame(matrix(nrow = length(rho.lst), ncol = 5))
names(RES) <- c('rho', 'true.FNA', 'true.beta', 'est.beta','ese.2000')

RES$rho <- rho.lst
RES$true.FNA <- true.FNA


n <- 2000 
dat.full <- gen.dat3(n) 
dat <- dat.full[,c('X1', 'X2', 'A', 'Y')]

for(i in 1:length(rho.lst)){
  
  rho <- rho.lst[i]
  
  # true bound for a given rho
  g.eta <- true.mu0*(1-true.mu1) - 
    rho * sqrt(true.mu0*(1-true.mu0)*true.mu1*(1-true.mu1))
  psi <- 1*(g.eta>=0) * g.eta    
  true.beta <- mean(psi)
  
  RES$true.beta[i] <- true.beta
  
  cat('========== rho =', rho, ' =========== \n')
  
  ## estimate the FNA for a given $rho$
  res <- EstFun(dat, rho = rho, K = 2)
  
  RES[i, 4] <- res$beta.hat 
  RES[i, 5] <- res$std
  
  ## generate data 
  
  print(RES)  
  
}  

RES$FH_low <- mean(pmax(true.mu0 - true.mu1, 0)) 
RES$FH_up <- mean(pmin(true.mu0, 1 - true.mu1)) 


write.csv(RES, file = 'case3.csv', row.names = F)


