
EstFun <- function(dat, rho, K){
  
  n <- nrow(dat)
  A <- dat$A
  Y <- dat$Y
  
  # step 1. estimate the nuisance parameters using K-folds cross-fitting
  ps <- rep(NA, n)
  mu1 <- rep(NA, n)
  mu0 <- rep(NA, n)
  m <- round(n / K)

  for(k in 1:K){
    
    sub.test <- (1:n)[(m*(k-1)+1):(m*k)]  # index of test data 
    sub.train <- setdiff(1:n, sub.test)
    dat.train <- dat[sub.train, ]
    dat.test <- dat[sub.test, ]
    
    # propensity score 
    mod.ps <- glm(A ~ . -Y, family = binomial(link = 'logit'), 
                  data = dat.train)
    ps[sub.test] <- predict.glm(mod.ps, newdata = dat.test, type = 'response')
   
    # outcome regression functions
    mod.mu1 <- glm(Y ~. -A, family = binomial(link = 'logit'),
                   data = dat.train, subset = A == 1)
    mod.mu0 <- glm(Y ~. -A, family = binomial(link = 'logit'),
                   data = dat.train, subset = A == 0)
    mu1[sub.test] <- predict.glm(mod.mu1, newdata = dat.test, type = 'response')
    mu0[sub.test] <- predict.glm(mod.mu0, newdata = dat.test, type = 'response')
    
  }
  

  # step 2. calculate the influence function
  
  phi.0 <- mu0*(1-mu1) + (1-A)*(Y-mu0)*(1-mu1)/(1-ps) - A*(Y-mu1)*mu0/ps
  phi.r <- (1-2*mu1)*sqrt(mu0*(1-mu0)/(mu1*(1-mu1)))*A*(Y-mu1)/(2*ps) +
    (1-2*mu0)*sqrt(mu1*(1-mu1)/(mu0*(1-mu0)))*(1-A)*(Y-mu0)/(2*(1-ps)) +
    sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  
  g.eta <- mu0*(1-mu1) - rho * sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  psi <- 1*(g.eta>=0) * (phi.0 - rho*phi.r)   
  
  # step 3. obtain the point estimate and confidence interval
  beta.hat <- mean(psi)
  std <- sd(psi)/sqrt(n)
  
  return(list(beta.hat = beta.hat, std = std))
}





# add variable selection for propensity score and outcome regression functions. 
EstFun2 <- function(dat, rho, K, p){
  
  n <- nrow(dat)
  A <- dat$A
  Y <- dat$Y
  X <- as.matrix(dat[, 1:p])
  
  # step 1. estimate the nuisance parameters using K-folds cross-fitting
  ps <- rep(NA, n)
  mu1 <- rep(NA, n)
  mu0 <- rep(NA, n)
  m <- round(n / K)
  
  for(k in 1:K){
    
    sub.test <- (1:n)[(m*(k-1)+1):(m*k)]  # index of test data 
    sub.train <- setdiff(1:n, sub.test)
    A.train <- A[sub.train]
    Y.train <- Y[sub.train]
    X.train <- X[sub.train,]
    X.test <- X[sub.test,]
    
    # propensity score 
    mod.ps = cv.glmnet(X.train, A.train, alpha = 1,  
                       family="binomial", nfolds = 5,  type.measure="class")
    ps.lasso <- predict(mod.ps, newx = X.test, type = "response", s = 'lambda.min')
    ps[sub.test] <- drop(ps.lasso)
    

    # outcome regression functions
    mod.mu1 = cv.glmnet(X.train, Y.train, alpha = 1,  
                       family="binomial", nfolds = 5,  type.measure="class",
                       subset = A.train == 1)
    mod.mu0 = cv.glmnet(X.train, Y.train, alpha = 1,  
                        family="binomial", nfolds = 5,  type.measure="class",
                        subset = A.train == 0)
    mu1[sub.test] <- predict(mod.mu1, newx = X.test, type = "response", s = 'lambda.min')
    mu0[sub.test] <-  predict(mod.mu0, newx = X.test, type = "response", s = 'lambda.min')
    

  }
  ps <- ifelse(ps <= 0.02, 0.02, ps)
  ps <- ifelse(ps >= 0.98, 0.98, ps)
  mu1 <- ifelse(mu1 <= 0.02, 0.02, mu1)
  mu1 <- ifelse(mu1 >= 0.98, 0.98, mu1)
  mu0 <- ifelse(mu0 <= 0.02, 0.02, mu0)
  mu0 <- ifelse(mu0 >= 0.98, 0.98, mu0)
  
  phi.0 <- mu0*(1-mu1) + (1-A)*(Y-mu0)*(1-mu1)/(1-ps) - A*(Y-mu1)*mu0/ps
  phi.r <- (1-2*mu1)*sqrt(mu0*(1-mu0)/(mu1*(1-mu1)))*A*(Y-mu1)/(2*ps) +
    (1-2*mu0)*sqrt(mu1*(1-mu1)/(mu0*(1-mu0)))*(1-A)*(Y-mu0)/(2*(1-ps)) +
    sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  
  g.eta <- mu0*(1-mu1) - rho * sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  psi <- 1*(g.eta>=0) * (phi.0 - rho*phi.r)   
  
  # step 3. obtain the point estimate and confidence interval
  beta.hat <- mean(psi)
  std <- sd(psi)/sqrt(n)
  
  return(list(beta.hat = beta.hat, std = std))
}

