

# --------------------- signal-noise ratio  ------------------
# gen.dat1, gen.dat2, gen.dat3, 
# from low signal-noise ratio to high signal-noise ratio;

gen.dat1 <- function(n){
  
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  A <-  rbinom(n, size = 1, prob = plogis(0.5*X1 - 0.5*X2))
  
  U <-  rnorm(n)  # runif(n, -0.5, 0.5)
  Y0 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + 3*U))
  Y1 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + 1 + 1.5*U))
  
  Y <- A*Y1 + (1-A)*Y0
  dat.full <- data.frame(X1 = X1, X2 = X2, A = A, Y = Y, Y1 = Y1, Y0 = Y0)
  
  return(dat.full)
}


gen.dat2 <- function(n){
  
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  A <-  rbinom(n, size = 1, prob = plogis(0.5*X1 - 0.5*X2))
  
  U <-  rnorm(n)  # runif(n, -0.5, 0.5)
  Y0 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + 2*U))
  Y1 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + 1 + U))
  
  Y <- A*Y1 + (1-A)*Y0
  dat.full <- data.frame(X1 = X1, X2 = X2, A = A, Y = Y, Y1 = Y1, Y0 = Y0)
  
  return(dat.full)
}




gen.dat3 <- function(n){
  
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  A <-  rbinom(n, size = 1, prob = plogis(0.5*X1 - 0.5*X2))
  
  U <-  rnorm(n)  # runif(n, -0.5, 0.5)
  Y0 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + U))
  Y1 <- rbinom(n, size = 1, prob = plogis(0.5*X1 + 0.5*X2 + 1 + 0.5*U))
  
  Y <- A*Y1 + (1-A)*Y0
  dat.full <- data.frame(X1 = X1, X2 = X2, A = A, Y = Y, Y1 = Y1, Y0 = Y0)
  
  return(dat.full)
}





# -------------------- dimension of covariates ------------------
# gen.dat4, gen.dat5, gen.dat6, 

gen.dat4 <- function(n, p = 20){
  
  X <-  matrix(rnorm(n*p), ncol = p)
  A <-  rbinom(n, size = 1, prob = plogis(0.5*X[,1] - 0.5*X[,2]))
  
  #alpha <- sapply(1:p, function(x) 1/(2*x*x) )  # 1/(x*x) is not well, but it is better than 1/(2*x*x)  
  alpha <- sapply(1:p, function(x) 1/(2^(x-1)) )
  
  
  U <-  rnorm(n)  # runif(n, -0.5, 0.5)
  Y0 <- rbinom(n, size = 1, prob = plogis( X %*% alpha + U))
  Y1 <- rbinom(n, size = 1, prob = plogis( X %*% alpha + 1 + U))
  
  Y <- A*Y1 + (1-A)*Y0
  dat.full <- data.frame(X = X, A = A, Y = Y, Y1 = Y1, Y0 = Y0)
  
  
  return(dat.full)
}


gen.dat5 <- function(n, p = 50){
  
  X <-  matrix(rnorm(n*p), ncol = p)
  A <-  rbinom(n, size = 1, prob = plogis(0.5*X[,1] - 0.5*X[,2]))
  
  alpha <- sapply(1:p, function(x) 1/(2^(x-1)) )
  
  U <-  rnorm(n)  # runif(n, -0.5, 0.5)
  Y0 <- rbinom(n, size = 1, prob = plogis( X %*% alpha + U))
  Y1 <- rbinom(n, size = 1, prob = plogis( X %*% alpha + 1 + U))
  
  Y <- A*Y1 + (1-A)*Y0
  dat.full <- data.frame(X = X, A = A, Y = Y, Y1 = Y1, Y0 = Y0)
  
  
  return(dat.full)
}


gen.dat6 <- function(n, p = 100){
  
  X <-  matrix(rnorm(n*p), ncol = p)
  A <-  rbinom(n, size = 1, prob = plogis(0.5*X[,1] - 0.5*X[,2]))
  
  # alpha <- sapply(1:p, function(x) 1/(2^x) ) 
  alpha <- sapply(1:p, function(x) 1/(2^(x-1)) )
  
  U <-  rnorm(n)  # runif(n, -0.5, 0.5)
  Y0 <- rbinom(n, size = 1, prob = plogis( X %*% alpha + U))
  Y1 <- rbinom(n, size = 1, prob = plogis( X %*% alpha + 1 + U))
  
  Y <- A*Y1 + (1-A)*Y0
  dat.full <- data.frame(X = X, A = A, Y = Y, Y1 = Y1, Y0 = Y0)
  
  
  return(dat.full)
}
