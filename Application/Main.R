# Encoding: UTF-8

rm(list = ls( ))

# ---------------------- read the data ----------------------
dat <- read.csv(file = 'rhc.csv')[, -1]

# --------------------- clean the data ----------------------
## treatment 
A <- dat$swang1  # whether or not a patient received a RHC  
table(A)         # 3551, control; 2184, treated
dat$swang1 <- NULL

## outcome
Y <- dat$dth30   # whether or not death at 30 days, Yes denotes death
dat$dth30 <- NULL

## some invalid variables
dat$ptid <- NULL     # Patient ID    
dat$death <- NULL    # Death at any time up to 180 Days
dat$dschdte <- NULL  # Hospital Discharge Date
dat$lstctdte <- NULL # Date of Last Contact
dat$dthdte <- NULL   # Date of Death
dat$sadmdte <- NULL  # Study Admission Date
dat$t3d30 <- NULL

# deal with missing values
dat$cat2[is.na(dat$cat2)] <- 'ref'   
dat$adld3p <- NULL # ADL
dat$urin1 <- NULL   # Urine output

##  encoding discrete variable with dummy variables
dis_names <- c('race', 'income', 'ninsclas', 'cat1', 'cat2', 'ca')
for(i in dis_names){
  dat[, i] <- as.factor(dat[, i])
}
dat.new <- model.matrix(~., dat)
dat.new <- data.frame(dat.new)   # a total of 72 covariates

dat.new <- dat.new[,-1]  # remove the intercept/constant term  

dat.new <- scale(dat.new)

rm(dat)


# --------------------- analysis the data ----------------------

## Step 1. estimate the nuisance parameters using K-folds cross-fitting
K <- 2                          # two-folds cross-fitting
n <- nrow(dat.new)  
A <- ifelse(A == 'RHC', 1, 0)   # treatment 
Y <- ifelse(Y == 'Yes', 0, 1)   # whether or not survival at 30 days, 1 denotes survival
# names(dat.new)

# (a). estimate propensity score
ps <- rep(NA, n)         # propensity score

# selecting the covariates in propensity score model
vars_lst <- c()
for(i in 1:ncol(dat.new)){
  Xi <- dat.new[,i]
  mod <- glm(A ~ Xi, family = binomial(link = 'logit'))
  if(abs(summary(mod)$coefficients[2,3]) >= 15){
    vars_lst <- c(vars_lst, i)
  }
}
dat.ps <- data.frame(dat.new[, vars_lst], A = A)

# estimate propensity score with K-fold cross-fitting 
K <- 2                  
m <- round(n / K)
set.seed(455)
ind.rand <- sample(1:n, m)
sub_test_list <- list(ind.rand, setdiff(1:n, ind.rand) ) 
for(k in 1:K){
  
  sub.test <- sub_test_list[[k]]  # index of test data 
  sub.train <- setdiff(1:n, sub.test)
  dat.train <- dat.ps[sub.train, ]
  dat.test <- dat.ps[sub.test, ]
  
  # propensity score 
  mod.ps <- glm(A ~ ., family = binomial(link = 'logit'), 
                data = dat.train)
  # mod.ps <- step(mod.ps, trace = 0)
  ps[sub.test] <- predict.glm(mod.ps, newdata = dat.test, type = 'response')
  
}

# (b) estimate outcome regression models 

# selecting the covariates in outcome regression model 
vars_lst <- c()
for(i in 1:ncol(dat.new)){
  Xi <- dat.new[,i]
  mod <- lm(Y ~ A + Xi)
  if(abs(summary(mod)$coefficients[3,3]) >= 8){
    vars_lst<- c(vars_lst, i)
  }
}


dat.or <- data.frame(dat.new[, vars_lst], A = A, Y = Y)

mu1 <- rep(NA, n)        # outcome regression function
mu0 <- rep(NA, n)

# estimate outcome regression models with K-fold cross-fitting
for(k in 1:K){
  
  sub.test <- sub_test_list[[k]]  # index of test data 
  sub.train <- setdiff(1:n, sub.test)
  
  dat.train <- dat.or[sub.train, ]
  dat.test <- dat.or[sub.test, ]

  # outcome regression functions
  mod.mu1 <- glm(Y ~. -A, family = binomial(link = 'logit'),
                 data = dat.train, subset = A == 1)  
  mod.mu0 <- glm(Y ~. -A, family = binomial(link = 'logit'),
                 data = dat.train, subset = A == 0)    
  mu1[sub.test] <- predict.glm(mod.mu1, newdata = dat.test, type = 'response')
  mu0[sub.test] <- predict.glm(mod.mu0, newdata = dat.test, type = 'response')
  
}


# truncation to obtain more stable results 
mu1 <- ifelse(mu1 <= 0.02, 0.02, mu1)
mu0 <- ifelse(mu0 <= 0.02, 0.02, mu0)


## Step 2. estimate the FNA

# ATE
psi.1 <- A*Y/ps - (A-ps)*mu1 / ps
psi.0 <- (1-A)*Y/(1-ps) + (A-ps)*mu0/(1-ps)
mean(psi.1 - psi.0)                  # DR
sd(psi.1 - psi.0) / sqrt(n)

# value ranges of rho, restricted by Frechet Hoeffiding Bounds
L_rho <- pmax(-(1-mu0)*(1-mu1),  -mu0*mu1)/sqrt(mu0*(1-mu0)*mu1*(1-mu1))
U_rho <- pmin( mu0*(1-mu1),  (1-mu0)*mu1)/sqrt(mu0*(1-mu0)*mu1*(1-mu1))


# Frechet Hoeffiding Bounds
FH_low <- mean(pmax(mu0-mu1, 0))
FH_up <- mean(pmin(mu0, 1-mu1))
FH_low
FH_up


# estimate the bound on FNA for a given set of rho
rho.lst <- seq(-0.45, 0.70, 0.05)    # specify the value of rho
rho.lst <- c(rho.lst, 0.706)

RES <- data.frame(matrix(nrow = length(rho.lst), ncol = 3))
names(RES) <- c('rho', 'est.beta',  'ese')


RES$rho <- rho.lst

for(i in 1:length(rho.lst)){
  
  rho <- rho.lst[i]
  cat('========== rho =', rho,  ' =========== \n')
  
  ## estimate the FNA for a given $rho$
  phi.0 <- mu0*(1-mu1) + (1-A)*(Y-mu0)*(1-mu1)/(1-ps) - A*(Y-mu1)*mu0/ps
  phi.r <- (1-2*mu1)*sqrt(mu0*(1-mu0)/(mu1*(1-mu1)))*A*(Y-mu1)/(2*ps) +
    (1-2*mu0)*sqrt(mu1*(1-mu1)/(mu0*(1-mu0)))*(1-A)*(Y-mu0)/(2*(1-ps)) +
    sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  
  g.eta <- mu0*(1-mu1) - rho * sqrt(mu0*(1-mu0)*mu1*(1-mu1))
  psi <- 1*(g.eta>=0) * (phi.0 - rho*phi.r)   
  
  # step 3. obtain the point estimate and confidence interval
  beta.hat <- mean(psi)
  ese <- sd(psi)/sqrt(n)
  
  RES[i, 2] <- beta.hat
  RES[i, 3] <- ese
  
  print(RES, digits = 3)
}  

# write.csv(RES, file = 'result.csv', row.names = F)

color <- c('#009E73','#009E73','#984ea3','#984ea3','#984ea3')
color_ribon <- c('#dceef2','#dceef2','#ebe6fa','#ebe6fa','#ebe6fa')

ggplot(data = RES) + 
  geom_ribbon(aes(x = rho, 
                  ymax = est.beta + 1.96*ese,
                  ymin = est.beta - 1.96*ese), 
              fill = color_ribon[1], color = NA, alpha = 0.8) +
  geom_hline(aes(yintercept = FH_low), linetype = 'dotted', linewidth=1, col = 'blue') + 
  geom_hline(aes(yintercept = FH_up), linetype = 'dotted',linewidth=1,  col = 'blue') + 
  geom_line(aes(x = rho, y = est.beta), linewidth = 0.6, col = 'red') + 
  theme_bw() +  
  geom_point(aes(x = 0.00, y = 0.2073), color = "red", size = 2.5) +
  geom_text(aes(x = 0.00, y = 0.2073, label = '(0.00, 0.207)'), size =3.5,
            vjust = -1.5, hjust=-0.05, color = "purple") + 
  geom_point(aes(x = -0.45, y = 0.2896), color = "red", size = 2.5) +
  geom_text(aes(x = -0.45, y = 0.2896, label = '(-0.45, 0.290)'), size =3.5, 
            vjust = -1.5, color = "purple") + 
  labs(x = expression(paste('rho (', rho,')')), y = 'Bounds on FNA') + 
  geom_point(aes(x = 0.706, y = 0.0742), color = "red", size = 2.5) +
  geom_text(aes(x = 0.706, y = 0.0742, label = '(0.706, 0.074)'), size =3.5,
            vjust = -1.5, hjust=0.8, color = "purple") + 
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  coord_cartesian(ylim = c(0.05, 0.32), xlim = c(-0.55, 0.75)) +
  theme(axis.title.x  = element_text(size = 14, family = "Times New Roman"),
        axis.title.y  = element_text(size = 14, family = "Times New Roman"),
        strip.text.x = element_text(size = 16))
