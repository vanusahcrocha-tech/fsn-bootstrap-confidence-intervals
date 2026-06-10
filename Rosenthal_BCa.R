
#BCa condidence intervals

# Install packages if necessary
# install.packages("boot")
# install.packages("sn")

library(boot)
library(sn)

##### Standard Normal distribution

# Number of Monte Carlo replications
R <- 10000

# Number of bootstrap resamples
B <- 1000

# One-sided 5% critical value of the standard normal distribution
za <- qnorm(0.95)

# Sample sizes considered
n <- c(5, 10, 15, 30, 50)

# Theoretical mean and variance of the normal
m <- 0
v <- 1

# Storage matrix for simulated FSN
rosent1 <- matrix(NA, nrow = R, ncol = length(n))  # 

# Storage array for bootstrap confidence interval limits
mat1 <- array(NA, dim = c(R, 2, length(n)))        

# Coverage probabilities 
coverage1 <- numeric(length(n))

# Expected FSN
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

fsn_normal <- (n^2 * m^2 + n * v) / za^2 - n + eps


# Function to compute the FSN estimator for bootstrap samples
fsn_estimator <- function(data, indices) {
  z_boot <- data[indices]
  S_boot <- sum(z_boot)
  fsn_boot <- (S_boot / za)^2 - length(z_boot)
  return(fsn_boot)
}

## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate standard normal observations
    z <- rnorm(n[k])
    S <- sum(z)
    
    # Compute Rosenthal's FSN
    if (S >= za * sqrt(n[k])) {
      
      # Significant meta-analysis
      fsn <- (S / za)^2 - n[k]  
      sig_count <- sig_count + 1
      
    } else {
      
      # Non-significant meta-analysis
      fsn <- 0
    }
    
    rosent1[i, k] <- fsn
    
    # Bootstrap distribution of the FSN estimator
    tryCatch({
      boot_obj <- boot(data = z, statistic = fsn_estimator, R = B)
      
      # BCa confidence interval
      bcaci <- boot.ci(boot_obj, type = "bca", conf = 0.95)
      
      # Extract the lower and upper BCa limits
      if (!is.null(bcaci$bca)) {
        mat1[i, , k] <- bcaci$bca[4:5]
      }
    }, error = function(e) {
      
      # Return NA if BCa computation fails
      mat1[i, , k] <- c(NA, NA)
    })
  }
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat1[, 1, k]) & (rosent1[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability for the theoretical FSN
    outside_fsn <- sum(mat1[valid, 1, k] > fsn_normal[k]) +
      sum(mat1[valid, 2, k] < fsn_normal[k])
    
    coverage1[k] <- 1 - outside_fsn / sum(valid)
  }
}


##### Half Normal distribution

# Number of Monte Carlo replications
R <- 10000

# Number of bootstrap resamples
B <- 1000

# One-sided 5% critical value of the standard normal distribution
za <- qnorm(0.95)

# Sample sizes considered
n <- c(5, 10, 15, 30, 50)

# Theoretical mean and variance of the half-normal
m <- sqrt(2/pi)
v <- 1-2/pi 

# Storage matrix for simulated FSN
rosent2 <- matrix(NA, nrow = R, ncol = length(n))  # 

# Storage array for bootstrap confidence interval limits
mat2 <- array(NA, dim = c(R, 2, length(n)))        

# Coverage probabilities 
coverage2 <- numeric(length(n))

# Expected FSN
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

fsn_hnormal <- (n^2 * m^2 + n * v) / za^2 - n + eps

# Function to compute the FSN estimator for bootstrap samples
fsn_estimator <- function(data, indices) {
  z_boot <- data[indices]
  S_boot <- sum(z_boot)
  fsn_boot <- (S_boot / za)^2 - length(z_boot)
  return(fsn_boot)  
}

## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
   # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate standard normal observations
    z <- abs(rnorm(n[k]))
    S <- sum(z)
    
    # Compute Rosenthal's FSN
    if (S >= za * sqrt(n[k])) {
      
      # Significant meta-analysis
      fsn <- (S / za)^2 - n[k]  
 
    } else {
      
      # Non-significant meta-analysis
      fsn <- 0
    }
    
    rosent2[i, k] <- fsn
    
    # Bootstrap distribution of the FSN estimator
    tryCatch({
      boot_obj <- boot(data = z, statistic = fsn_estimator, R = B)
      
      # BCa confidence interval
      bcaci <- boot.ci(boot_obj, type = "bca", conf = 0.95)
      
      # Extract the lower and upper BCa limits
      if (!is.null(bcaci$bca)) {
        mat2[i, , k] <- bcaci$bca[4:5]
      }
    }, error = function(e) {
      
      # Return NA if BCa computation fails
      mat2[i, , k] <- c(NA, NA)
    })
  }
  
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat2[, 1, k]) & (rosent2[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability for the theoretical FSN
    outside_fsn <- sum(mat2[valid, 1, k] > fsn_hnormal[k]) +
      sum(mat2[valid, 2, k] < fsn_hnormal[k])
    
    coverage2[k] <- 1 - outside_fsn / sum(valid)
  }
}

#### Skew-normal distribution (negative skewness)
# This distribution was excluded from the BCa confidence interval analysis,
# as the proportion of positive FSN values was zero or negligibly small.



#### Skew Normal distribution (positive skewness)

# Number of Monte Carlo replications
R <- 10000

# Number of bootstrap resamples
B <- 1000

# One-sided 5% critical value of the standard normal distribution
za <- qnorm(0.95)

# Sample sizes considered
n <- c(5, 10, 15, 30, 50)

# Theoretical mean and variance of the skew-normal(positive skewness)
# Set delta to one of the desired skewness levels:0.2, 0.5, 0.8
delta <- 0.8
m <- delta*sqrt(2/pi)
v <- 1-2*(delta^2)/pi

# Storage matrix for simulated FSN
rosent4 <- matrix(NA, nrow = R, ncol = length(n))  # 

# Storage array for bootstrap confidence interval limits
mat4 <- array(NA, dim = c(R, 2, length(n)))        

# Coverage probabilities 
coverage4 <- numeric(length(n))

# Expected FSN
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

fsn_sn <- (n^2 * m^2 + n * v) / za^2 - n + eps


# Function to compute the FSN estimator for bootstrap samples
fsn_estimator <- function(data, indices) {
  z_boot <- data[indices]
  S_boot <- sum(z_boot)
  fsn_boot <- (S_boot / za)^2 - length(z_boot)
  return(fsn_boot)  
}

## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate skew-normal observations
    # alpha values correspond to skewness levels delta = 0.2, 0.5, and 0.8
    #z <- rsn(n[k],xi=0,omega=1,alpha= 0.2041)
    #z <- rsn(n[k],xi=0,omega=1,alpha= 0.5774)
    z <- rsn(n[k],xi=0,omega=1,alpha= 1.3333)
    
    S <- sum(z)
    
    # Compute Rosenthal's FSN
    if (S >= za * sqrt(n[k])) {
      
      # Significant meta-analysis
      fsn <- (S / za)^2 - n[k]  
      sig_count <- sig_count + 1
      
    } else {
      
      # Non-significant meta-analysis
      fsn <- 0
    }
    
    rosent4[i, k] <- fsn
    
    # Bootstrap distribution of the FSN estimator
    tryCatch({
      boot_obj <- boot(data = z, statistic = fsn_estimator, R = B)
      
      # BCa confidence interval
      bcaci <- boot.ci(boot_obj, type = "bca", conf = 0.95)
      
      # Extract the lower and upper BCa limits
      if (!is.null(bcaci$bca)) {
        mat4[i, , k] <- bcaci$bca[4:5]
      }
    }, error = function(e) {
      
      # Return NA if BCa computation fails
      mat4[i, , k] <- c(NA, NA)
    })
  }
  
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat4[, 1, k]) & (rosent4[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability for the theoretical FSN
    outside_fsn <- sum(mat4[valid, 1, k] > fsn_sn[k]) +
      sum(mat4[valid, 2, k] < fsn_sn[k])
    
    coverage4[k] <- 1 - outside_fsn / sum(valid)
  }
}


# Skew-t distributions

# Number of Monte Carlo replications
R <- 10000

# Number of bootstrap resamples
B <- 1000

# One-sided 5% critical value of the standard normal distribution
za <- qnorm(0.95)

# Sample sizes considered
n <- c(5, 10, 15, 30, 50)

# Skewness parameter of the skew-t distribution
# Uncomment the desired value
# beta <- 0.2041
beta <- 0.5774
# beta <- 1.3333

# Degrees of freedom of the skew-t distribution
# Change to one of: 3, 5, 10, or 30
nu <- 30

# Corresponding skewness parameter
delta <- beta / sqrt(1 + beta^2)

# Skew-t moment constant
b_nu <- sqrt(nu) * gamma((nu - 1)/2) /
  (sqrt(pi) * gamma(nu/2))

# Theoretical mean and variance of the skew-t
m <- delta * b_nu
v <- nu/(nu - 2) - (delta^2 * b_nu^2)

# Storage matrix for simulated FSN
rosent5 <- matrix(NA, nrow = R, ncol = length(n))  # 

# Storage array for bootstrap confidence interval limits
mat5 <- array(NA, dim = c(R, 2, length(n)))        

# Coverage probabilities 
coverage5 <- numeric(length(n))

# Expected FSN
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

fsn_st <- (n^2 * m^2 + n * v) / za^2 - n + eps

# Function to compute the FSN estimator for bootstrap samples
fsn_estimator <- function(data, indices) {
  z_boot <- data[indices]
  S_boot <- sum(z_boot)
  fsn_boot <- (S_boot / za)^2 - length(z_boot)
  return(fsn_boot)  
}

## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate skew-t observations
    z <- rst(n[k],xi = 0,omega = 1,alpha = beta,nu = nu)
    
    S <- sum(z)
    
    # Compute Rosenthal's FSN
    if (S >= za * sqrt(n[k])) {
      
      # Significant meta-analysis
      fsn <- (S / za)^2 - n[k]  
 
    } else {
      
      # Non-significant meta-analysis
      fsn <- 0
    }
    
    rosent5[i, k] <- fsn
    
     # Bootstrap distribution of the FSN estimator
    tryCatch({
      boot_obj <- boot(data = z, statistic = fsn_estimator, R = B)
      
      # BCa confidence interval
      bcaci <- boot.ci(boot_obj, type = "bca", conf = 0.95)

      # Extract the lower and upper BCa limits
      if (!is.null(bcaci$bca)) {
        mat5[i, , k] <- bcaci$bca[4:5]
      }
    }, error = function(e) {
      
      # Return NA if BCa computation fails
      mat5[i, , k] <- c(NA, NA)
    })
  }
    
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat5[, 1, k]) & (rosent5[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability for the theoretical FSN
    outside_fsn <- sum(mat5[valid, 1, k] > fsn_st[k]) +
      sum(mat5[valid, 2, k] < fsn_st[k])
    
    coverage5[k] <- 1 - outside_fsn / sum(valid)
  }
}

