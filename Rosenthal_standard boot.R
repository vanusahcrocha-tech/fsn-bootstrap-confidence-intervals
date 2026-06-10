### Standard normal aproximation confidence intervals


# install.packages("sn")
# install.packages("ggplot2")
# install.packages("patchwork")

library(sn)
library(ggplot2)
library(patchwork)

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

# Storage matrix for:
rosent1 <- matrix(NA, nrow = R, ncol = length(n))  # simulated FSN
v1 <- matrix(NA, nrow = R, ncol = length(n))   # estimated variance
coverage1 <- matrix(NA, nrow = length(n), ncol = 2) # coverage probabilities 

# Storage array for bootstrap confidence interval limits
mat1 <- array(NA, dim = c(R, 2, length(n))) 

# Proportion of statistically significant meta-analyses
prop_sig <- numeric(length(n))   

# Expected FSN (uncorrected)
fsn_normal <- (n^2 * m^2 + n * v) / za^2 - n

# Epsilon correction
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

# Epsilon-corrected expected FSN
fsn_norm_eps <- fsn_normal + eps


## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Counter for significant meta-analyses
  sig_count <- 0  
  
  # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate standard normal observations
    z <- rnorm(n[k])
    S <- sum(z)
    
    # Compute Rosenthal's FSN
    if (S >= za * sqrt(n[k])) {
      
      # Significant meta-analysis
      fsn <- (S / za)^2 - n[k]  
      
    } else {
      
      # Non-significant meta-analysis
      fsn <- 0
    }
    
    rosent1[i, k] <- fsn
    
    # Bootstrap distribution of the FSN estimator
    boot_fsn <- numeric(B)
    
    for (j in 1:B) {
      
      # Resample indices
      indices <- sample(seq_len(n[k]), n[k], replace = TRUE)
      
      # Bootstrap FSN estimate
      boot_fsn_j <- (sum(z[indices]) / za)^2 - n[k]
      
      # Truncate at zero (FSN definition)
      boot_fsn[j] <- pmax(0,boot_fsn_j)  
      
    }
    
    # Bootstrap variance and standard error
    v1[i, k] <- var(boot_fsn)
    se <- sqrt(v1[i, k])
    
    # Standard bootstrap confidence interval
    mat1[i, , k] <- c(
      fsn - 1.96 * se,
      fsn + 1.96 * se
    )
  }

  # Proportion of significant meta-analyses
  prop_sig[k] <- sig_count/R
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat1[, 1, k]) & (rosent1[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability (without epsilon)
    outside_no_eps <- sum(mat1[valid, 1, k] > fsn_normal[k]) +
      sum(mat1[valid, 2, k] < fsn_normal[k])
    
    coverage1[k, 1] <- 1 - outside_no_eps / sum(valid)
    
    
    # Coverage probability (with epsilon)
    outside_eps <- sum(mat1[valid, 1, k] > fsn_norm_eps[k]) +
      sum(mat1[valid, 2, k] < fsn_norm_eps[k])
    
    coverage1[k, 2] <- 1 - outside_eps / sum(valid)
    
  }
  
}
colnames(coverage1) <- c("Cov_no_eps", "Cov_eps")
rownames(coverage1) <- paste0("n=", n)


#### Half-normal distribution

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

# Storage matrix for:
rosent2 <- matrix(NA, nrow = R, ncol = length(n))  # simulated FSN
v2 <- matrix(NA, nrow = R, ncol = length(n))   # estimated variance
coverage2 <- matrix(NA, nrow = length(n), ncol = 2) # coverage probabilities 

# Storage array for bootstrap confidence interval limits
mat2 <- array(NA, dim = c(R, 2, length(n))) 

# Proportion of statistically significant meta-analyses
prop_sig <- numeric(length(n))   

# Expected FSN (uncorrected)
fsn_hn <- (n^2 * m^2 + n * v) / za^2 - n

# Epsilon correction
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

# Epsilon-corrected expected FSN
fsn_hn_eps <- fsn_hn + eps


## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Counter for significant meta-analyses
  sig_count <- 0  
  
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
    boot_fsn <- numeric(B)
    
    for (j in 1:B) {
      
      # Resample indices
      indices <- sample(seq_len(n[k]), n[k], replace = TRUE)
      
      # Bootstrap FSN estimate
      boot_fsn_j <- (sum(z[indices]) / za)^2 - n[k]
      
      # Truncate at zero (FSN definition)
      boot_fsn[j] <- pmax(0,boot_fsn_j)  
      
    }
    
    # Bootstrap variance and standard error
    v2[i, k] <- var(boot_fsn)
    se <- sqrt(v2[i, k])
    
    # Standard bootstrap confidence interval
    mat2[i, , k] <- c(
      fsn - 1.96 * se,
      fsn + 1.96 * se
    )
  }
  
  # Proportion of significant meta-analyses
  prop_sig[k] <- sig_count/R
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat2[, 1, k]) & (rosent2[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability (no epsilon)
    outside_no_eps <- sum(mat2[valid, 1, k] > fsn_hn[k]) +
      sum(mat2[valid, 2, k] < fsn_hn[k])
    
    coverage2[k, 1] <- 1 - outside_no_eps / sum(valid)
    
    
    # Coverage probability (with epsilon)
    outside_eps <- sum(mat2[valid, 1, k] > fsn_hn_eps[k]) +
      sum(mat2[valid, 2, k] < fsn_hn_eps[k])
    
    coverage2[k, 2] <- 1 - outside_eps / sum(valid)
    
  }
  
}
colnames(coverage2) <- c("Cov_no_eps", "Cov_eps")
rownames(coverage2) <- paste0("n=", n)



#### Skew Normal distribution (negative skewness)

# Number of Monte Carlo replications
R <- 10000

# Number of bootstrap resamples
B <- 1000

# One-sided 5% critical value of the standard normal distribution
za <- qnorm(0.95)

# Sample sizes considered
n <- c(5, 10, 15, 30, 50)

# Theoretical mean and variance of the skew-normal
# Set delta to one of the desired skewness levels:-0.2, -0.5, -0.8
delta <- -0.5
m <- delta*sqrt(2/pi)
v <- 1-2*(delta^2)/pi

# Storage matrix for:
rosent3 <- matrix(NA, nrow = R, ncol = length(n))  # simulated FSN
v3 <- matrix(NA, nrow = R, ncol = length(n))   # estimated variance
coverage3 <- matrix(NA, nrow = length(n), ncol = 2) # coverage probabilities 

# Storage array for bootstrap confidence interval limits
mat3 <- array(NA, dim = c(R, 2, length(n))) 

# Proportion of statistically significant meta-analyses
prop_sig <- numeric(length(n))   

# Expected FSN (uncorrected)
fsn_sn1 <- (n^2 * m^2 + n * v) / za^2 - n

# Epsilon correction
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

# Epsilon-corrected expected FSN
fsn_sn_eps1 <- fsn_sn1 + eps


## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Counter for significant meta-analyses
  sig_count <- 0  
  
  # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate skew-normal observations
    # alpha values correspond to skewness levels delta = -0.2, -0.5, and -0.8
    #z <- rsn(n[k],xi=0,omega=1,alpha= -0.2041)
    z <- rsn(n[k], xi=0, omega=1, alpha= -0.5773503)
    #z <- rsn(n[k],xi=0,omega=1,alpha= -1.3333)

    S <- sum(z)
  
    # Compute Rosenthal's FSN
    if (S >= za * sqrt(n[k])) {
      
      # Significant meta-analysis
      fsn <- (S / za)^2 - n[k] 
      
    } else {
      
      # Non-significant meta-analysis
      fsn <- 0
    }
    
    rosent3[i, k] <- fsn
    
    # Bootstrap distribution of the FSN estimator
    boot_fsn <- numeric(B)
    
    for (j in 1:B) {
      
      # Resample indices
      indices <- sample(seq_len(n[k]), n[k], replace = TRUE)
      
      # Bootstrap FSN estimate
      boot_fsn_j <- (sum(z[indices]) / za)^2 - n[k]
      
      # Truncate at zero (FSN definition)
      boot_fsn[j] <- pmax(0,boot_fsn_j)  
      
    }
    
    # Bootstrap variance and standard error
    v3[i, k] <- var(boot_fsn)
    se <- sqrt(v3[i, k])
    
    # Standard bootstrap confidence interval
    mat3[i, , k] <- c(
      fsn - 1.96 * se,
      fsn + 1.96 * se
    )
  }
  
  # Proportion of significant meta-analyses
  prop_sig[k] <- sig_count/R
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat3[, 1, k]) & (rosent3[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability (no epsilon)
    outside_no_eps <- sum(mat3[valid, 1, k] > fsn_sn1[k]) +
      sum(mat3[valid, 2, k] < fsn_sn1[k])
    
    coverage3[k, 1] <- 1 - outside_no_eps / sum(valid)
    
    
    # Coverage probability (with epsilon)
    outside_eps <- sum(mat3[valid, 1, k] > fsn_sn_eps1[k]) +
      sum(mat3[valid, 2, k] < fsn_sn_eps1[k])
    
    coverage3[k, 2] <- 1 - outside_eps / sum(valid)
    
  }
  
}
colnames(coverage3) <- c("Cov_no_eps", "Cov_eps")
rownames(coverage3) <- paste0("n=", n)



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
delta <- 0.5
m <- delta*sqrt(2/pi)
v <- 1-2*(delta^2)/pi

# Storage matrix for:
rosent4 <- matrix(NA, nrow = R, ncol = length(n))  # simulated FSN
v4 <- matrix(NA, nrow = R, ncol = length(n))   # estimated variance
coverage4 <- matrix(NA, nrow = length(n), ncol = 2) # coverage probabilities 

# Storage array for bootstrap confidence interval limits
mat4 <- array(NA, dim = c(R, 2, length(n))) 

# Proportion of statistically significant meta-analyses
prop_sig <- numeric(length(n))   

# Expected FSN (uncorrected)
fsn_sn2 <- (n^2 * m^2 + n * v) / za^2 - n

# Epsilon correction
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

# Epsilon-corrected expected FSN
fsn_sn_eps2 <- fsn_sn2 + eps


## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Counter for significant meta-analyses
  sig_count <- 0  
  
  # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate skew-normal observations
    # alpha values correspond to skewness levels delta = 0.2, 0.5, and 0.8
    #z <- rsn(n[k],xi=0,omega=1,alpha= 0.2041)
    z <- rsn(n[k], xi=0, omega=1, alpha= 0.5773503)
    #z <- rsn(n[k],xi=0,omega=1,alpha= 1.3333)
    
    S <- sum(z)
    
    # Compute Rosenthal's FSN
    if (S >= za * sqrt(n[k])) {
      
      # Significant meta-analysis
      fsn <- (S / za)^2 - n[k]  
      
    } else {
      
      # Non-significant meta-analysis
      fsn <- 0
    }
    
    rosent4[i, k] <- fsn
    
    # Bootstrap distribution of the FSN estimator
    boot_fsn <- numeric(B)
    
    for (j in 1:B) {
      
      # Resample indices
      indices <- sample(seq_len(n[k]), n[k], replace = TRUE)
      
      # Bootstrap FSN estimate
      boot_fsn_j <- (sum(z[indices]) / za)^2 - n[k]
      
      # Truncate at zero (FSN definition)
      boot_fsn[j] <- pmax(0,boot_fsn_j)  
      
    }
    
    # Bootstrap variance and standard error
    v4[i, k] <- var(boot_fsn)
    se <- sqrt(v4[i, k])
    
    # Standard bootstrap confidence interval
    mat4[i, , k] <- c(
      fsn - 1.96 * se,
      fsn + 1.96 * se
    )
  }
  
  # Proportion of significant meta-analyses
  prop_sig[k] <- sig_count/R
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat4[, 1, k]) & (rosent4[, k] > 0)
  
  if (sum(valid) > 0) {
    
    # Coverage probability (no epsilon)
    outside_no_eps <- sum(mat4[valid, 1, k] > fsn_sn2[k]) +
      sum(mat4[valid, 2, k] < fsn_sn2[k])
    
    coverage4[k, 1] <- 1 - outside_no_eps / sum(valid)
    
    
    # Coverage probability (with epsilon)
    outside_eps <- sum(mat4[valid, 1, k] > fsn_sn_eps2[k]) +
      sum(mat4[valid, 2, k] < fsn_sn_eps2[k])
    
    coverage4[k, 2] <- 1 - outside_eps / sum(valid)
    
  }
  
}
colnames(coverage4) <- c("Cov_no_eps", "Cov_eps")
rownames(coverage4) <- paste0("n=", n)


#################################################
# Figure: Behavior of Phi(lambda*) and epsilon under normal, half-normal,
# and skew-normal distributions (delta = +/-0.5)

# Parameters
alpha <- 0.05
Z_alpha <- qnorm(1 - alpha)
k <- 1:50

mu_norm <- 0
sigma_norm <- 1

mu_half <- sqrt(2/pi)
sigma_half <- sqrt(1 - 2/pi)

mu_skew_neg <- -sqrt(1/(2*pi))
mu_skew_pos <-  sqrt(1/(2*pi))
sigma_skew <- sqrt(1 - 1/(2*pi))

# Lambda

lambda_norm <- (sqrt(k)*mu_norm - Z_alpha)/sigma_norm
lambda_half <- (sqrt(k)*mu_half - Z_alpha)/sigma_half
lambda_skew_neg <- (sqrt(k)*mu_skew_neg - Z_alpha)/sigma_skew
lambda_skew_pos <- (sqrt(k)*mu_skew_pos - Z_alpha)/sigma_skew

# Phi e phi

Phi_norm <- pnorm(lambda_norm)
Phi_half <- pnorm(lambda_half)
Phi_skew_neg <- pnorm(lambda_skew_neg)
Phi_skew_pos <- pnorm(lambda_skew_pos)

phi_norm <- dnorm(lambda_norm)
phi_half <- dnorm(lambda_half)
phi_skew_neg <- dnorm(lambda_skew_neg)
phi_skew_pos <- dnorm(lambda_skew_pos)

# Epsilon

eps_norm <- (phi_norm / Phi_norm) * (k * sigma_norm * (sqrt(k)*mu_norm + Z_alpha)) / (Z_alpha^2)
eps_half <- (phi_half / Phi_half) * (k * sigma_half * (sqrt(k)*mu_half + Z_alpha)) / (Z_alpha^2)
eps_skew_neg <- (phi_skew_neg / Phi_skew_neg) * (k * sigma_skew * (sqrt(k)*mu_skew_neg + Z_alpha)) / (Z_alpha^2)
eps_skew_pos <- (phi_skew_pos / Phi_skew_pos) * (k * sigma_skew * (sqrt(k)*mu_skew_pos + Z_alpha)) / (Z_alpha^2)

# Data frames

df_phi <- data.frame(
  k = rep(k, 4),
  value = c(Phi_norm, Phi_half, Phi_skew_neg, Phi_skew_pos),
  Distribution = factor(rep(c("N", "HN", "SN-", "SN+"), each = length(k)))
)

df_eps <- data.frame(
  k = rep(k, 4),
  value = c(eps_norm, eps_half, eps_skew_neg, eps_skew_pos),
  Distribution = factor(rep(c("N", "HN", "SN-", "SN+"), each = length(k)))
)

df_phi$Distribution <- factor(df_phi$Distribution, levels = c("N", "HN", "SN-", "SN+"))
df_eps$Distribution <- factor(df_eps$Distribution, levels = c("N", "HN", "SN-", "SN+"))

# Plot Phi

p_phi <- ggplot(df_phi, aes(x = k, y = value, color = Distribution)) +
  geom_line(linewidth = 0.9) +
  labs(
    x = "Number of studies (k)",
    y = expression(Phi(lambda^"*")),
    color = NULL
  ) +
  scale_color_manual(
    values = c("blue","purple", "red", "darkgreen"),
    labels = c(
      "N(0,1)",
      "HN(0,1)",
      expression(SN~(delta == -0.5)),
      expression(SN~(delta == 0.5))
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(size = 9),
    axis.text.y  = element_text(size = 9),
    legend.position = "bottom"
  )

# Plot Epsilon

p_eps <- ggplot(df_eps, aes(x = k, y = value, color = Distribution)) +
  geom_line(linewidth = 0.9) +
  labs(
    x = "Number of studies (k)",
    y = expression(epsilon),
    color = NULL
  ) +
  scale_color_manual(
    values = c("blue","purple", "red", "darkgreen"),
    labels = c(
      "N(0,1)",
      "HN(0,1)",
      expression(SN~(delta == -0.5)),
      expression(SN~(delta == 0.5))
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(size = 9),
    axis.text.y  = element_text(size = 9),
    legend.position = "bottom"
  )


# Combine plots
p_final <- (p_phi + p_eps) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")


# Save figure
ggsave("figure1.pdf", p_final, width = 10, height = 4)


#####################################

### Boxplots of positive FSN values
### For skew-normal delta=0.5

make_boxplot <- function(data_matrix, n, panel_label) {
  
  df <- data.frame(
    FSN = as.vector(data_matrix),
    SampleSize = factor(rep(n, each = nrow(data_matrix)),
                        levels = c(5, 10, 15, 30, 50))
  )
  
  # Keep only positive FSN values
  df <- subset(df, FSN > 0)
  
  ggplot(df, aes(x = SampleSize, y = FSN)) +
    geom_boxplot(fill = "gray70", colour = "black") +
    scale_x_discrete(drop = FALSE) +   
    labs(
      title = panel_label,
      x = "Number of studies (k)",
      y = "FSN"
    ) +
    theme_classic() +
    theme(plot.margin = margin(5, 5, 5, 5))
}


# Generate plots
p1 <- make_boxplot(rosent1, n, "A")
p2 <- make_boxplot(rosent2, n, "B")
p3 <- make_boxplot(rosent3, n, "C")
p4 <- make_boxplot(rosent4, n, "D")

p <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(plot.margin = margin(5, 5, 5, 5))
p

# Save figure
ggsave("figure3.pdf", plot = p, width = 6, height = 4)


#########################################
   
# Skew-t distributions
# we calculated the coverage only taking account 
## the inclusion of the adjustment epsilon in the expected FSN value.

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
s <- nu/(nu - 2) - (delta^2 * b_nu^2)

# Storage matrix for:
rosent5 <- matrix(NA, nrow = R, ncol = length(n))  # simulated FSN
v5 <- matrix(NA, nrow = R, ncol = length(n))   # estimated variance

# Storage array for bootstrap confidence interval limits
mat5 <- array(NA, dim = c(R, 2, length(n)))        

# Coverage probabilities 
coverage5 <- numeric(length(n))

# Proportion of statistically significant meta-analyses
prop_sig <- numeric(length(n))   

# Expected FSN
lam_ast <- (sqrt(n) * m - za) / sqrt(v)
eps <- (dnorm(lam_ast) / pnorm(lam_ast)) * 
  (n * sqrt(v) * (sqrt(n) * m + za) / za^2)

fsn_st <- (n^2 * m^2 + n * v) / za^2 - n + eps


## seed number
set.seed(123456)

# Loop over sample sizes
for (k in 1:length(n)) {
  
  # Counter for significant meta-analyses
  sig_count <- 0 
  
  # Monte Carlo replications
  for (i in 1:R) {
    
    # Generate skew-t observations
    z <- rst(n[k],xi = 0,omega = 1,alpha = beta,nu = nu)
    
    # Sum of study Z-statistics
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
    
    rosent5[i, k] <- fsn
    
    # Bootstrap distribution of the FSN estimator
    boot_fsn <- numeric(B)
    for (j in 1:B) {
      
      # Bootstrap resampling indices
      indices <- sample(seq_len(n[k]), n[k], replace = TRUE)
      
      # Bootstrap FSN estimate
      boot_fsn_j <- (sum(z[indices]) / za)^2 - n[k]
      
      # Truncate at zero, following the FSN definition
      boot_fsn[j] <- pmax(0,boot_fsn_j)
    }
    
    # Bootstrap variance and standard error
    v5[i, k] <- var(boot_fsn)
    se <- sqrt(v5[i, k])
    
    # Standard bootstrap confidence interval
    mat5[i, , k] <- c(
      fsn - 1.96 * se,
      fsn + 1.96 * se
    )
  }
  
  # Proportion of significant meta-analyses
  prop_sig[k] <- sig_count/R
  
  # Consider only significant meta-analyses
  valid <- !is.na(mat5[, 1, k]) & (rosent5[, k] > 0)
  
  if (sum(valid) > 0) {
    # Coverage probability for the theoretical FSN
    outside_fsn <- sum(mat5[valid, 1, k] > fsn_st[k]) +
      sum(mat5[valid, 2, k] < fsn_st[k])
    
    coverage5[k] <- 1 - outside_fsn / sum(valid)

  }
}

#########################################
# Density functions of the distributions considered in the simulation study

x <- seq(-4, 6, length.out = 2000)


# A: Normal + Half-normal

df_norm <- data.frame(x = x, y = dnorm(x), dist = "N(0,1)")

df_half <- data.frame(x = x[x >= 0], y = sqrt(2/pi) * exp(-x[x >= 0]^2 / 2),
  dist = "HN(0,1)")

df_vert <- data.frame(x = 0, y = 0, yend = sqrt(2/pi))

pA <- ggplot() +
  geom_line(data = df_norm,
            aes(x = x, y = y, color = dist),
            linewidth = 0.9) +
  geom_line(data = df_half,
            aes(x = x, y = y, color = dist),
            linewidth = 0.9) +
  geom_segment(data = df_vert,
               aes(x = x, xend = x, y = y, yend = yend),
               linetype = "dashed",
               color = "gray50",
               linewidth = 0.6) +
  annotate("text",x = -3, y = Inf, label = "A", size = 5, vjust = 1.5) +
  scale_color_manual(values = c("N(0,1)" = "red3","HN(0,1)" = "forestgreen")) +
  labs(x = "x", y = "Density") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 10))


# B: Skew-normal (negative)

alphas_neg <- c(-0.2041, -0.5774, -1.3333)
labels_neg <- c("δ = -0.2", "δ = -0.5", "δ = -0.8")

dfB <- do.call(rbind,lapply(seq_along(alphas_neg), function(i) {
  data.frame(x = x, density = dsn(x,xi = 0,omega = 1,
                                  alpha = alphas_neg[i]),dist = labels_neg[i])}))

pB <- ggplot(dfB,aes(x = x, y = density, color = dist,linetype = dist)) +
  geom_line(linewidth = 0.9) +
  annotate("text",x = -3,y = Inf,label = "B",size = 5,vjust = 1.5) +
  scale_color_manual(values = c("δ = -0.2" = "red3","δ = -0.5" = "red2",
      "δ = -0.8" = "red4")) +
  scale_linetype_manual(values = c("δ = -0.2" = 1,"δ = -0.5" = 2,
                                   "δ = -0.8" = 3)) +
  labs(x = "x", y = "Density") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 10))


# C: Skew-normal (positive)

alphas_pos <- c(0.2041, 0.5774, 1.3333)
labels_pos <- c("δ = 0.2", "δ = 0.5", "δ = 0.8")

dfC <- do.call(rbind,
               lapply(seq_along(alphas_pos), function(i) {
                 data.frame(x = x,
                   density = dsn(x,xi = 0,omega = 1,alpha = alphas_pos[i]),
                   dist = labels_pos[i])
                 }))

pC <- ggplot(dfC,aes(x = x,y = density,color = dist,linetype = dist)) +
  geom_line(linewidth = 0.9) +
  annotate("text",x = -3,y = Inf,label = "C",size = 5, vjust = 1.5) +
  scale_color_manual(values = c("δ = 0.2" = "purple3","δ = 0.5" = "purple2",
      "δ = 0.8" = "purple4")) +
  scale_linetype_manual(
    values = c("δ = 0.2" = 1,"δ = 0.5" = 2,"δ = 0.8" = 3)) +
  labs(x = "x", y = "Density") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 10))


# D: Skew-t

nu_values <- c(3, 5, 10, 30)
labels_D <- c("ν = 3", "ν = 5", "ν = 10", "ν = 30")

dfD <- do.call(rbind,lapply(seq_along(nu_values), function(i) {
                 data.frame(x = x,density = dst(x,xi = 0,omega = 1,
                     alpha = 0.5774,nu = nu_values[i]),dist = labels_D[i])
               }))

pD <- ggplot(dfD,aes(x = x,y = density,color = dist,linetype = dist)) +
  geom_line(linewidth = 0.9) +
  annotate("text",x = -3,y = Inf,label = "D",size = 5,vjust = 1.5) +
  scale_color_manual(values = c("ν = 3" = "blue","ν = 5" = "red",
      "ν = 10" = "green4","ν = 30" = "black")) +
  scale_linetype_manual(values = c("ν = 3" = 1,"ν = 5" = 2,"ν = 10" = 3,
      "ν = 30" = 4)) +
  labs(x = "x", y = "Density") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 10))


# Final figure
p <- (pA | pB) / (pC | pD)

p


ggsave("figure2.png", plot = p,width = 8,height = 6, dpi = 600)
