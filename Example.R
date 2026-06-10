
# Example: Meta-analysis and FSN Bootstrap
# Data source: Liu et al. (2015) meta-analysis on physical activity interventions
# and their effects on self-esteem and self-concept in children and adolescents.

library(readxl)
library(meta)
library(boot)

# Load dataset
data<- read_excel("Data_example.xlsx")
#View(data)

# Meta-analysis: RCT studies with PA intervention only

data1=data[data$Desenho=="RCT" & data$Intervencao=="AF",
           c("StudyID", "Outcome", "Total_Exp", "Mean_Exp", 
             "SD_Exp", "Total_Cont", "Mean_Cont", "SD_Cont")]

## Random-effects meta-analysis using Hedges' SMD
library(meta)
meta_analysis1=metacont(n.e = data1$Total_Exp, mean.e = data1$Mean_Exp, 
                        sd.e = data1$SD_Exp, n.c = data1$Total_Cont, 
                        mean.c = data1$Mean_Cont, sd.c = data1$SD_Cont,
                        data = data1, studlab = data1$StudyID, sm = "SMD",
                        method.smd = "Hedges", subgroup = data1$Outcome, 
                        random = T, common = F)

summary(meta_analysis1)
#forest(meta_analysis1)
#baujat(meta_analysis1)


# Prepare data for Rosenthal's Fail-Safe N
# (using metafor-style structure)

data_rosenthal <- data.frame(
  yi = meta_analysis1$TE,
  vi = meta_analysis1$seTE^2
)

# Compute Rosenthal's FSN
rosenthal_result <- fsn(x = data_rosenthal$yi, vi = data_rosenthal$vi,
                        type = "Rosenthal")
print(rosenthal_result)

fsn_obs <- rosenthal_result$fsnum


# Bootstrap function for FSN
fsn_boot <- function(data, indices) {
  d <- data[indices, ]
  
  tryCatch({
    fsn(x = d$yi,vi = d$vi,type = "Rosenthal" )$fsnum
    }, error = function(e) NA)
}

set.seed(123)

# Bootstrap resampling
boot_obj <- boot(data = data_rosenthal,statistic = fsn_boot,R = 10000)

# Remove failed bootstrap replications
boot_vals <- boot_obj$t[!is.na(boot_obj$t)]


# Confidence intervals for FSN

# Standard normal bootstrap CI
theta_hat <- boot_obj$t0
se_boot <- sd(boot_vals)

z <- qnorm(0.975)
ci_norm <- c(theta_hat - z * se_boot,theta_hat + z * se_boot)


### Percentile bootstrap CI

ci_perc <- quantile(boot_vals, probs = c(0.025, 0.975))

### BCa bootstrap CI

ci_bca <- boot.ci(boot_obj, type = "bca",conf = 0.95)$bca[4:5]


# Results 
results_ci <- data.frame(
  Method = c("Standard Normal", "Percentile", "BCa"),
  Lower = c(ci_norm[1], ci_perc[1], ci_bca[1]),
  Upper = c(ci_norm[2], ci_perc[2], ci_bca[2])
)

print(fsn_obs)
print(results_ci)
