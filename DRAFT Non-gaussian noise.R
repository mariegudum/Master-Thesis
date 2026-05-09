# =============================
# Setup
# =============================
set.seed(1)

library(ggplot2)
library(tidyr)

# -----------------------------
# Parameters
# -----------------------------
n_list <- c(25, 100, 400, 1600)
M <- 200
n_test <- 2000

alpha <- 0.1
beta  <- 2
a     <- 1
sigma <- 1

f <- function(x) beta * x + a

# t-distribution scaling (df = 3)
df <- 3
scale_t <- sigma / sqrt(df / (df - 2))

# "Gaussian oracle" (misspecified now)
z <- qnorm(1 - alpha/2)

# True oracle (optional)
q_true <- qt(1 - alpha/2, df = df) * scale_t

# -----------------------------
# Storage
# -----------------------------
results_all <- data.frame()

# -----------------------------
# Simulation loop
# -----------------------------
for (n in n_list) {
  
  cat("Running n =", n, "\n")
  
  n_train <- n/2
  n_cal   <- n/2
  
  cov_oracle <- cov_lin <- cov_qr <- cov_cp <- cov_naive <- numeric(M)
  len_oracle <- len_lin <- len_qr <- len_cp <- len_naive <- numeric(M)
  
  for (m in 1:M) {
    
    # -------------------------
    # Data (NON-GAUSSIAN)
    # -------------------------
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      eps <- rt(n, df = df) * scale_t
      y <- f(x) + eps
      data.frame(x = x, y = y)
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # -------------------------
    # (i) Oracle (Gaussian WRONG)
    # -------------------------
    oracle_lower <- f(test$x) - z * sigma
    oracle_upper <- f(test$x) + z * sigma
    
    # (optional) TRUE oracle:
    # oracle_lower <- f(test$x) - q_true
    # oracle_upper <- f(test$x) + q_true
    
    # -------------------------
    # (ii) Linear Gaussian
    # -------------------------
    lm_fit <- lm(y ~ x, data = train)
    
    pred_train <- predict(lm_fit, newdata = train)
    pred_cal   <- predict(lm_fit, newdata = cal)
    pred_test  <- predict(lm_fit, newdata = test)
    
    sigma_hat <- sqrt(mean((train$y - pred_train)^2))
    
    lin_lower <- pred_test - z * sigma_hat
    lin_upper <- pred_test + z * sigma_hat
    
    # -------------------------
    # (iii) QR (residual quantiles)
    # -------------------------
    res_train <- train$y - pred_train
    
    q_low  <- quantile(res_train, alpha/2)
    q_high <- quantile(res_train, 1 - alpha/2)
    
    qr_lower <- pred_test + q_low
    qr_upper <- pred_test + q_high
    
    # -------------------------
    # (iv) Conformal
    # -------------------------
    res_cal <- abs(cal$y - pred_cal)
    
    k <- ceiling((n_cal + 1) * (1 - alpha))
    k <- min(k, n_cal)
    
    q_cp <- sort(res_cal)[k]
    
    cp_lower <- pred_test - q_cp
    cp_upper <- pred_test + q_cp
    
    # -------------------------
    # (v) Naive (NO split!)
    # -------------------------
    q_low_naive  <- quantile(res_train, alpha/2)
    q_high_naive <- quantile(res_train, 1 - alpha/2)
    
    naive_lower <- pred_test + q_low_naive
    naive_upper <- pred_test + q_high_naive
    
    # -------------------------
    # Metrics
    # -------------------------
    coverage <- function(l, u, y) mean(y >= l & y <= u)
    length   <- function(l, u) mean(u - l)
    
    cov_oracle[m] <- coverage(oracle_lower, oracle_upper, test$y)
    cov_lin[m]    <- coverage(lin_lower, lin_upper, test$y)
    cov_qr[m]     <- coverage(qr_lower, qr_upper, test$y)
    cov_cp[m]     <- coverage(cp_lower, cp_upper, test$y)
    cov_naive[m]  <- coverage(naive_lower, naive_upper, test$y)
    
    len_oracle[m] <- length(oracle_lower, oracle_upper)
    len_lin[m]    <- length(lin_lower, lin_upper)
    len_qr[m]     <- length(qr_lower, qr_upper)
    len_cp[m]     <- length(cp_lower, cp_upper)
    len_naive[m]  <- length(naive_lower, naive_upper)
  }
  
  # -------------------------
  # Store results
  # -------------------------
  results_all <- rbind(results_all, data.frame(
    n = n,
    Method = rep(c("Oracle", "Linear", "QR", "CP", "Naive"), each = M),
    Coverage = c(cov_oracle, cov_lin, cov_qr, cov_cp, cov_naive),
    Length   = c(len_oracle, len_lin, len_qr, len_cp, len_naive)
  ))
}

# =============================
# Reshape
# =============================
results_long <- results_all |>
  pivot_longer(cols = c(Coverage, Length),
               names_to = "Metric",
               values_to = "Value")

results_long$Method <- factor(
  results_long$Method,
  levels = c("Oracle", "Linear", "Naive", "QR", "CP")
)

# =============================
# Plot: Methods vs n (farve = n)
# =============================
p_all <- ggplot(results_long,
                aes(x = Method, y = Value, fill = factor(n))) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.75)) +
  facet_wrap(~ Metric, scales = "free_y") +
  geom_hline(
    data = data.frame(Metric = "Coverage"),
    aes(yintercept = 1 - alpha),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  labs(fill = "n") +
  theme_bw()

print(p_all)
