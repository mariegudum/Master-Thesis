set.seed(1)

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

z <- qnorm(1 - alpha/2)

# -----------------------------
# Storage
# -----------------------------
results <- data.frame()

# -----------------------------
# Simulation loop
# -----------------------------
for (n in n_list) {
  
  cat("Running n =", n, "\n")
  
  n_train <- n/2
  n_cal   <- n/2
  
  cov_oracle <- cov_lin <- cov_qr <- cov_cp <- numeric(M)
  len_oracle <- len_lin <- len_qr <- len_cp <- numeric(M)
  
  for (m in 1:M) {
    
    # -------------------------
    # Data
    # -------------------------
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      eps <- rnorm(n, 0, sigma)
      y <- f(x) + eps
      list(x=x, y=y)
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # -------------------------
    # (i) Oracle
    # -------------------------
    oracle_lower <- f(test$x) - z * sigma
    oracle_upper <- f(test$x) + z * sigma
    
    # -------------------------
    # (ii) Linear Gaussian
    # -------------------------
    lm_fit <- lm(train$y ~ train$x)
    
    pred_train <- predict(lm_fit, data.frame(train$x=train$x))
    sigma_hat <- sqrt(mean((train$y - pred_train)^2))
    
    pred_test <- predict(lm_fit, data.frame(train$x=test$x))
    
    lin_lower <- pred_test - z * sigma_hat
    lin_upper <- pred_test + z * sigma_hat
    
    # -------------------------
    # (iii) "QR" (residual method)
    # -------------------------
    res_train <- train$y - pred_train
    
    q_low  <- quantile(res_train, alpha/2)
    q_high <- quantile(res_train, 1 - alpha/2)
    
    qr_lower <- pred_test + q_low
    qr_upper <- pred_test + q_high
    
    # -------------------------
    # (iv) Conformal
    # -------------------------
    pred_cal <- predict(lm_fit, data.frame(train$x=cal$x))
    res_cal <- abs(cal$y - pred_cal)
    
    q_cp <- quantile(res_cal, 1 - alpha)
    
    cp_lower <- pred_test - q_cp
    cp_upper <- pred_test + q_cp
    
    # -------------------------
    # Metrics
    # -------------------------
    coverage <- function(l, u, y) mean(y >= l & y <= u)
    length   <- function(l, u) mean(u - l)
    
    cov_oracle[m] <- coverage(oracle_lower, oracle_upper, test$y)
    cov_lin[m]    <- coverage(lin_lower, lin_upper, test$y)
    cov_qr[m]     <- coverage(qr_lower, qr_upper, test$y)
    cov_cp[m]     <- coverage(cp_lower, cp_upper, test$y)
    
    len_oracle[m] <- length(oracle_lower, oracle_upper)
    len_lin[m]    <- length(lin_lower, lin_upper)
    len_qr[m]     <- length(qr_lower, qr_upper)
    len_cp[m]     <- length(cp_lower, cp_upper)
  }
  
  # -----------------------------
  # Aggregate
  # -----------------------------
  results <- rbind(results, data.frame(
    n = n,
    Method = c("Oracle", "Linear", "QR", "CP"),
    Coverage = c(mean(cov_oracle),
                 mean(cov_lin),
                 mean(cov_qr),
                 mean(cov_cp)),
    Length = c(mean(len_oracle),
               mean(len_lin),
               mean(len_qr),
               mean(len_cp))
  ))
}

print(results)