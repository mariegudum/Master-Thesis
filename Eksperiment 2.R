set.seed(1)

library(ggplot2)
library(dplyr)
library(ranger)

# =============================
# Setup
# =============================
n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta * x + a

# Exponential noise (mean 0)
gen_eps <- function(n) rexp(n, rate = 1/sigma) - sigma

# Quantile regression (linear)
rho <- function(u, tau) u * (tau - (u < 0))
fit_qr <- function(X, y, tau) {
  optim(rep(0, ncol(X)),
        function(b) sum(rho(y - X %*% b, tau)),
        method = "BFGS")$par
}

# =============================
# Progress setup
# =============================
total_iters <- length(n_list) * M
counter <- 1
start_time <- Sys.time()

results_all_2 <- vector("list", total_iters)

# =============================
# Simulation
# =============================
for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    # ===== Progress =====
    if (counter %% 5 == 0 || counter == 1) {
      elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
      progress <- counter / total_iters
      
      eta <- if (counter > 1) {
        elapsed / (counter - 1) * (total_iters - counter + 1)
      } else NA
      
      cat(sprintf(
        "Iter %d/%d (%.1f%%) | n=%d, m=%d | ETA=%.1fs\n",
        counter, total_iters, 100*progress, n, m, eta
      ))
    }
    
    # =========================
    # Data
    # =========================
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      data.frame(x = x, y = f(x) + gen_eps(n))
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =========================
    # LM
    # =========================
    lm_fit <- lm(y ~ x, data = train)
    lm_pred_train <- predict(lm_fit, train)
    lm_pred_cal   <- predict(lm_fit, cal)
    lm_pred_test  <- predict(lm_fit, test)
    
    # =========================
    # RF (tuned)
    # =========================
    rf_fit <- ranger(
      y ~ x,
      data = train,
      num.trees = 500,
      min.node.size = 30,
      sample.fraction = 0.8,
      num.threads = parallel::detectCores()
    )
    
    rf_pred_cal  <- predict(rf_fit, cal)$predictions
    rf_pred_test <- predict(rf_fit, test)$predictions
    
    # =========================
    # Oracle
    # =========================
    q_low  <- qexp(alpha/2,     rate = 1/sigma) - sigma
    q_high <- qexp(1 - alpha/2, rate = 1/sigma) - sigma
    
    oracle_l <- f(test$x) + q_low
    oracle_u <- f(test$x) + q_high
    
    # =========================
    # Gaussian (LM)
    # =========================
    z <- qnorm(1 - alpha/2)
    sigma_hat <- sqrt(mean((train$y - lm_pred_train)^2))
    
    gauss_l <- lm_pred_test - z*sigma_hat
    gauss_u <- lm_pred_test + z*sigma_hat
    
    # =========================
    # RQ
    # =========================
    rq_l_lm <- lm_pred_test + quantile(cal$y - lm_pred_cal, alpha/2)
    rq_u_lm <- lm_pred_test + quantile(cal$y - lm_pred_cal, 1 - alpha/2)
    
    rq_l_rf <- rf_pred_test + quantile(cal$y - rf_pred_cal, alpha/2)
    rq_u_rf <- rf_pred_test + quantile(cal$y - rf_pred_cal, 1 - alpha/2)
    
    # =========================
    # QR
    # =========================
    X_train <- model.matrix(~x, train)
    X_test  <- model.matrix(~x, test)
    
    b_l <- fit_qr(X_train, train$y, alpha/2)
    b_u <- fit_qr(X_train, train$y, 1 - alpha/2)
    
    qr_l_lm <- X_test %*% b_l
    qr_u_lm <- X_test %*% b_u
    
    # =========================
    # QR_RF per-n tuning
    # =========================
    
    grid <- expand.grid(
      min.node.size = c(10, 20, 30, 50, 80, 120, 200),
      sample.fraction = c(0.6, 0.8)
    )
    
    best_score <- Inf
    best_model <- NULL
    
    for (g in 1:nrow(grid)) {
      
      rf_tmp <- ranger(
        y ~ x,
        data = train,
        quantreg = TRUE,
        num.trees = 300,
        min.node.size = grid$min.node.size[g],
        sample.fraction = grid$sample.fraction[g],
        num.threads = parallel::detectCores()
      )
      
      # predictions på calibration
      pred_cal <- predict(
        rf_tmp, cal,
        type = "quantiles",
        quantiles = c(alpha/2, 1-alpha/2)
      )$predictions
      
      # pinball loss
      u1 <- cal$y - pred_cal[,1]
      u2 <- cal$y - pred_cal[,2]
      
      loss <- mean(
        u1 * (alpha/2 - (u1 < 0)) +
          u2 * ((1 - alpha/2) - (u2 < 0))
      )
      
      if (loss < best_score) {
        best_score <- loss
        best_model <- rf_tmp
      }
    }
    
    # ===== final prediction =====
    qr_pred <- predict(
      best_model, test,
      type = "quantiles",
      quantiles = c(alpha/2, 1-alpha/2)
    )$predictions
    
    qr_l_rf <- qr_pred[,1]
    qr_u_rf <- qr_pred[,2]
    
    # =========================
    # CP
    # =========================
    k <- ceiling((n_cal + 1)*(1 - alpha))
    
    q_lm <- sort(abs(cal$y - lm_pred_cal))[k]
    cp_l_lm <- lm_pred_test - q_lm
    cp_u_lm <- lm_pred_test + q_lm
    
    q_rf <- sort(abs(cal$y - rf_pred_cal))[k]
    cp_l_rf <- rf_pred_test - q_rf
    cp_u_rf <- rf_pred_test + q_rf
    
    # =========================
    # Metrics
    # =========================
    cov <- function(l,u) mean(test$y >= l & test$y <= u)
    wid <- function(l,u) mean(u - l)
    
    results_all_2[[counter]] <- data.frame(
      n=n, m=m,
      method=c("Oracle",
               "Gaussian_LM",
               "RQ_LM","RQ_RF_tuned",
               "QR_LM","QR_RF_tuned",
               "CP_LM","CP_RF_tuned"),
      coverage=c(
        cov(oracle_l,oracle_u),
        cov(gauss_l,gauss_u),
        cov(rq_l_lm,rq_u_lm),
        cov(rq_l_rf,rq_u_rf),
        cov(qr_l_lm,qr_u_lm),
        cov(qr_l_rf,qr_u_rf),
        cov(cp_l_lm,cp_u_lm),
        cov(cp_l_rf,cp_u_rf)
      ),
      length=c(
        wid(oracle_l,oracle_u),
        wid(gauss_l,gauss_u),
        wid(rq_l_lm,rq_u_lm),
        wid(rq_l_rf,rq_u_rf),
        wid(qr_l_lm,qr_u_lm),
        wid(qr_l_rf,qr_u_rf),
        wid(cp_l_lm,cp_u_lm),
        wid(cp_l_rf,cp_u_rf)
      )
    )
    
    counter <- counter + 1
  }
}

results_all_2 <- bind_rows(results_all_2)

# =============================
# Runtime
# =============================
total_time <- Sys.time() - start_time
cat("\nTotal runtime:", round(as.numeric(total_time, units="secs"),1), "seconds\n")

results_all_2 %>%
  group_by(n, method) %>%
  summarise(
    mean_cov = mean(coverage),
    sd_cov   = sd(coverage),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size (Exponential noise)",
    y = "Coverage"
  ) +
  theme_minimal()

oracle_mean_2 <- results_all_2 %>%
  filter(method == "Oracle") %>%
  group_by(n) %>%
  summarise(
    oracle_length = mean(length),
    .groups = "drop"
  )

rel_results_2 <- results_all_2 %>%
  left_join(oracle_mean_2, by = "n") %>%
  mutate(rel_length = length / oracle_length)

rel_results_2 %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(
    mean_rel = mean(rel_length),
    sd_rel   = sd(rel_length),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = mean_rel, color = method)) +
  geom_line() +
  geom_point() +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  scale_x_log10() +
  labs(
    title = "Relative interval length (Exponential noise)",
    y = "Length / Oracle length"
  ) +
  theme_minimal()


# Gaussian metode (LM Gaussian) → dårlig
# RQ / QR → bedre (asymmetri håndteres)
# CP → stadig korrekt coverage