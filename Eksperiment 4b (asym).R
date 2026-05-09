# =============================
# Setup
# =============================
set.seed(1)
library(dplyr)
library(ggplot2)
library(ranger)

n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1

beta <- 2; a <- 1
f <- function(x) beta * x + a

# =============================
# Asymmetric heteroskedastic noise
# =============================
gen_eps <- function(n) rexp(n, rate = 1) - 1

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  sigma_x <- 1 + abs(x)
  data.frame(x = x, y = f(x) + sigma_x * gen_eps(n))
}

# =============================
# CQR function
# =============================
cqr_rf <- function(train, cal, test, alpha) {
  
  rf <- ranger(
    y ~ x,
    data = train,
    quantreg = TRUE,
    num.trees = 500,
    min.node.size = 30,
    sample.fraction = 0.8
  )
  
  pred_cal <- predict(
    rf, cal,
    type = "quantiles",
    quantiles = c(alpha/2, 1 - alpha/2)
  )$predictions
  
  E <- pmax(
    pred_cal[,1] - cal$y,
    cal$y - pred_cal[,2]
  )
  
  k <- ceiling((length(E) + 1)*(1 - alpha))
  q_hat <- sort(E)[k]
  
  pred_test <- predict(
    rf, test,
    type = "quantiles",
    quantiles = c(alpha/2, 1 - alpha/2)
  )$predictions
  
  l <- pred_test[,1] - q_hat
  u <- pred_test[,2] + q_hat
  
  list(l = l, u = u)
}

# =============================
# Storage
# =============================
results_all_4b <- data.frame()

start_time <- Sys.time()
counter <- 1
total_iters <- length(n_list) * M

# =============================
# Simulation
# =============================
for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    # Progress
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    progress <- counter / total_iters
    
    if (counter > 1) {
      eta <- elapsed / (counter - 1) * (total_iters - counter + 1)
    } else {
      eta <- NA
    }
    
    if (counter %% 5 == 0) {
      cat(sprintf(
        "Iter %d/%d (%.1f%%) | n=%d, m=%d | ETA=%.1fs\n",
        counter, total_iters, 100*progress, n, m, eta
      ))
    }
    
    # Data
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
    # RF (mean)
    # =========================
    rf_fit <- ranger(
      y ~ x,
      data = train,
      num.trees = 500,
      min.node.size = 30,
      sample.fraction = 0.8
    )
    
    rf_pred_cal  <- predict(rf_fit, cal)$predictions
    rf_pred_test <- predict(rf_fit, test)$predictions
    
    # =========================
    # Oracle (correct!)
    # =========================
    q_low  <- qexp(alpha/2, rate = 1) - 1
    q_high <- qexp(1 - alpha/2, rate = 1) - 1
    
    sigma_x <- 1 + abs(test$x)
    
    oracle_l <- f(test$x) + sigma_x * q_low
    oracle_u <- f(test$x) + sigma_x * q_high
    
    # =========================
    # RQ
    # =========================
    rq_l_lm <- lm_pred_test + quantile(cal$y - lm_pred_cal, alpha/2)
    rq_u_lm <- lm_pred_test + quantile(cal$y - lm_pred_cal, 1 - alpha/2)
    
    rq_l_rf <- rf_pred_test + quantile(cal$y - rf_pred_cal, alpha/2)
    rq_u_rf <- rf_pred_test + quantile(cal$y - rf_pred_cal, 1 - alpha/2)
    
    # =========================
    # QR (RF)
    # =========================
    qr_fit <- ranger(
      y ~ x,
      data = train,
      quantreg = TRUE,
      num.trees = 500,
      min.node.size = 30,
      sample.fraction = 0.8
    )
    
    qr_pred <- predict(
      qr_fit, test,
      type = "quantiles",
      quantiles = c(alpha/2, 1 - alpha/2)
    )$predictions
    
    qr_l <- qr_pred[,1]
    qr_u <- qr_pred[,2]
    
    # =========================
    # CP (LM + RF)
    # =========================
    k <- ceiling((n_cal + 1)*(1 - alpha))
    
    q_lm <- sort(abs(cal$y - lm_pred_cal))[k]
    cp_l_lm <- lm_pred_test - q_lm
    cp_u_lm <- lm_pred_test + q_lm
    
    q_rf <- sort(abs(cal$y - rf_pred_cal))[k]
    cp_l_rf <- rf_pred_test - q_rf
    cp_u_rf <- rf_pred_test + q_rf
    
    # =========================
    # CQR
    # =========================
    cqr <- cqr_rf(train, cal, test, alpha)
    
    # =========================
    # Metrics
    # =========================
    cov <- function(l,u) mean(test$y >= l & test$y <= u)
    wid <- function(l,u) mean(u-l)
    
    results_all_4b <- rbind(results_all_4b, data.frame(
      n = n, m = m,
      method = c("Oracle",
                 "RQ_LM","RQ_RF",
                 "QR_RF",
                 "CP_LM","CP_RF",
                 "CQR_RF"),
      coverage = c(
        cov(oracle_l,oracle_u),
        cov(rq_l_lm,rq_u_lm),
        cov(rq_l_rf,rq_u_rf),
        cov(qr_l,qr_u),
        cov(cp_l_lm,cp_u_lm),
        cov(cp_l_rf,cp_u_rf),
        cov(cqr$l,cqr$u)
      ),
      length = c(
        wid(oracle_l,oracle_u),
        wid(rq_l_lm,rq_u_lm),
        wid(rq_l_rf,rq_u_rf),
        wid(qr_l,qr_u),
        wid(cp_l_lm,cp_u_lm),
        wid(cp_l_rf,cp_u_rf),
        wid(cqr$l,cqr$u)
      )
    ))
    
    counter <- counter + 1
  }
}