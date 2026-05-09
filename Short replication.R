# ============================================================
# FAST Simulation (~10 minutes)
# ============================================================

library(MASS)
library(ranger)
library(ggplot2)

set.seed(1)

# ============================================================
# 1. Data generating process
# ============================================================

generate_X <- function(n, p, correlated = FALSE, rho = 0.6) {
  if (!correlated) {
    matrix(rnorm(n * p), n, p)
  } else {
    Sigma <- rho^abs(outer(1:p, 1:p, "-"))
    mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  }
}

m_function <- function(X, type) {
  x1 <- X[,1]; x2 <- X[,2]
  
  if (type == "linear") {
    x1 + x2
  } else if (type == "nonlinear") {
    2 * exp(-abs(x1) - abs(x2))
  } else {
    2 * exp(-abs(x1) - abs(x2)) + x1 * x2
  }
}

generate_epsilon <- function(n, type, m) {
  if (type == "homoscedastic") {
    rnorm(n)
  } else if (type == "heavy") {
    rt(n, df = 3) / sqrt(3)
  } else {
    sigma <- sqrt(0.5 + 0.5 * abs(m))
    rnorm(n, sd = sigma)
  }
}

generate_data <- function(n, p, m_type, eps_type, correlated) {
  X <- generate_X(n, p, correlated)
  m <- m_function(X, m_type)
  eps <- generate_epsilon(n, eps_type, m)
  Y <- m + eps
  list(X = X, Y = Y)
}

# ============================================================
# 2. Methods (FAST: no tuning)
# ============================================================

rf_fit <- function(X, Y) {
  ranger(
    Y ~ .,
    data = data.frame(Y = Y, X),
    num.trees = 500,
    mtry = floor(sqrt(ncol(X))),
    min.node.size = 5,
    keep.inbag = TRUE
  )
}

# ---------------- OOB ----------------
oob_intervals <- function(rf, X_train, Y_train, X_test, alpha = 0.1) {
  res <- Y_train - rf$predictions
  
  q_lo <- quantile(res, alpha/2)
  q_hi <- quantile(res, 1 - alpha/2)
  
  pred <- predict(rf, data.frame(X_test))$predictions
  
  list(lower = pred + q_lo, upper = pred + q_hi)
}

# ---------------- Split conformal ----------------
split_conformal <- function(X, Y, X_test, alpha = 0.1) {
  n <- nrow(X)
  idx <- sample(1:n, n/2)
  
  X_train <- X[idx,]; Y_train <- Y[idx]
  X_cal   <- X[-idx,]; Y_cal <- Y[-idx]
  
  rf <- ranger(Y ~ ., data = data.frame(Y = Y_train, X_train), num.trees = 500)
  
  pred_cal <- predict(rf, data.frame(X_cal))$predictions
  scores <- abs(Y_cal - pred_cal)
  
  q <- quantile(scores, 1 - alpha)
  
  pred_test <- predict(rf, data.frame(X_test))$predictions
  
  list(lower = pred_test - q, upper = pred_test + q)
}

# ---------------- QRF approximation ----------------
qrf_approx <- function(X, Y, X_test, alpha = 0.1) {
  
  rf <- ranger(
    Y ~ .,
    data = data.frame(Y = Y, X),
    num.trees = 500,
    write.forest = TRUE
  )
  
  preds_all <- predict(rf, data.frame(X_test), predict.all = TRUE)$predictions
  
  lower <- apply(preds_all, 1, quantile, probs = alpha/2)
  upper <- apply(preds_all, 1, quantile, probs = 1 - alpha/2)
  
  list(lower = lower, upper = upper)
}

# ============================================================
# 3. Evaluation
# ============================================================

evaluate <- function(Y, lower, upper) {
  coverage <- mean(Y >= lower & Y <= upper)
  width <- mean(upper - lower)
  c(coverage = coverage, width = width)
}

# ============================================================
# 4. Simulation (FAST VERSION)
# ============================================================

run_simulation <- function() {
  
  n_vals <- c(200, 500)
  m_types <- c("linear", "nonlinear", "interaction")
  eps_types <- c("homoscedastic", "heteroscedastic")
  corr_types <- c(FALSE)
  
  results <- data.frame()
  
  for (n in n_vals) {
    for (m_type in m_types) {
      for (eps_type in eps_types) {
        for (corr in corr_types) {
          
          cat("Running:", n, m_type, eps_type, "\n")
          
          for (rep in 1:10) {
            
            train <- generate_data(n, 10, m_type, eps_type, corr)
            test  <- generate_data(1000, 10, m_type, eps_type, corr)
            
            rf <- rf_fit(train$X, train$Y)
            
            oob <- oob_intervals(rf, train$X, train$Y, test$X)
            sc  <- split_conformal(train$X, train$Y, test$X)
            qrf <- qrf_approx(train$X, train$Y, test$X)
            
            res_oob <- evaluate(test$Y, oob$lower, oob$upper)
            res_sc  <- evaluate(test$Y, sc$lower, sc$upper)
            res_qrf <- evaluate(test$Y, qrf$lower, qrf$upper)
            
            results <- rbind(results,
                             data.frame(
                               n = n,
                               m = m_type,
                               eps = eps_type,
                               method = c("OOB","SC","QRF_approx"),
                               coverage = c(res_oob[1], res_sc[1], res_qrf[1]),
                               width = c(res_oob[2], res_sc[2], res_qrf[2])
                             )
            )
          }
        }
      }
    }
  }
  
  results
}

# ============================================================
# 5. Run simulation
# ============================================================

results <- run_simulation()

# ============================================================
# 6. Plot: Coverage (Figure 1 style)
# ============================================================

ggplot(results, aes(x = method, y = coverage)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  facet_grid(m ~ eps) +
  theme_bw()

# ============================================================
# 7. Plot: Width ratios (Figure 2 style)
# ============================================================

df_ratio <- do.call(rbind, lapply(split(results, interaction(results$n, results$m, results$eps)), function(df) {
  oob_width <- df$width[df$method == "OOB"]
  df$ratio <- NA
  df$ratio[df$method != "OOB"] <- log(df$width[df$method != "OOB"] / oob_width)
  df
}))

ggplot(na.omit(df_ratio), aes(x = method, y = ratio)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(m ~ eps) +
  theme_bw()
