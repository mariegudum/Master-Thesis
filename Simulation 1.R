# =============================
# Setup
# =============================
set.seed(1)
library(ggplot2)
library(dplyr)

# Parameters
n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1
beta   <- 2; a <- 1; sigma <- 1
f <- function(x) beta * x + a

# Quantile regression loss function
rho <- function(u, tau) u * (tau - (u < 0))

fit_qr <- function(X, y, tau) {
  p <- ncol(X)
  optim(
    par    = rep(0, p),
    fn     = function(beta) sum(rho(y - X %*% beta, tau)),
    method = "BFGS"
  )$par
}


# =============================
# Simulation
# =============================
results_all <- data.frame()

for (n in n_list) {
  cat("Running n =", n, "\n")
  n_train <- n / 2
  n_cal   <- n / 2
  
  for (m in 1:M) {
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    lm_fit     <- lm(y ~ x, data = train)
    pred_train <- predict(lm_fit, train)
    pred_cal   <- predict(lm_fit, cal)
    pred_test  <- predict(lm_fit, test)
    
    # Oracle
    z            <- qnorm(1 - alpha / 2)
    oracle_lower <- f(test$x) - z * sigma
    oracle_upper <- f(test$x) + z * sigma
    
    # Gaussian
    sigma_hat  <- sqrt(sum((train$y - pred_train)^2) / (n_train - 2))
    XtX_inv    <- solve(crossprod(model.matrix(lm_fit)))
    x_test_mat <- model.matrix(~x, data = test)
    se_pred    <- sqrt(sigma_hat^2 * (1 + rowSums((x_test_mat %*% XtX_inv) * x_test_mat)))
    t_val      <- qt(1 - alpha / 2, df = n_train - 2)
    gauss_lower <- pred_test - t_val * se_pred
    gauss_upper <- pred_test + t_val * se_pred
    
    # QR
    X_cal <- model.matrix(~ x, data = cal)
    X_test <- model.matrix(~ x, data = test)
    
    beta_low  <- fit_qr(X_cal, cal$y, alpha / 2)
    beta_high <- fit_qr(X_cal, cal$y, 1 - alpha / 2)
    
    qr_lower <- X_test %*% beta_low
    qr_upper <- X_test %*% beta_high
    
    # Conformal
    k      <- ceiling((n_cal + 1) * (1 - alpha))
    q_cp   <- sort(abs(cal$y - pred_cal))[k]
    cp_lower <- pred_test - q_cp
    cp_upper <- pred_test + q_cp
    
    # Metrics
    cov_fn <- function(l, u) mean(test$y >= l & test$y <= u)
    len_fn <- function(l, u) mean(u - l)
    
    results_all <- rbind(results_all, data.frame(
      n      = n,
      m = m,
      method = c("Oracle", "Gaussian", "QR", "CP"),
      coverage = c(cov_fn(oracle_lower, oracle_upper),
                   cov_fn(gauss_lower, gauss_upper),
                   cov_fn(qr_lower, qr_upper),
                   cov_fn(cp_lower, cp_upper)),
      length = c(len_fn(oracle_lower, oracle_upper),
                 len_fn(gauss_lower, gauss_upper),
                 len_fn(qr_lower, qr_upper),
                 len_fn(cp_lower, cp_upper))
    ))
  }
}

# =============================
# Post-processing
# =============================
summary_df <- results_all |>
  mutate(length_scaled = ifelse(
    method == "Oracle",
    length,
    sqrt(n) * length
  )) |>
  group_by(n, method) |>
  summarise(
    mean_cov    = mean(coverage),
    se_cov      = sd(coverage) / sqrt(M),
    mean_len    = mean(length_scaled),
    se_len      = sd(length_scaled) / sqrt(M),
    .groups = "drop"
  )

# =============================
# Plot: Coverage
# =============================
ggplot(summary_df, aes(x = n, y = mean_cov, color = method, group = method)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = mean_cov - 2 * se_cov,
                    ymax = mean_cov + 2 * se_cov), width = 0.05) +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed", color = "black") +
  scale_x_log10() +
  labs(y = "Coverage", x = "Sample size n") +
  theme_bw()


mean_cov <- results_all |> 
  filter(method != "Oracle") |>
  group_by(n, method) |>
  summarise(mean_cov = mean(coverage), .groups = "drop")

method_colors <- c("CP" = "#F8766D", "Gaussian" = "#00BA38", "QR" = "#619CFF")

method_colors <- c(
  "CP"       = "darkorange1",  # orange-rød
  "Gaussian" = "green3",  # grøn
  "QR"       = "blue"
)

ggplot(results_all |> filter(method != "Oracle"),
       aes(x = as.numeric(as.character(n)), y = coverage, color = method)) +
  geom_jitter(shape = 20, size = 1, width = 0.20) +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed", color = "black") +
  geom_text(data = mean_cov,
            aes(x = as.numeric(as.character(n)), y = 1.02,
                label = scales::percent(mean_cov, accuracy = 1)),
            size = 3) +
  coord_cartesian(ylim = c(NA, 1.05)) +
  scale_x_log10(breaks = n_list, labels = n_list) +
  facet_wrap(~ method) +
  labs(y = "Coverage", x = "Sample size n") +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  theme(legend.position = "none")


# =============================
# Plot: Length (scaled)
# =============================
results_all <- left_join(results_all,
                         results_all |> filter(method == "Oracle") |>
                           group_by(n) |> summarise(L_oracle = mean(length)),
                         by = "n")


summary_diff <- results_all |>
  filter(method != "Oracle") |>
  mutate(length_diff = sqrt(n) * (length - L_oracle)) |>
  group_by(n, method) |>
  summarise(
    mean_len = mean(length_diff),
    se_len   = sd(length_diff) / sqrt(M),
    .groups  = "drop"
  )


ggplot(summary_diff, aes(x = n, y = mean_len, color = method, group = method)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = mean_len - 2 * se_len,
                    ymax = mean_len + 2 * se_len), width = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10(breaks = n_list) +
  labs(y = expression(sqrt(n) * (L[n] - L[oracle])), x = "Sample size n") +
  scale_color_manual(values = method_colors) +
  theme_bw()


