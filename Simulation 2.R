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
error_dist <- "exp"

gen_eps <- function(n) {
  switch(error_dist,
         "normal" = rnorm(n, 0, sigma),
         "exp"    = rexp(n, rate = 1/sigma) - sigma
  )
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
      data.frame(x = x, y = f(x) + gen_eps(n))
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    lm_fit     <- lm(y ~ x, data = train)
    pred_train <- predict(lm_fit, train)
    pred_cal   <- predict(lm_fit, cal)
    pred_test  <- predict(lm_fit, test)
    
    # Oracle
    q_oracle_low  <- qexp(alpha/2,     rate = 1/sigma) - sigma
    q_oracle_high <- qexp(1 - alpha/2, rate = 1/sigma) - sigma
    oracle_lower  <- f(test$x) + q_oracle_low
    oracle_upper  <- f(test$x) + q_oracle_high
    
    # Gaussian
    sigma_hat  <- sqrt(sum((train$y - pred_train)^2) / (n_train - 2))
    XtX_inv    <- solve(crossprod(model.matrix(lm_fit)))
    x_test_mat <- model.matrix(~x, data = test)
    se_pred    <- sqrt(sigma_hat^2 * (1 + rowSums((x_test_mat %*% XtX_inv) * x_test_mat)))
    t_val      <- qt(1 - alpha / 2, df = n_train - 2)
    gauss_lower <- pred_test - t_val * se_pred
    gauss_upper <- pred_test + t_val * se_pred
    
    # Residual quantile
    res_cal <- cal$y - pred_cal
    RQ_lower <- pred_test + quantile(res_cal, alpha / 2)
    RQ_upper <- pred_test + quantile(res_cal, 1 - alpha / 2)
    
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
      method = c("Oracle", "Gaussian", "RQ", "QR", "CP"),
      coverage = c(cov_fn(oracle_lower, oracle_upper),
                   cov_fn(gauss_lower, gauss_upper),
                   cov_fn(RQ_lower, RQ_upper),
                   cov_fn(qr_lower, qr_upper),
                   cov_fn(cp_lower, cp_upper)),
      length = c(len_fn(oracle_lower, oracle_upper),
                 len_fn(gauss_lower, gauss_upper),
                 len_fn(RQ_lower, RQ_upper),
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



# =============================
# non-scaled length
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
    mean_len    = mean(length),
    se_len      = sd(length) / sqrt(M),
    .groups = "drop"
  )

# =============================
# Plot: Length 
# =============================
ggplot(summary_df, aes(x = n, y = mean_len, color = method, group = method)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = mean_len - 2 * se_len,
                    ymax = mean_len + 2 * se_len), width = 0.05) +
  scale_x_log10(breaks = n_list) +
  labs(y = expression(L[n]), x = "Sample size n") +
  theme_bw()

results_all |> 
  filter(method == "Oracle") |> 
  group_by(n) |> 
  summarise(mean_length = mean(length))

method_colors <- c(
  "CP"       = "darkorange1",  # orange-rød
  "Gaussian" = "green3",  # grøn
  "RQ"       = "red",
  "QR"       = "blue"
)


# scaled difference 
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


# Plot
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













# Generer ét datasæt for et givet n
n <- 80
n_train <- n/2
n_cal   <- n/2

train <- gen_data(n_train)
cal   <- gen_data(n_cal)

lm_fit     <- lm(y ~ x, data = train)
pred_train <- predict(lm_fit, train)
pred_cal   <- predict(lm_fit, cal)

# Grid af x-værdier
x_grid     <- seq(-1, 1, length.out = 300)
x_grid_mat <- model.matrix(~ x, data = data.frame(x = x_grid))
pred_grid  <- x_grid_mat %*% coef(lm_fit)

# Gaussian
sigma_hat <- sqrt(sum((train$y - pred_train)^2) / (n_train - 2))
XtX_inv   <- solve(crossprod(model.matrix(lm_fit)))
se_grid   <- sqrt(sigma_hat^2 * (1 + rowSums((x_grid_mat %*% XtX_inv) * x_grid_mat)))
t_val     <- qt(1 - alpha/2, df = n_train - 2)

# CP
q_cp <- sort(abs(cal$y - pred_cal))[ceiling((n_cal + 1) * (1 - alpha))]

# Oracle
q_low  <- qexp(alpha/2,     rate = 1/sigma) - sigma
q_high <- qexp(1 - alpha/2, rate = 1/sigma) - sigma

# Saml i data frame
band_df <- data.frame(
  x       = x_grid,
  pred    = as.vector(pred_grid),
  # Oracle
  oracle_low  = f(x_grid) + q_low,
  oracle_high = f(x_grid) + q_high,
  # Gaussian
  gauss_low   = as.vector(pred_grid) - t_val * se_grid,
  gauss_high  = as.vector(pred_grid) + t_val * se_grid,
  # CP
  cp_low      = as.vector(pred_grid) - q_cp,
  cp_high     = as.vector(pred_grid) + q_cp
)

# Plot - ét panel per metode
band_long <- tidyr::pivot_longer(band_df, 
                                 cols = c(oracle_low, oracle_high, gauss_low, gauss_high, cp_low, cp_high),
                                 names_to = c("method", ".value"),
                                 names_pattern = "(.*)_(low|high)"
) |> rename(low = low, high = high)

ggplot(band_long, aes(x = x)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = method), alpha = 0.3) +
  geom_line(aes(y = pred), color = "black", linetype = "dashed") +
  geom_point(data = train, aes(x = x, y = y), size = 0.8, alpha = 0.5) +
  facet_wrap(~ method) +
  labs(x = expression(chi), y = expression(gamma)) +
  theme_bw() +
  theme(legend.position = "none")

# Generer testdata med eksponentiel støj
set.seed(42)
n_show <- 1000
x_show <- runif(n_show, -1, 1)
eps_show <- rexp(n_show, rate = 1/sigma) - sigma
y_show <- f(x_show) + eps_show

# Træn på separat data
train <- gen_data(50)
cal   <- gen_data(50)

lm_fit     <- lm(y ~ x, data = train)
pred_train <- predict(lm_fit, train)
pred_cal   <- predict(lm_fit, cal)

# Grid
x_grid     <- seq(-1, 1, length.out = 300)
x_grid_mat <- model.matrix(~ x, data = data.frame(x = x_grid))
pred_grid  <- as.vector(x_grid_mat %*% coef(lm_fit))

# Oracle bånd (asymmetrisk)
q_low  <- qexp(alpha/2,     rate = 1/sigma) - sigma
q_high <- qexp(1 - alpha/2, rate = 1/sigma) - sigma

# Gaussian bånd
sigma_hat <- sqrt(sum((train$y - pred_train)^2) / (n_train - 2))
XtX_inv   <- solve(crossprod(model.matrix(lm_fit)))
se_grid   <- sqrt(sigma_hat^2 * (1 + rowSums((x_grid_mat %*% XtX_inv) * x_grid_mat)))
t_val     <- qt(1 - alpha/2, df = n_train - 2)

# CP bånd
q_cp <- sort(abs(cal$y - pred_cal))[ceiling((n_cal + 1) * (1 - alpha))]

band_df <- data.frame(
  x            = x_grid,
  pred         = pred_grid,
  oracle_low   = f(x_grid) + q_low,
  oracle_high  = f(x_grid) + q_high,
  gauss_low    = pred_grid - t_val * se_grid,
  gauss_high   = pred_grid + t_val * se_grid,
  cp_low       = pred_grid - q_cp,
  cp_high      = pred_grid + q_cp
)

point_df <- data.frame(x = x_show, y = y_show,
                       covered_oracle  = y_show >= f(x_show) + q_low  & y_show <= f(x_show) + q_high,
                       covered_gauss   = y_show >= pred_grid[findInterval(x_show, x_grid)] - t_val &
                         y_show <= pred_grid[findInterval(x_show, x_grid)] + t_val,
                       covered_cp      = y_show >= f(x_show) - q_cp & y_show <= f(x_show) + q_cp
)

# Plot
p_oracle <- ggplot(band_df, aes(x = x)) +
  geom_ribbon(aes(ymin = oracle_low, ymax = oracle_high), fill = "turquoise", alpha = 0.3) +
  geom_line(aes(y = f(x)), color = "black", linetype = "dashed") +
  geom_point(data = data.frame(x = x_show, y = y_show),
             aes(x = x, y = y), size = 0.8, alpha = 0.5) +
  labs(title = "Oracle", x = expression(x), y = expression(y)) +
  theme_bw()

p_gauss <- ggplot(band_df, aes(x = x)) +
  geom_ribbon(aes(ymin = gauss_low, ymax = gauss_high), fill = "green4", alpha = 0.3) +
  geom_line(aes(y = pred), color = "black", linetype = "dashed") +
  geom_point(data = data.frame(x = x_show, y = y_show),
             aes(x = x, y = y), size = 0.8, alpha = 0.5) +
  labs(title = "Gaussian", x = expression(x), y = expression(y)) +
  theme_bw()

p_cp <- ggplot(band_df, aes(x = x)) +
  geom_ribbon(aes(ymin = cp_low, ymax = cp_high), fill = "red", alpha = 0.3) +
  geom_line(aes(y = pred), color = "black", linetype = "dashed") +
  geom_point(data = data.frame(x = x_show, y = y_show),
             aes(x = x, y = y), size = 0.8, alpha = 0.5) +
  labs(title = "CP", x = expression(x), y = expression(y)) +
  theme_bw()

library(patchwork)
p_oracle + p_gauss + p_cp

# RQ kvantiler
res_cal   <- cal$y - pred_cal
q_rq_low  <- quantile(res_cal, alpha/2)
q_rq_high <- quantile(res_cal, 1 - alpha/2)

# QR kvantiler
X_cal <- model.matrix(~ x, data = cal)

beta_low  <- fit_qr(X_cal, cal$y, alpha/2)
beta_high <- fit_qr(X_cal, cal$y, 1 - alpha/2)

# Tilføj til band_df
band_df <- band_df |>
  mutate(
    rq_low  = pred + q_rq_low,
    rq_high = pred + q_rq_high,
    qr_low  = as.vector(x_grid_mat %*% beta_low),
    qr_high = as.vector(x_grid_mat %*% beta_high)
  )

# Plot funktion så vi undgår gentagelse
plot_band <- function(low_col, high_col, title, fill_col) {
  ggplot(band_df, aes(x = x)) +
    geom_ribbon(aes(ymin = .data[[low_col]], ymax = .data[[high_col]]),
                fill = fill_col, alpha = 0.3) +
    geom_line(aes(y = pred), color = "black", linetype = "dashed") +
    geom_point(data = data.frame(x = x_show, y = y_show),
               aes(x = x, y = y), size = 0.8, alpha = 0.5) +
    labs(title = title, x = expression(x), y = expression(y)) +
    theme_bw()
}

p_oracle <- plot_band("oracle_low", "oracle_high", "Oracle",   "turquoise")
p_gauss  <- plot_band("gauss_low",  "gauss_high",  "Gaussian", "green4")
p_cp     <- plot_band("cp_low",     "cp_high",     "CP",       "red")
p_rq     <- plot_band("rq_low",     "rq_high",     "RQ",       "purple")
p_qr     <- plot_band("qr_low",     "qr_high",     "QR",       "blue")

library(patchwork)
p_oracle + p_gauss + p_cp + p_rq + p_qr
