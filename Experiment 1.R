# =============================
# Setup
# =============================
set.seed(1)
library(ggplot2)
library(dplyr)
library(ranger)

n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta * x + a

# Quantile regression (linear)
rho <- function(u, tau) u * (tau - (u < 0))
fit_qr <- function(X, y, tau) {
  optim(rep(0, ncol(X)),
        function(b) sum(rho(y - X %*% b, tau)),
        method = "BFGS")$par
}

results_all_1 <- data.frame()

# =============================
# Simulation
# =============================
for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
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
    # RF
    # =========================
    rf_fit <- ranger(y ~ x, data = train)
    
    rf_pred_train <- predict(rf_fit, train)$predictions
    rf_pred_cal   <- predict(rf_fit, cal)$predictions
    rf_pred_test  <- predict(rf_fit, test)$predictions
    
    # =========================
    # Oracle
    # =========================
    z <- qnorm(1 - alpha/2)
    oracle_l <- f(test$x) - z*sigma
    oracle_u <- f(test$x) + z*sigma
    
    # =========================
    # Gaussian (LM only)
    # =========================
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
    X_cal  <- model.matrix(~x, cal)
    X_test <- model.matrix(~x, test)
    
    b_l <- fit_qr(X_cal, cal$y, alpha/2)
    b_u <- fit_qr(X_cal, cal$y, 1 - alpha/2)
    
    qr_l_lm <- X_test %*% b_l
    qr_u_lm <- X_test %*% b_u
    
    qr_rf <- ranger(y ~ x, data = cal, quantreg = TRUE)
    qr_pred <- predict(qr_rf, test,
                       type = "quantiles",
                       quantiles = c(alpha/2, 1-alpha/2))$predictions
    
    qr_l_rf <- qr_pred[,1]
    qr_u_rf <- qr_pred[,2]
    
    # =========================
    # CP
    # =========================
    k <- ceiling((n_cal + 1)*(1-alpha))
    
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
    
    results_all_1 <- rbind(results_all_1, data.frame(
      n=n, m=m,
      method=c("Oracle",
               "Gaussian_LM",
               "RQ_LM","RQ_RF",
               "QR_LM","QR_RF",
               "CP_LM","CP_RF"),
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
    ))
  }
}

saveRDS(results_all_1, "partial_results.rds")

rel_results_1 <- results_all_1 %>%
  group_by(n, m) %>%
  mutate(
    oracle_length = length[method == "Oracle"],
    rel_length = length / oracle_length
  ) %>%
  ungroup()

rel_results_1 %>%
  group_by(n, method) %>%
  summarise(mean_rel = mean(rel_length), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_rel, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "Relative interval length",
       y = "Length / Oracle length")



results_all_1 %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "Coverage vs sample size",
       y = "Coverage")

ggplot(results_all_1, aes(x = coverage, fill = method)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 1 - alpha, linetype = "dashed") +
  labs(title = "Coverage distribution")



# without oracle


results_no_oracle_1 <- results_all_1 %>%
  filter(method != "Oracle")


results_no_oracle_1 %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  )

results_no_oracle %>%
  group_by(n, method) %>%
  summarise(mean_len = mean(length), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_len, color = method)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(
    title = "Interval length vs sample size",
    y = "Length" 
  )

oracle_len <- results_all_1 %>%
  filter(method == "Oracle") %>%
  summarise(mean_len = mean(length)) %>%
  pull(mean_len)


# final
rel_results_1 %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_rel = mean(rel_length), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_rel, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Relative interval length vs sample size",
    y = "Length / Oracle length"
  )

results_all_1 %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  )

ggplot(
  results_all_1 %>% filter(method != "Oracle"),
  aes(x = coverage, fill = method)
) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 1 - alpha, linetype = "dashed") +
  labs(title = "Coverage distribution")




# residualer
rf_residuals <- cal$y - rf_pred_cal
lm_residuals <- cal$y - lm_pred_cal

# samlet data
df_res <- data.frame(
  x = cal$x,
  rf_res = rf_residuals,
  lm_res = lm_residuals
)

# =========================
# Plot RF residualer
# =========================
p_rf <- ggplot(df_res, aes(x = x, y = rf_res)) +
  geom_point(alpha = 0.3) +
  geom_smooth(se = FALSE, color = "red") +
  labs(
    title = "RF residuals vs X",
    x = "X",
    y = "Residual"
  ) +
  theme_minimal()

# =========================
# Plot LM residualer
# =========================
p_lm <- ggplot(df_res, aes(x = x, y = lm_res)) +
  geom_point(alpha = 0.3) +
  geom_smooth(se = FALSE, color = "blue") +
  labs(
    title = "LM residuals vs X",
    x = "X",
    y = "Residual"
  ) +
  theme_minimal()

print(p_rf)
print(p_lm)





# Tuning of QR with RF


library(ranger)

set.seed(1)

# =========================
# Data (samme som dit setup)
# =========================
n <- 500
n_train <- 250
n_cal   <- 250
n_test  <- 20000

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta*x + a

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
}

train <- gen_data(n_train)
cal   <- gen_data(n_cal)
test  <- gen_data(n_test)

alpha <- 0.1

# =========================
# 1. Default QR_RF
# =========================
qr_rf_default <- ranger(
  y ~ x,
  data = cal,
  quantreg = TRUE,
  num.trees = 300
)

pred_default <- predict(
  qr_rf_default, test,
  type = "quantiles",
  quantiles = c(alpha/2, 1-alpha/2)
)$predictions

l_def <- pred_default[,1]
u_def <- pred_default[,2]

# =========================
# 2. Tuned QR_RF
# =========================
qr_rf_tuned <- ranger(
  y ~ x,
  data = cal,
  quantreg = TRUE,
  num.trees = 800,
  min.node.size = 50,
  sample.fraction = 0.8
)

pred_tuned <- predict(
  qr_rf_tuned, test,
  type = "quantiles",
  quantiles = c(alpha/2, 1-alpha/2)
)$predictions

l_tuned <- pred_tuned[,1]
u_tuned <- pred_tuned[,2]

# =========================
# 3. Evaluering
# =========================
cov <- function(l,u) mean(test$y >= l & test$y <= u)
wid <- function(l,u) mean(u - l)

results <- data.frame(
  method = c("QR_RF_default", "QR_RF_tuned"),
  coverage = c(
    cov(l_def, u_def),
    cov(l_tuned, u_tuned)
  ),
  length = c(
    wid(l_def, u_def),
    wid(l_tuned, u_tuned)
  )
)

print(results)



library(ranger)
library(dplyr)

set.seed(1)

# =========================
# Setup (optimeret)
# =========================
n_list <- c(20, 80, 320, 1280)
M      <- 30          # ↓ fra 200
n_test <- 10000       # ↓ fra 100000
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta*x + a

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
}

# =========================
# Pre-allocation (MEGET vigtigt)
# =========================
results <- vector("list", length(n_list) * M)
counter <- 1

# =========================
# Loop
# =========================
for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =========================
    # Default QR_RF
    # =========================
    rf_def <- ranger(
      y ~ x,
      data = cal,
      quantreg = TRUE,
      num.trees = 300,
      num.threads = parallel::detectCores()
    )
    
    pred_def <- predict(
      rf_def, test,
      type = "quantiles",
      quantiles = c(alpha/2, 1-alpha/2)
    )$predictions
    
    # =========================
    # Tuned QR_RF (optimeret)
    # =========================
    rf_tuned <- ranger(
      y ~ x,
      data = cal,
      quantreg = TRUE,
      num.trees = 500,          # ↓ fra 800
      min.node.size = 30,
      sample.fraction = 0.8,
      num.threads = parallel::detectCores()
    )
    
    pred_tuned <- predict(
      rf_tuned, test,
      type = "quantiles",
      quantiles = c(alpha/2, 1-alpha/2)
    )$predictions
    
    # =========================
    # Metrics (vectorized)
    # =========================
    cov_def   <- mean(test$y >= pred_def[,1] & test$y <= pred_def[,2])
    cov_tuned <- mean(test$y >= pred_tuned[,1] & test$y <= pred_tuned[,2])
    
    wid_def   <- mean(pred_def[,2] - pred_def[,1])
    wid_tuned <- mean(pred_tuned[,2] - pred_tuned[,1])
    
    results[[counter]] <- data.frame(
      n = n,
      m = m,
      method = c("QR_RF_default", "QR_RF_tuned"),
      coverage = c(cov_def, cov_tuned),
      length   = c(wid_def, wid_tuned)
    )
    
    counter <- counter + 1
  }
}

# =========================
# Combine (én gang!)
# =========================
qr_results <- bind_rows(results)

names(results_all_1)
names(qr_results)

results_extended <- bind_rows(results_all_1, qr_results)

results_extended <- results_extended %>%
  filter(method != "QR_RF")

results_extended %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  theme_minimal()

rel_results_ext <- results_extended %>%
  group_by(n, m) %>%
  mutate(
    oracle_length = length[method == "Oracle"][1],
    rel_length = length / oracle_length
  ) %>%
  ungroup()

rel_results_ext %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_rel = mean(rel_length), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_rel, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Relative interval length",
    y = "Length / Oracle length"
  ) +
  theme_minimal()

# 100 times 

library(ranger)
library(dplyr)

set.seed(1)

n_list <- c(20, 80, 320, 1280)
M      <- 100
n_test <- 20000
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta*x + a

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
}

results_list <- list()
counter <- 1

for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =========================
    # Default QR_RF
    # =========================
    rf_def <- ranger(y ~ x, data = cal, quantreg = TRUE, num.trees = 300)
    
    pred_def <- predict(rf_def, test,
                        type = "quantiles",
                        quantiles = c(alpha/2, 1-alpha/2))$predictions
    
    # =========================
    # Tuned QR_RF
    # =========================
    rf_tuned <- ranger(
      y ~ x,
      data = cal,
      quantreg = TRUE,
      num.trees = 800,
      min.node.size = 50,
      sample.fraction = 0.8
    )
    
    pred_tuned <- predict(rf_tuned, test,
                          type = "quantiles",
                          quantiles = c(alpha/2, 1-alpha/2))$predictions
    
    # metrics
    cov <- function(l,u) mean(test$y >= l & test$y <= u)
    wid <- function(l,u) mean(u - l)
    
    results_list[[counter]] <- data.frame(
      n = n,
      method = c("QR_RF_default", "QR_RF_tuned"),
      coverage = c(
        cov(pred_def[,1], pred_def[,2]),
        cov(pred_tuned[,1], pred_tuned[,2])
      ),
      length = c(
        wid(pred_def[,1], pred_def[,2]),
        wid(pred_tuned[,1], pred_tuned[,2])
      )
    )
    
    counter <- counter + 1
  }
}

qr_results <- bind_rows(results_list)

qr_results$method <- as.character(qr_results$method)

results_extended_1 <- bind_rows(
  results_all_1,
  qr_results
)

results_extended_1 <- results_extended_1 %>%
  filter(method != "QR_RF")


# samlet plot

results_extended_1 %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  theme_minimal()


rel_results_ext <- results_extended_1 %>%
  group_by(n, m) %>%
  mutate(
    oracle_length = length[method == "Oracle"],
    rel_length = length / oracle_length
  ) %>%
  ungroup()

rel_results_ext %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_rel = mean(rel_length), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_rel, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Relative interval length",
    y = "Length / Oracle length"
  ) +
  theme_minimal()












# tuned CP



library(ranger)
library(dplyr)

set.seed(1)

n_list <- c(20, 80, 320, 1280)
M      <- 100        # mindre for speed
n_test <- 20000      # mindre for speed
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta*x + a

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
}

results_list <- list()
counter <- 1

for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =========================
    # CP_RF default
    # =========================
    rf_def <- ranger(y ~ x, data = train, num.trees = 300)
    
    pred_cal_def  <- predict(rf_def, cal)$predictions
    pred_test_def <- predict(rf_def, test)$predictions
    
    scores_def <- abs(cal$y - pred_cal_def)
    k <- ceiling((length(scores_def) + 1)*(1-alpha))
    q_def <- sort(scores_def)[k]
    
    l_def <- pred_test_def - q_def
    u_def <- pred_test_def + q_def
    
    # =========================
    # CP_RF tuned
    # =========================
    rf_tuned <- ranger(
      y ~ x,
      data = train,
      num.trees = 800,
      min.node.size = 30,
      sample.fraction = 0.8
    )
    
    pred_cal_tuned  <- predict(rf_tuned, cal)$predictions
    pred_test_tuned <- predict(rf_tuned, test)$predictions
    
    scores_tuned <- abs(cal$y - pred_cal_tuned)
    q_tuned <- sort(scores_tuned)[k]
    
    l_tuned <- pred_test_tuned - q_tuned
    u_tuned <- pred_test_tuned + q_tuned
    
    # metrics
    cov <- function(l,u) mean(test$y >= l & test$y <= u)
    wid <- function(l,u) mean(u - l)
    
    results_list[[counter]] <- data.frame(
      n = n,
      method = c("CP_RF_default", "CP_RF_tuned"),
      coverage = c(
        cov(l_def, u_def),
        cov(l_tuned, u_tuned)
      ),
      length = c(
        wid(l_def, u_def),
        wid(l_tuned, u_tuned)
      )
    )
    
    counter <- counter + 1
  }
}




# =========================
# 1. Saml alle resultater
# =========================
results_extended <- bind_rows(
  results_all_1,   # dine gamle (med Oracle)
  cp_rf_extra,     # CP tuned
  qr_results       # QR tuned
)

# (valgfrit: fjern gammel QR_RF hvis den findes)
results_extended <- results_extended %>%
  filter(method != "RQ_RF_default")

# =========================
# 2. COVERAGE PLOT
# =========================
coverage_plot <- results_extended %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(
    mean_cov = mean(coverage),
    sd_cov   = sd(coverage),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  theme_minimal()

print(coverage_plot)

# =========================
# 3. ORACLE (mean per n)
# =========================
oracle_mean <- results_extended %>%
  filter(method == "Oracle") %>%
  group_by(n) %>%
  summarise(
    oracle_length = mean(length),
    .groups = "drop"
  )

# =========================
# 4. RELATIVE LENGTH
# =========================
rel_results <- results_extended %>%
  left_join(oracle_mean, by = "n") %>%
  mutate(rel_length = length / oracle_length)

# =========================
# 5. WIDTH PLOT
# =========================
width_plot <- rel_results %>%
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
    title = "Relative interval length",
    y = "Length / Oracle length"
  ) +
  theme_minimal()

print(width_plot)




# tuned rq

library(ranger)
library(dplyr)

set.seed(1)

n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta*x + a

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
}

results <- vector("list", length(n_list)*M)
counter <- 1

for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =========================
    # Default RF
    # =========================
    rf_def <- ranger(
      y ~ x,
      data = train,
      num.trees = 300,
      num.threads = parallel::detectCores()
    )
    
    pred_cal_def  <- predict(rf_def, cal)$predictions
    pred_test_def <- predict(rf_def, test)$predictions
    
    # =========================
    # Tuned RF
    # =========================
    rf_tuned <- ranger(
      y ~ x,
      data = train,
      num.trees = 500,
      min.node.size = 30,
      sample.fraction = 0.8,
      num.threads = parallel::detectCores()
    )
    
    pred_cal_tuned  <- predict(rf_tuned, cal)$predictions
    pred_test_tuned <- predict(rf_tuned, test)$predictions
    
    # =========================
    # RQ intervals
    # =========================
    q_low_def  <- quantile(cal$y - pred_cal_def, alpha/2)
    q_high_def <- quantile(cal$y - pred_cal_def, 1 - alpha/2)
    
    q_low_tuned  <- quantile(cal$y - pred_cal_tuned, alpha/2)
    q_high_tuned <- quantile(cal$y - pred_cal_tuned, 1 - alpha/2)
    
    l_def <- pred_test_def + q_low_def
    u_def <- pred_test_def + q_high_def
    
    l_tuned <- pred_test_tuned + q_low_tuned
    u_tuned <- pred_test_tuned + q_high_tuned
    
    # =========================
    # Metrics
    # =========================
    cov_def   <- mean(test$y >= l_def & test$y <= u_def)
    cov_tuned <- mean(test$y >= l_tuned & test$y <= u_tuned)
    
    wid_def   <- mean(u_def - l_def)
    wid_tuned <- mean(u_tuned - l_tuned)
    
    results[[counter]] <- data.frame(
      n = n,
      m = m,
      method = c("RQ_RF_default", "RQ_RF_tuned"),
      coverage = c(cov_def, cov_tuned),
      length   = c(wid_def, wid_tuned)
    )
    
    counter <- counter + 1
  }
}

rq_results <- bind_rows(results)


results_extended <- bind_rows(
  results_extended,   # fra før (CP + QR osv.)
  rq_results
)


coverage_plot <- results_extended %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(
    mean_cov = mean(coverage),
    sd_cov   = sd(coverage),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  theme_minimal()

print(coverage_plot)

# =========================
# 3. ORACLE (mean per n)
# =========================
oracle_mean <- results_extended %>%
  filter(method == "Oracle") %>%
  group_by(n) %>%
  summarise(
    oracle_length = mean(length),
    .groups = "drop"
  )

# =========================
# 4. RELATIVE LENGTH
# =========================
rel_results <- results_extended %>%
  left_join(oracle_mean, by = "n") %>%
  mutate(rel_length = length / oracle_length)

# =========================
# 5. WIDTH PLOT
# =========================
width_plot <- rel_results %>%
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
    title = "Relative interval length",
    y = "Length / Oracle length"
  ) +
  theme_minimal()

print(width_plot)





# samlet tunet QR og CP


library(ranger)
library(dplyr)

set.seed(1)

# =========================
# Setup (optimeret)
# =========================
n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta*x + a

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
}

# =========================
# Pre-allocation
# =========================
total_iters <- length(n_list) * M
results <- vector("list", total_iters)
counter <- 1

start_time <- Sys.time()

# =========================
# Loop
# =========================
for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    # =========================
    # Progress + ETA
    # =========================
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
    
    # =========================
    # Data
    # =========================
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =========================
    # QR_RF tuned
    # =========================
    rf_qr <- ranger(
      y ~ x,
      data = train,
      quantreg = TRUE,
      num.trees = 500,
      min.node.size = 30,
      sample.fraction = 0.8,
      num.threads = parallel::detectCores()
    )
    
    pred_qr <- predict(
      rf_qr, test,
      type = "quantiles",
      quantiles = c(alpha/2, 1-alpha/2)
    )$predictions
    
    l_qr <- pred_qr[,1]
    u_qr <- pred_qr[,2]
    
    # =========================
    # CP_RF tuned
    # =========================
    rf_cp <- ranger(
      y ~ x,
      data = train,
      num.trees = 500,
      min.node.size = 30,
      sample.fraction = 0.8,
      num.threads = parallel::detectCores()
    )
    
    pred_cal_cp  <- predict(rf_cp, cal)$predictions
    pred_test_cp <- predict(rf_cp, test)$predictions
    
    scores <- abs(cal$y - pred_cal_cp)
    k <- ceiling((length(scores) + 1)*(1 - alpha))
    q_hat <- sort(scores)[k]
    
    l_cp <- pred_test_cp - q_hat
    u_cp <- pred_test_cp + q_hat
    
    # =========================
    # Metrics
    # =========================
    cov_qr <- mean(test$y >= l_qr & test$y <= u_qr)
    cov_cp <- mean(test$y >= l_cp & test$y <= u_cp)
    
    wid_qr <- mean(u_qr - l_qr)
    wid_cp <- mean(u_cp - l_cp)
    
    results[[counter]] <- data.frame(
      n = n,
      m = m,
      method = c("QR_RF_tuned", "CP_RF_tuned"),
      coverage = c(cov_qr, cov_cp),
      length   = c(wid_qr, wid_cp)
    )
    
    counter <- counter + 1
  }
}

# =========================
# Combine
# =========================
tuned_results <- bind_rows(results)


# =========================
# Runtime
# =========================
total_time <- Sys.time() - start_time
cat("\nTotal runtime:", round(as.numeric(total_time, units="secs"),1), "seconds\n")


results_extended <- bind_rows(
  results_all_1,   # dine oprindelige (incl. Oracle)
  tuned_results,   # QR_RF_tuned + CP_RF_tuned
  rq_results       # RQ_RF_tuned
)

results_extended_new <- results_extended_new %>% 
  filter(method != "RQ_RF_default") 

unique(results_extended$method)

coverage_summary <- results_extended_new %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(
    mean_cov = mean(coverage),
    sd_cov   = sd(coverage),
    .groups = "drop"
  )
# with errorbars
coverage_plot <- results_extended_new %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(
    mean_cov = mean(coverage),
    sd_cov   = sd(coverage),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    data = coverage_summary %>% filter(method == "QR_RF_tuned"),
    aes(ymin = mean_cov - sd_cov, ymax = mean_cov + sd_cov),
    width = 0.05
  ) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  theme_minimal()


# without errorbars
coverage_plot <- results_extended_new %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(
    mean_cov = mean(coverage),
    sd_cov   = sd(coverage),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  theme_minimal()



print(coverage_plot)


oracle_mean_new <- results_extended_new %>%
  filter(method == "Oracle") %>%
  group_by(n) %>%
  summarise(
    oracle_length = mean(length),
    .groups = "drop"
  )

rel_results_new <- results_extended_new %>%
  left_join(oracle_mean_new, by = "n") %>%
  mutate(rel_length = length / oracle_length)

width_plot <- rel_results_new %>%
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
    title = "Relative interval length",
    y = "Length / Oracle length"
  ) +
  theme_minimal()

print(width_plot)




# test af fejl

library(ranger)

set.seed(1)

n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1

beta <- 2; a <- 1; sigma <- 1
f <- function(x) beta*x + a

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
}

# kandidater
grid <- c(10, 30, 80)

results <- list()
counter <- 1

for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    best_score <- Inf
    best_model <- NULL
    
    # ===== tuning =====
    for (node_size in grid) {
      
      rf_tmp <- ranger(
        y ~ x,
        data = train,
        quantreg = TRUE,
        num.trees = 300,
        min.node.size = node_size,
        num.threads = parallel::detectCores()
      )
      
      pred_cal <- predict(
        rf_tmp, cal,
        type = "quantiles",
        quantiles = c(alpha/2, 1-alpha/2)
      )$predictions
      
      l <- pred_cal[,1]
      u <- pred_cal[,2]
      
      # score: coverage deviation + længde
      cov <- mean(cal$y >= l & cal$y <= u)
      len <- mean(u - l)
      
      score <- abs(cov - (1-alpha)) + 0.01*len
      
      if (score < best_score) {
        best_score <- score
        best_model <- rf_tmp
      }
    }
    
    # ===== evaluering =====
    pred_test <- predict(
      best_model, test,
      type = "quantiles",
      quantiles = c(alpha/2, 1-alpha/2)
    )$predictions
    
    l <- pred_test[,1]
    u <- pred_test[,2]
    
    cov <- mean(test$y >= l & test$y <= u)
    len <- mean(u - l)
    
    results[[counter]] <- data.frame(
      n = n,
      m = m,
      method = "QR_RF_tuned_per_n",
      coverage = cov,
      length = len
    )
    
    counter <- counter + 1
    print(counter)
  }
}

qr_tuned_per_n <- dplyr::bind_rows(results)

qr_tuned_per_n %>%
  group_by(n) %>%
  summarise(
    mean_cov = mean(coverage),
    sd_cov   = sd(coverage),
    mean_len = mean(length),
    .groups = "drop"
  )

combined <- bind_rows(
  qr_results,          # din gamle tuned, måske tuned_results
  qr_tuned_per_n       # ny per-n tuning
)


combined %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "QR: fixed vs per-n tuning",
       y = "Coverage")


# kun n=1280
n <- 1280
M <- 200
n_test <- 100000
alpha <- 0.1

grid <- c(50, 80, 120, 200)

results <- list()
counter <- 1

for (m in 1:M) {
  
  train <- gen_data(n/2)
  cal   <- gen_data(n/2)
  test  <- gen_data(n_test)
  
  best_score <- Inf
  best_model <- NULL
  
  for (node_size in grid) {
    
    rf_tmp <- ranger(
      y ~ x,
      data = train,
      quantreg = TRUE,
      num.trees = 500,
      min.node.size = node_size,
      sample.fraction = 0.7,
      num.threads = parallel::detectCores()
    )
    
    pred_cal <- predict(
      rf_tmp, cal,
      type = "quantiles",
      quantiles = c(alpha/2, 1-alpha/2)
    )$predictions
    
    l <- pred_cal[,1]
    u <- pred_cal[,2]
    
    cov <- mean(cal$y >= l & cal$y <= u)
    len <- mean(u - l)

    score <- abs(cov - 0.9) + 0.005*len
    
    if (score < best_score) {
      best_score <- score
      best_model <- rf_tmp
    }
  }
  
  # test evaluation
  pred_test <- predict(
    best_model, test,
    type = "quantiles",
    quantiles = c(alpha/2, 1-alpha/2)
  )$predictions
  
  l <- pred_test[,1]
  u <- pred_test[,2]
  
  results[[counter]] <- data.frame(
    coverage = mean(test$y >= l & test$y <= u),
    length   = mean(u - l)
  )
  
  counter <- counter + 1
}

df <- dplyr::bind_rows(results)

cat("Mean coverage:", mean(df$coverage), "\n")
cat("SD:", sd(df$coverage), "\n")
cat("Mean length:", mean(df$length), "\n")


# Combine all

# =========================
# 1. Klargør QR (n ≤ 320)
# =========================
qr_small <- qr_tuned_per_n %>%
  filter(n <= 320)

# =========================
# 2. Klargør strong tuning (n = 1280)
# =========================
qr_large <- df %>%   # df = din strong tuning output
  mutate(
    n = 1280,
    m = row_number(),
    method = "QR_RF_tuned_per_n"
  )

# =========================
# 3. Kombinér QR
# =========================
qr_combined <- bind_rows(qr_small, qr_large)

# =========================
# 4. Fjern gammel QR og indsæt ny
# =========================
results_extended_new_1 <- results_extended_new %>%
  filter(method != "QR_RF_tuned") %>%
  bind_rows(qr_combined)

# =========================
# 5. Coverage plot
# =========================
results_extended_new_1 %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups = "drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  theme_minimal()



oracle_mean_new_1 <- results_extended_new_1 %>%
  filter(method == "Oracle") %>%
  group_by(n) %>%
  summarise(
    oracle_length = mean(length),
    .groups = "drop"
  )

rel_results_new_1 <- results_extended_new_1 %>%
  left_join(oracle_mean_new_1, by = "n") %>%
  mutate(rel_length = length / oracle_length)

rel_results_new_1 %>%
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
    title = "Relative interval length",
    y = "Length / Oracle length"
  ) +
  theme_minimal()





