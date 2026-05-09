# =============================
# SETUP
# =============================
library(dplyr)
library(ggplot2)
library(ranger)

# (antager disse allerede findes i din kode)
# n_list, M, n_test, alpha
# results_all_1, ..., results_all_4

# =============================
# CQR FUNCTION
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
  
  # calibration predictions
  pred_cal <- predict(
    rf, cal,
    type = "quantiles",
    quantiles = c(alpha/2, 1 - alpha/2)
  )$predictions
  
  # nonconformity scores
  E <- pmax(
    pred_cal[,1] - cal$y,
    cal$y - pred_cal[,2]
  )
  
  k <- ceiling((length(E) + 1)*(1 - alpha))
  q_hat <- sort(E)[k]
  
  # test predictions
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
# GENERIC RUN FUNCTION
# =============================
run_cqr_experiment <- function(gen_data, n_list, M, n_test, alpha) {
  
  results <- vector("list", length(n_list)*M)
  counter <- 1
  
  for (n in n_list) {
    
    n_train <- n/2
    n_cal   <- n/2
    
    for (m in 1:M) {
      
      train <- gen_data(n_train)
      cal   <- gen_data(n_cal)
      test  <- gen_data(n_test)
      
      cqr <- cqr_rf(train, cal, test, alpha)
      
      cov <- mean(test$y >= cqr$l & test$y <= cqr$u)
      wid <- mean(cqr$u - cqr$l)
      
      results[[counter]] <- data.frame(
        n = n,
        m = m,
        method = "CQR_RF",
        coverage = cov,
        length = wid
      )
      
      counter <- counter + 1
    }
  }
  
  bind_rows(results)
}

# =============================
# DATA GENERATORS
# =============================

# Exp 1 (Gaussian)
gen_data_1 <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = 2*x + 1 + rnorm(n))
}

# Exp 2 (asymmetric)
gen_data_2 <- function(n) {
  x <- runif(n, -1, 1)
  eps <- rexp(n, rate = 1) - 1
  data.frame(x = x, y = 2*x + 1 + eps)
}

# Exp 3 (nonlinear)
gen_data_3 <- function(n) {
  x <- runif(n, -1, 1)
  data.frame(x = x, y = sin(2*pi*x) + rnorm(n))
}

# Exp 4 (heteroskedastic)
gen_data_4 <- function(n) {
  x <- runif(n, -1, 1)
  sigma_x <- 1 + abs(x)
  data.frame(x = x, y = 2*x + 1 + sigma_x * rnorm(n))
}

# =============================
# RUN CQR
# =============================
cqr_1 <- run_cqr_experiment(gen_data_1, n_list, M, n_test, alpha)
cqr_2 <- run_cqr_experiment(gen_data_2, n_list, M, n_test, alpha)
cqr_3 <- run_cqr_experiment(gen_data_3, n_list, M, n_test, alpha)
# cqr_4 <- run_cqr_experiment(gen_data_4, n_list, M, n_test, alpha)

# =============================
# MERGE WITH EXISTING RESULTS
# =============================
results_all_1_ext <- bind_rows(results_all_1, cqr_1)
results_all_2_ext <- bind_rows(results_all_2, cqr_2)
results_all_3_ext <- bind_rows(results_all_3, cqr_3)
# results_all_4_ext <- bind_rows(results_all_4, cqr_4)
head(results_all_1)

# =============================
# PLOTTING FUNCTIONS
# =============================
plot_cov <- function(data) {
  data %>%
    group_by(n, method) %>%
    summarise(mean_cov = mean(coverage), .groups="drop") %>%
    ggplot(aes(x = n, y = mean_cov, color = method)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
    scale_x_log10() +
    theme_minimal()
}

plot_len <- function(data) {
  
  oracle <- data %>%
    filter(method == "Oracle") %>%
    group_by(n) %>%
    summarise(oracle_length = mean(length), .groups="drop")
  
  data %>%
    left_join(oracle, by="n") %>%
    mutate(rel_length = length / oracle_length) %>%
    filter(method != "Oracle") %>%
    group_by(n, method) %>%
    summarise(mean_rel = mean(rel_length), .groups="drop") %>%
    ggplot(aes(x = n, y = mean_rel, color = method)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    theme_minimal()
}


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






# =============================
# EXAMPLE PLOTS
# =============================
plot_cov(results_all_1_ext)
plot_len(results_all_1_ext)
