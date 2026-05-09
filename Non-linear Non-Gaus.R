# =============================
# Setup
# =============================
set.seed(1)
library(ggplot2)
library(patchwork)

# -----------------------------
# Parameters
# -----------------------------
n <- 200
alpha <- 0.1
beta <- 2
a <- 1
sigma <- 1

df <- 3
scale_t <- sigma / sqrt(df/(df-2))

# Nonlinear function
f <- function(x) sin(2*pi*x)

# True quantile
q_true <- qt(1 - alpha/2, df = df) * scale_t

# -----------------------------
# Generate data
# -----------------------------
gen_data <- function(n) {
  x <- runif(n, -1, 1)
  eps <- rt(n, df = df) * scale_t
  y <- f(x) + eps
  data.frame(x, y)
}

train <- gen_data(n/2)
cal   <- gen_data(n/2)
test  <- gen_data(2000)

# Smooth grid for plotting
grid <- data.frame(x = seq(-1, 1, length.out = 200))
grid$f <- f(grid$x)

# -----------------------------
# Fit linear model
# -----------------------------
lm_fit <- lm(y ~ x, data = train)

b0 <- coef(lm_fit)[1]
b1 <- coef(lm_fit)[2]

pred_train <- b0 + b1 * train$x
pred_cal   <- b0 + b1 * cal$x
pred_grid  <- b0 + b1 * grid$x

# -----------------------------
# Oracle
# -----------------------------
grid$oracle_l <- grid$f - q_true
grid$oracle_u <- grid$f + q_true

# -----------------------------
# Linear
# -----------------------------
z <- qnorm(1 - alpha/2)
sigma_hat <- sqrt(mean((train$y - pred_train)^2))

grid$lin_l <- pred_grid - z * sigma_hat
grid$lin_u <- pred_grid + z * sigma_hat

# -----------------------------
# Naive
# -----------------------------
res_train <- train$y - pred_train
ql <- quantile(res_train, alpha/2)
qu <- quantile(res_train, 1 - alpha/2)

grid$naive_l <- pred_grid + ql
grid$naive_u <- pred_grid + qu

# -----------------------------
# QR (optim)
# -----------------------------
qr_fit <- function(x, y, tau) {
  X <- cbind(1, x)
  loss <- function(beta) {
    r <- y - X %*% beta
    sum((tau - (r < 0)) * r)
  }
  beta_init <- coef(lm(y ~ x))
  optim(beta_init, loss, method="BFGS")$par
}

qr_predict <- function(beta, x) {
  as.vector(cbind(1, x) %*% beta)
}

b_low  <- qr_fit(train$x, train$y, alpha/2)
b_high <- qr_fit(train$x, train$y, 1 - alpha/2)

grid$qr_l <- qr_predict(b_low, grid$x)
grid$qr_u <- qr_predict(b_high, grid$x)

# -----------------------------
# CP
# -----------------------------
res_cal <- abs(cal$y - pred_cal)
q_cp <- quantile(res_cal, 1 - alpha)

grid$cp_l <- pred_grid - q_cp
grid$cp_u <- pred_grid + q_cp

# -----------------------------
# Plot function
# -----------------------------
plot_method <- function(lower, upper, title) {
  ggplot(grid, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "skyblue", alpha = 0.4) +
    geom_line(aes(y = f), color = "black", size = 1) +
    geom_point(data = test,
               aes(x = x, y = y),
               alpha = 0.1, size = 0.5) +
    labs(title = title, y = "", x = "") +
    theme_bw()
}

# -----------------------------
# Create plots
# -----------------------------
p_oracle <- plot_method(grid$oracle_l, grid$oracle_u, "Oracle")
p_lin    <- plot_method(grid$lin_l, grid$lin_u, "Linear")
p_naive  <- plot_method(grid$naive_l, grid$naive_u, "Naive")
p_qr     <- plot_method(grid$qr_l, grid$qr_u, "QR")
p_cp     <- plot_method(grid$cp_l, grid$cp_u, "CP")

# Combine
(p_oracle | p_lin | p_naive) /
  (p_qr | p_cp)
