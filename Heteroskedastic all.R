library(ggplot2)
library(gridExtra)

set.seed(1)

# -----------------------------
# Data
# -----------------------------
n <- 1000
alpha <- 0.1
z <- qnorm(1 - alpha/2)

f <- function(x) x   # or sin(2*pi*x)

# Heteroskedastic variance function
sigma_fun <- function(x) 0.1 + 0.3 * x

X <- runif(n, 0, 5)
Y <- f(X) + rnorm(n, 0, sigma_fun(X))

x_grid <- seq(0, 5, length.out = 200)

# -----------------------------
# Split
# -----------------------------
idx <- sample(1:n, n/2)
X_train <- X[idx]; Y_train <- Y[idx]
X_cal <- X[-idx]; Y_cal <- Y[-idx]

# -----------------------------
# CP
# -----------------------------
model_cp <- lm(Y_train ~ X_train)

y_hat_cal <- predict(model_cp, newdata=data.frame(X_train=X_cal))
R <- abs(Y_cal - y_hat_cal)
q_hat <- quantile(R, 1-alpha)

y_hat_grid <- predict(model_cp, newdata=data.frame(X_train=x_grid))

lower_cp <- y_hat_grid - q_hat
upper_cp <- y_hat_grid + q_hat

# -----------------------------
# Quantile regression (approx)
# -----------------------------
quantile_lm <- function(X, Y, tau, n_iter=10) {
  w <- rep(1, length(Y))
  for (k in 1:n_iter) {
    fit <- lm(Y ~ X, weights = w)
    res <- Y - predict(fit)
    w <- ifelse(res >= 0, tau, 1 - tau)
  }
  return(fit)
}

model_l <- quantile_lm(X_train, Y_train, alpha/2)
model_u <- quantile_lm(X_train, Y_train, 1-alpha/2)

ql_grid <- predict(model_l, newdata=data.frame(X=x_grid))
qu_grid <- predict(model_u, newdata=data.frame(X=x_grid))

# -----------------------------
# QR
# -----------------------------
lower_qr <- ql_grid
upper_qr <- qu_grid

# -----------------------------
# CQR
# -----------------------------
ql_cal <- predict(model_l, newdata=data.frame(X=X_cal))
qu_cal <- predict(model_u, newdata=data.frame(X=X_cal))

R_cqr <- pmax(ql_cal - Y_cal, Y_cal - qu_cal)
q_hat_cqr <- quantile(R_cqr, 1-alpha)

lower_cqr <- ql_grid - q_hat_cqr
upper_cqr <- qu_grid + q_hat_cqr

# -----------------------------
# Oracle (UPDATED!)
# -----------------------------
lower_oracle <- f(x_grid) - z * sigma_fun(x_grid)
upper_oracle <- f(x_grid) + z * sigma_fun(x_grid)

# -----------------------------
# Lengths
# -----------------------------
df_len <- data.frame(
  x = rep(x_grid, 4),
  length = c(upper_cp - lower_cp,
             upper_qr - lower_qr,
             upper_cqr - lower_cqr,
             upper_oracle - lower_oracle),
  method = rep(c("CP","QR","CQR","Oracle"),
               each = length(x_grid))
)

# -----------------------------
# Plots
# -----------------------------

p1 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, y_hat_grid)) +
  geom_ribbon(aes(x=x_grid, ymin=lower_cp, ymax=upper_cp), alpha=0.3) +
  ggtitle("(a) Split Conformal")

p2 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, f(x_grid))) +
  geom_ribbon(aes(x=x_grid, ymin=lower_oracle, ymax=upper_oracle), alpha=0.3) +
  ggtitle("(b) Oracle")

p3 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, ql_grid), linetype="dashed") +
  geom_line(aes(x_grid, qu_grid), linetype="dashed") +
  geom_ribbon(aes(x=x_grid, ymin=lower_qr, ymax=upper_qr), alpha=0.3) +
  ggtitle("(c) QR")

p4 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, ql_grid), linetype="dashed") +
  geom_line(aes(x_grid, qu_grid), linetype="dashed") +
  geom_ribbon(aes(x=x_grid, ymin=lower_cqr, ymax=upper_cqr), alpha=0.3) +
  ggtitle("(d) CQR")

grid.arrange(p2, p1, p3, p4, ncol=2)

# -----------------------------
# Length plot (heteroskedastic effect!)
# -----------------------------
p5 <- ggplot(df_len, aes(x, length, color = method, linetype = method)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(
    "CP" = "#1b9e77",
    "QR" = "#d95f02",
    "CQR" = "#7570b3",
    "Oracle" = "#000000"
  )) +
  scale_linetype_manual(values = c(
    "CP" = "solid",
    "QR" = "dashed",
    "CQR" = "dotted",
    "Oracle" = "longdash"
  )) +
  ggtitle("(e) Interval Length") +
  xlab("X") +
  ylab("Length") +
  theme_minimal()

p5