set.seed(1)

# Parameters
n <- 1000
n_test <- 100
sigma <- 0.2
alpha <- 0.1
z <- qnorm(1 - alpha/2)

# True function
f <- function(x) sin(2*pi*x)

# Training data
X <- runif(n)
Y <- f(X) + rnorm(n, 0, sigma)

# Test data
X_test <- runif(n_test)
Y_test <- f(X_test) + rnorm(n_test, 0, sigma)

# True prediction intervals
lower <- f(X_test) - z * sigma
upper <- f(X_test) + z * sigma

# Coverage
coverage <- mean(Y_test >= lower & Y_test <= upper)

# Width
width <- mean(upper - lower)

coverage # 0.9
width # approx 2*z*sigma = 0.658


# Compare to conformal


set.seed(1)

# Parameters
n <- 1000
n_test <- 100
sigma <- 0.2
alpha <- 0.1
z <- qnorm(1 - alpha/2)
B <- 200  # repetitions

# True function
f <- function(x) sin(2*pi*x)

# Storage
cov_oracle <- numeric(B)
width_oracle <- numeric(B)

cov_cp <- numeric(B)
width_cp <- numeric(B)

for (b in 1:B) {
  
  # Generate data
  X <- runif(n)
  Y <- f(X) + rnorm(n, 0, sigma)
  
  X_test <- runif(n_test)
  Y_test <- f(X_test) + rnorm(n_test, 0, sigma)
  
  # ----------------------
  # ORACLE
  # ----------------------
  lower_o <- f(X_test) - z * sigma
  upper_o <- f(X_test) + z * sigma
  
  cov_oracle[b] <- mean(Y_test >= lower_o & Y_test <= upper_o)
  width_oracle[b] <- mean(upper_o - lower_o)
  
  # ----------------------
  # SPLIT CONFORMAL
  # ----------------------
  idx <- sample(1:n, n/2)
  X_train <- X[idx]
  Y_train <- Y[idx]
  
  X_cal <- X[-idx]
  Y_cal <- Y[-idx]
  
  # Misspecificeret model (lineær regression)
  model <- lm(Y_train ~ X_train)
  
  # Calibration scores
  y_hat_cal <- predict(model, newdata = data.frame(X_train = X_cal))
  R <- abs(Y_cal - y_hat_cal)
  
  q_hat <- quantile(R, 1 - alpha)
  
  # Test predictions
  y_hat_test <- predict(model, newdata = data.frame(X_train = X_test))
  
  lower_cp <- y_hat_test - q_hat
  upper_cp <- y_hat_test + q_hat
  
  cov_cp[b] <- mean(Y_test >= lower_cp & Y_test <= upper_cp)
  width_cp[b] <- mean(upper_cp - lower_cp)
}

# Summary
mean(cov_oracle)
mean(width_oracle)

mean(cov_cp)
mean(width_cp)


plot(width_cp, cov_cp,
     xlab="Width", ylab="Coverage",
     main="CP: Coverage vs Width")
abline(h=0.9, col="red")

boxplot(width_oracle, width_cp,
        names=c("Oracle","CP"),
        main="Interval width comparison")


# Now with CQR


set.seed(1)

# Parameters
n <- 1000
n_test <- 100
sigma <- 0.2
alpha <- 0.1
z <- qnorm(1 - alpha/2)
B <- 200

f <- function(x) sin(2*pi*x)

# Storage
cov_oracle <- width_oracle <- numeric(B)
cov_cp <- width_cp <- numeric(B)
cov_cqr <- width_cqr <- numeric(B)

# Simple quantile approximation (iterative reweighting)
quantile_lm <- function(X, Y, tau, n_iter=10) {
  w <- rep(1, length(Y))
  for (k in 1:n_iter) {
    fit <- lm(Y ~ X, weights = w)
    res <- Y - predict(fit)
    w <- ifelse(res >= 0, tau, 1 - tau)
  }
  return(fit)
}

for (b in 1:B) {
  
  # Data
  X <- runif(n)
  Y <- f(X) + rnorm(n, 0, sigma)
  
  X_test <- runif(n_test)
  Y_test <- f(X_test) + rnorm(n_test, 0, sigma)
  
  # ---------------- ORACLE ----------------
  lower_o <- f(X_test) - z*sigma
  upper_o <- f(X_test) + z*sigma
  
  cov_oracle[b] <- mean(Y_test >= lower_o & Y_test <= upper_o)
  width_oracle[b] <- mean(upper_o - lower_o)
  
  # ---------------- SPLIT CP ----------------
  idx <- sample(1:n, n/2)
  
  X_train <- X[idx]; Y_train <- Y[idx]
  X_cal <- X[-idx]; Y_cal <- Y[-idx]
  
  model <- lm(Y_train ~ X_train)
  
  y_hat_cal <- predict(model, newdata=data.frame(X_train=X_cal))
  R <- abs(Y_cal - y_hat_cal)
  q_hat <- quantile(R, 1-alpha)
  
  y_hat_test <- predict(model, newdata=data.frame(X_train=X_test))
  
  lower_cp <- y_hat_test - q_hat
  upper_cp <- y_hat_test + q_hat
  
  cov_cp[b] <- mean(Y_test >= lower_cp & Y_test <= upper_cp)
  width_cp[b] <- mean(upper_cp - lower_cp)
  
  # ---------------- CQR ----------------
  
  # Fit quantiles
  model_l <- quantile_lm(X_train, Y_train, alpha/2)
  model_u <- quantile_lm(X_train, Y_train, 1 - alpha/2)
  
  ql_cal <- predict(model_l, newdata=data.frame(X=X_cal))
  qu_cal <- predict(model_u, newdata=data.frame(X=X_cal))
  
  R_cqr <- pmax(ql_cal - Y_cal, Y_cal - qu_cal)
  q_hat_cqr <- quantile(R_cqr, 1 - alpha)
  
  ql_test <- predict(model_l, newdata=data.frame(X=X_test))
  qu_test <- predict(model_u, newdata=data.frame(X=X_test))
  
  lower_cqr <- ql_test - q_hat_cqr
  upper_cqr <- qu_test + q_hat_cqr
  
  cov_cqr[b] <- mean(Y_test >= lower_cqr & Y_test <= upper_cqr)
  width_cqr[b] <- mean(upper_cqr - lower_cqr)
}

# Results
mean(cov_oracle); mean(width_oracle)
mean(cov_cp); mean(width_cp)
mean(cov_cqr); mean(width_cqr)


boxplot(width_oracle, width_cp, width_cqr,
        names = c("Oracle", "CP", "CQR"),
        main = "Interval Width Comparison",
        xlab = "Method",
        ylab = "Width")


plot(width_cp, cov_cp,
     xlab = "Width",
     ylab = "Coverage",
     main = "Coverage vs Width",
     pch = 16)

points(width_cqr, cov_cqr, pch = 16)

abline(h = 0.9)
legend("bottomright",
       legend = c("CP", "CQR"),
       pch = 16)



library(ggplot2)

df <- data.frame(
  width = c(width_oracle, width_cp, width_cqr),
  coverage = c(cov_oracle, cov_cp, cov_cqr),
  method = factor(rep(c("Oracle", "CP", "CQR"),
                      each = length(width_oracle)))
)

# Boxplot
ggplot(df, aes(x = method, y = width)) +
  geom_boxplot() +
  ggtitle("Interval Width Comparison")

# Coverage vs width
ggplot(df[df$method != "Oracle",], 
       aes(x = width, y = coverage, shape = method)) +
  geom_point() +
  geom_hline(yintercept = 0.9) +
  ggtitle("Coverage vs Width")


aggregate(width ~ method, df, mean)
aggregate(coverage ~ method, df, mean)









# Combine results
df <- data.frame(
  width = c(width_cp, width_cqr),
  coverage = c(cov_cp, cov_cqr),
  method = factor(rep(c("CP", "CQR"),
                      each = length(width_cp)))
)

# Romano-style plot
ggplot(df, aes(x = width, y = coverage, shape = method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(
    title = "Coverage vs Interval Width",
    x = "Average Interval Width",
    y = "Empirical Coverage"
  ) +
  theme_minimal()

oracle_df <- data.frame(
  width = mean(width_oracle),
  coverage = mean(cov_oracle)
)

ggplot(df, aes(x = width, y = coverage, shape = method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_point(data = oracle_df, aes(x = width, y = coverage),
             shape = 4, size = 4) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(
    title = "Coverage vs Interval Width (Romano-style)",
    x = "Average Interval Width",
    y = "Empirical Coverage"
  ) +
  theme_minimal()










set.seed(1)
library(ggplot2)
library(gridExtra)

# Data
n <- 1000
sigma <- 0.2
alpha <- 0.1

f <- function(x) sin(2*pi*x)

X <- runif(n, 0, 5)
Y <- f(X) + rnorm(n, 0, sigma)

# Grid for visualization
x_grid <- seq(0, 5, length.out = 200)

# ---------------------------
# SPLIT CONFORMAL
# ---------------------------
idx <- sample(1:n, n/2)

X_train <- X[idx]; Y_train <- Y[idx]
X_cal <- X[-idx]; Y_cal <- Y[-idx]

model <- lm(Y_train ~ X_train)

y_hat_cal <- predict(model, newdata=data.frame(X_train=X_cal))
R <- abs(Y_cal - y_hat_cal)
q_hat <- quantile(R, 1-alpha)

y_hat_grid <- predict(model, newdata=data.frame(X_train=x_grid))

lower_cp <- y_hat_grid - q_hat
upper_cp <- y_hat_grid + q_hat



# ---------------------------
# CQR (simple approximation)
# ---------------------------
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

ql_cal <- predict(model_l, newdata=data.frame(X=X_cal))
qu_cal <- predict(model_u, newdata=data.frame(X=X_cal))

R_cqr <- pmax(ql_cal - Y_cal, Y_cal - qu_cal)
q_hat_cqr <- quantile(R_cqr, 1-alpha)

ql_grid <- predict(model_l, newdata=data.frame(X=x_grid))
qu_grid <- predict(model_u, newdata=data.frame(X=x_grid))

lower_cqr <- ql_grid - q_hat_cqr
upper_cqr <- qu_grid + q_hat_cqr

# ---------------------------
# ORACLE
# ---------------------------
z <- qnorm(1 - alpha/2)
lower_oracle <- f(x_grid) - z*sigma
upper_oracle <- f(x_grid) + z*sigma

# ---------------------------
# PLOTS
# ---------------------------

# (a) Split CP
p1 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, y_hat_grid)) +
  geom_ribbon(aes(x=x_grid, ymin=lower_cp, ymax=upper_cp), alpha=0.3) +
  ggtitle("(a) Split Conformal") +
  xlab("X") + ylab("Y")

# (b) Oracle
p2 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, f(x_grid))) +
  geom_ribbon(aes(x=x_grid, ymin=lower_oracle, ymax=upper_oracle), alpha=0.3) +
  ggtitle("(b) Oracle") +
  xlab("X") + ylab("Y")

# (c) CQR
p3 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, ql_grid), linetype="dashed") +
  geom_line(aes(x_grid, qu_grid), linetype="dashed") +
  geom_ribbon(aes(x=x_grid, ymin=lower_cqr, ymax=upper_cqr), alpha=0.3) +
  ggtitle("(c) CQR") +
  xlab("X") + ylab("Y")

# (d) Interval length
length_cp <- upper_cp - lower_cp
length_cqr <- upper_cqr - lower_cqr
length_oracle <- upper_oracle - lower_oracle

df_len <- data.frame(
  x = rep(x_grid, 3),
  length = c(length_cp, length_cqr, length_oracle),
  method = rep(c("CP", "CQR", "Oracle"), each=length(x_grid))
)

p4 <- ggplot(df_len, aes(x, length, linetype=method)) +
  geom_line() +
  ggtitle("(d) Interval Length") +
  xlab("X") + ylab("Length")

# Combine
grid.arrange(p1, p2, p3, p4, ncol=2)





# With QR




set.seed(1)

# -----------------------------
# Parameters
# -----------------------------
n <- 1000
n_test <- 100
sigma <- 0.2
alpha <- 0.1
z <- qnorm(1 - alpha/2)
B <- 200

f <- function(x) sin(2*pi*x)

# -----------------------------
# Storage
# -----------------------------
cov_oracle <- width_oracle <- numeric(B)
cov_cp <- width_cp <- numeric(B)
cov_cqr <- width_cqr <- numeric(B)
cov_qr <- width_qr <- numeric(B)

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

# -----------------------------
# Simulation loop
# -----------------------------
for (b in 1:B) {
  
  # Data
  X <- runif(n)
  Y <- f(X) + rnorm(n, 0, sigma)
  
  X_test <- runif(n_test)
  Y_test <- f(X_test) + rnorm(n_test, 0, sigma)
  
  # ---------------- ORACLE ----------------
  lower_o <- f(X_test) - z*sigma
  upper_o <- f(X_test) + z*sigma
  
  cov_oracle[b] <- mean(Y_test >= lower_o & Y_test <= upper_o)
  width_oracle[b] <- mean(upper_o - lower_o)
  
  # ---------------- SPLIT ----------------
  idx <- sample(1:n, n/2)
  
  X_train <- X[idx]; Y_train <- Y[idx]
  X_cal <- X[-idx]; Y_cal <- Y[-idx]
  
  model <- lm(Y_train ~ X_train)
  
  y_hat_cal <- predict(model, newdata=data.frame(X_train=X_cal))
  R <- abs(Y_cal - y_hat_cal)
  q_hat <- quantile(R, 1-alpha)
  
  y_hat_test <- predict(model, newdata=data.frame(X_train=X_test))
  
  lower_cp <- y_hat_test - q_hat
  upper_cp <- y_hat_test + q_hat
  
  cov_cp[b] <- mean(Y_test >= lower_cp & Y_test <= upper_cp)
  width_cp[b] <- mean(upper_cp - lower_cp)
  
  # ---------------- QR ----------------
  model_l_qr <- quantile_lm(X_train, Y_train, alpha/2)
  model_u_qr <- quantile_lm(X_train, Y_train, 1-alpha/2)
  
  ql_test_qr <- predict(model_l_qr, newdata=data.frame(X=X_test))
  qu_test_qr <- predict(model_u_qr, newdata=data.frame(X=X_test))
  
  lower_qr <- ql_test_qr
  upper_qr <- qu_test_qr
  
  cov_qr[b] <- mean(Y_test >= lower_qr & Y_test <= upper_qr)
  width_qr[b] <- mean(upper_qr - lower_qr)
  
  # ---------------- CQR ----------------
  ql_cal <- predict(model_l_qr, newdata=data.frame(X=X_cal))
  qu_cal <- predict(model_u_qr, newdata=data.frame(X=X_cal))
  
  R_cqr <- pmax(ql_cal - Y_cal, Y_cal - qu_cal)
  q_hat_cqr <- quantile(R_cqr, 1-alpha)
  
  ql_test <- predict(model_l_qr, newdata=data.frame(X=X_test))
  qu_test <- predict(model_u_qr, newdata=data.frame(X=X_test))
  
  lower_cqr <- ql_test - q_hat_cqr
  upper_cqr <- qu_test + q_hat_cqr
  
  cov_cqr[b] <- mean(Y_test >= lower_cqr & Y_test <= upper_cqr)
  width_cqr[b] <- mean(upper_cqr - lower_cqr)
}

# -----------------------------
# Summary tables
# -----------------------------
df <- data.frame(
  width = c(width_oracle, width_qr, width_cp, width_cqr),
  coverage = c(cov_oracle, cov_qr, cov_cp, cov_cqr),
  method = factor(rep(c("Oracle","QR","CP","CQR"),
                      each = B))
)

aggregate(width ~ method, df, mean)
aggregate(coverage ~ method, df, mean)

# -----------------------------
# Boxplot
# -----------------------------
boxplot(width ~ method, df,
        main = "Interval Width Comparison",
        xlab = "Method",
        ylab = "Width")

# -----------------------------
# Romano-style scatter
# -----------------------------
plot(df$width[df$method=="QR"],
     df$coverage[df$method=="QR"],
     pch=16, xlab="Width", ylab="Coverage",
     main="Coverage vs Width")

points(df$width[df$method=="CP"],
       df$coverage[df$method=="CP"], pch=16)

points(df$width[df$method=="CQR"],
       df$coverage[df$method=="CQR"], pch=16)

abline(h=0.9)

legend("bottomright",
       legend=c("QR","CP","CQR"),
       pch=16)



# Romano plot with all


library(ggplot2)
library(gridExtra)

set.seed(1)

# -----------------------------
# Data
# -----------------------------
n <- 1000
sigma <- 0.2
alpha <- 0.1
z <- qnorm(1 - alpha/2)

f <- function(x) x^2 # or sin(2*pi*x) 

X <- runif(n, 0, 5)
Y <- f(X) + rnorm(n, 0, sigma)

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
# QR (no conformal)
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
# Oracle
# -----------------------------
lower_oracle <- f(x_grid) - z*sigma
upper_oracle <- f(x_grid) + z*sigma

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

# (a) CP
p1 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, y_hat_grid)) +
  geom_ribbon(aes(x=x_grid, ymin=lower_cp, ymax=upper_cp), alpha=0.3) +
  ggtitle("(a) Split Conformal")

# (b) Oracle
p2 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, f(x_grid))) +
  geom_ribbon(aes(x=x_grid, ymin=lower_oracle, ymax=upper_oracle), alpha=0.3) +
  ggtitle("(b) Oracle")

# (c) QR
p3 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, ql_grid), linetype="dashed") +
  geom_line(aes(x_grid, qu_grid), linetype="dashed") +
  geom_ribbon(aes(x=x_grid, ymin=lower_qr, ymax=upper_qr), alpha=0.3) +
  ggtitle("(c) QR")

# (d) CQR
p4 <- ggplot() +
  geom_point(aes(X, Y), alpha=0.3) +
  geom_line(aes(x_grid, ql_grid), linetype="dashed") +
  geom_line(aes(x_grid, qu_grid), linetype="dashed") +
  geom_ribbon(aes(x=x_grid, ymin=lower_cqr, ymax=upper_cqr), alpha=0.3) +
  ggtitle("(d) CQR")

# (e) Length
p5 <- ggplot(df_len, aes(x, length, linetype=method)) +
  geom_line() +
  ggtitle("(e) Interval Length")

# Layout (4 panels)
grid.arrange(p2, p1, p3, p4, ncol=2)

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

