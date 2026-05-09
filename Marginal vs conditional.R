n <- 100
p <- 200
X <- matrix(rnorm(n*p), n, p)
beta <- c(rep(2,10), rep(0,p-10))
Y <- X %*% beta + rnorm(n)

# marginal
mean(Y_test >= lower & Y_test <= upper)

# conditional 
bins <- cut(X_test[,1], breaks=quantile(X_test[,1], probs=seq(0,1,.2)))
tapply(covered, bins, mean)