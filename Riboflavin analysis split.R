# Setup

library(glmnet)
library(hdi)
library(ggplot2)

set.seed(123)

data("riboflavin", package="hdi")

x <- as.matrix(riboflavin[,-1])
y <- riboflavin[,1]

x_scaled <- scale(x)

n <- nrow(x_scaled)

# Test / train split

train_idx <- sample(1:n, floor(0.7*n))

X_train <- x_scaled[train_idx,]
Y_train <- y[train_idx]

X_test <- x_scaled[-train_idx,]
Y_test <- y[-train_idx]

# Ridge prediction intervals

cv_ridge <- cv.glmnet(X_train, Y_train, alpha=0) # maybe use nfolds = 5 or 3 instead

ridge_model <- glmnet(
  X_train,
  Y_train,
  alpha=0,
  lambda=cv_ridge$lambda.min
)

ridge_pred <- predict(ridge_model, X_test)

ridge_train_pred <- predict(ridge_model, X_train)

sigma_ridge <- sqrt(mean((Y_train - ridge_train_pred)^2))

z <- qnorm(.975)

ridge_lower <- ridge_pred - z*sigma_ridge
ridge_upper <- ridge_pred + z*sigma_ridge

# lasso prediction intervals

cv_lasso <- cv.glmnet(X_train, Y_train, alpha=1) # maybe use nfolds = 5 or 3 instead

lasso_model <- glmnet(
  X_train,
  Y_train,
  alpha=1,
  lambda=cv_lasso$lambda.min
)

lasso_pred <- predict(lasso_model, X_test)

lasso_train_pred <- predict(lasso_model, X_train)

sigma_lasso <- sqrt(mean((Y_train - lasso_train_pred)^2))

lasso_lower <- lasso_pred - z*sigma_lasso
lasso_upper <- lasso_pred + z*sigma_lasso


# Split conformal function

split_conformal <- function(X,Y,X_new,alpha=0.05,penalty=1){
  
  n <- nrow(X)
  
  idx <- sample(1:n,floor(n/2))
  
  X_train <- X[idx,]
  Y_train <- Y[idx]
  
  X_cal <- X[-idx,]
  Y_cal <- Y[-idx]
  
  m <- length(Y_cal)
  
  cvfit <- cv.glmnet(X_train,Y_train,alpha=penalty)
  
  model <- glmnet(
    X_train,
    Y_train,
    alpha=penalty,
    lambda=cvfit$lambda.min
  )
  
  pred_cal <- predict(model,X_cal)
  
  R <- abs(Y_cal - pred_cal)
  
  k <- ceiling((m+1)*(1-alpha))
  
  q_hat <- sort(R)[k]
  
  pred_test <- as.vector(predict(model,X_new))
  
  lower <- pred_test - q_hat
  upper <- pred_test + q_hat
  
  list(lower=lower,upper=upper)
}

# Conformal ridge and lasso

conf_ridge <- split_conformal(X_train,Y_train,X_test,penalty=0)

conf_lasso <- split_conformal(X_train,Y_train,X_test,penalty=1)

# Coverage and width

ridge_cov <- mean(Y_test >= ridge_lower & Y_test <= ridge_upper)
lasso_cov <- mean(Y_test >= lasso_lower & Y_test <= lasso_upper)

conf_ridge_cov <- mean(Y_test >= conf_ridge$lower & Y_test <= conf_ridge$upper)
conf_lasso_cov <- mean(Y_test >= conf_lasso$lower & Y_test <= conf_lasso$upper)

ridge_width <- mean(ridge_upper - ridge_lower)
lasso_width <- mean(lasso_upper - lasso_lower)

conf_ridge_width <- mean(conf_ridge$upper - conf_ridge$lower)
conf_lasso_width <- mean(conf_lasso$upper - conf_lasso$lower)

comparison <- data.frame(
  Method=c("Ridge","Lasso","Conformal Ridge","Conformal Lasso"),
  Coverage=c(ridge_cov,lasso_cov,conf_ridge_cov,conf_lasso_cov),
  Width=c(ridge_width,lasso_width,conf_ridge_width,conf_lasso_width)
)
comparison

# Coverage-width tradeoff plot 

coverage_plot <- ggplot(comparison,
                        aes(x=Width,y=Coverage,color=Method))+
  geom_point(size=4)+
  geom_text(aes(label=Method),vjust=-1)+
  geom_hline(yintercept=.95,linetype="dashed",color="red")+
  labs(
    title="Coverage–Width Tradeoff",
    x="Average Interval Width",
    y="Empirical Coverage"
  )+
  theme_minimal()

ggsave("coverage_width_tradeoff.png",
       coverage_plot,
       width=6,
       height=4)

# Conditional coverage 

covered_conf <- (Y_test >= conf_lasso$lower &
                   Y_test <= conf_lasso$upper)

bins <- cut(
  X_test[,1],
  breaks=quantile(X_test[,1],probs=seq(0,1,0.2)),
  include.lowest=TRUE
)

conditional_coverage <- tapply(covered_conf,bins,mean)

cond_df <- data.frame(
  bin=names(conditional_coverage),
  coverage=as.numeric(conditional_coverage)
)

cond_plot <- ggplot(cond_df,
                    aes(x=bin,y=coverage))+
  geom_col(fill="steelblue")+
  geom_hline(yintercept=.95,
             linetype="dashed",
             color="red")+
  ylim(0,1)+
  labs(
    title="Conditional Coverage",
    x="Feature Bin",
    y="Coverage"
  )+
  theme_minimal()

ggsave("conditional_coverage.png",
       cond_plot,
       width=6,
       height=4)


