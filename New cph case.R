# =========================
# 0. Pakker
# =========================
library(ranger)
library(dplyr)

# =========================
# 1. Preprocessing funktion
# =========================
prep_cph_data <- function(data, use_time = FALSE) {
  
  data <- data %>%
    mutate(
      cph_log_price = log(KONTANT_KOEBESUM),
      BYG_ALDER_AAR = pmax(0, BYG_ALDER_AAR),
      KOMMUNE_NAVN = as.factor(KOMMUNE_NAVN),
      LANDSDEL_NAVN = as.factor(LANDSDEL_NAVN)
    )
  
  if (use_time) {
    data <- data %>%
      mutate(cph_time = as.numeric(as.Date(SALG_KVARTAL)))
  }
  
  # simpel imputering
  num_vars <- sapply(data, is.numeric)
  data[num_vars] <- lapply(data[num_vars], function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  
  cat_vars <- sapply(data, is.character)
  data[cat_vars] <- lapply(data[cat_vars], function(x) {
    x[is.na(x)] <- "missing"
    as.factor(x)
  })
  
  data <- data %>% select(-KONTANT_KOEBESUM)
  
  return(data)
}

# =========================
# 2. Split funktion
# =========================
split_data <- function(data, split_type = "random") {
  
  n <- nrow(data)
  
  if (split_type == "random") {
    idx <- sample(1:n)
  } else {
    data <- data %>% arrange(SALG_KVARTAL)
    idx <- 1:n
  }
  
  train <- data[idx[1:floor(0.6*n)], ]
  cal   <- data[idx[(floor(0.6*n)+1):floor(0.8*n)], ]
  test  <- data[idx[(floor(0.8*n)+1):n], ]
  
  list(train=train, cal=cal, test=test)
}

# =========================
# 3. CP funktion
# =========================
run_cp <- function(train, cal, test, alpha=0.1) {
  
  # model
  model <- ranger(
    cph_log_price ~ .,
    data = train,
    num.trees = 300,
    seed = 1
  )
  
  # predictions
  pred_cal <- predict(model, cal)$predictions
  pred_test <- predict(model, test)$predictions
  
  # scores
  scores <- abs(cal$cph_log_price - pred_cal)
  n_cal <- length(scores)
  
  q_hat <- quantile(
    scores,
    probs = ceiling((n_cal + 1) * (1 - alpha)) / n_cal
  )
  
  lower <- pred_test - q_hat
  upper <- pred_test + q_hat
  
  coverage <- mean(test$cph_log_price >= lower &
                     test$cph_log_price <= upper)
  
  width <- mean(upper - lower)
  
  return(list(coverage=coverage, width=width))
}

# =========================
# 4. Adaptive CP funktion
# =========================
run_adaptive_cp <- function(train, cal, test, alpha=0.1) {
  
  # base model
  model <- ranger(cph_log_price ~ ., data = train, num.trees=300)
  
  pred_train <- predict(model, train)$predictions
  pred_cal <- predict(model, cal)$predictions
  pred_test <- predict(model, test)$predictions
  
  # residual model
  res_train <- abs(train$cph_log_price - pred_train)
  train$residual <- res_train
  
  sigma_model <- ranger(residual ~ ., data = train, num.trees=300)
  
  sigma_cal <- predict(sigma_model, cal)$predictions
  sigma_test <- predict(sigma_model, test)$predictions
  
  scores <- abs(cal$cph_log_price - pred_cal) / sigma_cal
  n_cal <- length(scores)
  
  q_hat <- quantile(
    scores,
    probs = ceiling((n_cal + 1) * (1 - alpha)) / n_cal
  )
  
  lower <- pred_test - q_hat * sigma_test
  upper <- pred_test + q_hat * sigma_test
  
  coverage <- mean(test$cph_log_price >= lower &
                     test$cph_log_price <= upper)
  
  width <- mean(upper - lower)
  
  return(list(coverage=coverage, width=width))
}

# =========================
# 5. Eksperiment
# =========================
run_experiment <- function(data, B=20) {
  
  results <- data.frame()
  
  for (b in 1:B) {
    
    for (split in c("random", "time")) {
      
      for (use_time in c(FALSE, TRUE)) {
        
        data_prep <- prep_cph_data(data, use_time)
        splits <- split_data(data_prep, split)
        splits$train <- splits$train %>% select(-SALG_KVARTAL)
        splits$cal   <- splits$cal   %>% select(-SALG_KVARTAL)
        splits$test  <- splits$test  %>% select(-SALG_KVARTAL)
        res_cp <- run_cp(splits$train, splits$cal, splits$test)
        res_adapt <- run_adaptive_cp(splits$train, splits$cal, splits$test)
        
        results <- rbind(results, data.frame(
          repetition = b,
          split = split,
          use_time = use_time,
          method = "CP",
          coverage = res_cp$coverage,
          width = res_cp$width
        ))
        
        results <- rbind(results, data.frame(
          repetition = b,
          split = split,
          use_time = use_time,
          method = "Adaptive",
          coverage = res_adapt$coverage,
          width = res_adapt$width
        ))
      }
    }
  }
  
  return(results)
}

# =========================
# 6. Kør eksperiment
# =========================
set.seed(1)

cph_results <- run_experiment(CphHousingPrices, B = 20)

# =========================
# 7. Opsummering
# =========================
summary_results <- cph_results %>%
  group_by(split, use_time, method) %>%
  summarise(
    coverage_mean = mean(coverage),
    coverage_sd = sd(coverage),
    width_mean = mean(width)
  )

print(summary_results)