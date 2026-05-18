# =========================================================
# Conditional coverage + width analysis for Split CP
# =========================================================

library(dplyr)
library(ggplot2)

# =========================================================
# 1. Create analysis dataframe
# =========================================================

cph_analysis_df <- cph_intervals_r %>%
  mutate(
    
    # Coverage indicator
    covered = faktisk >= lower & faktisk <= upper,
    
    # Absolute interval width
    width = upper - lower,
    
    # Relative width
    rel_width = width / faktisk,
    
    # Log width
    log_width = log(width),
    
    # Price in millions
    price_million = faktisk / 1e6
  ) %>%
  
  bind_cols(
    cph_test_r %>%
      select(
        cph_salgsaar,
        KOMMUNE_NAVN,
        LANDSDEL_NAVN,
        ENH_BEBO_ARL
      )
  )

# =========================================================
# 2. Coverage by sale year
# =========================================================

coverage_year <- cph_analysis_df %>%
  group_by(cph_salgsaar) %>%
  summarise(
    coverage = mean(covered),
    n = n(),
    .groups = "drop"
  )

p_cov_year <- ggplot(
  coverage_year,
  aes(x = cph_salgsaar,
      y = coverage)
) +
  
  geom_line() +
  geom_point(size = 2) +
  
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    color = "red"
  ) +
  
  ylim(0,1) +
  
  theme_minimal() +
  
  labs(
    title = "Coverage by sale year",
    x = "Sale year",
    y = "Coverage"
  )

print(p_cov_year)

# =========================================================
# 3. Width by sale year
# =========================================================

width_year <- cph_analysis_df %>%
  group_by(cph_salgsaar) %>%
  summarise(
    mean_width = mean(width),
    mean_rel_width = mean(rel_width),
    .groups = "drop"
  )

p_width_year <- ggplot(
  width_year,
  aes(x = cph_salgsaar,
      y = mean_width)
) +
  
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen", size = 2) +
  
  theme_minimal() +
  
  labs(
    title = "Interval width by sale year",
    x = "Sale year",
    y = "Mean width (DKK)"
  )

print(p_width_year)

# =========================================================
# 4. Create price bins
# =========================================================

breaks_price <- quantile(
  cph_analysis_df$faktisk,
  probs = seq(0,1,length.out = 11),
  na.rm = TRUE
)

cph_analysis_df <- cph_analysis_df %>%
  mutate(
    price_bin = cut(
      faktisk,
      breaks = breaks_price,
      include.lowest = TRUE
    )
  )

# =========================================================
# 5. Coverage by price
# =========================================================

coverage_price <- cph_analysis_df %>%
  group_by(price_bin) %>%
  summarise(
    coverage = mean(covered),
    n = n(),
    .groups = "drop"
  )

p_cov_price <- ggplot(
  coverage_price,
  aes(x = price_bin,
      y = coverage,
      group = 1)
) +
  
  geom_line() +
  geom_point(size = 2) +
  
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    color = "red"
  ) +
  
  ylim(0,1) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "Coverage by price bin",
    x = "Price bin",
    y = "Coverage"
  )

print(p_cov_price)

# =========================================================
# 6. Width by price
# =========================================================

width_price <- cph_analysis_df %>%
  group_by(price_bin) %>%
  summarise(
    mean_width = mean(width),
    mean_rel_width = mean(rel_width),
    .groups = "drop"
  )

# Absolute width

p_width_price <- ggplot(
  width_price,
  aes(x = price_bin,
      y = mean_width,
      group = 1)
) +
  
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen", size = 2) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "Interval width by price",
    x = "Price bin",
    y = "Mean width (DKK)"
  )

print(p_width_price)

# =========================================================
# 7. Relative width by price
# =========================================================

p_rel_width <- ggplot(
  width_price,
  aes(x = price_bin,
      y = mean_rel_width,
      group = 1)
) +
  
  geom_line(color = "firebrick") +
  geom_point(color = "firebrick", size = 2) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "Relative interval width by price",
    x = "Price bin",
    y = "Relative width"
  )

print(p_rel_width)

# =========================================================
# 8. Coverage by municipality
# =========================================================

coverage_area <- cph_analysis_df %>%
  group_by(KOMMUNE_NAVN) %>%
  summarise(
    coverage = mean(covered),
    width = mean(width),
    n = n(),
    .groups = "drop"
  ) %>%
  
  filter(n > 50)

p_cov_area <- ggplot(
  coverage_area,
  aes(x = reorder(KOMMUNE_NAVN, coverage),
      y = coverage)
) +
  
  geom_point(size = 2) +
  
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    color = "red"
  ) +
  
  coord_flip() +
  
  ylim(0,1) +
  
  theme_minimal() +
  
  labs(
    title = "Coverage by municipality",
    x = "",
    y = "Coverage"
  )

print(p_cov_area)

# =========================================================
# 9. Width by municipality
# =========================================================

p_width_area <- ggplot(
  coverage_area,
  aes(x = reorder(KOMMUNE_NAVN, width),
      y = width)
) +
  
  geom_point(
    size = 2,
    color = "darkgreen"
  ) +
  
  coord_flip() +
  
  theme_minimal() +
  
  labs(
    title = "Interval width by municipality",
    x = "",
    y = "Mean width"
  )

print(p_width_area)

# =========================================================
# 10. Scatterplot:
# Relative width vs price
# =========================================================

p_scatter <- ggplot(
  cph_analysis_df,
  aes(x = price_million,
      y = rel_width)
) +
  
  geom_point(alpha = 0.15) +
  
  geom_smooth() +
  
  theme_minimal() +
  
  labs(
    title = "Relative width vs price",
    x = "Price (million DKK)",
    y = "Relative width"
  )

print(p_scatter)

# =========================================================
# 11. Optional summaries
# =========================================================

cat("\n============================\n")
cat("Overall coverage:\n")
print(mean(cph_analysis_df$covered))

cat("\n============================\n")
cat("Overall mean width:\n")
print(mean(cph_analysis_df$width))

cat("\n============================\n")
cat("Overall mean relative width:\n")
print(mean(cph_analysis_df$rel_width))


# =========================================================
# Localized / grouped CP analysis by municipality
# =========================================================

library(dplyr)
library(ggplot2)

# =========================================================
# 1. Create analysis dataframe
# =========================================================

cph_analysis_df <- cph_intervals_r %>%
  mutate(
    
    # Coverage indicator
    covered = faktisk >= lower & faktisk <= upper,
    
    # Absolute width
    width = upper - lower,
    
    # Relative width
    rel_width = width / faktisk,
    
    # Price in millions
    price_million = faktisk / 1e6
  ) %>%
  
  bind_cols(
    cph_test_r %>%
      select(
        cph_salgsaar,
        KOMMUNE_NAVN,
        LANDSDEL_NAVN,
        ENH_BEBO_ARL
      )
  )

# =========================================================
# 2. Compute municipality average prices
# =========================================================

kommune_price <- cph_analysis_df %>%
  group_by(KOMMUNE_NAVN) %>%
  summarise(
    mean_price = mean(faktisk),
    median_price = median(faktisk),
    n = n(),
    .groups = "drop"
  )

# =========================================================
# 3. Create municipality bins
# =========================================================

kommune_price <- kommune_price %>%
  mutate(
    
    kommune_bin = cut(
      mean_price,
      
      breaks = quantile(
        mean_price,
        probs = seq(0,1,length.out = 4),
        na.rm = TRUE
      ),
      
      labels = c(
        "Low-price municipalities",
        "Medium-price municipalities",
        "High-price municipalities"
      ),
      
      include.lowest = TRUE
    )
  )

kommune_price %>%
  arrange(mean_price) 
print(kommune_price %>%
        arrange(mean_price) , n = 29)
# =========================================================
# 4. Merge bins back into analysis dataframe
# =========================================================

cph_analysis_df <- cph_analysis_df %>%
  left_join(
    kommune_price %>%
      select(
        KOMMUNE_NAVN,
        kommune_bin,
        mean_price,
        median_price
      ),
    
    by = "KOMMUNE_NAVN"
  )

# =========================================================
# 5. Coverage + width by municipality bin
# =========================================================

coverage_local <- cph_analysis_df %>%
  group_by(kommune_bin) %>%
  summarise(
    
    coverage = mean(covered),
    
    mean_width = mean(width),
    
    mean_rel_width = mean(rel_width),
    
    mean_price = mean(faktisk),
    
    n = n(),
    
    .groups = "drop"
  )

print(coverage_local)

# =========================================================
# 6. Coverage plot
# =========================================================

p_cov_local <- ggplot(
  coverage_local,
  
  aes(
    x = kommune_bin,
    y = coverage,
    group = 1
  )
) +
  
  geom_point(size = 3) +
  
  geom_line() +
  
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    color = "red"
  ) +
  
  ylim(0,1) +
  
  theme_minimal() +
  
  labs(
    title = "Conditional coverage by municipality price group",
    x = "",
    y = "Coverage"
  )

print(p_cov_local)

# =========================================================
# 7. Absolute width plot
# =========================================================

p_width_local <- ggplot(
  coverage_local,
  
  aes(
    x = kommune_bin,
    y = mean_width,
    group = 1
  )
) +
  
  geom_point(
    size = 3,
    color = "darkgreen"
  ) +
  
  geom_line(
    color = "darkgreen"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Interval width by municipality price group",
    x = "",
    y = "Mean width (DKK)"
  )

print(p_width_local)

# =========================================================
# 8. Relative width plot
# =========================================================

p_rel_local <- ggplot(
  coverage_local,
  
  aes(
    x = kommune_bin,
    y = mean_rel_width,
    group = 1
  )
) +
  
  geom_point(
    size = 3,
    color = "firebrick"
  ) +
  
  geom_line(
    color = "firebrick"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Relative interval width by municipality price group",
    x = "",
    y = "Relative width"
  )

print(p_rel_local)

# =========================================================
# 9. Municipality-level scatterplot
# =========================================================

p_scatter <- ggplot(
  kommune_price,
  
  aes(
    x = mean_price / 1e6,
    y = n
  )
) +
  
  geom_point(size = 3) +
  
  theme_minimal() +
  
  labs(
    title = "Municipality mean price vs sample size",
    x = "Mean municipality price (million DKK)",
    y = "Number of observations"
  )

print(p_scatter)

# =========================================================
# 10. Coverage by individual municipality
# =========================================================

coverage_kommune <- cph_analysis_df %>%
  group_by(KOMMUNE_NAVN) %>%
  summarise(
    
    coverage = mean(covered),
    
    mean_width = mean(width),
    
    mean_rel_width = mean(rel_width),
    
    n = n(),
    
    .groups = "drop"
  ) %>%
  
  filter(n > 50)

# =========================================================
# 11. Municipality coverage plot
# =========================================================

p_cov_kommune <- ggplot(
  coverage_kommune,
  
  aes(
    x = reorder(KOMMUNE_NAVN, coverage),
    y = coverage
  )
) +
  
  geom_point(size = 2) +
  
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    color = "red"
  ) +
  
  coord_flip() +
  
  ylim(0,1) +
  
  theme_minimal() +
  
  labs(
    title = "Coverage by municipality",
    x = "",
    y = "Coverage"
  )

print(p_cov_kommune)

# =========================================================
# 12. Municipality relative width plot
# =========================================================

p_rel_kommune <- ggplot(
  coverage_kommune,
  
  aes(
    x = reorder(KOMMUNE_NAVN, mean_rel_width),
    y = mean_rel_width
  )
) +
  
  geom_point(
    size = 2,
    color = "firebrick"
  ) +
  
  coord_flip() +
  
  theme_minimal() +
  
  labs(
    title = "Relative interval width by municipality",
    x = "",
    y = "Relative width"
  )

print(p_rel_kommune)

# =========================================================
# 13. Summary
# =========================================================

cat("\n=====================================\n")
cat("Localized CP analysis summary\n")
cat("=====================================\n\n")

print(coverage_local)


