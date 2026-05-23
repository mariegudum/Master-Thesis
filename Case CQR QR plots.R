# =====================================================
# Libraries
# =====================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# =====================================================
# Helper function
# =====================================================

make_coverage_plot <- function(df,
                               xvar,
                               title,
                               xlab = NULL,
                               reorder_var = FALSE,
                               rotate_x = FALSE,
                               ylim = c(0.75, 1.0)) {
  
  plot_df <- df %>%
    group_by(.data[[xvar]]) %>%
    summarise(
      coverage = mean(covered),
      .groups = "drop"
    )
  
  if(reorder_var){
    
    plot_df[[xvar]] <- reorder(
      plot_df[[xvar]],
      plot_df$coverage
    )
  }
  
  p <- ggplot(
    plot_df,
    
    aes(
      x = .data[[xvar]],
      y = coverage,
      group = 1
    )
  ) +
    
    geom_line(
      color = "black",
      linewidth = 0.8
    ) +
    
    geom_point(
      color = "black",
      size = 1.8
    ) +
    
    geom_hline(
      yintercept = 0.9,
      color = "red",
      linetype = "dashed"
    ) +
    
    scale_y_continuous(
      limits = ylim
    ) +
    
    labs(
      title = title,
      x = xlab,
      y = "Coverage"
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(
        size = 11,
        face = "bold"
      ),
      
      axis.title = element_text(size = 10),
      
      axis.text = element_text(size = 8)
    )
  
  if(rotate_x){
    
    p <- p +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      )
  }
  
  return(p)
}

# =====================================================
# PRICE BINS
# =====================================================

# If not already created

breaks_price <- quantile(
  cph_analysis_df$faktisk / 1e6,
  probs = seq(0,1,length.out = 11),
  na.rm = TRUE
)

qr_analysis_df <- qr_analysis_df %>%
  mutate(
    price_bin = cut(
      faktisk / 1e6,
      breaks = breaks_price,
      include.lowest = TRUE
    )
  )

cqr_analysis_df <- cqr_analysis_df %>%
  mutate(
    price_bin = cut(
      faktisk / 1e6,
      breaks = breaks_price,
      include.lowest = TRUE
    )
  )

# =====================================================
# QR plots
# =====================================================

p_qr_year <-
  make_coverage_plot(
    qr_analysis_df,
    "cph_salgsaar",
    "QR: Sale year",
    xlab = "Sale year"
  )

p_qr_price <-
  make_coverage_plot(
    qr_analysis_df,
    "price_bin",
    "QR: Price bin",
    xlab = "Price interval (million DKK)",
    rotate_x = TRUE,
    ylim = c(0.6,1.0)
  )

p_qr_area <-
  make_coverage_plot(
    qr_analysis_df,
    "KOMMUNE_NAVN",
    "QR: Municipality",
    reorder_var = TRUE,
    ylim = c(0.6,1.0)
  ) +
  coord_flip()

# =====================================================
# CQR plots
# =====================================================

p_cqr_year <-
  make_coverage_plot(
    cqr_analysis_df,
    "cph_salgsaar",
    "CQR: Sale year",
    xlab = "Sale year"
  )

p_cqr_price <-
  make_coverage_plot(
    cqr_analysis_df,
    "price_bin",
    "CQR: Price bin",
    xlab = "Price interval (million DKK)",
    rotate_x = TRUE,
    ylim = c(0.6,1.0)
  )

p_cqr_area <-
  make_coverage_plot(
    cqr_analysis_df,
    "KOMMUNE_NAVN",
    "CQR: Municipality",
    reorder_var = TRUE,
    ylim = c(0.6,1.0)
  ) +
  coord_flip()

# =====================================================
# Combine all plots
# =====================================================

(p_qr_year | p_qr_price | p_qr_area) /
  (p_cqr_year | p_cqr_price | p_cqr_area)