library(dplyr)
library(gridExtra)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

# Read the data
#setwd("~/Nextcloud/Recherche/1_Methods/PBD_analog/")
long_df <- read.csv("simulations_output/4-Truncated_PBD/combined_long_df.csv")

# Custom labeling functions
value_labels <- function(value, param_var) {
  # Labels for the columns
  labels <- list(bquote(lambda[1]), bquote(lambda[2]), bquote(lambda[3]), bquote(mu[1]), bquote(mu[2]))
  column <- labels[[as.numeric(param_var)]]
  
  # Unique values for the corresponding column
  columns_original <- c("PBD.l1", "PBD.l2", "PBD.l3", "PBD.mu1", "PBD.mu2")
  column_original <- columns_original[as.numeric(param_var)]
  
  # Return the formatted expression
  bquote(.(column) == .(round(subset_df[[column_original]][1], 2)))
}

# Create plots for speciation rates (lambda)
create_plot <- function(data, y_var, x_label, y_label, eq_var) {
  ggplot(data, aes(x = split_times, y = !!sym(y_var))) +
    geom_jitter(alpha = 0.1, col=brewer.pal(n = 7, name = "YlOrBr")[i_var + 2]) +
    geom_hline(aes(yintercept = !!sym(eq_var))) +
    scale_x_continuous(trans = "pseudo_log") +
    scale_y_continuous(trans = "sqrt", limits = c(0.0, 1.5)) +
    labs(title = plot_title, x = x_label, y = y_label) + 
    theme_tufte() +
    theme(
      plot.margin = margin(0, 5, 0, 5),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.line = element_line(size = 0.5, colour = "grey70"),
      axis.ticks = element_line(size = 0.5, colour = "grey70")
    )
}

create_plot <- function(data, y_var, x_label, y_label, eq_var) {
  # Calculate the median of the filtered data for each split time
  median_data <- data %>%
    group_by(split_times) %>%
    summarise(median_value = median(!!sym(y_var), na.rm = TRUE))
  
  # Filter out high values (some inferred rates explode)
  filtered_data <- data[data[y_var] < 3,]
  
  ggplot(filtered_data, aes(x = split_times, y = !!sym(y_var))) +
    geom_jitter(alpha = 0.1, col = brewer.pal(n = 7, name = "YlOrBr")[i_var + 2]) +
    #geom_line(alpha = 0.1, col = brewer.pal(n = 7, name = "YlOrBr")[i_var + 2], linetype = replicate) +
    geom_hline(aes(yintercept = !!sym(eq_var))) +
    geom_line(data = median_data, aes(x = split_times, y = median_value), linetype = "dotdash", color = "darkblue", lwd=0.8) +  # Add the median trend line
    scale_x_continuous(trans = "pseudo_log") +
    scale_y_continuous(trans = "sqrt", limits = c(0.0, NA)) +
    labs(title = plot_title, x = x_label, y = y_label) + 
    theme_tufte() +
    theme(
      plot.margin = margin(0, 5, 0, 5),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.line = element_line(size = 0.5, colour = "grey70"),
      axis.ticks = element_line(size = 0.5, colour = "grey70")
    )
}

# Lists for plots
plots_lambda <- list()
plots_mu <- list()

# Unique values for param_vary and i_param_var
unique_i_param_var <- unique(long_df$i_param_var)
unique_param_vary <- unique(long_df$param_vary)

for (i_var in unique_i_param_var) {
  for (param_var in unique_param_vary) {
    subset_df <- long_df %>% filter(param_vary == param_var, i_param_var == i_var)
    
    x_label <- if (i_var == 5) "Time of truncation" else ""
    y_label_lambda <- if (param_var == 1 & i_var == 3) "Inferred Speciation Rate" else ""
    y_label_mu <- if (param_var == 1 & i_var == 3) "Inferred Extinction Rate" else ""
    plot_title <- value_labels(i_var, param_var)
    
    # Create the plot for lamb_par
    p_lambda <- create_plot(subset_df, "lamb_par", x_label, y_label_lambda, "eq_lamb_par")
    plots_lambda <- c(plots_lambda, list(p_lambda))
    
    # Create the plot for mu_par
    p_mu <- create_plot(subset_df, "mu_par", x_label, y_label_mu, "eq_mu_par")
    plots_mu <- c(plots_mu, list(p_mu))
  }
}

# Arrange the plots in a grid and save the plots
plot_lambda_grid <- grid.arrange(grobs = plots_lambda, ncol = length(unique_param_vary))
plot_mu_grid <- grid.arrange(grobs = plots_mu, ncol = length(unique_param_vary))

# Save the figures
ggsave("fig/SM_truncation_lambda.pdf", plot_lambda_grid, width = 12, height = 7)
ggsave("fig/SM_truncation_mu.pdf", plot_mu_grid, width = 12, height = 7)
