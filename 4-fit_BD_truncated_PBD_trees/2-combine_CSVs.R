library(dplyr)
library(tidyr)

# List all CSV files in the specified directory
csv_files <- list.files(path = "simulations_output/4-Truncated_PBD/Fits", pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSV files, extracting information from filenames
combined_df <- csv_files %>%
  lapply(function(file) {
    df <- read.csv(file)
    file_info <- strsplit(basename(file), "_|\\.")[[1]]
    df <- df %>%
      mutate(param_vary = as.integer(gsub("par", "", file_info[3])),
             i_param_var = as.integer(gsub("var", "", file_info[4])),
             replicate = as.integer(gsub("rep", "", file_info[5])))
    return(df)
  }) %>%
  bind_rows()
combined_df <- combined_df[-1]

# Display the first few rows of the combined dataframe
head(combined_df)

# Read the parameters dataframe
pars_df <- read.csv("simulations_output/1-PBD/all_simulations_inference.csv")

# Join the combined dataframe with the parameters dataframe
combined_df <- combined_df %>%
  left_join(pars_df, by = c("param_vary", "i_param_var", "replicate"))

# Create new columns for "split_times" and "model" from the "X.x" column
combined_df <- combined_df %>%
  pivot_longer(cols = matches("^[VX][1-9]?[0-9]$"), names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = X.x, values_from = value) %>%
  rename(split_times = split_times,
         model = model,
         lamb_par = lamb_par,
         mu_par = mu_par)

combined_df$lamb_par <- abs(as.numeric(combined_df$lamb_par))
combined_df$mu_par <- abs(as.numeric(combined_df$mu_par))
combined_df$split_times <- as.numeric(combined_df$split_times)

# Compute equivalent constant rates
equivalent_bd_rates <- function(param) {
  l1 <- param[1]
  l2 <- param[2]
  l3 <- param[3]
  m1 <- param[4]
  m2 <- param[5]
  p <- 0.5 * (l2 + l3 + m2) / l3 * (1 - sqrt(1 - 4 * l3 * m2 / ((l2 + l3 + m2)^2)))
  l <- (1 - p) * l1
  m <- m1
  rates <- c(l, m)
  names(rates) <- NULL
  return(rates)
}
eq_rates <- apply(combined_df, 1, function(X)equivalent_bd_rates(as.numeric(c(X["PBD.l1"], X["PBD.l2"], X["PBD.l3"], X["PBD.mu1"], X["PBD.mu2"]))))
combined_df$eq_lamb_par <- eq_rates[1,]
combined_df$eq_mu_par <- eq_rates[2,]

# Select relevant columns for long dataframe
long_df <- combined_df %>%
  select(split_times, model, param_vary, i_param_var, replicate, variable, lamb_par, mu_par, eq_lamb_par, eq_mu_par, PBD.l1, PBD.l2, PBD.l3, PBD.mu1, PBD.mu2)

# Display the first few rows of the long dataframe
head(long_df)

# Save the long dataframe to a CSV file
write.csv(long_df, "simulations_output/4-Truncated_PBD/combined_long_df.csv", row.names = FALSE)
