#' DiagnLamp.R
#'
#' @importFrom ggplot2 geom_histogram geom_density geom_vline scale_color_manual labs theme
#' @importFrom stats coef cor dbinom density dexp dnorm dpois dunif influence.measures predict qnorm sd var
#'
#' @name simulate_estimator_mse
#' @title Simulation of Estimator MSE
#' @description This function performs a Monte Carlo simulation to evaluate the Mean Squared Error (MSE) of a given estimator.
# It generates data using a specified data generation function, applies the estimator function to the generated data, and computes the MSE, bias, and variance of the estimator. Additionally, it checks for empirical consistency of the estimator by evaluating how the MSE changes with increasing sample sizes.
#' @param n Integer. The sample size for each simulation.
#' @param true_param Numeric. The true parameter value that the estimator is trying to estimate.
#' @param data_gen_func Function. A function that generates data given a sample size n.
#' @param estimator_func Function. A function that computes the estimator given a dataset.
#' @param nsim Integer. The number of Monte Carlo simulations to perform (default is 1000).
#' @return A list containing the MSE, bias, variance, consistency check, and the estimates from the simulations.
#' @examples
#' # Define a data generation function for normal data
#' generate_data <- function(n) { rnorm(n, mean = 0, sd = 2) }
#'
#' # Define an estimator function (sample variance)
#' correct_estimator <- function(data) { var(data) }
#'
#' # Run the simulation
#' results <- simulate_estimator_mse(n = 10,
#' true_param = 4,
#' data_gen_func = generate_data,
#' estimator_func = correct_estimator, nsim = 1000)
#'
#' @note
#' The consistency check is based on the Spearman correlation between sample sizes and their corresponding MSEs. A strongly negative correlation (<= -0.9) indicates empirical consistency of the estimator.
#' @export

globalVariables(c("Dataset", "MSE", "Estimate", "x", "y", "Observation", "Leverage", "Influential_Leverage", "Cooks_Distance", "Influential_Cooks", "generate_data", "true_param_vector", "mean_estimate_vector"))

simulate_estimator_mse <- function(n, true_param, data_gen_func, estimator_func, nsim = 1000) {

  # initialize a vector to store estimations
  estimates <- numeric(nsim)

  # Perform Monte Carlo simulations
  estimates <- sapply(1:nsim, function(i) {
    data <- generate_data(n)
    estimator_func(data)
  })

  # Calculate metrics: MSE = Bias^2 + Variance
  mean_estimate <- mean(estimates)
  bias <- mean_estimate - true_param
  variance <- var(estimates)
  mse <- mean((estimates - true_param)^2) # Equal to: bias^2 + variance

  # We use a sequence of logarithmic n but restrained to don't weight on RAM
  n_seq <- c(10, 50, 100, 500, 1000)
  mse_seq <- numeric(length(n_seq))

  mse_seq <- sapply(seq_along(n_seq), function(k) {
    current_n <- n_seq[k]
    estimates2 <- numeric(nsim)

    estimates2 <- sapply(1:nsim, function(i) {
      data <- generate_data(current_n)
      estimator_func(data)
    })

    mean((estimates2 - true_param)^2)

  })

  # Verify if the MSE monotonously decreases as n rises. We use the Spearman correlation that only evaluates the ranking: if n increases, the MSE decreases.
  correlation <- cor(n_seq, mse_seq, method = "spearman")

  # If the correlation is strongly negative (es. <= -0.9), the estimator is empirically consistent
  is_consistent <- !is.na(correlation) && correlation <= -0.9

  # Return result as a list
  res <- list(MSE = mse, Bias = bias, Variance = variance, consistency = is_consistent, estimates = estimates, nsim = nsim, mse = mse, bias = bias, variance = variance, true_param = true_param, mean_estimate = mean_estimate)
  class(res) <- "diagnlamp_mse"
  return(res)
}

#' @name plot_diagnlamp_mse
#' @title Plot of the output about 'simulate_estimator_mse'
#' @description This function realize a density plot of the estimator distribution, with vertical lines for the true parameter and the mean of estimates.
#' @param res List. The output of 'simulate_estimator_mse'
#' @return A plot about the distribution of the estimator, with Monte Carlo simulation.
#' @examples
#' # Given a population with mean = 400 and variance = 3:
#'
#' # True parameters
#' num_samples <- 3
#' true_variance <- 350
#'
#' generate_data <- function(n) { rnorm(n, mean = 400, sd = sqrt(350)) }
#'
#' # Calculate the expected value of the
#' # following estimator for the variance
#' # and say if it is distorted, using three samples
#' estimator_func <- function(data) {
#' (0.5 * (data[1]^2 * data[2]^2 * data[3]^3) - (data[1]*data[2]*data[3])) }
#'
#' # Run the simulation
#' results <- simulate_estimator_mse(n = 10,
#' true_param = true_variance,
#' data_gen_func = generate_data,
#' estimator_func = estimator_func, nsim = 1000)
#'
#' # Realize the plot
#' distorted_estimator <- plot_diagnlamp_mse(results)
#'
#' print(distorted_estimator)
#'
#' @note
#' The consistency check is based on the Spearman correlation between sample sizes and their corresponding MSEs. A strongly negative correlation (<= -0.9) indicates empirical consistency of the estimator.
#' @export
plot_diagnlamp_mse <- function(res) {
  nsim <- res[['nsim']]
  mse <- res[['mse']]
  bias <- res[['bias']]
  variance <- res[['variance']]
  true_param <- res[['true_param']]
  mean_estimate <- res[['mean_estimate']]

  df_plot <- data.frame(
    Estimate = res[["estimates"]],
    true_param_vector = true_param,
    mean_estimate_vector = mean_estimate
  )

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Estimate)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    ggplot2::geom_density(color = "darkblue", linewidth = 1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = true_param_vector, color = "True parameter"), linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean_estimate_vector, color = "Mean of estimations"), linetype = "dotted", linewidth = 1) +
    ggplot2::scale_color_manual(name = "Legend",
                       values = c("True parameter" = "red", "Mean of estimations" = "darkgreen")) +
    labs(title = "Distribution of the estimator (Monte Carlo simulation)",
         subtitle = sprintf("Simulations: %d | MSE: %.4f | Bias: %.4f | Variance: %.4f",
                            nsim, mse, bias, variance),
         x = "Estimated value",
         y = "Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  print(p)
}

#' @name find_method_mm
#' @title Method of Moments Estimation and Diagnostics
#' @description This function estimates a specified moment (mean, variance, skewness, or kurtosis) from the provided data using the Method of Moments (MoM). It also checks for bias and computes confidence intervals based on the selected method (Normal or Hoeffding). Additionally, it can generate a plot comparing the empirical distribution of the data with the theoretical distribution implied by the estimated moment.
#' @param data Numeric vector. The dataset from which to estimate the moment.
#' @param moment Character. The moment to estimate. Options are "mean", "variance", "skewness" or "kurtosis".
#' @param confidence_intervals Character. The method for computing confidence intervals. Options are "Normal" or "Hoeffding".
#' @param level Numeric. The confidence level for the intervals (default is 0.90).
#' @param assumption_about_data Character. Normal as default.
#' @return A list containing the estimated moment, bias status, confidence intervals, and other relevant information.
#' @examples
#'
#' # Example usage of find_method_mm
#'
#' data <- rnorm(100, mean = 5, sd = 2)
#'
#' result <- find_method_mm(data, moment = "mean", confidence_intervals = "Normal", level = 0.95)
#'
#' @note
#' The Hoeffding confidence interval requires that the data be bounded; the function uses the empirical minimum and maximum of the data as bounds.
#' @export
find_method_mm <- function(data,
                           moment = c("mean", "variance", "skewness", "kurtosis"),
                           confidence_intervals = c("Normal", "Hoeffding"),
                           level = 0.90, assumption_about_data = c("Normal", "Binomial", "Poisson", "Exponential", "Uniform")) {

  moment <- match.arg(moment)
  ci_type <- match.arg(confidence_intervals)
  assumption_about_data <- match.arg(assumption_about_data)
  n <- length(data)

  # Calculation of the sample momento
  sample_moment <- switch(moment,
                          "mean"     = mean(data),
                          "variance" = var(data),
                          "skewness" = {
                            m3 <- sum((data - mean(data))^3) / n
                            m2 <- sum((data - mean(data))^2) / n
                            m3 / (m2^(1.5))
                          },
                          "kurtosis" = {
                            m4 <- sum((data - mean(data))^4) / n
                            m2 <- sum((data - mean(data))^2) / n
                            m4 / (m2^2)
                          }
  )

  # Theorical verification of the biasedness (general asymptotic properties)
  is_biased <- switch(moment,
                      "mean"     = FALSE,
                      "variance" = FALSE,
                      "skewness" = TRUE,
                      "kurtosis" = TRUE
  )

  # Calculation of the confidence interval (based on mean)
  ci <- c(NA, NA)
  alpha <- 1 - level

  if (moment == "mean") {
    if (ci_type == "Normal") {
      se <- sd(data) / sqrt(n)
      margin <- qnorm(1 - alpha/2) * se
      ci <- c(sample_moment - margin, sample_moment + margin)
    } else if (ci_type == "Hoeffding") {
      # Hoeffding requires data limited in [a, b]
      a <- min(data)
      b <- max(data)
      R <- b - a
      margin <- R * sqrt(log(2/alpha) / (2 * n))
      ci <- c(sample_moment - margin, sample_moment + margin)
    }
  } else if (moment == "variance" && ci_type == "Normal") {
    # Asymptotic approximation for the variance based on kurtosi
    kurt <- sum((data - mean(data))^4) / (n * var(data)^2)
    se_var <- var(data) * sqrt((kurt - 1) / n)
    margin <- qnorm(1 - alpha/2) * se_var
    ci <- c(sample_moment - margin, sample_moment + margin)
  }

  # Return result
  return(list(
    data = data,
    moment = moment,
    estimate = sample_moment,
    biased = is_biased,
    ci = ci,
    ci_level = level,
    ci_type = ci_type,
    assumption_about_data = assumption_about_data
  ))
}

#' @name plot_find_method_mm
#' @title Plot of the output about 'find_method_mm'
#' @description This function realize an histogram plot about the data of the estimator distribution, with vertical lines for the true parameter and the mean of estimates.
#' @param res List. The output of 'simulate_estimator_mse'
#' @return A plot about the distribution of the estimator, with Monte Carlo simulation.
#' @examples
#'
#' data <- rnorm(100, mean = 5, sd = 2)
#'
#' result <- find_method_mm(data, moment = "mean", confidence_intervals = "Normal", level = 0.95)
#'
#' plot_find_method_mm(result)
#'
#' @export
plot_find_method_mm <- function(res) {
  data <- res[['data']]
  moment <- res[['moment']]
  sample_moment <- res[['estimate']]
  assumption_about_data <- res[['assumption_about_data']]

  if (assumption_about_data == "Normal") {
    df_plot <- data.frame(x = data)
    mu_mom <- mean(data)
    sd_mom <- sd(data)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 0.5, fill = "cadetblue2", color = "white", alpha = 0.6) +
      ggplot2::stat_function(fun = dnorm, args = list(mean = mu_mom, sd = sd_mom), color = "darkblue", linewidth = 1.2) +
      ggplot2::labs(
        title = paste("DiagnLamp: MoM adaptment (Normal assumption)"),
        subtitle = sprintf("Analyzed Moment: %s | Estimated value: %.4f", toupper(moment), sample_moment),
        x = "Value",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
    print(p)
  }

  else if (assumption_about_data == "Binomial") {
    cat("Ensure you are passing integers. However, in case of double numbers, casting is performed to the closest integer.")
    df_plot <- data.frame(x = as.integer(data))
    prob_mom <- 1- (var(as.integer(data)) / mean(as.integer(data)))
    size_mom <- (mean(as.integer(data)) / (prob_mom))
    size_mom <- as.integer(size_mom)

    x_vals <- 0:size_mom
    df_theo <- data.frame(
      x = x_vals,
      y = dbinom(x_vals, size = size_mom, prob = prob_mom)
    )

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 0.5, fill = "cadetblue2", color = "white", alpha = 0.6) +
      ggplot2::geom_point(data = df_theo, ggplot2::aes(x = x, y = y), color = 'darkblue', size = 2) +
      ggplot2::geom_line(data = df_theo, ggplot2::aes(x = x, y = y), color = 'darkblue', linewidth = 1) +
      ggplot2::labs(
        title = paste("DiagnLamp: MoM adaptment (Binomial assumption)"),
        subtitle = sprintf("Analyzed Moment: %s | Estimated value: %.4f", toupper(moment), sample_moment),
        x = "Value",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
    print(p)
  }
  else if (assumption_about_data == "Poisson") {
    cat("Ensure you are passing integers. However, in case of double numbers, casting is performed to the closest integer.")
    df_plot <- data.frame(x = as.integer(data))
    lambda_mom <- mean(data)

    df_theo <- data.frame(
      x = as.integer(data),
      y = dpois(as.integer(data), lambda = lambda_mom)
    )

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 0.5, fill = "cadetblue2", color = "white", alpha = 0.6) +
      ggplot2::geom_point(data = df_theo, ggplot2::aes(x = x, y = y), color = 'darkblue', size = 2) +
      ggplot2::geom_line(data = df_theo, ggplot2::aes(x = x, y = y), color = 'darkblue', linewidth = 1) +
      ggplot2::labs(
        title = paste("DiagnLamp: MoM adaptment (Poisson assumption)"),
        subtitle = sprintf("Analyzed Moment: %s | Estimated value: %.4f", toupper(moment), sample_moment),
        x = "Value",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
    print(p)
  }
  else if (assumption_about_data == "Exponential") {
    df_plot <- data.frame(x = data)
    rate_mom <- 1 / mean(data)

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 0.5, fill = "cadetblue2", color = "white", alpha = 0.6) +
      ggplot2::stat_function(fun = dexp, args = list(rate = rate_mom), color = "darkblue", linewidth = 1.2) +
      ggplot2::labs(
        title = paste("DiagnLamp: MoM adaptment (Exponential assumption)"),
        subtitle = sprintf("Analyzed Moment: %s | Estimated value: %.4f", toupper(moment), sample_moment),
        x = "Value",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
    print(p)
  }
  else if (assumption_about_data == "Uniform") {
    df_plot <- data.frame(x = data)
    mean_data <- mean(data)
    a_unif <- (mean_data * 2) - max(data)
    b_unif <- (mean_data * 2) - min(data)

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), binwidth = 0.5, fill = "cadetblue2", color = "white", alpha = 0.6) +
      ggplot2::stat_function(fun = dunif, args = list(min = a_unif, max = b_unif), color = "darkblue", linewidth = 1.2) +
      ggplot2::labs(
        title = paste("DiagnLamp: MoM adaptment (Uniform assumption)"),
        subtitle = sprintf("Analyzed Moment: %s | Estimated value: %.4f", toupper(moment), sample_moment),
        x = "Value",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
    print(p)
  }
}

#' @name Errf
#' @title Model Diagnostics: Empirical vs Analytical Optimism
#' @description This function evaluates the performance of a statistical model by comparing the Mean Squared Error (MSE) on training and testing datasets. It computes both empirical and analytical optimism, as well as Mallows' Cp for both datasets. Additionally, it generates a dumbell plot to visualize the generalization gap between training and testing errors.
#' @param model A fitted model object (e.g., from lm()).
#' @param train_data A data frame containing the training dataset.
#' @param test_data A data frame containing the testing dataset.
#' @param response_var A string specifying the name of the response variable in the datasets.
#' @return A list containing MSE for training and testing datasets, empirical and analytical optimism, and Mallows' Cp for both datasets, along with a dumbell plot visualizing the generalization gap.
#' @examples
#'
#' # Prepare data
#' df <- na.omit(airquality)
#' idx <- sample(1:nrow(df), size = 0.7 * nrow(df))
#' train <- df[idx, ]
#' test <- df[-idx, ]
#'
#' # Fit a linear model on the training data
#' model <- lm(Ozone ~ Wind + Temp, data = train)
#'
#' # Evaluate the model diagnostics
#' diagnostics <- Errf(model, train, test, response_var = "Ozone")
#'
#' @export
Errf <- function(model, train_data, test_data, response_var) {

  # 1. Extraction of real values of the response
  y_train <- train_data[[response_var]]
  y_test <- test_data[[response_var]]

  # 2. Prediction of the model
  pred_train <- predict(model, newdata = train_data)
  pred_test <- predict(model, newdata = test_data)

  # 3. MSE calculation
  mse_train <- mean((y_train - pred_train)^2)
  mse_test <- mean((y_test - pred_test)^2)

  # 4. Analytical and empirical optimism
  opt_empirical <- mse_test - mse_train

  n <- nrow(train_data)
  p <- length(coef(model)) # Number of estimated parameters (with intercept)
  sigma2_hat <- summary(model)$sigma^2 # Variance of residuals in training
  opt_analytical <- 2 * sigma2_hat * (p / n)

  # 5. Mallows' Cp (training and test)
  sse_train <- sum((y_train - pred_train)^2)
  cp_train <- (sse_train / sigma2_hat) - n + 2*p

  sse_test <- sum((y_test - pred_test)^2)
  cp_test <- (sse_test / sigma2_hat) - nrow(test_data) + 2*p

  # Return diagnostic data
  return(list(
    MSE_train = mse_train,
    MSE_test = mse_test,
    optimism_empirical = opt_empirical,
    optimism_analytical = opt_analytical,
    mallows_cp_train = cp_train,
    mallows_cp_test = cp_test
  ))
}

#' @name plot_Errf
#' @title Plot of the output about 'Errf'
#' @description This function realize an histogram plot about the data of the estimator distribution, with vertical lines for the true parameter and the mean of estimates.
#' @param res List. The output of 'Errf'
#' @return A plot about the distribution of the estimator, with Monte Carlo simulation.
#' @examples
#'
#' # Prepare data
#' df <- na.omit(airquality)
#' idx <- sample(1:nrow(df), size = 0.7 * nrow(df))
#' train <- df[idx, ]
#' test <- df[-idx, ]
#'
#' # Fit a linear model on the training data
#' model <- lm(Ozone ~ Wind + Temp, data = train)
#'
#' # Evaluate the model diagnostics
#' diagnostics <- Errf(model, train, test, response_var = "Ozone")
#'
#' plot_Errf(diagnostics)
#'@export
plot_Errf <- function(res) {
  MSE_train <- res[['MSE_train']]
  MSE_test <- res[['MSE_test']]
  opt_empirical <- res[['optimism_empirical']]
  opt_analytical <- res[['optimism_analytical']]

  # Dumbbell Plot
  df_dumbbell <- data.frame(
    Dataset = c("Train", "Test"),
    MSE = c(MSE_train, MSE_test)
  )

  p_opt <- ggplot2::ggplot(df_dumbbell, ggplot2::aes(x = Dataset, y = MSE)) +
    ggplot2::geom_line(ggplot2::aes(group = 1), color = "coral1", linetype = "solid", linewidth = 1.2) +
    ggplot2::geom_point(ggplot2::aes(color = Dataset), size = 6) +
    ggplot2::scale_color_manual(values = c("Train" = "chartreuse4", "Test" = "firebrick3")) +
    ggplot2::labs(
      title = "DiagnLamp: Optimism plot (Generalization Gap)",
      subtitle = sprintf("Empirical optimism: %.4f | analytical (2*sigma^2*p/n): %.4f",
                         opt_empirical, opt_analytical),
      x = "Dataset Split",
      y = "Mean Squared Error (MSE)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  print(p_opt)
}

#' @name influence_analysis
#' @title Influence Analysis of Model Points
#' @description This function analyzes the influence of individual data points on a fitted statistical model. It computes influence statistics such as leverage, Cook's distance, and DFBETAS, and identifies points that have a disproportionate impact on the parameter estimates. The function generates plots to visualize influential points and returns a detailed report.
#' @param model A fitted model object (e.g., from lm()).
#' @param data A data frame containing the dataset used for model evaluation.
#' @param response_var A string specifying the name of the response variable in the dataset.
#' @return A list containing a data frame of influence statistics, thresholds for leverage and Cook's distance, and plots visualizing influential points.
#' @examples
#'
#' # Prepare data
#' df <- na.omit(airquality)
#'
#' # Fit a linear model on the training data
#' model <- lm(Ozone ~ Wind + Temp, data = df)
#'
#' # Evaluate the model diagnostics
#' diagnostics <- influence_analysis(model, df, response_var = "Ozone")
#'
#' @note
#' Influence statistics help identify data points that may unduly affect the model's parameter estimates, which is crucial for model diagnostics and validation.
#' @export
influence_analysis <- function(model, data, response_var) {
  # Calculation of influence statistics
  influence_stats <- influence.measures(model)

  # Extraction of main statistics
  leverage <- influence_stats$infmat[, "hat"]
  cooks_distance <- influence_stats$infmat[, "cook.d"]
  dfbetas <- influence_stats$infmat[, grep("dfb", colnames(influence_stats$infmat))]

  # Dataframe with results
  results_df <- data.frame(
    Observation = seq_len(nrow(data)),
    Leverage = leverage,
    Cooks_Distance = cooks_distance,
    DFBETAS = dfbetas
  )

  # Influential points based on common points
  leverage_threshold <- 2 * mean(leverage)
  cooks_threshold <- 4 / nrow(data)

  results_df$Influential_Leverage <- leverage > leverage_threshold
  results_df$Influential_Cooks <- cooks_distance > cooks_threshold

  # Return of results as a list
  return(list(
    results_df = results_df,
    leverage_threshold = leverage_threshold,
    cooks_threshold = cooks_threshold
  ))
}

#' @name plot_influence_analysis
#' @title Plot of the output about 'influence_analysis'
#' @description This function realize two plots about the Cook's distance and leverage's distance, specifying the influential points.
#' @param res List. The output of 'influence_analysis'
#' @return Two plots about the influential points.
#' @examples
#'
#' # Prepare data
#' df <- na.omit(airquality)
#'
#' # Fit a linear model
#' model <- lm(Ozone ~ Wind + Temp, data = df)
#'
#' # Evaluate the model diagnostics
#' diagnostics <- influence_analysis(model, df, response_var = "Ozone")
#'
#' plot_influence_analysis(diagnostics)
#' @export
plot_influence_analysis <- function(res) {
  results_df <- res[['results_df']]
  leverage_threshold = res[['leverage_threshold']]
  cooks_threshold = res[['cooks_threshold']]

  # Plot to visualize the influential points
  p_leverage <- ggplot2::ggplot(results_df, ggplot2::aes(x = Observation, y = Leverage)) +
    ggplot2::geom_point(ggplot2::aes(color = Influential_Leverage), size = 2) +
    ggplot2::geom_hline(yintercept = leverage_threshold, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "Leverage Plot", x = "Observation", y = "Leverage") +
    ggplot2::theme_minimal()

  p_cooks <- ggplot2::ggplot(results_df, ggplot2::aes(x = Observation, y = Cooks_Distance)) +
    ggplot2::geom_point(ggplot2::aes(color = Influential_Cooks), size = 2) +
    ggplot2::geom_hline(yintercept = cooks_threshold, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "Cook's Distance Plot", x = "Observation", y = "Cook's Distance") +
    ggplot2::theme_minimal()

  print(p_leverage)
  print(p_cooks)
}

#' @name train_val_test_split
#' @title Dataset splitting.
#' @description Split a dataset into three popular sub-elements: training, validation and test datasets.
#' @param data Dataframe. The starting dataframe.
#' @param train_size Dataframe. The size of the training set.
#' @param val_size Dataframe. The size of the validation set.
#' @param test_size Dataframe. The size of the test set.
#' @return A list with the three sets.
#' @examples
#'
#' # Prepare data
#' df <- na.omit(airquality)
#'
#' # Split the dataset
#' res <-  train_val_test_split(df, train_size = 0.7, val_size = 0.15, test_size = 0.15)
#'
#' @export
train_val_test_split <- function(df, train_size = 0.7, val_size = 0.15, test_size = 0.15) {
  if (train_size + val_size + test_size != 1) {
    stop("The sum of sizes must be 1")
  }
  else {
  set.seed(123)
  idx <- sample(1:nrow(df), size = train_size * nrow(df))
  train <- df[idx, ]
  idx2 <- sample(1:nrow(train), size = val_size * nrow(train))
  val <- train[idx2, ]
  test <- df[-idx,]
  return(list(
    train = train,
    validation = val,
    test = test
  ))
  }
}



