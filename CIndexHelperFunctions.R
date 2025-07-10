set.seed(123)

# Lifelines wrapper for R
lifelinesR <- function(time, predicted, censoring){
  # Load python libraries
  np <- import("numpy", convert = TRUE)
  concordance_index <- import("lifelines.utils", convert = TRUE)$concordance_index
  # Calculate C-index
  concordance_index(event_times      = time,
                    predicted_scores = 1 - predicted,
                    event_observed   =  censoring)

}


# Pycox wrapper for R (Antolini)
pycoxR_SurvMatrix_Ant <- function(time, surv_matrix, censoring){
  # Load python libraries
  np       <- import("numpy", convert = FALSE)  
  pd       <- import("pandas", convert = FALSE)
  EvalSurv <- import("pycox.evaluation", convert = FALSE)$EvalSurv
  # Library to calculate separately all pairs
  pc_eval  <- import("pycox.evaluation.concordance")
  
  # Convert R objects to python with correct type
  time_py <- np$array(as.numeric(time), dtype = "float64")  
  censoring_py <- np$array(as.integer(censoring), dtype = "int32")
  # Transpose the survival matrix and convert to python and correct type
  surv_mat <- t(as.matrix(surv_matrix))  # now: rows = times, cols = patients
  surv_np <- np$array(surv_mat, dtype = "float64")
  time_index <- np$array(as.numeric(colnames(surv_matrix)), dtype = "float64") # time grid
  
  # Create pandas dataframe
  surv_df <- pd$DataFrame(data = surv_np, index = time_index)
  
  # Very important and quite problematic!! 
  # censor_surv = "km" does not adjust for censoring. 
  # censor_surv is made for other performance metrics, but does NOT apply to c-index
  # Leaving it as NULL, as it makes no difference than default or "km"
  # Build the evaluation object
  ev <- EvalSurv(surv_df, time_py, censoring_py, censor_surv = NULL)
  
  # Calculate C-index and convert back to R objects
  cindex <- py_to_r(ev$concordance_td(method = "antolini"))
  
  # Calculate the concordant and comparable pairs separately to analyse differences
  surv_idx <- np$searchsorted(time_index, time_py)
  num_concordant <- pc_eval$`_sum_concordant_disc`(surv_np, time_py, censoring_py, 
                                                   surv_idx, pc_eval$`_is_concordant_antolini`)
  num_comparable <- pc_eval$`_sum_comparable`(time_py, censoring_py, 
                                              pc_eval$`_is_comparable_antolini`)
  # Convert them back to R objects
  num_concordant <- py_to_r(num_concordant)
  num_comparable <- py_to_r(num_comparable)
  
  return(list(
    cindex = cindex,
    concordant = num_concordant,
    comparable = num_comparable
  ))
}

# Pycox wrapper for R (Antolini modified version)
pycoxR_SurvMatrix_AdjAnt <- function(time, surv_matrix, censoring){
  # Load python libraries
  np <- import("numpy", convert = FALSE)  
  pd <- import("pandas", convert = FALSE)
  EvalSurv <- import("pycox.evaluation", convert = FALSE)$EvalSurv
  # Library to calculate separately all pairs
  pc_eval <- import("pycox.evaluation.concordance")
  
  # Convert R objects to python with correct type
  time_py <- np$array(as.numeric(time), dtype = "float64")  
  censoring_py <- np$array(as.integer(censoring), dtype = "int32")
  # Transpose the survival matrix and convert to python and correct type
  surv_mat <- t(as.matrix(surv_matrix))  # now: rows = times, cols = patients
  surv_np <- np$array(surv_mat, dtype = "float64")
  
  # Transpose the survival matrix and convert to python and correct type
  time_index <- np$array(as.numeric(colnames(surv_matrix)), dtype = "float64") # time grid
  
  # Create pandas dataframe
  surv_df <- pd$DataFrame(data = surv_np, index = time_index)
  
  # Very important and quite problematic!! 
  # censor_surv = "km" does not adjust for censoring. 
  # censor_surv is made for other performance metrics, but does NOT apply to c-index
  # Leaving it as NULL, as it makes no difference than default or "km"
  # Build the evaluation object
  ev <- EvalSurv(surv_df, time_py, censoring_py, censor_surv = NULL)
  
  # Calculate C-index and convert back to R objects
  cindex <- py_to_r(ev$concordance_td(method = "adj_antolini"))
  
  # Calculate the concordant and comparable pairs separately to analyse differences
  surv_idx <- np$searchsorted(time_index, time_py)
  num_concordant <- pc_eval$`_sum_concordant_disc`(surv_np, time_py, censoring_py, 
                                                   surv_idx, pc_eval$`_is_concordant`)
  num_comparable <- pc_eval$`_sum_comparable`(time_py, censoring_py, 
                                              pc_eval$`_is_comparable`)
  
  # Convert them back to R objects
  num_concordant <- py_to_r(num_concordant)
  num_comparable <- py_to_r(num_comparable)
  
  return(list(
    cindex = cindex,
    concordant = num_concordant,
    comparable = num_comparable
  ))
}



# Python wrapper for sksurv censored version
sksurv.censoredR <- function(time, predicted, censoring) {
  # Load python libraries
  np     <- import("numpy", convert = FALSE)
  sksurv <- import("sksurv.metrics", convert = FALSE)
  # Convert R objects into R
  time_py      <- np$array(as.numeric(time), dtype = "float64")  
  censoring_py <- np$array(as.logical(censoring), dtype = "bool")  
  predicted_py <- np$array(as.numeric(predicted), dtype = "float64")
  
  # Calculate C-index
  c_index <- sksurv$concordance_index_censored(
    event_indicator = censoring_py,  
    event_time = time_py,
    estimate = predicted_py
  )
  # Return as R
  py_to_r(c_index)
  
}

# Python wrapper for sksurv ipcw version with default settings
sksurv.ipcwR <- function(time, predicted, censoring, 
                         eval.times, 
                         sksurv_train_time, 
                         sksurv_train_status, 
                         sksurv_tied_tol = 1e-8) { # Default

  # Load libraries
  np     <- import("numpy", convert = FALSE)
  sksurv <- import("sksurv.metrics", convert = FALSE)
  # Convert time to R
  time_train_py <- np$array(as.numeric(sksurv_train_time), dtype = "float64")
  time_test_py  <- np$array(as.numeric(time), dtype = "float64")
  # Convert status to R (event = 1, censored = 0)
  censoring_train_py <- np$array(as.logical(sksurv_train_status), dtype = "bool")
  censoring_test_py  <- np$array(as.logical(censoring), dtype = "bool")  
  # Convert predictions to R
  predicted_py <- np$array(as.numeric(predicted), dtype = "float64")
  
  # Create survival objects
  survival_dtype <- np$dtype(list(
    tuple("event", np$bool_),
    tuple("time", np$float64)
  ))
  
  train_survival_py <- np$empty(length(time_train_py), dtype = survival_dtype)
  test_survival_py <- np$empty(length(time_test_py), dtype = survival_dtype)
  
  train_survival_py["event"] <- censoring_train_py 
  train_survival_py["time"] <- time_train_py
  
  test_survival_py["event"] <- censoring_test_py
  test_survival_py["time"] <- time_test_py
  
  # Compute IPCW-adjusted C-index
  c_index_ipcw <- sksurv$concordance_index_ipcw(
    survival_train = train_survival_py,
    survival_test = test_survival_py,
    estimate = predicted_py,
    tau = eval.times,
    tied_tol = sksurv_tied_tol
  )
  
  py_to_r(c_index_ipcw)
}

# ---- Bootstrapping C-index -----
# The bootstrap is done based on a data subset (test set)
# 100 bootstrap re-samples are generated from it
bootstrap.metric <- function(metrics.wrapper, 
                             dataset, 
                             implementation, 
                             eval.times = NULL, 
                             sampled_data = NULL,
                             additional = NULL) { 
  # Empty lists length of number of boostraps
  iterations <- vector("list", length(sampled_data))
  # Number of boostraps
  N_bootstraps <- length(sampled_data)
  # Loop through the sampled_data
  for (i in seq_along(sampled_data)) {
    # Subset each sampled data from list
    resampled_data <- sampled_data[[i]]
    # If the implementation is Ctd (antolinis) the input is the survival matrix
    if (("pycox.Ant" %in% implementation) | ("pycox.Adj.Ant" %in% implementation)) {
      # Get status and time
      data <- data.frame(
        censoring = resampled_data[[dataset$censoring]],
        time = resampled_data[[dataset$time]]
      )
      # Subset the correct model matrix
      surv_matrix = resampled_data[, grep(paste0("^", dataset$predicted, "\\."), 
                                          names(resampled_data), value = TRUE), drop = FALSE]
      # Subtract only the time points as column names
      colnames(surv_matrix) <- as.numeric(sub(".*\\.", "", colnames(surv_matrix)))
      # Run the metrics wrapper to calculate C-index
      iterations[[i]] <- do.call(metrics.wrapper, c(data, list(surv_matrix = surv_matrix), 
                                                    list(implementation = implementation), 
                                                    eval.times, additional = NULL)) 
    # Other C and C tau implementations
    } else {
      # Create data with status, time and predictions (RMST or EM)
      data <- data.frame(
        predicted = resampled_data[[dataset$predicted]],
        censoring = resampled_data[[dataset$censoring]],
        time = resampled_data[[dataset$time]]
      )
      # Run the metrics wrapper to calculate C-index
      iterations[[i]] <- do.call(metrics.wrapper, 
                                 c(data, list(implementation = implementation), 
                                   eval.times, additional)) 
    }
  }
  # Aggregate metrics
  metrics_df <- do.call(rbind, lapply(iterations, function(x) {
    as.data.frame(as.list(x), check.names = FALSE)  # Keep "::" in column names
  }))
  
  # Calculate mean and 95% confidence intervals for each metric
  subset_metrics <- metrics_df[, colnames(metrics_df) %in% implementation, drop = FALSE] 
  mean_metrics <- colMeans(subset_metrics, na.rm = TRUE)
  # conf_intervals <- t(apply(subset_metrics, 2, function(metric) {
  #   mean_metric <- mean(metric, na.rm = TRUE)
  #   stderr <- sd(metric, na.rm = TRUE) / sqrt(N_bootstraps)
  #   mean_metric + qt(c(0.025, 0.975), df = N_bootstraps - 1) * stderr
  # }))
  conf_intervals <- t(apply(subset_metrics, 2, function(metric) {
    quantile(metric, probs = c(0.025, 0.975), na.rm = TRUE)
  }))
  
  return(list(
    mean = mean_metrics,
    confidence.intervals = conf_intervals,
    eval.times = eval.times,
    batch.metrics = metrics_df
  ))
}

# ---- Bootstrapping C-index in Parallel-----
# The bootstrap is done based on a data subset (test set)
# 100 bootstrap re-samples are generated from it
# Parallel computation
bootstrap.metric.parallel <- function(metrics.wrapper, 
                             dataset, 
                             implementation, 
                             eval.times = NULL,
                             sampled_data = NULL, 
                             additional = NULL) {
  # Empty lists length of number of boostraps
  iterations <- vector("list", length(sampled_data))
  # Calculate the number of bootstraps
  N_bootstraps <- length(sampled_data)
  # Create empty lists
  resampled_data <- vector("list", length = N_bootstraps)
  surv_matrixes <- vector("list", length = N_bootstraps)
  # Loop through the sampled_data
  for (i in seq_along(sampled_data)) {
    # Subset each sampled data from list
    sampled_data_ <- sampled_data[[i]]
    # If the implementation is Ctd (antolinis) the input is the survival matrix
    if (("pycox.Ant" %in% implementation) | ("pycox.Adj.Ant" %in% implementation)) {
      resampled_data[[i]] <- data.frame(
        censoring = sampled_data_[[dataset$censoring]],
        time = sampled_data_[[dataset$time]]
      )
      # Just take the survival matrix of the model
      surv_mat = sampled_data_[, grep(dataset$predicted, names(sampled_data_), value = TRUE), drop = FALSE]
      # Columns are reduce to time grid
      colnames(surv_mat) <- as.numeric(sub(".*\\.", "", colnames(surv_mat)))
      # Save matrixes
      surv_matrixes[[i]] <- surv_mat
    # Other C and C tau implementations
    } else {
      # Create data with status, time and predictions (RMST or EM)
      resampled_data[[i]] <- data.frame(
        predicted = sampled_data_[[dataset$predicted]],
        censoring = sampled_data_[[dataset$censoring]],
        time = sampled_data_[[dataset$time]]
      )
    }
  }
  # Parallel function
  parallel_func <- function(resampled_data, implementation, eval.times, surv_matrixes, additional) { 
    # build progress bar
    pro_bar <- progressor(along = resampled_data)
    # Run foreach
    results <- suppressWarnings(foreach(i = seq_along(resampled_data),
                       data = resampled_data,
                       .combine = rbind,  ## Ensures all results are collected
                       .options.future = list(seed = TRUE)) %dofuture% {
                         # Update progress bar
                         pro_bar(sprintf("Iteration %d of %d", i, length(resampled_data)))
                         ## Call function with proper arguments
                         do.call(metrics.wrapper, c(data, list(surv_matrix = surv_matrixes[[i]]),
                                                    list(implementation = implementation),
                                                    eval.times, additional))
                       })
    
  }
  
  # Set the parallel backend to future package
  registerDoFuture() 
  
  plan(multisession) # Linux/Mac/Windows
  #plan(multicore) # Linux/Mac
  #plan(sequential)

  # Set progress bar style
  handlers("cli")
  # Run the parallel function
  metrics_df <- with_progress(parallel_func(resampled_data,
                                            implementation, 
                                            eval.times, 
                                            surv_matrixes, additional))

  # Remove parallel index 
  rownames(metrics_df) <- NULL

  # Calculate mean and 95% confidence intervals for each metric
  subset_metrics <- metrics_df[, colnames(metrics_df) %in% implementation, drop = FALSE] 
  mean_metrics <- apply(subset_metrics, 2, function(x) mean(as.numeric(x), na.rm = TRUE))
  # conf_intervals <- t(apply(subset_metrics, 2, function(metric) {
  #   mean_metric <- mean(as.numeric(metric), na.rm = TRUE)
  #   stderr <- sd(as.numeric(metric), na.rm = TRUE) / sqrt(N_bootstraps)
  #   mean_metric + qt(c(0.025, 0.975), df = N_bootstraps - 1) * stderr
  conf_intervals <- t(apply(subset_metrics, 2, function(metric) {
    quantile(as.numeric(metric), probs = c(0.025, 0.975), na.rm = TRUE)
  }))
  
  return(list(
    mean = mean_metrics,
    confidence.intervals = conf_intervals,
    eval.times = eval.times,
    batch.metrics = metrics_df
  ))
}

# ---- C-index implementations-----
# All the implementations in python and R
metrics.wrapper <- function(predicted, surv_matrix = NULL, 
                            censoring, time, 
                            implementation, eval.times, 
                            sksurv_train_time = NULL, # sksurv.ipcw
                            sksurv_train_status = NULL, # sksurv.ipcw
                            sksurv_tied_tol = NULL) { # sksurv.ipcw

  # Initialize an empty array for the results
  results <- list()
  # Initialize eval.times list
  batch.eval.times <- c()
  # Flatten the implementation list into a single vector
  implementation <- unlist(implementation)
  
  # If eval.times is missing, set to the maximum uncensoring time as in pec
  # The maximum uncensored time is when the C-index stops updating
  if (missing(eval.times) || is.null(eval.times) || (eval.times == "resample_max_uncensored_time")) {
    eval.times <- max(time[censoring == 1], na.rm=TRUE) 
  }
  # Store the eval.times for the results
  results$batch.eval.times <- c(results$batch.eval.times, eval.times)
  
  # For pec::cindex
  if ("pec::cindex" %in% implementation) {
    # Create a temporary dataframe
    tmp <- data.frame(predicted, 
                      censoring, 
                      time)

    pec <- pec::cindex(
      list("res" = as.matrix(1 - predicted)),
      formula = Hist(time, censoring) ~ 1, 
      data = tmp, 
      eval.times = eval.times,
      cens.model = "marginal", # default marginal
    )
    # Store Cindex
    results[["pec::cindex"]] <- pec$AppCindex$res
    # Store elegible pairs
    results[["pec::cindex.EP"]] <- pec$Pairs$res
    # Store concordant pairs
    results[["pec::cindex.CP"]] <- pec$Concordant$res
    
  }
  

  # Uno's real implementation: 
  # https://github.com/cran/survC1/blob/master/R/FUN-cstat-ver003b.R
  # https://github.com/cran/survC1/blob/master/src/husurvC1v1.f
  if ("survC1::Est.Cval" %in% implementation) {
    # Need to create this survival object: 
    results[["survC1::Est.Cval"]] <- 
      survC1::Est.Cval(data.frame(time, 
                                  censoring, 
                                  predicted), # Expects risk scores
                       tau=eval.times, 
                       nofit=TRUE)$Dhat 
  }
  
  # For Harrel's real implementation
  if ("Hmisc::rcorr.cens" %in% implementation) {

    hmisc <- Hmisc::rcorr.cens(x = 1 - predicted, 
                               S = Surv(time, censoring))
    
    results[["Hmisc::rcorr.cens"]] <- hmisc[["C Index"]]
    results[["Hmisc::rcorr.cens.EP"]] <- hmisc[["Relevant Pairs"]]
    results[["Hmisc::rcorr.cens.CP"]] <- hmisc[["Concordant"]]
    #results[["Hmisc::rcorr.cens.CP"]] <- hmisc[["Uncertain"]]

  }
  
  # For SurvMetrics::Cindex
  if ("SurvMetrics::Cindex" %in% implementation) {
    # Need to create this survival object: 
    surv_obj = Surv(time, censoring)
    # Calculate the metric
    # t_star is only used to extract surv at a specific time point
    # not to calculate 
    results[["SurvMetrics::Cindex"]] <-
      SurvMetrics::Cindex(surv_obj, predicted=1-predicted, 
                          t_star = eval.times)[[1]]

  }
  
  # For randomForestSRC::get.cindex 
  if ("randomForestSRC::get.cindex" %in% implementation) {
    results[["randomForestSRC::get.cindex"]] <- 
      randomForestSRC::get.cindex(predicted = 1 - predicted,
                                  censoring = censoring,
                                  time = time) 
  }

  # lifelines::concordance_index implementation
  if ("lifelines" %in% implementation) {
    # Handle Python implementations
    results[["lifelines"]] <- lifelinesR(time, predicted, censoring)
    
  }
  # https://github.com/havakv/pycox/blob/master/pycox/evaluation/concordance.py#L12
  # using the eval.times only to subset the survival probabilities at a time point
  if ("pycox.Ant" %in% implementation) {
    # Handle Python implementations
    #results[["pycox"]] <- pycoxR(time, predicted, censoring, eval.times)
    pycox_ant <- pycoxR_SurvMatrix_Ant(time=time, 
                                       surv_matrix=surv_matrix, 
                                       censoring=censoring)
    
    results[["pycox.Ant"]] <- pycox_ant$cindex
    results[["pycox.Ant.EP"]] <- pycox_ant$comparable
    results[["pycox.Ant.CP"]] <- pycox_ant$concordant
  }
  if ("pycox.Adj.Ant" %in% implementation) {
    # Handle Python implementations
    #results[["pycox"]] <- pycoxR(time, predicted, censoring, eval.times)
    pycox_adj <- pycoxR_SurvMatrix_AdjAnt(time=time, 
                                          surv_matrix=surv_matrix, 
                                          censoring=censoring)
    
    results[["pycox.Adj.Ant"]] <- pycox_adj$cindex
    results[["pycox.Adj.Ant.EP"]] <- pycox_adj$comparable
    results[["pycox.Adj.Ant.CP"]] <- pycox_adj$concordant
  }
  
  # For pysurvival we have extracted the function of interest in file pysurvivalR.tar.gz
  # https://github.com/square/pysurvival/blob/master/pysurvival/cpp_extensions/metrics.cpp
  if ("pysurvival" %in% implementation){
    # R wrapped version of python/c++ library pysurvival

    pysurv <- pysurvivalR::concordance_index(risk=predicted,
                                             T=time,
                                             E=censoring,
                                             include_ties = TRUE) # risk ties default\
    
    results[["pysurvival"]] <- pysurv$`0` 
    results[["pysulvival.EP"]] <- pysurv$`1` 
    results[["pysulvival.CP"]] <- pysurv$`2` 
    
  }
  
  # For sksurv without ipcw
  if ("sksurv.censored" %in% implementation) { # Like Harrels
    sksurv <- sksurv.censoredR(time,
                               predicted, 
                               censoring)
    
    results[["sksurv.censored"]] <- sksurv[[1]] # C index
    results[["sksurv.censored.CP"]] <- sksurv[[2]] # Concordant
    results[["sksurv.censored.DP"]] <- sksurv[[3]] # Disconcordant
    results[["sksurv.censored.TR"]] <- sksurv[[4]] # tied risk
    results[["sksurv.censored.TT"]] <- sksurv[[5]] # tied time
  }
  
  # For sksurv with ipcw
  if ("sksurv.ipcw" %in% implementation) { # Like Uno's
    sksurv <- sksurv.ipcwR(time,
                           predicted,
                           censoring,
                           eval.times,
                           sksurv_train_time, 
                           sksurv_train_status, 
                           sksurv_tied_tol)

    results[["sksurv.ipcw"]] <- sksurv[[1]] # C index
    results[["sksurv.ipcw.CP"]] <- sksurv[[2]] # Concordant
    results[["sksurv.ipcw.DP"]] <- sksurv[[3]] # Disconcordant
    results[["sksurv.ipcw.TR"]] <- sksurv[[4]] # tied risk
    results[["sksurv.ipcw.TT"]] <- sksurv[[5]] # tied time
  }
  
  # For survival with IPCW
  if ("survival.n/G2" %in% implementation) {
    
    predicted = 1- predicted 
    surv.conc <- survival::concordance(Surv(time, censoring) ~ predicted, 
                                                        timewt ="n/G2", ymax = eval.times)
    results[["survival.n/G2"]] <- surv.conc$concordance
    results[["survival.n/G2.CP"]]  <- surv.conc$count[["concordant"]] # Concordanct
    results[["survival.n/G2.DP"]]  <- surv.conc$count[["discordant"]] # Discordant
    results[["survival.n/G2.TR"]]  <- surv.conc$count[["tied.x"]] # Tied risk
    results[["survival.n/G2.TT"]]  <- surv.conc$count[["tied.y"]] # Tied time
    results[["survival.n/G2.TB"]]  <- surv.conc$count[["tied.xy"]] # Tied both time and risk
    
  }
  
  # For survival without IPCW
  if ("survival.n" %in% implementation) {
    
    predicted = predicted # expects survival? confusing....
    surv.conc <- survival::concordance(Surv(time, censoring) ~ predicted, 
                                       timewt = "n", ymax = eval.times)
    results[["survival.n"]] <- surv.conc$concordance
    results[["survival.n.CP"]]  <- surv.conc$count[["concordant"]] # Concordanct
    results[["survival.n.DP"]]  <- surv.conc$count[["discordant"]] # Discordant
    results[["survival.n.TR"]]  <- surv.conc$count[["tied.x"]] # Tied risk
    results[["survival.n.TT"]]  <- surv.conc$count[["tied.y"]] # Tied time
    results[["survival.n.TB"]]  <- surv.conc$count[["tied.xy"]] # Tied both time and risk
    
  }
  
  if (length(results) == 0) {
    stop("Missing implementation list.
         It is required a list with all the c-index 
         implementations that need to be computed")
  }
  
  return(results)
  
}


# Function to extract confidence intervals 
# from lambda in Weibull PH models
get_lambda_CI <- function(model, level = 0.95) {
  sumry <- summary(model)
  mu <- sumry$coefficients[1]
  sigma <- sumry$scale
  vcov <- sumry$var
  
  k <- length(sumry$coefficients) - 1
  idx_mu <- 1
  idx_logsigma <- k + 2
  se_logsigma <- sumry$table["Log(scale)", "Std. Error"]
  
  var_mu <- vcov[idx_mu, idx_mu]
  var_sigma <- (sigma^2) * (se_logsigma^2)
  cov_mu_logsigma <- vcov[idx_mu, idx_logsigma]
  cov_mu_sigma <- cov_mu_logsigma * sigma
  
  lambda <- exp(-mu / sigma)
  var_lambda <- lambda^2 * (
    (var_mu / sigma^2) -
      (2 * mu / sigma^3) * cov_mu_sigma +
      (mu^2 / sigma^4) * var_sigma
  )
  se_lambda <- sqrt(var_lambda)
  
  z <- qnorm(1 - (1 - level)/2)
  ci <- c(lambda - z * se_lambda, lambda + z * se_lambda)
  return(list(lambda = lambda, SE = se_lambda, CI = unname(ci)))
}


# ---- Generate synthetic data -----
# Syntehtic event times follow Weibull PH
generate_synthetic_event_times <- function(survreg_model, 
                                           n = 1000, 
                                           covariates, 
                                           seed = 123, 
                                           verbose = FALSE) {
  # Set seed
  set.seed(seed)
  # Number of patients
  total_n <- nrow(covariates)
  # If the user given number of patients is less that the total
  # Then shuffle and sample from the the covariates
  if (n < total_n) {
    covariates <- covariates[sample(1:total_n, n, replace = FALSE), ]
  # Otherwise no shuffling or sampling needed
  } else {
    n <- total_n 
  }

  # Get parameters from the original model (this is AFT scale)
  mu <- survreg_model$coefficients[1]  # Intercept
  sigma <- survreg_model$scale  # Scale
  alpha <- survreg_model$coefficients[-1] # Regression coefficients
  
  # PARAMETER TRANSFORMATION
  # Covert into WeibullPH:
  lambda <- exp(-mu/sigma)
  gamma <- 1/sigma
  beta <- -alpha/sigma
  
  # Calculate the linear predictor
  betas_cov <- as.matrix(covariates) %*% beta
  
  # Simulate survival times
  U <- runif(n)
  T_s <- (-log(U) / (lambda * exp(betas_cov)))^(1/gamma)
  
  if (verbose) {
    cat("Lambda T", lambda, "\n")
    cat("Gamma T", gamma, "\n")
    cat("Betas", beta, "\n")
  }
  
  # Generate status based on censoring
  status_obs <- rep(1, times = n) # <= or < ?
  
  # Calculate percentage of censoring
  censoring_percentage <- mean(status_obs == 0) * 100
  
  # Create output data.frame
  synth_data <- cbind(
    data.frame("time" = T_s, "censoring_time" = rep(Inf, n), 
               "observed_time" = T_s, "status" = status_obs), covariates)
  
  # Set the attributes
  attr(synth_data, "lambda0") <- lambda
  attr(synth_data, "gamma")   <- gamma
  attr(synth_data, "beta")    <- beta
  attr(synth_data, "betas_cov") <- betas_cov # linear predictor
  attr(synth_data, "seed")    <- seed
  attr(synth_data, "mu")      <- mu
  attr(synth_data, "sigma")   <- sigma
  attr(synth_data, "alpha")   <- alpha
  
  # Censoring percentage as attribute
  attr(synth_data, "censoring_percentage") <- censoring_percentage
  
  return(synth_data)
}

# ---- Adding synthetic censoring -----
# Censoring can be: administrative, 
# weibull PH (informative or not), 
# and uniform
add_censoring <- function(my_surv_data,
                          cens.type = "administrative",
                          cens.params = list(
                            cens_limit_admin = NULL, # used only for admin cens
                            min_cens_unif = NULL, # used only for uniform cens
                            max_cens_unif = NULL, # used only for uniform cens
                            cens_increase_unif = NULL,  # used only for uniform cens
                            weibull_cens_model = NULL, # used only for weibull cens
                            lambda_c_factor = NULL, # used only for weibull cens
                            covariates = NULL ## when informative
                          ),
                          seed = 123) {
  # Set seed
  set.seed(seed)
  # Get sample size
  n <- nrow(my_surv_data)
  
  # Administrative censoring
  if(cens.type == "administrative") {
    censoring_times <- rep(cens.param$cens_limit_admin, times = n)
  }
  # Uniform censroing
  if(cens.type == "uniform") {
    # Lower bound of Uniform
    T_s <- my_surv_data$time
    # Lower bound is the minimum of event times
    C_min <- min(T_s)
    # Upper bound can be modified
    C_max <- quantile(T_s, 1 - cens.params$cens_increase_unif)
    # Generate censoring times
    censoring_times  <- runif(n, min = C_min, max = C_max)
    # It can also have administrative
    if (!is.null(cens.params$cens_limit_admin)) {
      censoring_times <- pmin(censoring_times, 
                              cens.params$cens_limit_admin)
    }
  }
  # Weibull censoring
  if(cens.type == "weibull") {
    # Get parameter values
    weibull_cens_model <- cens.params$weibull_cens_model
    lambda_c_factor <- cens.params$lambda_c_factor
    
    mu_c    <-  weibull_cens_model$coefficients[1]  # Intercept
    sigma_c <-  weibull_cens_model$scale  # Scale
    #alpha_c <-  weibull_cens_model$coefficients[-1] # Regression coefficients
    
    # Transform to WeibullPH
    lambda_c <- exp(-mu_c/sigma_c) * lambda_c_factor
    gamma_c <- 1/sigma_c
    #beta_c <- -alpha_c/sigma_c

    # Generate censoring times
    U_c <- runif(n)
    censoring_times <- (-log(U_c) / lambda_c)^(1/gamma_c)
    
    # add administrative censoring also 
    if (!is.null(cens.params$cens_limit_admin)) {
      censoring_times <- pmin(censoring_times, 
                              cens.params$cens_limit_admin)
    }
    
  } 
  # If Weibull informative
  if(cens.type == "informative") {
    # Get parameter values
    weibull_cens_model <- cens.params$weibull_cens_model
    lambda_c_factor <- cens.params$lambda_c_factor
    
    mu_c    <-  weibull_cens_model$coefficients[1]  # Intercept
    sigma_c <-  weibull_cens_model$scale  # Scale
    alpha_c <-  weibull_cens_model$coefficients[-1] # Regression coefficients
    
    # Transform to WeibullPH
    lambda_c <- exp(-mu_c/sigma_c) * lambda_c_factor
    gamma_c <- 1/sigma_c
    beta_c <- -alpha_c/sigma_c
    
    # Calculate the linear predictor
    betas_cov_c <- as.matrix(cens.params$covariates) %*% beta_c

    # Generate censoring times
    U_c <- runif(n)
    censoring_times <- (-log(U_c) /(lambda_c * exp(betas_cov_c)))^(1/gamma_c)
    
    # add administrative censoring also 
    if (!is.null(cens.params$cens_limit_admin)) {
      censoring_times <- pmin(censoring_times, 
                              cens.params$cens_limit_admin)
    }
    
  }

  # Create a new data frame with observed times and status
  my_surv_data$censoring_time <- censoring_times
  my_surv_data$observed_time <- pmin(my_surv_data$time, my_surv_data$censoring_time)
  my_surv_data$status <- ifelse(my_surv_data$time <= my_surv_data$censoring_time, 1, 0)
  
  return(my_surv_data)
  
}

# Easy load of synthetic datasets
load_synthetic_datasets <- function(file_path) {
  # Load the structured list
  loaded_list <- readRDS(file_path)
  
  # Extract datasets and attributes
  datasets <- loaded_list$data
  attributes_list <- loaded_list$attributes
  
  # Store the attributes 
  for (i in seq_along(datasets)) {
    for (j in seq_along(datasets[[i]])) {
      attributes(datasets[[i]][[j]]) <- attributes_list[[i]][[j]]
    }
  }
  
  return(datasets)
}

# ---- Extraction of model predictions -----
get_model_preds2 <- function(stacked_predictions, 
                             model_names = "all", 
                             input_type = c("Distribution", "ExpectedMortality", "RiskAtT", "RMST"),
                             specific_time = NULL, 
                             bootstrap_patient_ids = NULL) {
  # Extract the input transformation
  input_type <- match.arg(input_type)
  base_cols <- c("test_time", "test_status", "patients_ids")
  df <- stacked_predictions  # work on a copy

  # Subset rows by bootstrap_patient_ids
  if (!is.null(bootstrap_patient_ids)) {
    match_ids <- match(bootstrap_patient_ids, df$patients_ids)
    df <- df[match_ids, , drop = FALSE]
  }
  
  # Detect model names if "all" (but we might only want to do this for 1 model)
  if (identical(model_names, "all")) {
    all_cols <- names(df)
    if (input_type %in% c("Distribution", "RiskAtT")) {
      model_names <- unique(sub("\\..*", "", grep("^[A-Za-z]+\\.\\d+$", all_cols, value = TRUE)))
    } else if (input_type == "ExpectedMortality") {
      model_names <- unique(sub("Exp\\Mort\\.", "", grep("^Exp\\Mort\\.", all_cols, value = TRUE)))
    } else if (input_type == "RMST") {
      model_names <- unique(sub("RMST\\.", "", grep("^RMST\\.", all_cols, value = TRUE)))
    }
  }
  
  # Store final prediction columns
  all_model_cols <- list()
  
  # Go through each model and input type
  for (model in model_names) {
    # Survival probability distribution as input
    if (input_type == "Distribution") {
      pattern <- paste0("^", model, "\\.")
      cols <- grep(pattern, names(df), value = TRUE)
      all_model_cols[[length(all_model_cols) + 1]] <- df[, cols, drop = FALSE]
    
    # Expected mortality
    } else if (input_type == "ExpectedMortality") {
      col_name <- paste0("ExpMort.", model)
      if (!col_name %in% names(df)) {
        stop("Missing column: ", col_name)
      }
      all_model_cols[[model]] <- df[, col_name, drop = FALSE]
      
    # Risk (1 - S(t)) at specific time t
    } else if (input_type == "RiskAtT") {
      if (is.null(specific_time)) {
        stop("specific_time must be provided for input_type = 'RiskAtT'")
      }
      cols <- paste0(model, ".", specific_time)
      missing_cols <- setdiff(cols, names(df))
      if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "))
      }
      # Transform survival to risk (1-S(t))
      transformed <- lapply(cols, function(col) 1 - df[[col]])
      names(transformed) <- cols
      all_model_cols[[length(all_model_cols) + 1]] <- as.data.frame(transformed)
    
    # RMST
    } else if (input_type == "RMST") {
      col_name <- paste0("RMST.", model)
      if (!col_name %in% names(df)) {
        stop("Missing column: ", col_name)
      }
      all_model_cols[[model]] <- df[, col_name, drop = FALSE]
    }
  }
  # Combine everything to return
  predictions <- do.call(cbind, all_model_cols)
  final_df <- cbind(df[, base_cols, drop = FALSE], predictions)
  
  return(final_df)
}

# ---- Extract info from bootstrap and point estimates -----
# When calculated with risk at t 
make_risk_plot_entries <- function(results_list, point_estimate) {
  # Extract each result and convert to data.frame rows
  risk_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    model <- sub("\\..*", "", name)
    t <- sub(".*?\\.", "", name)
    tau <- result$eval.times
    metrics <- names(result$mean)
    
    entries <- lapply(seq_along(metrics), function(i) {
      metric_name <- metrics[i]
      mean_c <- point_estimate[[name]][[metric_name]]
      ci <- result$confidence.intervals[i, ]
      
      data.frame(
        Time = t,
        Metric = metric_name,
        Model = model,
        InputType = paste0("Risk(t=", t, ", tau=", tau, ")"),
        cindex = mean_c,
        lower = ci[1],
        upper = ci[2],
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, entries)
  })
  
  # Combine all data frames into one
  do.call(rbind, risk_entries)
}

# ---- Extract info from bootstrap and point estimates -----
# When calculated with Expected mortality
make_expm_plot_entries <- function(results_list, point_estimate) {
  expm_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    model <- sub("ExpMort\\.", "", name)
    metrics <- names(result$mean)
    
    entries <- lapply(seq_along(metrics), function(i) {
      metric_name <- metrics[i]
      #mean_c <- result$mean[i] # mean of bootstrap
      mean_c <- point_estimate[[name]][[metric_name]]
      ci <- result$confidence.intervals[i, ]

      data.frame(
        Metric = metric_name,
        Model = model,
        cindex = mean_c,
        lower = ci[1],
        upper = ci[2],
        InputType = "Exp.Mort",
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, entries)
  })
  
  do.call(rbind, expm_entries)
}
# ---- Extract info from bootstrap and point estimates -----
# When calculated with survival probabilities
make_surv_plot_entries <- function(results_list, point_estimate) {
  surv_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    metrics <- names(result$mean)
    
    entries <- lapply(seq_along(metrics), function(i) {
      metric_name <- metrics[i]
      #mean_c <- result$mean[i] # mean of bootstrap
      mean_c <- point_estimate[[name]][[metric_name]]
      ci <- result$confidence.intervals[i, ]
      
      data.frame(
        Metric = metric_name,
        Model = name,
        cindex = mean_c,
        lower = ci[1],
        upper = ci[2],
        InputType = "Distrib.",
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, entries)
  })
  
  do.call(rbind, surv_entries)
}

# ---- Extract info from bootstrap and point estimates -----
# When calculated with RMST
make_rmst_plot_entries <- function(results_list, point_estimate) {
  expm_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    model <- sub("RMST\\.", "", name)
    metrics <- names(result$mean)
    
    entries <- lapply(seq_along(metrics), function(i) {
      metric_name <- metrics[i]
      #mean_c <- result$mean[i] # mean of bootstrap
      mean_c <- point_estimate[[name]][[metric_name]]
      ci <- result$confidence.intervals[i, ]
      
      data.frame(
        Metric = metric_name,
        Model = model,
        cindex = mean_c,
        lower = ci[1],
        upper = ci[2],
        InputType = "RMST",
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, entries)
  })
  
  do.call(rbind, expm_entries)
}

# ---- Survival curves plotting -----
plot_survival_curves <- function(surv_matrix, 
                                 patient_ids = NULL, 
                                 seed = 42, 
                                 n_patients = 5,
                                 title = "Survival Curves") {
  if (is.null(patient_ids)) {
    set.seed(seed)
    if (nrow(surv_matrix) <= n_patients){
      n_patients = nrow(surv_matrix)
    }
    patient_ids <- sample(1:nrow(surv_matrix), n_patients)
  }
  
  subset_curves <- surv_matrix[patient_ids, , drop = FALSE]
  df <- as.data.frame(subset_curves)
  df$patient_id <- rownames(subset_curves)
  df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
  
  df_long <- reshape2::melt(df, id.vars = "patient_id", variable.name = "time", value.name = "surv_prob")
  df_long$time <- as.numeric(as.character(df_long$time))
  
  # Plot
  ggplot(df_long, aes(x = time, y = surv_prob, color = patient_id)) +
    geom_line(linewidth = 1) +
    labs(title = title, x = "Time", y = "Survival Probability",
         color = "Patient ID") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
}

# ---- Extract info from bootstrap and point estimates synthetic data-----
# When calculated with Expected mortality
extract_expm_entries_synthetic <- function(list_exp, list_exp_pe, lambda) {
  lapply(seq_along(list_exp), function(dataset_idx) {
    dataset_results <- list_exp[[dataset_idx]]
    dataset_results_pe <- list_exp_pe[[dataset_idx]]
    
    lapply(names(dataset_results), function(name) {
      if (name == "batch.metrics") return(NULL)
      
      result <- dataset_results[[name]]
      metrics <- names(result$mean)
      
      entries <- lapply(seq_along(metrics), function(i) {
        metric_name <- metrics[i]
        mean_c <- dataset_results_pe[[name]][[metric_name]]
        ci <- result$confidence.intervals[i, ]
        
        data.frame(
          Lambda = lambda,
          Dataset = dataset_idx,
          Metric = metric_name,
          Model = "CoxPH",
          Cindex = mean_c,
          Lower = ci[1],
          Upper = ci[2],
          InputType = "Exp.Mort",
          stringsAsFactors = FALSE
        )
      })
      
      do.call(rbind, entries)
    }) |> 
      bind_rows()  
  })
}

# ---- Extract info from bootstrap and point estimates synthetic data-----
# When calculated with RMST
extract_rmst_entries_synthetic <- function(list_rmst, list_rmst_pe, lambda) {

  lapply(seq_along(list_rmst), function(dataset_idx) {
    dataset_results <- list_rmst[[dataset_idx]]
    dataset_results_pe <- list_rmst_pe[[dataset_idx]]

    lapply(names(dataset_results), function(name) {
      if (name == "batch.metrics") return(NULL)
      
      result <- dataset_results[[name]]
      metrics <- names(result$mean)
      
      entries <- lapply(seq_along(metrics), function(i) {
        metric_name <- metrics[i]
        mean_c <- dataset_results_pe[[name]][[metric_name]]
        ci <- result$confidence.intervals[i, ]
       
        data.frame(
          Lambda = lambda,
          Dataset = dataset_idx,
          Metric = metric_name,
          Model = "CoxPH",
          Cindex = mean_c,
          Lower = ci[1],
          Upper = ci[2],
          InputType = "RMST",
          stringsAsFactors = FALSE
        )
        
      })
      
      do.call(rbind, entries)
    }) |> 
      bind_rows()  
  })
}

# ---- Extract info from bootstrap and point estimates synthetic data version 2 -----
# When calculated with Expected mortality
extract_rmst_entries_synthetic2 <- function(list_rmst_pe, lambda, model) {
  
  lapply(seq_along(list_rmst_pe), function(dataset_idx) {
    dataset_results_pe <- list_rmst_pe[[dataset_idx]]
    
    lapply(names(dataset_results_pe), function(name) {
      if (name == "batch.metrics") return(NULL)

      result <- dataset_results_pe[[name]]
      metrics <- names(result)
      metrics <- metrics[metrics %in% c("Hmisc::rcorr.cens",
                                        "pysurvival",
                                        "SurvMetrics::Cindex", 
                                        "lifelines", 
                                        "sksurv.censored",
                                        "survC1::Est.Cval", 
                                        "pec::cindex", 
                                        "survival.n", 
                                        "survival.n/G2")]
      
      entries <- lapply(seq_along(metrics), function(i) {
        metric_name <- metrics[i]
        mean_c <- dataset_results_pe[[name]][[metric_name]]

        data.frame(
          Lambda = lambda,
          Dataset = dataset_idx,
          Metric = metric_name,
          Model = model,
          Cindex = mean_c,
          InputType = "RMST",
          stringsAsFactors = FALSE
        )
        
      })
      
      do.call(rbind, entries)
    }) |> 
      bind_rows()  
  })
}
# ---- Extract info from bootstrap and point estimates synthetic data -----
# When calculated with Survival probs
extract_surv_entries_synthetic <- function(list_surv, list_surv_pe, lambda) {
  lapply(seq_along(list_surv), function(dataset_idx) {
    dataset_results <- list_surv[[dataset_idx]]
    dataset_results_pe <- list_surv_pe[[dataset_idx]]
    
    lapply(names(dataset_results), function(name) {
      if (name == "batch.metrics") return(NULL)
      
      result <- dataset_results[[name]]
      metrics <- names(result$mean)
      
      entries <- lapply(1:2, function(i) {
        metric_name <- metrics[i]
        mean_c <- dataset_results_pe[[name]][[metric_name]]
        ci <- result$confidence.intervals[i, ]
        
        data.frame(
          Lambda = lambda,
          Dataset = dataset_idx,
          Metric = metric_name,
          Model = "CoxPH",
          Cindex = mean_c,
          Lower = ci[1],
          Upper = ci[2],
          InputType = "Distrib.",
          stringsAsFactors = FALSE
        )
      })
      
      do.call(rbind, entries)
    }) |>
      bind_rows()
  })
}
# ---- Extract info from bootstrap and point estimates synthetic data version 2 -----
# When calculated with Survival probs
extract_surv_entries_synthetic2 <- function(list_surv_pe, lambda, model) {
  lapply(seq_along(list_surv_pe), function(dataset_idx) {
    dataset_results_pe <- list_surv_pe[[dataset_idx]]

    lapply(names(dataset_results_pe), function(name) {
      if (name == "batch.metrics") return(NULL)
      
      result <- dataset_results_pe[[name]]
      metrics <- names(result)
      metrics <- metrics[metrics %in% c("pycox.Ant", "pycox.Adj.Ant")]

      entries <- lapply(seq_along(metrics), function(i) {
        metric_name <- metrics[i]
        mean_c <- dataset_results_pe[[name]][[metric_name]]
        data.frame(
          Lambda = lambda,
          Dataset = dataset_idx,
          Metric = metric_name,
          Model = model,
          Cindex = mean_c,
          InputType = "Distrib.",
          stringsAsFactors = FALSE
        )
      })
      
      do.call(rbind, entries)
    }) |>
      bind_rows()
  })
}


# extract_rms_pairs <- function(list_rms, lambda, extract) {
#   all_batches <- lapply(seq_along(list_rms), function(dataset_idx) {
#     dataset_results <- list_rms[[dataset_idx]]
# 
#     filtered_list <- lapply(names(dataset_results), function(name) {
#       
#       filtered_batch <- dataset_results[[name]]
#       filtered_batch$lambda <- lambda
#       filtered_batch$dataset <- dataset_idx
#       return(filtered_batch)
#       })
#     }) 
#     
#     bind_rows(unlist(all_batches, recursive = FALSE))
#   }
# 
# extract_batch_metrics <- function(results_list) {
#   models <- names(results_list)
#   
#   # Loop over models
#   metrics_list <- lapply(models, function(model_name) {
#     model_data <- results_list[[model_name]]
#     
#     if (!is.null(model_data$batch.metrics)) {
#       df <- model_data$batch.metrics
#       df$Model <- model_name
#       return(df)
#     } else {
#       return(NULL)
#     }
#   })
#   
#   # Combine all into a single data.frame
#   metrics_df <- do.call(rbind, metrics_list)
#   return(metrics_df)
# }

# ---- Compute RMST and Expected mortality -----
# During 5 fold cross validation
compute_measures <- function(df, model_name, time_step = 1) {
  model_cols <- grep(paste0("^", model_name, "\\."), names(df), value = TRUE)
  surv_mat <- as.matrix(df[, model_cols])
  
  rmst <- (-rowSums(surv_mat) * time_step)
  
  if (any(as.vector(surv_mat) == 0))  {
    cat("There is a 0 ", model_name," survival preds \n")
    second_smallest_survival_prob = unique(sort(as.vector(surv_mat), partial = 2))[2]
    surv_mat <- surv_mat + second_smallest_survival_prob
  } else {
    surv_mat <- surv_mat
  }
  exp_mort <- rowSums(-log(surv_mat))
  
  return(list(rmst = rmst,
              exp_mort = exp_mort,
              surv_mat = surv_mat)
  )
}

# ---- Make alluvial -----
make_alluvial_plot <- function(df, metric, notation, 
                               title, custom_colors) {
  data <- df %>%
    filter(Metric == metric, Notation == notation) %>%
    group_by(fold_n) %>%
    mutate(Rank = rank(-cindex, ties.method = "first")) %>%
    ungroup() %>%
    mutate(fold_n = factor(fold_n),
           Rank = as.factor(Rank))
  
  ggplot(data,
         aes(x = fold_n, stratum = Rank, alluvium = Model, y = 1,
             fill = Model)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.9, color = "grey") +
    geom_stratum(width = 0.25, color = "grey30") +
    scale_fill_manual(values = custom_colors, labels = scales::parse_format()) +
    labs(x = "Cross-validation Fold", y = NULL, title = NULL) +
    ggtitle(parse(text = title)) +
    theme_minimal(base_size = 17) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(size = 17),  
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_text(size = 17),
      axis.title.y = element_text(size = 17),
      legend.position = "bottom"
    )
}

extract_rmst_pairs <- function(list_rmst, lambda) {
  all_batches <- lapply(seq_along(list_rmst), function(dataset_idx) {
    dataset_results <- list_rmst[[dataset_idx]]
    filtered_list <- lapply(names(dataset_results), function(name) {
      results <- dataset_results[[name]]
      results$lambda <- lambda
      results$dataset <- dataset_idx
      return(results)
    })
  }) 
  
  bind_rows(unlist(all_batches, recursive = FALSE))
}


### ---- Same as metrics.wrapper but modification of ties ----
metrics.wrapper2 <- function(predicted, surv_matrix = NULL, 
                            censoring, time, 
                            implementation, eval.times, 
                            sksurv_train_time = NULL, # sksurv.ipcw
                            sksurv_train_status = NULL, # sksurv.ipcw
                            sksurv_tied_tol = NULL,
                            pec_incl_ties_preds = TRUE, # pec TRUE default 
                            pec_incl_ties_times = TRUE, # pec TRUE default 
                            pec_incl_both_tied = TRUE, # pec TRUE default 
                            hmisc_excl_ties_preds = FALSE, # default is False 
                            pysurvival_incl_ties_preds = TRUE # default is TRUE
                            
) { # sksurv.ipcw
  
  # Initialize an empty array for the results
  results <- list()
  # Initialize eval.times list
  batch.eval.times <- c()
  # Flatten the implementation list into a single vector
  implementation <- unlist(implementation)
  
  # If eval.times is missing, set to the maximum uncensoring time as in pec
  # The maximum uncensored time is when the C-index stops updating
  if (missing(eval.times) || is.null(eval.times) || (eval.times == "resample_max_uncensored_time")) {
    eval.times <- max(time[censoring == 1], na.rm=TRUE) 
  }
  # Store the eval.times for the results
  results$batch.eval.times <- c(results$batch.eval.times, eval.times)
  
  # For pec::cindex
  if ("pec::cindex" %in% implementation) {
    # Create a temporary dataframe
    tmp <- data.frame(predicted, 
                      censoring, 
                      time)
    
    pec <- pec::cindex(
      list("res" = as.matrix(1 - predicted)),
      formula = Hist(time, censoring) ~ 1, 
      data = tmp, 
      eval.times = eval.times,
      cens.model = "marginal", # default marginal
      tiedPredictionsIn = pec_incl_ties_preds, # by default TRUE including them
      tiedOutcomeIn = pec_incl_ties_times, # by default TRUE including tied times (their docu is wrong)
      tiedMatchIn = pec_incl_both_tied # by default TRUE including pairs with identical times and preds.
    )
    # Store Cindex
    results[["pec::cindex"]] <- pec$AppCindex$res
    # Store elegible pairs
    results[["pec::cindex.EP"]] <- pec$Pairs$res
    # Store concordant pairs
    results[["pec::cindex.CP"]] <- pec$Concordant$res
    
  }
  
  # Uno's real implementation: 
  # https://github.com/cran/survC1/blob/master/R/FUN-cstat-ver003b.R
  # https://github.com/cran/survC1/blob/master/src/husurvC1v1.f
  if ("survC1::Est.Cval" %in% implementation) {
    # Need to create this survival object: 
    results[["survC1::Est.Cval"]] <- 
      survC1::Est.Cval(data.frame(time, 
                                  censoring, 
                                  predicted), # Expects risk scores
                       tau=eval.times, 
                       nofit=TRUE)$Dhat 
  }
  
  # For Harrel's real implementation
  if ("Hmisc::rcorr.cens" %in% implementation) {
    
    hmisc <- Hmisc::rcorr.cens(x = 1 - predicted, 
                               S = Surv(time, censoring), 
                               outx = hmisc_excl_ties_preds) # to control for it
    
    results[["Hmisc::rcorr.cens"]] <- hmisc[["C Index"]]
    results[["Hmisc::rcorr.cens.EP"]] <- hmisc[["Relevant Pairs"]]
    results[["Hmisc::rcorr.cens.CP"]] <- hmisc[["Concordant"]]
    #results[["Hmisc::rcorr.cens.CP"]] <- hmisc[["Uncertain"]]
    
  }
  
  # For SurvMetrics::Cindex
  if ("SurvMetrics::Cindex" %in% implementation) {
    # Need to create this survival object: 
    surv_obj = Surv(time, censoring)
    # Calculate the metric
    # t_star is only used to extract surv at a specific time point
    # not to calculate 
    results[["SurvMetrics::Cindex"]] <-
      SurvMetrics::Cindex(surv_obj, predicted=1-predicted, 
                          t_star = eval.times)[[1]]
    
  }
  
  # For randomForestSRC::get.cindex 
  if ("randomForestSRC::get.cindex" %in% implementation) {
    results[["randomForestSRC::get.cindex"]] <- 
      randomForestSRC::get.cindex(predicted = 1 - predicted,
                                  censoring = censoring,
                                  time = time) 
  }
  
  # lifelines::concordance_index implementation
  if ("lifelines" %in% implementation) {
    # Handle Python implementations
    results[["lifelines"]] <- lifelinesR(time, predicted, censoring)
    
  }
  # https://github.com/havakv/pycox/blob/master/pycox/evaluation/concordance.py#L12
  # using the eval.times only to subset the survival probabilities at a time point
  if ("pycox.Ant" %in% implementation) {
    # Handle Python implementations
    #results[["pycox"]] <- pycoxR(time, predicted, censoring, eval.times)
    pycox_ant <- pycoxR_SurvMatrix_Ant(time=time, 
                                       surv_matrix=surv_matrix, 
                                       censoring=censoring)
    
    results[["pycox.Ant"]] <- pycox_ant$cindex
    results[["pycox.Ant.EP"]] <- pycox_ant$comparable
    results[["pycox.Ant.CP"]] <- pycox_ant$concordant
  }
  if ("pycox.Adj.Ant" %in% implementation) {
    # Handle Python implementations
    #results[["pycox"]] <- pycoxR(time, predicted, censoring, eval.times)
    pycox_adj <- pycoxR_SurvMatrix_AdjAnt(time=time, 
                                          surv_matrix=surv_matrix, 
                                          censoring=censoring)
    
    results[["pycox.Adj.Ant"]] <- pycox_adj$cindex
    results[["pycox.Adj.Ant.EP"]] <- pycox_adj$comparable
    results[["pycox.Adj.Ant.CP"]] <- pycox_adj$concordant
  }
  
  # For pysurvival we have extracted the function of interest in file pysurvivalR.tar.gz
  # https://github.com/square/pysurvival/blob/master/pysurvival/cpp_extensions/metrics.cpp
  if ("pysurvival" %in% implementation){
    # R wrapped version of python/c++ library pysurvival
    
    pysurv <- pysurvivalR::concordance_index(risk=predicted,
                                             T=time,
                                             E=censoring,
                                             include_ties = pysurvival_incl_ties_preds) # risk ties default\
    
    results[["pysurvival"]] <- pysurv$`0` 
    results[["pysulvival.EP"]] <- pysurv$`1` 
    results[["pysulvival.CP"]] <- pysurv$`2` 
    
  }
  
  # For sksurv without ipcw
  if ("sksurv.censored" %in% implementation) { # Like Harrels
    sksurv <- sksurv.censoredR(time,
                               predicted, 
                               censoring)
    
    results[["sksurv.censored"]] <- sksurv[[1]] # C index
    results[["sksurv.censored.CP"]] <- sksurv[[2]] # Concordant
    results[["sksurv.censored.DP"]] <- sksurv[[3]] # Disconcordant
    results[["sksurv.censored.TR"]] <- sksurv[[4]] # tied risk
    results[["sksurv.censored.TT"]] <- sksurv[[5]] # tied time
  }
  
  # For sksurv with ipcw
  if ("sksurv.ipcw" %in% implementation) { # Like Uno's
    sksurv <- sksurv.ipcwR(time,
                           predicted,
                           censoring,
                           eval.times,
                           sksurv_train_time, 
                           sksurv_train_status, 
                           sksurv_tied_tol)
    
    results[["sksurv.ipcw"]] <- sksurv[[1]] # C index
    results[["sksurv.ipcw.CP"]] <- sksurv[[2]] # Concordant
    results[["sksurv.ipcw.DP"]] <- sksurv[[3]] # Disconcordant
    results[["sksurv.ipcw.TR"]] <- sksurv[[4]] # tied risk
    results[["sksurv.ipcw.TT"]] <- sksurv[[5]] # tied time
  }
  
  # For survival with IPCW
  if ("survival.n/G2" %in% implementation) {
    
    predicted = 1- predicted 
    surv.conc <- survival::concordance(Surv(time, censoring) ~ predicted, 
                                       timewt ="n/G2", ymax = eval.times)
    results[["survival.n/G2"]] <- surv.conc$concordance
    results[["survival.n/G2.CP"]]  <- surv.conc$count[["concordant"]] # Concordanct
    results[["survival.n/G2.DP"]]  <- surv.conc$count[["discordant"]] # Discordant
    results[["survival.n/G2.TR"]]  <- surv.conc$count[["tied.x"]] # Tied risk
    results[["survival.n/G2.TT"]]  <- surv.conc$count[["tied.y"]] # Tied time
    results[["survival.n/G2.TB"]]  <- surv.conc$count[["tied.xy"]] # Tied both time and risk
    
  }
  
  # For survival without IPCW
  if ("survival.n" %in% implementation) {
    
    predicted = predicted # expects survival? confusing....
    surv.conc <- survival::concordance(Surv(time, censoring) ~ predicted, 
                                       timewt = "n", ymax = eval.times)
    results[["survival.n"]] <- surv.conc$concordance
    results[["survival.n.CP"]]  <- surv.conc$count[["concordant"]] # Concordanct
    results[["survival.n.DP"]]  <- surv.conc$count[["discordant"]] # Discordant
    results[["survival.n.TR"]]  <- surv.conc$count[["tied.x"]] # Tied risk
    results[["survival.n.TT"]]  <- surv.conc$count[["tied.y"]] # Tied time
    results[["survival.n.TB"]]  <- surv.conc$count[["tied.xy"]] # Tied both time and risk
    
  }
  
  if (length(results) == 0) {
    stop("Missing implementation list.
         It is required a list with all the c-index 
         implementations that need to be computed")
  }
  
  return(results)
  
}