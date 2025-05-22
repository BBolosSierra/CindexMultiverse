set.seed(123)


# Lifelines wrapper for R
lifelinesR <- function(time, predicted, censoring){
  # Function to convert python into R wrapping lifelines
  np   <- import("numpy", convert = TRUE)
  concordance_index <- import("lifelines.utils", convert = TRUE)$concordance_index
  
  concordance_index(event_times = time,
                    predicted_scores = 1 - predicted,
                    event_observed =  censoring)

}

# pycox wrapper for R
# If this is done, there is no antolinis... 
# pycoxR <- function(time, predicted, censoring, eval.times) {  
#   
#   # Import numpy and pandas with False to use the python objects
#   np <- import("numpy", convert = FALSE)  
#   pd <- import("pandas", convert = FALSE)
#   EvalSurv <- import("pycox.evaluation", convert = FALSE)$EvalSurv
# 
#   # Convert R to python with correct float or int
#   time_py <- np$array(as.numeric(time), dtype = "float64")  
#   censoring_py <- np$array(as.integer(censoring), dtype = "int32")  
#   predicted_py <- np$array(as.numeric(predicted), dtype = "float64")
#   
#   St <- np$subtract(1, predicted_py)
#   
#   St_reshaped <- St$reshape(1L, -1L)
#   
#   eval_times_py <-  np$array(list(as.numeric(eval.times)), dtype = "float64")
#   
#   surv_df <- pd$DataFrame(data = St_reshaped, index = eval_times_py)  
#   
#   ev <- EvalSurv(surv_df, time_py, censoring_py, censor_surv = 'km')
# 
#   # Two versions of C-index: here antolini used in DeepHit
#   # Alternatively "adj_antolini"
#   py_to_r(ev$concordance_td("antolini"))
# 
# }

pycoxR_SurvMatrix_Ant <- function(time, surv_matrix, censoring){
  ### As difference from the previous one, this function calculates
  ### c-index across all times 

  np <- import("numpy", convert = FALSE)  
  pd <- import("pandas", convert = FALSE)
  EvalSurv <- import("pycox.evaluation", convert = FALSE)$EvalSurv
  
  # Convert R to python with correct float or int
  time_py <- np$array(as.numeric(time), dtype = "float64")  
  #time_py <- np$round(time_py, as.integer(2))
  #time_py <- np$array(as.numeric(time_py), dtype = "float64")  
  censoring_py <- np$array(as.integer(censoring), dtype = "int32")
  # Transpose the survival matrix
  surv_mat <- t(as.matrix(surv_matrix))  # now: rows = times, cols = patients
  
  # Convert to numpy
  surv_np <- np$array(surv_mat, dtype = "float64")

  # Time index 
  time_index <- np$array(as.numeric(colnames(surv_matrix)), dtype = "float64")
  
  # Create pandas dataframe
  surv_df <- pd$DataFrame(data = surv_np, index = time_index)

  ev <- EvalSurv(surv_df, time_py, censoring_py, censor_surv = "km")
  
  py_to_r(ev$concordance_td(method = "antolini"))
}


pycoxR_SurvMatrix_AdjAnt <- function(time, surv_matrix, censoring){
  ### As difference from the previous one, this function calculates
  ### c-index across all times 
  np <- import("numpy", convert = FALSE)  
  pd <- import("pandas", convert = FALSE)
  EvalSurv <- import("pycox.evaluation", convert = FALSE)$EvalSurv
  
  # Convert R to python with correct float or int
  time_py <- np$array(as.numeric(time), dtype = "float64")  
  censoring_py <- np$array(as.integer(censoring), dtype = "int32")
  # Transpose the survival matrix
  surv_mat <- t(as.matrix(surv_matrix))  # now: rows = times, cols = patients
  
  # Convert to numpy
  surv_np <- np$array(surv_mat, dtype = "float64")
  
  # Time index 
  time_index <- np$array(as.numeric(colnames(surv_matrix)), dtype = "float64")
  
  # Create pandas dataframe
  surv_df <- pd$DataFrame(data = surv_np, index = time_index)
  
  ev <- EvalSurv(surv_df, time_py, censoring_py, censor_surv = "km")
  
  py_to_r(ev$concordance_td(method = "adj_antolini"))
}




sksurv.censoredR <- function(time, predicted, censoring) {
  # like sksurvR::concordance_index_censored
  np <- import("numpy", convert = FALSE)
  sksurv <- import("sksurv.metrics", convert = FALSE)
  
  time_py <- np$array(as.numeric(time), dtype = "float64")  
  censoring_py <- np$array(as.logical(censoring), dtype = "bool")  
  predicted_py <- np$array(as.numeric(predicted), dtype = "float64")
  
  c_index <- sksurv$concordance_index_censored(
    event_indicator = censoring_py,  
    event_time = time_py,
    estimate = predicted_py
  )
  
  py_to_r(c_index)
  
}

  
sksurv.ipcwR <- function(time, predicted, censoring, 
                         eval.times, 
                         sksurv_train_time, 
                         sksurv_train_status, 
                         sksurv_tied_tol = 1e-8) { # Default

  # Like Uno's
  np <- import("numpy", convert = FALSE)
  sksurv <- import("sksurv.metrics", convert = FALSE)
  
  time_train_py <- np$array(as.numeric(sksurv_train_time), dtype = "float64")
  censoring_train_py <- np$array(as.logical(sksurv_train_status), dtype = "bool")
  
  time_test_py <- np$array(as.numeric(time), dtype = "float64")
  censoring_test_py <- np$array(as.logical(censoring), dtype = "bool")  
  
  predicted_py <- np$array(as.numeric(predicted), dtype = "float64")
  
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


bootstrap.metric <- function(metrics.wrapper, 
                             dataset, 
                             implementation, 
                             eval.times = NULL, 
                             sampled_data = NULL,
                             resample_indices = NULL, 
                             additional = NULL, 
                             seed = NULL) { 
  
  if (is.null(sampled_data)) {
    
    # Helper function to sample dataset based on indices
    sample_dataset <- function(dataset, sample_idx) {
      sampled_dataset <- lapply(dataset, function(value) value[sample_idx])
      return(sampled_dataset)
    }
    
    N <- 1:length(resample_indices)
    # If resample indices are not provided, generate them
    if (is.null(resample_indices)) {
      set.seed(seed)
      size <- length(dataset[[1]])  # Assuming all elements in the dataset are of the same length
      resample_indices <- replicate(N, sample(seq_len(size), size = size, replace = TRUE), simplify = FALSE)
    }
    # Calculate the number of bootstraps
    N_bootstraps <- length(resample_indices)
    # Initialize iterations
    iterations <- vector("list", length(resample_indices))  # Match the number of resample indices
    
    for (i in seq_along(resample_indices)) {
      # Use provided or internally generated resample indices
      resample_idx <- resample_indices[[i]]
      
      # Sample the dataset and calculate the metric
      sampled_data <- sample_dataset(dataset, resample_idx)
  
      # For the sampled data what happens if we sort
      iterations[[i]] <- do.call(metrics.wrapper, c(sampled_data, list(implementation = implementation), eval.times, additional)) # sksurv.ipcw
  
    }
  # if the data has already been sampled from indexes:  
  } else {
    iterations <- vector("list", length(sampled_data))
    # Calculate the number of bootstraps
    N_bootstraps <- length(sampled_data)
    # Loop through the sampled_datar
    for (i in seq_along(sampled_data)) {
      # Subset each sampled data from list
      resampled_data <- sampled_data[[i]]
      
      if (("pycox.Ant" %in% implementation) | ("pycox.Adj.Ant" %in% implementation)) {
        data <- data.frame(
          censoring = resampled_data[[dataset$censoring]],
          time = resampled_data[[dataset$time]]
        )
        surv_matrix = resampled_data[, grep(dataset$predicted, names(resampled_data), value = TRUE), drop = FALSE]
        # Run the implementations

        colnames(surv_matrix) <- as.numeric(sub(".*\\.", "", colnames(surv_matrix)))
        iterations[[i]] <- do.call(metrics.wrapper, c(data, list(surv_matrix = surv_matrix), list(implementation = implementation), eval.times, additional)) # sksurv.ipcw
        
      } else {
        # Create dataset for calculation with minimal columns
        data <- data.frame(
          predicted = resampled_data[[dataset$predicted]],
          censoring = resampled_data[[dataset$censoring]],
          time = resampled_data[[dataset$time]]
        )
        iterations[[i]] <- do.call(metrics.wrapper, 
                                   c(data, list(implementation = implementation), 
                                     eval.times, additional)) # sksurv.ipcw
      }
    }
  }
  # Aggregate metrics
  metrics_df <- do.call(rbind, lapply(iterations, function(x) {
    as.data.frame(as.list(x), check.names = FALSE)  # Keep "::" in column names
  }))
  
  # Calculate mean and 95% confidence intervals for each metric
  subset_metrics <- metrics_df[, colnames(metrics_df) %in% implementation, drop = FALSE] 
  mean_metrics <- colMeans(subset_metrics, na.rm = TRUE)
  conf_intervals <- t(apply(subset_metrics, 2, function(metric) {
    mean_metric <- mean(metric, na.rm = TRUE)
    stderr <- sd(metric, na.rm = TRUE) / sqrt(N_bootstraps)
    mean_metric + qt(c(0.025, 0.975), df = N_bootstraps - 1) * stderr
  }))
  
  return(list(
    mean = mean_metrics,
    confidence.intervals = conf_intervals,
    eval.times = eval.times,
    batch.metrics = metrics_df
  ))
}

# With parallelisation
bootstrap.metric.parallel <- function(metrics.wrapper, 
                             dataset, 
                             implementation, 
                             eval.times = NULL,
                             sampled_data = NULL, 
                             resample_indices = NULL, 
                             additional = NULL, 
                             seed_resample = NULL) {
  
  # If the data has been already sampled from the bootstrap index,
  # use the sampled data, otherwise, based on the resample index 
  # resample the dataset
  if (is.null(sampled_data)) {
   
    # Helper function to sample dataset based on indices
    sample_dataset <- function(dataset, sample_idx) {
      sampled_dataset <- lapply(dataset, function(value) value[sample_idx])
      return(sampled_dataset)
    }
    
    N <- 1:length(resample_indices)
    
    # If sample indices are not provided, generate them
    if (is.null(resample_indices)) {
      # generate resample indeces with a seed
      set.seed(seed_resample)
      size <- length(dataset[[1]])  # Assuming all elements in the dataset are of the same length
      resample_indices <- replicate(N, sample(seq_len(size), 
                                              size = size, replace = TRUE), 
                                    simplify = FALSE)
    }

    N_bootstraps <- length(resample_indices)
    # Initialize iterations
    resampled_data <- vector("list", length = N_bootstraps) 
    for (i in seq_along(resample_indices)) {
      # Use provided or internally generated resample indices
      resample_idx <- resample_indices[[i]]
  
      # Sample the dataset and calculate the metric
      resampled_data[[i]] <- sample_dataset(dataset, resample_idx)
    
    }
  } else {
    iterations <- vector("list", length(sampled_data))
    # Calculate the number of bootstraps
    N_bootstraps <- length(sampled_data)
    # Loop through the sampled_datar
    resampled_data <- vector("list", length = N_bootstraps)
    surv_matrixes <- vector("list", length = N_bootstraps)
    for (i in seq_along(sampled_data)) {
      # Subset each sampled data from list
      sampled_data_ <- sampled_data[[i]]
      
      if (("pycox.Ant" %in% implementation) | ("pycox.Adj.Ant" %in% implementation)) {
        resampled_data[[i]] <- data.frame(
          censoring = sampled_data_[[dataset$censoring]],
          time = sampled_data_[[dataset$time]]
        )
        # Just take the survival matrix for the selected 
        surv_mat = sampled_data_[, grep(dataset$predicted, names(sampled_data_), value = TRUE), drop = FALSE]
        # Run the implementations
        #browser()
        # iterations[[i]] <- do.call(metrics.wrapper, c(data, list(surv_matrix = surv_matrix), list(implementation = implementation), eval.times, additional)) # sksurv.ipcw
        colnames(surv_mat) <- as.numeric(sub(".*\\.", "", colnames(surv_mat)))
        surv_matrixes[[i]] <- surv_mat
      } else {
        # Create dataset for calculation with minimal columns
        resampled_data[[i]] <- data.frame(
          predicted = sampled_data_[[dataset$predicted]],
          censoring = sampled_data_[[dataset$censoring]],
          time = sampled_data_[[dataset$time]]
        )
      }
    }
  }
  #parallel_func <- function(resampled_data, implementation, eval.times, surv_matrixes) { 
  parallel_func <- function(resampled_data, implementation, eval.times, surv_matrixes, additional) { 
    
    ## build progress bar
    pro_bar <- progressor(along = resampled_data)
  
    #Run foreach
    results <- suppressWarnings(foreach(i = seq_along(resampled_data),
                       data = resampled_data,
                       .combine = rbind,  ## Ensures all results are collected
                       .options.future = list(seed = TRUE)) %dofuture% {

                         #py_config()
                         ## Update progress bar
                         pro_bar(sprintf("Iteration %d of %d", i, length(resampled_data)))

                         #print(sapply(ls(), function(x) format(object.size(get(x)), units="auto")))
                         #print(py_config())
                         #surv_matrix <- surv_matrixes[[i]]
                         
                         ## Call function with proper arguments
                         do.call(metrics.wrapper, c(data, list(surv_matrix = surv_matrixes[[i]]),
                                                    list(implementation = implementation),
                                                    eval.times, additional))
                       })
    
  }
  
  ## Set the parallel backend to future package
  ## tell foreach to use futures
  registerDoFuture() 
  ## Plan doesn't work with reticulate as multi-session is 

  plan(multisession) # Linux/Mac/Windows
  #plan(multicore) # Linux/Mac
  #plan(sequential)

  ## Set progress bar style
  handlers("cli")

  metrics_df <- with_progress(parallel_func(resampled_data,
                                            implementation, 
                                            eval.times, 
                                            surv_matrixes, additional))

  # Need to remove index from paralelization
  rownames(metrics_df) <- NULL
  # Calculate mean and 95% confidence intervals for each metric
  subset_metrics <- metrics_df[, colnames(metrics_df) %in% implementation, drop = FALSE] 
  mean_metrics <- apply(subset_metrics, 2, function(x) mean(as.numeric(x), na.rm = TRUE))
  conf_intervals <- t(apply(subset_metrics, 2, function(metric) {
    mean_metric <- mean(as.numeric(metric), na.rm = TRUE)
    stderr <- sd(as.numeric(metric), na.rm = TRUE) / sqrt(N_bootstraps)
    mean_metric + qt(c(0.025, 0.975), df = N_bootstraps - 1) * stderr
  }))
  
  return(list(
    mean = mean_metrics,
    confidence.intervals = conf_intervals,
    eval.times = eval.times,
    batch.metrics = metrics_df
  ))
}

# Metrics wrapper:
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
  if (missing(eval.times) || is.null(eval.times)) {
    eval.times <- max(time[censoring == 1], na.rm=TRUE) # ?
  }
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
      cens.model = "marginal" # default
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
    results[["pycox.Ant"]] <- pycoxR_SurvMatrix_Ant(time=time, 
                                                    surv_matrix=surv_matrix, 
                                                    censoring=censoring)
    
  }
  if ("pycox.Adj.Ant" %in% implementation) {
    # Handle Python implementations
    #results[["pycox"]] <- pycoxR(time, predicted, censoring, eval.times)
    results[["pycox.Adj.Ant"]] <- pycoxR_SurvMatrix_AdjAnt(time=time, 
                                                           surv_matrix=surv_matrix, 
                                                           censoring=censoring)
    
  }
  #https://github.com/square/pysurvival/blob/master/pysurvival/cpp_extensions/metrics.cpp
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
  
  # If no implementations match, return NA for all metrics
  if (length(results) == 0) {
    stop("Missing implementation list. 
         It is required a list with all the c-index 
         implementations that need to be computed")
  }
  
  return(results)
  
}



plot.cindex.distribution <- function(implementation, 
                                     bootstrap.object, 
                                     title, 
                                     y_limits = NULL) {
  
  # Extract metrics from bootstrap object
  batch.metrics <- bootstrap.object$batch.metrics 
  
  # Convert batch metrics to long format for ggplot2
  long_batch_metrics <- reshape2::melt(batch.metrics, 
                                       variable.name = "Metric", 
                                       value.name = "Value")
  long_batch_metrics <- long_batch_metrics[long_batch_metrics$Metric %in% implementation, ]
  
  # Create the plot
  ggplot() +
    # Add boxplots for each metric
    geom_violin(
      data = long_batch_metrics,
      aes(x = Metric, y = Value),
      width = 0.4,
      fill = "lightblue",
      alpha = 0.6
    ) +
    # Add jittered points for batch metrics
    geom_jitter(
      data = long_batch_metrics,
      aes(x = Metric, y = Value),
      width = 0.2,
      color = "darkgray",
      alpha = 0.5,
      size = 1
    ) +
    labs(
      title = title,
      x = "Metric",
      y = "Value"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    if (!is.null(y_limits)) ylim(y_limits[1], y_limits[2]) else NULL
  
}

plot.cindex.ci <- function(implementation, bootstrap.object, title, y_limits = NULL) {
  
  # Extract metrics from bootstrap object
  means <- bootstrap.object$mean
  confidence.intervals <- bootstrap.object$confidence.intervals
  
  # Convert the data into dataframe
  metrics_data <- data.frame(
    Metric = names(means),
    Mean = as.numeric(means),
    Lower = as.numeric(confidence.intervals[, 1]),
    Upper = as.numeric(confidence.intervals[, 2])
  )
  
  # Filter by implementation list
  metrics_data <- metrics_data[metrics_data$Metric %in% implementation, ]
  
  # Create boxplots with confidence intervals
  ggplot(metrics_data, aes(x = Metric, y = Mean)) +
    geom_boxplot(width = 0.4, fill = "lightblue", alpha = 0.6) +  # Boxplot
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "red") +  # Confidence intervals
    geom_point(size = 3, color = "darkblue") +  # Mean points
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray", size = 0.8) +
    theme_minimal() +
    labs(
      title = title,
      x = "Metric",
      y = "Mean Value"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + 
    if (!is.null(y_limits)) ylim(y_limits[1], y_limits[2]) else NULL
  
}

plot.cindex.ci.longitudinal <- function(implementation, bootstrap.object, title, y_limits = NULL, time.points = NULL) {
  
  metrics_data <- do.call(rbind, lapply(time.points, function(t) {
    
    # Access bootstrap results at this time point
    bootstrap_t <- bootstrap.object[[as.character(t)]]
    
    # Create dataframe for ggplot
    data.frame(
      Time = t,  # Store time point
      Metric = names(bootstrap_t$mean),  # Extract metric names
      Mean = as.numeric(bootstrap_t$mean),  # Extract mean values
      Lower = as.numeric(bootstrap_t$confidence.intervals[, 1]),  # Lower CI
      Upper = as.numeric(bootstrap_t$confidence.intervals[, 2])   # Upper CI
    )
  }))
  
  
  # Filter by implementation
  metrics_data <- metrics_data[metrics_data$Metric %in% implementation, ]
  # Plot
  ggplot(metrics_data, aes(x = Time, y = Mean, color = Metric, group = Metric)) +
    geom_line(size = 1) + 
    geom_point(size = 3) +  
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +  
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray", size = 0.8) +  
    theme_minimal() +
    scale_x_continuous(breaks = time.points) +
    labs(
      title = title,
      x = "Time Point",
      y = "Mean C-Index",
      color = "Implementation"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    if (!is.null(y_limits)) ylim(y_limits[1], y_limits[2]) else NULL
}  
  
  

plot.violin.metrics <- function(selected_values, 
                                column_name, 
                                bootstrap_object, 
                                stacked_predictions, 
                                y_limits, 
                                title=NULL,
                                implementation) {
  
  selected_values <- as.factor(selected_values)
  
  all_metrics_data <- data.frame()
  
  for (i in seq_along(selected_values)) {
    value <- selected_values[[i]]
    
    # Subset bootstrap object and stacked predictions based on column_name
    bootstrap_value <- bootstrap_object[[i]]
    num_datasets <- 1:length(bootstrap_value)
    cp_value <- stacked_predictions[stacked_predictions[[column_name]] == value, ]
    
    # Collect all data
    metrics_data <- do.call(rbind, lapply(num_datasets, function(d) {
      bootstrap_d <- bootstrap_value[[d]]
      cp_version <- unique(cp_value[cp_value$version == d,]$cp)
      
      metric_std_dev <- sapply(implementation, function(m) {
        sd(unlist(bootstrap_d$batch.metrics[, m]))
      })
      
      data.frame(
        Dataset = d,
        Group = factor(value),  
        Cp = cp_version,
        Metric = names(bootstrap_d$mean),
        Mean = as.numeric(bootstrap_d$mean),
        StdDev = as.numeric(metric_std_dev)
      )
    }))
    
    # Append to full dataset
    all_metrics_data <- rbind(all_metrics_data, metrics_data)
  }
  
  p <- ggplot(all_metrics_data, aes(x = Metric, y = Mean, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.6, position = position_dodge(0.8)) +  
    geom_jitter(aes(color = Group), alpha = 0.4, size = 1.8, position = position_dodge(0.8)) +  
    stat_summary(fun = "mean", geom = "point", color = "black", size = 2, position = position_dodge(0.8)) +
    stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(0.8)) + 
    theme_minimal() +
    labs(
      title = title,
      x = "Implementation",
      y = "C-index Value",
      fill = column_name 
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    scale_color_discrete(guide = "none") +
    ylim(y_limits[1], y_limits[2])
  
  print(p)
}

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

generate_synthetic_event_times <- function(survreg_model, 
                                           n = 1000, 
                                           covariates, 
                                           seed = 123, 
                                           verbose = FALSE) {
  set.seed(seed)
  
  total_n <- nrow(covariates)
  
  if (n < total_n) {
    covariates <- covariates[sample(1:total_n, n, replace = FALSE), ]
  } else {
    n <- total_n  # If n == total_n, no shuffling needed
  }
  # Maybe add replacement true for when >original size. Or Model covariates?
  
  ### Get parameters from the original model (this is AFT scale)
  mu <- survreg_model$coefficients[1]  # Intercept
  sigma <- survreg_model$scale  # Scale
  alpha <- survreg_model$coefficients[-1] # Regression coefficients
  
  ### PARAMETER TRANSFORMATION
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
  attr(synth_data, "lambda") <- lambda
  attr(synth_data, "gamma") <- gamma
  attr(synth_data, "beta") <- beta
  
  # Censoring percentage as attribute
  attr(synth_data, "censoring_percentage") <- censoring_percentage
  
  return(synth_data)
}

# Add censoring to the survival data
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
  
  set.seed(seed)
  
  # Get sample size
  n <- nrow(my_surv_data)
  
  # Generate censoring times
  
  ## Administrative censoring
  if(cens.type == "administrative") {
    censoring_times <- rep(cens.param$cens_limit_admin, times = n)
  }
  if(cens.type == "uniform") {
    # Lower bound of Uniform
    T_s <- my_surv_data$time
    
    C_min <- min(T_s)
    C_max <- quantile(T_s, 1 - cens.params$cens_increase_unif)
    
    censoring_times  <- runif(n, min = C_min, max = C_max)
    
    if (!is.null(cens.params$cens_limit_admin)) {
      censoring_times <- pmin(censoring_times, 
                              cens.params$cens_limit_admin)
    }
    # censoring_times <- runif(n, 
    #                          min = cens.params$min_cens_unif, 
    #                          max = cens.params$max_cens_unif)
  }
  if(cens.type == "weibull") {
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
    
    # Generate censoring times
    U_c <- runif(n)
    censoring_times <- (-log(U_c) / lambda_c)^(1/gamma_c)
    
    # add administrative censoring also 
    if (!is.null(cens.params$cens_limit_admin)) {
      censoring_times <- pmin(censoring_times, 
                              cens.params$cens_limit_admin)
    }
    
  } 
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

generate_modified_weibull_times <- function(survreg_model, 
                                            survreg_censoring_model,
                                            n = 1000, 
                                            covariates, 
                                            seed=123, 
                                            lambda_c_factor = 1,
                                            admin_censoring_times = NULL,
                                            verbose=FALSE) {
  set.seed(seed)
  
  total_n <- nrow(covariates)
  
  if (n < total_n) {
    covariates <- covariates[sample(1:total_n, n, replace = FALSE), ]
  } else {
    n <- total_n  # If n == total_n, no shuffling needed
  }
  # Maybe add replacement true for when >original size. Or Model covariates?
  
  ### EVENT
  # Fit the uncensoring data
  mu <- survreg_model$coefficients[1]  # Inntercept
  sigma <- survreg_model$scale  # Scale
  
  # Get confidence intervals
  conf_survreg <- confint(survreg_model) 
  mu_CI <- conf_survreg[1, ] 
  # Scale does not have confint:
  se_logsigma <- summary(survreg_model)$table["Log(scale)", "Std. Error"]
  # Estimate variance delta method
  var_sigma <- (sigma^2) * (se_logsigma^2)
  se_sigma <- sqrt(var_sigma)
  
  
  ### CENSORING
  # Fit censoring times
  mu_c <- survreg_censoring_model$coefficients[1]  # Intercept
  sigma_c <- survreg_censoring_model$scale  # Scale
  
  # Get confidence intervals
  conf_survreg_c <- confint(survreg_censoring_model)  
  mu_CI_c <- conf_survreg_c[1, ]  
  # Scale does not have confint:
  se_logsigma <- summary(survreg_censoring_model)$table["Log(scale)", "Std. Error"]
  # Estimate variance delta method
  var_sigma_c <- (sigma_c^2) * (se_logsigma^2)
  se_sigma_c <- sqrt(var_sigma_c)
  
  
  ### PARAMETER TRANSFORMATION
  # Covert into WeibullPH:
  lambda <- exp(-mu/sigma)
  gamma <- 1/sigma
  lambda_c <- exp(-mu_c/sigma_c) * lambda_c_factor
  gamma_c <- 1/sigma_c  
  
  
  # Confidence interval for gamma
  z <- qnorm(0.975)
  se_gamma <- (1 / sigma^2) * se_sigma
  gamma_CI <- c(gamma - z * se_gamma, gamma + z * se_gamma)
  
  se_gamma_c <- (1 / sigma_c^2) * se_sigma_c
  gamma_CI_c <- c(gamma_c - z * se_gamma_c, gamma_c + z * se_gamma_c)
  
  # Confidence interval for lambda
  #lambda_CI <- exp(-mu_CI / sigma)
  # Transform confidence intervals for WeibullPH
  #lambda_CI_c <- exp(-mu_CI_c / sigma_c)
  # Transform confidence intervals for WeibullPH
  lambda_CI <- get_lambda_CI(survreg_model)
  lambda_CI_c <- get_lambda_CI(survreg_censoring_model)
  
  n_coef <- dim(covariates)[2] + 1
  # Calculate survival times
  U <- runif(n)
  beta <-  -survreg_model$coefficients[2:n_coef]/sigma
  betas_cov <- as.matrix(covariates) %*% beta
  
  T_s <- (-log(U) / (lambda * exp(betas_cov)))^(1/gamma)
  
  # Calculate censoring times
  U_c <- runif(n)
  
  T_c <- (-log(U_c) / (lambda_c))^(1/gamma_c)
  
  if (verbose) {
    cat("Lambda C", lambda_c, "\n")
    cat("Gamma C", gamma_c, "\n")
    cat("Lambda T", lambda, "\n")
    cat("Gamma T", gamma, "\n")
    cat("Betas", beta, "\n")
  }
  
  # Get conficence intervals for beta
  conf_beta <- confint(survreg_model)[2:n_coef, ]
  beta_CI <- -conf_beta / sigma  
  beta_CI <- beta_CI[, c(2,1)]
  
  # Administrative censoring
  if (!is.null(admin_censoring_times)) {
    # for instance max(T) of the original dataset:
    T_c <- pmin(T_c, admin_censoring_times)
  } #?
  
  # Select minimum times
  T_obs      <- pmin(T_s, T_c)
  
  # Generate status based on censoring
  status_obs <- ifelse(T_s <= T_c, 1, 0) # <= or < ?
  
  # Calculate percentage of censoring
  censoring_percentage <- mean(status_obs == 0) * 100
  
  synth_data <- cbind(data.frame("time" = T_s, 
                                 "censoring_time" = T_c, 
                                 "observed_time" = T_obs,
                                 "status" = status_obs), covariates)
  # Set the attributes
  attr(synth_data, "lambda_T") <- lambda
  attr(synth_data, "lambda_Tc") <- lambda_c
  attr(synth_data, "gamma_T") <- gamma
  attr(synth_data, "gamma_Tc") <- gamma_c
  attr(synth_data, "beta") <- beta
  
  # Also the Ci
  attr(synth_data, "lambda_T_CI") <- lambda_CI$CI
  attr(synth_data, "lambda_Tc_CI") <- lambda_CI_c$CI
  attr(synth_data, "gamma_T_CI") <- gamma_CI
  attr(synth_data, "gamma_Tc_CI") <- gamma_CI_c
  attr(synth_data, "beta_CI") <- beta_CI
  
  # Censoring percentage as attribute
  attr(synth_data, "censoring_percentage") <- censoring_percentage
  
  return(synth_data)
  
}






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
  
get_model_preds2 <- function(stacked_predictions, 
                             model_names = "all", 
                             input_type = c("Distribution", "ExpectedMortality", "RiskAtT", "PartialHazAtT"),
                             specific_time = NULL, 
                             bootstrap_patient_ids = NULL) {
  
  input_type <- match.arg(input_type)
  base_cols <- c("test_time", "test_status", "patients_ids")
  df <- stacked_predictions  # work on a copy

  # Subset rows by bootstrap_patient_ids
  if (!is.null(bootstrap_patient_ids)) {
    match_ids <- match(bootstrap_patient_ids, df$patients_ids)
    df <- df[match_ids, , drop = FALSE]
  }
  
  # Detect model names if "all"
  if (identical(model_names, "all")) {
    all_cols <- names(df)
    if (input_type %in% c("Distribution", "RiskAtT")) {
      model_names <- unique(sub("\\..*", "", grep("^[A-Za-z]+\\.\\d+$", all_cols, value = TRUE)))
    } else if (input_type == "ExpectedMortality") {
      model_names <- unique(sub("Exp\\Mort\\.", "", grep("^Exp\\Mort\\.", all_cols, value = TRUE)))
    }
  }
  
  # Store final prediction columns
  all_model_cols <- list()
  
  for (model in model_names) {
    if (input_type == "Distribution") {
      pattern <- paste0("^", model, "\\.")
      cols <- grep(pattern, names(df), value = TRUE)
      all_model_cols[[length(all_model_cols) + 1]] <- df[, cols, drop = FALSE]
      
    } else if (input_type == "ExpectedMortality") {
      col_name <- paste0("ExpMort.", model)
      if (!col_name %in% names(df)) {
        stop("Missing column: ", col_name)
      }
      all_model_cols[[model]] <- df[, col_name, drop = FALSE]
      
    } else if (input_type == "RiskAtT") {
      if (is.null(specific_time)) {
        stop("specific_time must be provided for input_type = 'RiskAtT'")
      }

      cols <- paste0(model, ".", specific_time)
      
      missing_cols <- setdiff(cols, names(df))
      if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "))
      }
      
      # Transform survival to risk
      transformed <- lapply(cols, function(col) 1 - df[[col]])
      names(transformed) <- cols
      all_model_cols[[length(all_model_cols) + 1]] <- as.data.frame(transformed)
    } else if (input_type == "PartialHazAtT") {
      if (is.null(specific_time)) {
        stop("specific_time must be provided for input_type = 'PartialHazAtT'")
      }
      cols <- paste0(model, ".", specific_time)
      missing_cols <- setdiff(cols, names(df))
      if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "))
      }
      
      # Transform survival to risk
      transformed <- lapply(cols, function(col)  - log(df[[col]]))
      names(transformed) <- cols
      all_model_cols[[length(all_model_cols) + 1]] <- as.data.frame(transformed)
    }
  }
  
  # Combine everything to return
  predictions <- do.call(cbind, all_model_cols)
  final_df <- cbind(df[, base_cols, drop = FALSE], predictions)
  
  return(final_df)
}

make_risk_plot_entries <- function(results_list) {
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
      mean_c <- result$mean[i]
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

make_risk_table <- function(results_list) {
  risk_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    model <- sub("\\..*", "", name)
    t <- sub(".*\\.", "", name)
    tau <- result$eval.times
    metrics <- names(result$mean)
    
    entries <- lapply(seq_along(metrics), function(i) {
      metric_name <- metrics[i]
      mean_c <- result$mean[i]
      ci <- result$confidence.intervals[i, ]
      value <- sprintf("%.3f [%.3f–%.3f]", mean_c, ci[1], ci[2])
      
      data.frame(
        Metric = metric_name,
        Model = model,
        InputType = paste0("Risk(t=", t, ", tau=", tau, ")"),
        Value = value,
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, entries)
  })
  
  risk_df <- do.call(rbind, risk_entries)
  
  # Reshape for table display
  df_wide <- tidyr::pivot_wider(
    risk_df,
    names_from = Model,
    values_from = Value
  )
  
  return(df_wide)
}


make_expm_plot_entries <- function(results_list) {
  expm_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    model <- sub("ExpMort\\.", "", name)
    metrics <- names(result$mean)
    
    entries <- lapply(seq_along(metrics), function(i) {
      metric_name <- metrics[i]
      mean_c <- result$mean[i]
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

make_expm_table <- function(results_list) {
  expm_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    model <- sub("ExpMort\\.", "", name)
    metrics <- names(result$mean)
    
    entries <- lapply(seq_along(metrics), function(i) {
      metric_name <- metrics[i]
      mean_c <- result$mean[i]
      ci <- result$confidence.intervals[i, ]
      value <- sprintf("%.3f [%.3f–%.3f]", mean_c, ci[1], ci[2])
      
      data.frame(
        Metric = metric_name,
        Model = model,
        Value = value,
        InputType = "Exp.Mort",
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, entries)
  })
  
  expm_df <- do.call(rbind, expm_entries)
  
  # Reshape for table display
  df_wide <- tidyr::pivot_wider(
    expm_df,
    names_from = Model,
    values_from = Value
  )
  
  return(df_wide)
}

make_surv_table <- function(results_list) {
  surv_entries <- lapply(names(results_list), function(name) {
    if (name == "batch.metrics") return(NULL)
    
    result <- results_list[[name]]
    model <- name
    metrics <- names(result$mean)
    
    entries <- lapply(1:2, function(i) {
      mean_c <- result$mean[i]
      ci <- result$confidence.intervals[i, ]
      value <- sprintf("%.3f [%.3f–%.3f]", mean_c, ci[1], ci[2])
      
      data.frame(
        Metric = metrics[i],
        InputType = "Distrib.",
        Model = model,
        Value = value,
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, entries)
  })
  
  surv_df <- do.call(rbind, surv_entries)
  
  df_wide <- tidyr::pivot_wider(
    surv_df,
    names_from = Model,
    values_from = Value
  )
  
  return(df_wide)
}

plot_survival_curves <- function(surv_matrix, patient_ids = NULL, seed = 42, 
                                 title = "Survival Curves") {
  if (is.null(patient_ids)) {
    set.seed(seed)
    patient_ids <- sample(1:nrow(surv_matrix), 5)
  }
  
  subset_curves <- surv_matrix[patient_ids, , drop = FALSE]
  df <- as.data.frame(subset_curves)
  df$patient_id <- paste0("Patient_", patient_ids)
  df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
  
  df_long <- reshape2::melt(df, id.vars = "patient_id", variable.name = "time", value.name = "surv_prob")
  df_long$time <- as.numeric(as.character(df_long$time))
  
  # Plot
  ggplot(df_long, aes(x = time, y = surv_prob, color = patient_id)) +
    geom_line(linewidth = 1) +
    labs(title = title, x = "Time", y = "Survival Probability") +
    theme_minimal()
  
}

