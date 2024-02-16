#this imports dependencies using roxygen
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import collapse
#' @import mgcv
#' @import Matrix
#' @import sandwich

##############################################################
# Estimate Propensity Score
##############################################################
estimate_propensity_score <- function(train_data, predict_data, return_graph = TRUE, k = -1, weights = NULL) {
  
  ##############################################################
  # Estimate First Stage 
  #--> Use cubic kernel given importance of cutoff 
  ##############################################################
  #p_hat_below <- gam(treat ~ s(rv, k = k, bs= "bs", m = c(3,2)), data = train_data %>% filter(above_cut == FALSE), weights = weights)
  #data_below <- predict_data %>% filter(above_cut == FALSE) %>% mutate(p_hat = predict(p_hat_below, newdata = predict_data %>% filter(above_cut == FALSE)))
  
  #p_hat_above <- gam(treat ~ s(rv, k = k, bs= "bs", m = c(3,2)), data = train_data %>% filter(above_cut == TRUE), weights = weights)
  #data_above<- predict_data %>% filter(above_cut == TRUE) %>% mutate(p_hat = predict(p_hat_above, newdata = predict_data %>% filter(above_cut == TRUE)))
  
  #predict_data <- rbind(data_below, data_above)
  p_hat_results <- gam(treat ~ s(rv, k = -1, bs = "bs", m = c(3, 2), by = sections), data = train_data)
  predict_data$p_hat <- predict(p_hat_results, newdata = predict_data)
  
  ##############################################################
  # Return Graph
  ##############################################################
  if (return_graph) {
    #below_predictions <- predict(p_hat_below, newdata = predict_data %>% filter(above_cut == FALSE), se.fit = TRUE)
    #data_below <- cbind(data_below, below_predictions$se.fit) 
    #data_below <- data_below %>% rename(se = "below_predictions$se.fit")

    #above_predictions <- predict(p_hat_above, newdata = predict_data %>% filter(above_cut == TRUE), se.fit = TRUE)
    #data_above <- cbind(data_above, above_predictions$se.fit) 
    #data_above <- data_above %>% rename(se = "above_predictions$se.fit")
    
    #graph_data <- rbind(data_below, data_above)
    predictions <- predict(p_hat_results, newdata = predict_data, se.fit = TRUE)
    graph_data <- cbind(predict_data, predictions$se.fit) 
    graph_data <- graph_data %>% rename(se = "predictions$se.fit")
    
    graph_data <- graph_data %>% mutate(p95 = p_hat + 1.96*se, p5 = p_hat - 1.96*se)
    
    graph <- ggplot(graph_data) + theme_bw() + geom_line(aes(x = rv, y = p_hat, color = sections), linewidth = 1) + 
    geom_line(aes(x = rv, y = p95, group = sections), color = 'black', linetype = 'dashed') +
      geom_line(aes(x = rv, y = p5, group = sections), color = 'black', linetype = 'dashed') +
      xlab("Running Variable") + ylab("Propensity Score") + theme(legend.position="none")  

  }
  
  ##############################################################
  # Return
  ############################################################## 
  if (return_graph) {
  return(list(predict_data, graph))
  }
  else{
    return(predict_data)
  }
}


##############################################################
# Program with Known Propensity Score
##############################################################
global_rdd_known_pscore <- function(data, estimator = "ATE", covs = NULL, c= 0, return_graph = TRUE, linear_penalty = FALSE, sp = NULL, linear_sp = -1,  vcv_cluster = FALSE, weights = NULL) {
  
  ##############################################################
  # Replace c with min(c) in case there are multiple discontinuities
  ##############################################################
  full_c = c
  c = min(c)
  
  ##############################################################
  # Estimate Spline
  ##############################################################
  # Create formula, potentially with covariates
  if (is.null(covs) == FALSE) {
    # Create formula
    formula <- as.formula(paste("Y ~ s(rv) + treat + ti(rv, by = treat) + p_hat + p_hat*treat +", paste(colnames(covs), collapse = "+")))
  }
  else {
    formula <- as.formula("Y ~ s(rv) + treat + ti(rv, by = treat) + p_hat + p_hat*treat")
  }
  
  # Estimate
  if (linear_penalty == FALSE) {
    m <- gam(formula, data = data, sp = sp, weights = weights)
  }
  if (linear_penalty) {
    m <- gam(formula, data = data, paraPen = list(p_hat=list(diag(1), sp = linear_sp), "treat:p_hat" = list(diag(1), sp = linear_sp)), sp = sp, weights = weights)
  }
  
  # Save Model Matrix
  m_cov_matrix <- predict(m, newdata = data %>% filter(treat == 1), type = "lpmatrix")
  m_cov_matrix <- as_tibble(m_cov_matrix)

  # Counts
  #rv_smoother_num <- ncol(m_cov_matrix %>% select(starts_with('s(rv).')))
  #treat_tv_smoother_num <- ncol(m_cov_matrix %>% select(starts_with('ti(rv):treat')))
  #stopifnot("The number of coefficients is wrong." = (4 + rv_smoother_num + treat_tv_smoother_num == nrow(as.matrix(m$coefficients))))
  
  # Replace intercept and s(rv) variables with 0 b/c they're not the treatments
  m_cov_matrix <- m_cov_matrix %>% mutate('(Intercept)' = 0)
  m_cov_matrix <- m_cov_matrix %>% mutate(across(starts_with('s(rv).'), ~if_else(is.na(.x),0,0)))
  
  # Adjust p_hat for selection
  if (estimator == "ATE") {
    m_cov_matrix <- m_cov_matrix %>% mutate(p_hat = 1, 'treat:p_hat' = 1)
  }
  if (estimator == "ATT") {
    m_cov_matrix <- m_cov_matrix %>% mutate(p_hat = 1)
  }
  # if (estimator == "ATE") {
  #    m_cov_matrix <- m_cov_matrix %>% mutate(p_hat = 0, 'treat:p_hat' = 1)
  # }
  # if (estimator == "ATT") {
  #    m_cov_matrix <- m_cov_matrix %>% mutate(p_hat = 1-p_hat)
  # }
  
  # If covs, replace covariates with their mean value
  if (is.null(covs) == FALSE) {
    for (cov_name in colnames(covs)) {
      m_cov_matrix <- m_cov_matrix %>% mutate(!!cov_name := mean(!!sym(cov_name)))
    }
  }
  
  ##############################################################
  # Estimate variance matrix if cluster
  ##############################################################
  if (vcv_cluster == FALSE) {
    V <- m$Ve
  }
  if (vcv_cluster == TRUE) {
    # Full m_cov 
    full_m_cov_matrix <- predict(m, type = "lpmatrix")
    
    # Penalty matrix
    if (is.null(covs)) {
      if (linear_penalty == FALSE) {
        if (is.null(sp)) {
        p_matrix = bdiag(matrix(0, 4,4), m$sp[[1]]*as.matrix(m$smooth[[1]]$S[[1]]), m$sp[[2]]*as.matrix(m$smooth[[2]]$S[[1]]))
        }
        if (!is.null(sp)) {
          p_matrix = bdiag(matrix(0, 4,4), sp[1]*as.matrix(m$smooth[[1]]$S[[1]]), sp[2]*as.matrix(m$smooth[[2]]$S[[1]]))
        }
      }
      if (linear_penalty == TRUE) {
        if (linear_sp == -1) {
          if (is.null(sp)) {
          p_matrix = bdiag(matrix(0, 2,2), as.matrix(m$sp[[1]]), as.matrix(m$sp[[2]]), m$sp[[3]]*as.matrix(m$smooth[[1]]$S[[1]]), m$sp[[4]]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
          if (!is.null(sp)) {
            p_matrix = bdiag(matrix(0, 2,2), as.matrix(m$sp[[1]]), as.matrix(m$sp[[2]]), sp[1]*as.matrix(m$smooth[[1]]$S[[1]]),sp[2]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
        }
        else {
          if (is.null(sp)) {
          p_matrix = bdiag(matrix(0, 2,2), as.matrix(linear_sp), as.matrix(linear_sp), m$sp[[1]]*as.matrix(m$smooth[[1]]$S[[1]]), m$sp[[2]]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
          if (!is.null(sp)) {
            p_matrix = bdiag(matrix(0, 2,2), as.matrix(linear_sp), as.matrix(linear_sp), sp[1]*as.matrix(m$smooth[[1]]$S[[1]]), sp[2]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
        }
      }
    }
    if (is.null(covs) == FALSE) {
      if (linear_penalty == FALSE) {
        if (is.null(sp)) {
          p_matrix = bdiag(matrix(0, 4 + ncol(covs),4 + ncol(covs)), m$sp[[1]]*as.matrix(m$smooth[[1]]$S[[1]]), m$sp[[2]]*as.matrix(m$smooth[[2]]$S[[1]]))
        }
        if (!is.null(sp)) {
          p_matrix = bdiag(matrix(0, 4 + ncol(covs),4 + ncol(covs)), sp[1]*as.matrix(m$smooth[[1]]$S[[1]]), sp[2]*as.matrix(m$smooth[[2]]$S[[1]]))
        }
      }
      if (linear_penalty == TRUE) {
        if (linear_sp == -1) {
          if (is.null(sp)) {
            p_matrix = bdiag(matrix(0, 2,2), as.matrix(m$sp[[1]]), matrix(0, ncol(covs),ncol(covs)), as.matrix(m$sp[[2]]), m$sp[[3]]*as.matrix(m$smooth[[1]]$S[[1]]), m$sp[[4]]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
          if (!is.null(sp)) {
            p_matrix = bdiag(matrix(0, 2,2), as.matrix(m$sp[[1]]), matrix(0, ncol(covs),ncol(covs)), as.matrix(m$sp[[2]]), sp[1]*as.matrix(m$smooth[[1]]$S[[1]]), sp[2]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
        }
        else {
          if (is.null(sp)) {
            p_matrix = bdiag(matrix(0, 2,2), as.matrix(linear_sp), matrix(0, ncol(covs),ncol(covs)), as.matrix(linear_sp), m$sp[[1]]*as.matrix(m$smooth[[1]]$S[[1]]), m$sp[[2]]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
          if (!is.null(sp)) {
            p_matrix = bdiag(matrix(0, 2,2), as.matrix(linear_sp), matrix(0, ncol(covs),ncol(covs)), as.matrix(linear_sp), sp[1]*as.matrix(m$smooth[[1]]$S[[1]]), sp[2]*as.matrix(m$smooth[[2]]$S[[1]]))
          }
        }
      }
    }
    
    # Calculate "Bread"
    inv <- solve(t(as.matrix(full_m_cov_matrix)) %*% as.matrix(full_m_cov_matrix) + p_matrix)
    
    # Run linear model (unpenalized) & calculate "meat" matrix
    reg <- lm(data$Y ~ 0 + full_m_cov_matrix)
    meat_mat = nrow(data)*vcovCL(reg, cluster = data$cluster_id, sandwich = FALSE)
    
    # Estimate Variance
    V <- inv %*% meat_mat %*% inv
    
  }
  
  ##############################################################
  # Test Restricted Models
  ############################################################## 
  # No Endogenous Selection
  R <- matrix(0, nrow = 2, ncol = nrow(as.matrix(m$coefficients)))
  R[1, 3] <- 1
  R[2, 4] <- 1
  wald_stat <- t(R %*% as.matrix(m$coefficients)) %*% solve(R %*% (V + diag(x = 1e-7, nrow(V), ncol(V)))%*% t(R)) %*% R %*% as.matrix(m$coefficients)
  no_endog_selection <- 1 - pchisq(wald_stat[1,1], df = 2)
  
  # No Treatment Effects
  R <- matrix(0, nrow = 2 + m$smooth[[2]]$df, ncol = nrow(as.matrix(m$coefficients)))
  R[1, 2] <- 1
  R[2, 4] <- 1
  for (i in seq(3, 3 + m$smooth[[2]]$df - 1, 1)) {
    j <- (i-3) + m$smooth[[2]]$first.para
    R[i, j] <- 1    
  }
  wald_stat <- t(R %*% as.matrix(m$coefficients)) %*% solve(R %*% (V + diag(x = 1e-7, nrow(V), ncol(V)))%*% t(R)) %*% R %*% as.matrix(m$coefficients)
  no_treatment_effects <- 1 - pchisq(wald_stat[1,1], df = nrow(R))
  
  # No Treatment Effect Heterogeneity
  R <- matrix(0, nrow = 1 + m$smooth[[2]]$df, ncol = nrow(as.matrix(m$coefficients)))
  R[1, 4] <- 1
  for (i in seq(2, 2 + m$smooth[[2]]$df - 1, 1)) {
    j <- (i-2) + m$smooth[[2]]$first.para
    R[i, j] <- 1    
  }
  wald_stat <- t(R %*% as.matrix(m$coefficients)) %*% solve(R %*% (V + diag(x = 1e-7, nrow(V), ncol(V)))%*% t(R)) %*% R %*% as.matrix(m$coefficients)
  no_te_heterogeneity <- 1 - pchisq(wald_stat[1,1], df = nrow(R))
  
  # Save
  p_test_results <- c(no_endogenous_selection = no_endog_selection, no_tau = no_treatment_effects, constant_tau = no_te_heterogeneity, tau_indep_of_eta = summary(m)$p.pv[4],tau_indep_of_rv = summary(m)$s.pv[2])
  
  ##############################################################
  # Calculate Conditional ATEs
  ##############################################################
  # CATE estimates
  cate_hat <- as.matrix(m_cov_matrix) %*% m$coefficients
  
  ##############################################################
  # Calculate Overall Avg. Effect
  ##############################################################
  # Calculate weights 
  if (estimator == "ATE") {
    out <- data %>% filter(treat == 1) %>% select(rv, p_hat)
    aggregation_weights <- out$p_hat^(-1)
  }
  if (estimator == "ATT") {
    aggregation_weights <- matrix(1, nrow = nrow(m_cov_matrix))
  }
  aggregation_weights <- aggregation_weights/sum(aggregation_weights)
  
  # Overall ATE
  avg_effect <- t(as.matrix(aggregation_weights)) %*% as.matrix(cate_hat)
  
  # Overall ATE Standard Error
  avg_effect_se <- t(as.matrix(aggregation_weights)) %*% as.matrix(m_cov_matrix) %*% V %*% t(as.matrix(m_cov_matrix)) %*% as.matrix(aggregation_weights)
  avg_effect_se <- avg_effect_se^.5
  
  ##############################################################
  # Calculate LATE via RD
  ##############################################################
  # Calculate weights -- equal to zero for all those not
  out <- data %>% filter(treat == 1) %>% select(rv, p_hat)
  late_output <- cbind(out, cate_hat) %>%
    mutate(agg_weights_above = (rv > c & rv == min(out %>% filter(rv > c) %>% select(rv)))) %>%
    mutate(agg_weights_below = (rv < c & rv == max(out %>% filter(rv < c) %>% select(rv))))
  late_output <- late_output %>% mutate(agg_weights_above = agg_weights_above/sum(late_output$agg_weights_above)) %>%
    mutate(agg_weights_below = agg_weights_below/sum(late_output$agg_weights_below))
  late_output <- late_output %>% mutate(late_aggregation_weights = pmax(agg_weights_below, agg_weights_above)*p_hat*(1 - 2*(rv < c)))

  aggregation_weights = late_output$late_aggregation_weights

  # Calculate p_hat jump
  just_above <- out %>% filter(rv > c & rv == min(out %>% filter(rv > c) %>% select(rv))) %>% select(p_hat)
  just_below <- out %>% filter(rv < c & rv == max(out %>% filter(rv < c) %>% select(rv))) %>% select(p_hat)
  p_jump <- just_above[1,] - just_below[1,]

  # Overall LATE
  late_effect <- 1/p_jump*t(as.matrix(aggregation_weights)) %*% as.matrix(cate_hat)

  # Overall LATE Standard Error
  late_effect_se <- (1/p_jump)^2*t(as.matrix(aggregation_weights)) %*% as.matrix(m_cov_matrix) %*% V %*% t(as.matrix(m_cov_matrix)) %*% as.matrix(aggregation_weights)
  late_effect_se <- avg_effect_se^.5

  ##############################################################
  # Adjust outputted dataset
  #--> Mostly to reduce size
  ##############################################################
  # Filter to reduce size
  n_samples <- min(1000, nrow(m_cov_matrix))
  subsample <- cbind(data %>% filter(treat == 1) %>% select(rv), m_cov_matrix) %>% arrange(rv) %>% 
    slice(seq(1,nrow(m_cov_matrix), round(nrow(m_cov_matrix)/n_samples)))
  m_cov_subsample <- subsample %>% select(!c(rv))
  rvphat_subsample <- data %>% filter(treat == 1) %>% select(rv, p_hat) %>% arrange(rv) %>% 
    slice(seq(1,nrow(m_cov_matrix), round(nrow(m_cov_matrix)/n_samples)))
  
  # CATE short
  cate_hat <- as.matrix(m_cov_subsample) %*% m$coefficients
  
  # Standard Error
  cate_var <- as.matrix(m_cov_subsample) %*% V %*% t(as.matrix(m_cov_subsample))
  cate_hat_se <- diag(cate_var)
  cate_hat_se <- cate_hat_se^.5
  
  # Combine results into a tibble
  output <- cbind(rvphat_subsample, cate_hat, cate_hat_se)
  
  # Add 95 ptiles
  output <- output %>% mutate(p5 = cate_hat - cate_hat_se*1.65, p95 = cate_hat + cate_hat_se*1.65)

  ##############################################################
  # Potentially Graph Results
  ##############################################################  
  if (return_graph) {
    # Return Predictions
    predicted_outcomes <- predict(m, se.fit = TRUE)
    predicted_outcomes <- tibble(data, predicted_outcomes$fit, predicted_outcomes$se.fit) %>% rename(p_y = 'predicted_outcomes$fit', p_y_se = 'predicted_outcomes$se.fit')
    predicted_outcomes <- predicted_outcomes %>% mutate(treated = (treat > .5), p5 = p_y - p_y_se*1.65, p95 = p_y + p_y_se*1.65)
    
    # Graph Predicted Outcomes
    outcome_graph <- ggplot(predicted_outcomes) + theme_bw() + geom_line(aes(x = rv, y = p_y, group = interaction(sections, treated), colour = treated), size = 1) + 
      geom_line(aes(x = rv, y = p95, group = interaction(sections, treated)), color = 'black', linetype = 'dashed') +
      geom_line(aes(x = rv, y = p5, group = interaction(sections, treated)), color = 'black', linetype = 'dashed') +
      xlab("Running Variable") + ylab("Estimated Conditional Outcome") + theme(legend.position="none")
    
    if (estimator == "ATT") {
      # Allow for discontinuity at cutoff due to a changing treatment thresholds
      #output <- output %>% mutate(above_cut = (rv > c))
      output$sections <- cut(output$rv, sort(append(c(-Inf, Inf), full_c)))
      
      graph <- ggplot(output) + theme_bw() + geom_line(aes(x = rv, y = cate_hat, group = sections), linewidth = 1) + 
        geom_line(aes(x = rv, y = p95, group = sections), color = 'black', linetype = 'dashed') +
        geom_line(aes(x = rv, y = p5, group = sections), color = 'black', linetype = 'dashed') +
        xlab("Running Variable") + ylab("Estimated Conditional Average Effect") + theme(legend.position="none")
      
      # Variation in average effect holding treatment thresholds fixed
      adjusted_variation <- predict(m, newdata = data %>% filter(treat == 1), type = "terms", se.fit = TRUE)
      fixed_threshold_data <- tibble(data %>% filter(treat == 1) %>% select(rv), 
                                    as_tibble(adjusted_variation$fit) %>% select("ti(rv):treat"), 
                                    as_tibble(adjusted_variation$se.fit) %>% select("ti(rv):treat"),
                                    .name_repair = ~c('rv','estimate','se'))
      # Set mean to avg_effect & add 95% CI
      fixed_threshold_data <- fixed_threshold_data %>% mutate(estimate = estimate + matrix(avg_effect, nrow = nrow(fixed_threshold_data)))
      fixed_threshold_data <- fixed_threshold_data %>% mutate(p5 = estimate - se*1.65, p95 = estimate + se*1.65)
      graph_fixed_threshold <- ggplot(fixed_threshold_data) + theme_bw() + geom_line(aes(x = rv, y = estimate), size = 1) + 
        geom_line(aes(x = rv, y = p95), color = 'black', linetype = 'dashed') +
        geom_line(aes(x = rv, y = p5), color = 'black', linetype = 'dashed') +
        xlab("Running Variable") + ylab("Estimated Conditional Average Effect") + theme(legend.position="none")
      
    }
    if (estimator == "ATE") {
      graph <- ggplot(output) + theme_bw() + geom_line(aes(x = rv, y = cate_hat), size = 1) + 
        geom_line(aes(x = rv, y = p95), color = 'black', linetype = 'dashed') +
        geom_line(aes(x = rv, y = p5), color = 'black', linetype = 'dashed') +
        xlab("Running Variable") + ylab("Estimated Conditional Average Effect") + theme(legend.position="none") 
    }
  }
  
  ##############################################################
  # Return
  ############################################################## 
  if (return_graph == FALSE) {
  return(list(avg_effect, avg_effect_se, late_effect, late_effect_se, tibble(output), cate_var, m,  p_test_results))
  }
  if (return_graph) {
    if (estimator == "ATE") {
     return(list(avg_effect, avg_effect_se, late_effect, late_effect_se, tibble(output), cate_var, m, p_test_results, outcome_graph, graph))
    }
    if (estimator == "ATT") {
      return(list(avg_effect, avg_effect_se, late_effect, late_effect_se, tibble(output), cate_var, m,  p_test_results, outcome_graph, graph, graph_fixed_threshold))
    }
  }
}


##############################################################
# Full Program
##############################################################
#' Estimate A Global RDD 
#' 
#' This function implements the Global RDD method, as described in Opper and Ozek (2024). 
#' 
#' @param outcome A vector indicating the value of to the outcome of interest for each observation
#' @param treatment A vector indicating whether someone is treated (with a value of one) or not (with a value of zero) for each observation
#' @param running_variable A vector indicating the value of the running variable for each observation
#' @param c Where is the discontinuity in the running varaible?. Defaults to 0
#' @param split_esitmation Do we randomly split the sample and estimate the propensity score and Global RDD on different samples?
#' @param return_graph Does the function return graphs plotting the results? Defaults to TRUE
#' @param cluster_id An optional argument indicating which cluster the different observations are in, for use in calculating the standard errors. Defaults to Null
#' @param covs Optional argument with list of column names that correspond to the covariates that are included. Defaults to NULL
#' @param p_score_k Optimal argument indicating the number of knots in the propensity score estimation. Default is -1, which means it's automatically selected.
#' @param estimator Do we return Conditional ATE or Conditional ATTs? Defaults to ATE
#' @param linear_penalty Do we shrink the linear terms as well? Defaults to FALSE, which means they are not shrunken. 
#' @param sp What is the smoothing parameter? Defaults to Null, which means the smoothing parameters are automatically selected
#' @param linear_sp Do we shrink the linear terms as well? Defaults to FALSE, which means they are not shrunken. 
#' @param weights An optional argument indicating the weights. Defaults to NULL which means each observation receives the same weight.
#' @returns The code returns a variety of graphs illustrating the results (if return_graph == TRUE), measures of the ATE or ATT and its standard error as well as a table indicating the CATE or CATT and the standard errors. It also returns the output from teh estimated GAM and the the size of the first stage jump.  
#' @keywords GlobalRDD
#' @export
global_rdd <- function(outcome, treatment, running_variable, c = 0, split_estimation = FALSE, return_graph = TRUE, cluster_id = NULL, covs = NULL, p_score_k = -1, estimator = "ATE", linear_penalty = FALSE, sp = NULL, linear_sp = -1, weights = NULL) {
  
  ##############################################################
  # Set-Up Data
  ##############################################################
  # Above or below threshold
  #above_cut <- (running_variable > c)
    
  if (is.null(cluster_id)) {
    # Combine into new dataframe
    data <- data.frame(
      Y = outcome,
      treat = treatment,
      rv = running_variable,
      #above_cut = above_cut
    )
    
    # VCV_cluster variable
    vcv_cluster = FALSE
  }
  if (is.null(cluster_id) == FALSE) {
    # Combine into new dataframe
    data <- data.frame(
      Y = outcome,
      treat = treatment,
      rv = running_variable,
      #above_cut = above_cut,
      cluster_id = cluster_id
    )
    
    # VCV_cluster
    vcv_cluster = TRUE
  }
  
  # Subset data based on the cutoffs
  data$sections <- cut(data$rv, sort(append(c(-Inf, Inf), c)))
  
  # Add weights, if any
  if (is.null(weights) == FALSE) {
    data$weights = weights
  }
  
  # Add covs, if any
  if (is.null(covs) == FALSE) {
    data <- cbind(data, covs)
  }
  
  # If linear sp penalty is zero, hard code no penalty
  if (linear_sp == 0) {
    linear_penalty = FALSE
  }
  
  # Make sure treatment status is numeric
  data$treat <- as.numeric(data$treat)
  
  # Drop missing data
  data <- data %>% drop_na()
  
  # Replace covs matrix with version that drops missing values
  if (is.null(covs) == FALSE) {
    covs <- data %>% select(colnames(covs))
  }
  
  if (split_estimation == FALSE) {
    ##############################################################
    # Estimate propensity score
    ##############################################################  
    # Estimate PScore
    propensity_score_estimation <- estimate_propensity_score(data, data, return_graph = return_graph, k = p_score_k, weights = weights)
    
    # Save data & potentially graph
    if (return_graph) {
      data <- propensity_score_estimation[[1]]
      propensity_score_graph <- propensity_score_estimation[[2]]
    }
    else {
      data <- propensity_score_estimation
    }
    
    # Calculate first-stage jump
    #first_stage_jump = abs(min(data %>% filter(above_cut == FALSE) %>% top_n(1, rv) %>% select(p_hat)) - min(data %>% filter(above_cut == TRUE) %>% top_n(-1, rv) %>% select(p_hat)))
    
    ##############################################################
    # Run with known propensity score
    ##############################################################
    estimates <- global_rdd_known_pscore(data, c= c, covs = covs, return_graph = return_graph, estimator = estimator, linear_penalty = linear_penalty, linear_sp = linear_sp, sp = sp, vcv_cluster = vcv_cluster, weights = weights)
    
    ##############################################################
    # Return
    ############################################################## 
    avg_effect <- estimates[[1]]
    avg_effect_se <- estimates[[2]]
    late <- estimates[[3]]
    late_se <- estimates[[4]]
    output <- estimates[[5]]
    output_var <- estimates[[6]]
    m <- estimates[[7]]
    p_test_results <- estimates[[8]]
    if (return_graph) {
      outcome_graph <- estimates[[9]]
      effect_graph <- estimates[[10]]
      if (estimator == "ATT") {
        effect_graph_fixed_threshold <- estimates[[11]]
      }
    }
  }
  else {
      ##############################################################
      # Split into two
      ##############################################################
      sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.5,0.5))
      data1 <- data[sample,]
      data2 <- data[!sample,]

      ##############################################################
      # Estimate propensity score
      ##############################################################
      propensity_score_estimation1 <- estimate_propensity_score(data2, data1, return_graph = return_graph, k = p_score_k, weights = weights)
      propensity_score_estimation2 <- estimate_propensity_score(data1, data2, return_graph = return_graph, k = p_score_k, weights = weights)
      if (return_graph) {
        data1 <- propensity_score_estimation1[[1]]
        data2 <- propensity_score_estimation2[[1]]
        # Only return one of the graphs
        propensity_score_graph <- propensity_score_estimation1[[2]]
      }
      else {
        data1 <- propensity_score_estimation1
        data2 <- propensity_score_estimation2
      }
      
      # Calculate first-stage jump as average of two jumps
      #first_stage_jump1 = abs(min(data1 %>% filter(above_cut == FALSE) %>% top_n(1, rv) %>% select(p_hat)) - min(data1 %>% filter(above_cut == TRUE) %>% top_n(-1, rv) %>% select(p_hat)))
      #first_stage_jump2 = abs(min(data2 %>% filter(above_cut == FALSE) %>% top_n(1, rv) %>% select(p_hat)) - min(data2 %>% filter(above_cut == TRUE) %>% top_n(-1, rv) %>% select(p_hat)))
      #first_stage_jump = (first_stage_jump1 + first_stage_jump2)/2
      
      ##############################################################
      # Run with known propensity score
      ##############################################################
      estimates1 <- global_rdd_known_pscore(data1, c= c, covs = covs, return_graph = FALSE, estimator = estimator, linear_penalty = linear_penalty, sp = sp, weights = weights)
      estimates2 <- global_rdd_known_pscore(data2, c= c, covs = covs, return_graph = FALSE, estimator = estimator, linear_penalty = linear_penalty, sp = sp, weights = weights)
      
      # Just return one of the gam results
      m <- estimates1[[7]]
      
      ##############################################################
      # Combine Results
      ##############################################################
      avg_effect <- .5*(estimates1[[1]] + estimates2[[1]])
      avg_effect_se <- .25*(estimates1[[2]]^2 + estimates2[[2]]^2)^.5
      late <- .5*(estimates1[[3]] + estimates2[[3]])
      late_se <- .25*(estimates1[[4]]^2 + estimates2[[4]]^2)^.5
      output <- data.frame(
        rv = (estimates1[[5]]$rv + estimates2[[5]]$rv)*.5,
        p_hat = (estimates1[[5]]$p_hat + estimates2[[5]]$p_hat)*.5,
        cate_hat = (estimates1[[5]]$cate_hat + estimates2[[5]]$cate_hat)*.5,
        cate_hat_se = .25*(estimates1[[5]]$cate_hat^2 + estimates2[[5]]$cate_hat^2)^.5
      )
      output <- output %>% mutate(p5 = cate_hat - cate_hat_se*1.65, p95 = cate_hat + cate_hat_se*1.65)
      p_test_results <- .5*(estimates1[[8]] + estimates2[[8]])
      
      ##############################################################
      # Potentially Graph Results
      ##############################################################
      if (return_graph) {
        if (estimator == "ATT") {
          # Allow for discontinuity at cutoff due to a changing treatment thresholds
          output$sections <- cut(data$rv, sort(append(c(-Inf, Inf), c)))
          #output <- output %>% mutate(above_cut = (rv > c))
          effect_graph <- ggplot(output) + theme_bw() + geom_line(aes(x = rv, y = cate_hat, group = sections), linewidth = 1) + 
            geom_line(aes(x = rv, y = p95, group = sections), color = 'black', linetype = 'dashed') +
            geom_line(aes(x = rv, y = p5, group = sections), color = 'black', linetype = 'dashed') +
            xlab("Running Variable") + ylab("Estimated Conditional Average Effect") + theme(legend.position="none")
          
        }
        if (estimator == "ATE") {
          effect_graph <- ggplot(output) + theme_bw() + geom_line(aes(x = rv, y = cate_hat), linewidth = 1) + 
            geom_line(aes(x = rv, y = p95), color = 'black', linetype = 'dashed') +
            geom_line(aes(x = rv, y = p5), color = 'black', linetype = 'dashed') +
            xlab("Running Variable") + ylab("Estimated Conditional Average Effect") + theme(legend.position="none") 
        }
      }
  }

  ##############################################################
  # Return
  ############################################################## 
  if (return_graph == FALSE) {
    #return(list(first_stage_jump = first_stage_jump, avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m))
    return(list(avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m, p_test_results = p_test_results))
  }
  if (return_graph == TRUE & estimator == "ATE") {
    #return(list(first_stage_jump = first_stage_jump, avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m, first_stage_graph = propensity_score_graph, outcome_graph = outcome_graph, effect_graph=effect_graph))
    return(list(avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m, p_test_results = p_test_results, first_stage_graph = propensity_score_graph, outcome_graph = outcome_graph, effect_graph=effect_graph))
  } 
  if (return_graph == TRUE & estimator == "ATT" & split_estimation == TRUE) {
    #return(list(first_stage_jump = first_stage_jump, avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m, first_stage_graph = propensity_score_graph, outcome_graph = outcome_graph, effect_graph=effect_graphh))
    return(list(avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m, p_test_results = p_test_results, first_stage_graph = propensity_score_graph, outcome_graph = outcome_graph, effect_graph=effect_graph))
  } 
  if (return_graph == TRUE & estimator == "ATT" & split_estimation == FALSE) {
    #return(list(first_stage_jump = first_stage_jump, avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m, first_stage_graph = propensity_score_graph, outcome_graph = outcome_graph, effect_graph=effect_graph, effect_graph_fixed_threshold = effect_graph_fixed_threshold))
    return(list(avg_effect = avg_effect, avg_effect_se = avg_effect_se, late = late, late_se = late_se, CATE_tibble = tibble(output), CATE_var = output_var, gam_output = m, p_test_results = p_test_results, first_stage_graph = propensity_score_graph, outcome_graph = outcome_graph, effect_graph=effect_graph, effect_graph_fixed_threshold = effect_graph_fixed_threshold))
  } 
} 



