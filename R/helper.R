# ============================================================================ #
#                              DATA GENERATION
# ============================================================================ #

#' Simulate score and truth values
#'
#' The simulation randomly generates scores and true values of size n which is
#' randomly drawn between 100 and 2500. The generation is based on randomly
#' flipping values to get AUC values on the full range between 0.5 and 1 when
#' using the simulated data.
#'
#' @param i (`integer(1)`) The repetition (> 0) or basically an integer number to
#'   shift the seed.
#' @param seed (`integer(1)`) The base seed (> 0) used for the simulation. Note that
#'   the actual seed used for simulation is `seed + i`.
#' @return (`data.frame()`) Data frame containing columns `score` and `truth`. The
#'   score represents the predicted scores (e.g. from a statistical or prediction model).
generateROCData = function(i, seed) {
  checkmate::assertCount(i)
  checkmate::assertCount(seed)

  # Used seed:
  seed_k = seed + i

  # Simulate data:

  set.seed(seed_k)
  nsim = sample(x = seq(from = 100L, to = 2500L), size = 1L)
  npos = nsim

  set.seed(seed_k)
  scores = runif(n = nsim, min = 0, max = 1)
  truth = ifelse(scores > 0.5, 1, 0)

  # Sometimes, the simulation produces a too unbalanced data situation.
  # We catch this by checking if the number of positives divided by the
  # number of negatives is in [0.1, 0.9]. But, when a simulation falls
  # outside of this range, we have to increase the seed by one to not
  # get the same values again.
  seed_add_gamma = 0L
  while ((npos / nsim > 0.9) || (npos / nsim < 0.1)) {

    set.seed(seed_k + seed_add_gamma)
    shuffle = runif(n = 1L, min = 0, max = 1)
    nshuffle = trunc(nsim * shuffle)

    set.seed(seed_k + seed_add_gamma)
    idx_shuffle = sample(x = seq_len(nsim), size = nshuffle)

    set.seed(seed_k + seed_add_gamma)
    truth[idx_shuffle] = sample(x = c(0,1), size = nshuffle, replace = TRUE)
    npos = sum(truth)

    seed_add_gamma = seed_add_gamma + 1L
  }
  return(data.frame(score = scores, truth = truth))
}

#' Generate ROC data (tpr and fpr) from scores and labels
#'
#' @param score (`numeric()`) Scores of the respective label.
#' @param truth (`integer()`) True labels coded a 0 and 1.
#' @return (`data.frame()`) Data frame containing columns `tpr` and `fpr`.
rocData = function(score, truth) {

  checkmate::assertNumeric(x = score, any.missing = FALSE, len = length(truth))
  checkmate::assertIntegerish(x = truth, lower = 0, upper = 1, any.missing = FALSE)

  label = truth[order(score, decreasing = TRUE)]

  n_pos = sum(label == 1)
  n_neg = sum(label == 0)

  tpr = c(0, cumsum(label == 1) / n_pos)
  fpr = c(0, cumsum(label == 0) / n_neg)

  return(list(tpr = tpr, fpr = fpr))
}


# ============================================================================ #
#                                  AUC + CI
# ============================================================================ #

#' Calculate empirical AUC
#'
#' @param data (`data.frame()`) Data frame containing columns `score` and `truth`.
#' @param ind (`integer()`) Indices for subsetting the data.
#' @param unlogit (`logical(1)`) If `TRUE` the AUC is given as it is on the [0,1]
#'   scale. If `FALSE`, the AUC is transformed with log(AUC / (1 - AUC)).
#' @return (`numeric(1)`) The value of AUC either with or without logit transformation.
logitAUC = function(data, ind = NULL, unlogit = FALSE) {
  checkmate::assertDataFrame(x = data)
  nuisance = lapply(colnames(data), function(nm)
    checkmate::assertChoice(x = nm, choices = c("score", "truth")))
  checkmate::assertIntegerish(ind, lower = 1, null.ok = TRUE)

  if (is.null(ind[1])) ind = seq_len(nrow(data))
  scores = data$score[ind]
  truth = data$truth[ind]
  emp_auc = as.numeric(pROC::auc(truth, scores))
  #emp_auc = mlr::measureAUC(probabilities = scores, truth = truth, negative = 0, positive = 1)
  if (unlogit) return(emp_auc)
  return(toLogit(emp_auc))
}

#' Transform a logit value to original scale
#'
#' @param x (`numeric()`) The logit value x.
#' @return (`numeric()`) The transformed value (1 + exp(-x))^(-1).
logitToAUC = function(x) 1 / (1 + exp(-x))

#' Logit transformation
#'
#' @param x (`numeric()`) The value x.
#' @return (`numeric()`) The transformed value log(x / (1 - x)).
toLogit = function(x) log(x / (1 - x))

#' Calculate the variance of the logit AUC based on DeLong
#'
#' @param score (`numeric()`) Scores of the respective label.
#' @param truth (`integer()`) True labels coded a 0 and 1.
#' @return (`numeric(1)`) Variance of the AUC based on DeLong.
deLongVar = function(scores, truth) {
  checkmate::assertNumeric(x = scores, any.missing = FALSE, len = length(truth))
  checkmate::assertIntegerish(x = truth, lower = 0, upper = 1, any.missing = FALSE)

  # survivor functions for diseased and non diseased:
  s_d = function(x) 1 - ecdf(scores[truth == 1])(x)
  s_nond = function(x) 1 - ecdf(scores[truth == 0])(x)

  # Variance of empirical AUC after DeLong:
  var_auc = var(s_d(scores[truth == 0])) / sum(truth == 0) +
    var(s_nond(scores[truth == 1])) / sum(truth == 1)

  return(var_auc)
}

#' Calculate the confidence interval (CI) after Pepe
#'
#' @param logit_auc (`numeric(1)`) The AUC on the logit scale.
#' @param alpha (`numeric(1)`) The significance level.
#' @param var_auc (`numeric(1)`) The variance of the logit AUC.
#' @return (`numeric(2)`) Lower and upper CI of the logit AUC.
pepeCI = function(logit_auc, alpha, var_auc) {
  checkmate::assertNumeric(x = logit_auc, len = 1L)
  checkmate::assertNumeric(x = alpha, len = 1L, lower = 0, upper = 1)
  checkmate::assertNumeric(x = var_auc, len = 1L, lower = 0)
  quant = qnorm(1 - alpha / 2) * sqrt(var_auc) /
    (logitToAUC(logit_auc) * (1 - logitToAUC(logit_auc)))
  return(logit_auc + c(-1, 1) * quant)
}
                    
               

#' Helper function to calculate the discrepancy between ROC GLM and empirical ROC curve
discrepancy_function <- function(x, a, b, roc_curve_function){
  abs(pnorm(a + b*qnorm(x)) - roc_curve_function(x))
}
                    


#' Calculate discrepancy between ROC GLM and empirical ROC curve
#'
#' @param data (`data.frame()`) Data frame containing columns `score` and `truth`.
#' @param a (`numeric()`) first param from ROC GLM
#' @param b (`numeric()`) first param from ROC GLM
#' @param ind (`integer()`) Indices for subsetting the data. -- no tused
#' @return (`numeric(1)`) The value of AUC either with or without logit transformation.
calcDiscrepancy  = function(data, a, b, ind = NULL) {
  checkmate::assertNumeric(a, len = 1L)
  checkmate::assertNumeric(b, len = 1L)
  # no other checks here as everything is already tested by other functions that are called earlier
  
  if (is.null(ind[1])) ind = seq_len(nrow(data))
  
  # calc roc curve to get tpr and fpr
  curve_tmp = pROC::roc(response=data$truth[ind], predictor=data$score[ind])
  tpr_tmp <- rev(curve_tmp$sensitivities)
  fpr_tmp <- rev(1-curve_tmp$specificities)
  
  # cast tpr and fpr into a function
  roc_curve_function_tmp <- stepfun(fpr_tmp+.Machine$double.eps, c(0,tpr_tmp), f = 0, right=F)
  
  return(integrate(function(x) discrepancy_function(x, a = a, b = b, roc_curve_function = roc_curve_function_tmp), lower = 0, upper = 1, subdivisions = 500)$value)
}
       
                    
                    

# ============================================================================ #
#                              PROBIT REGRESSION
# ============================================================================ #

#' Calculate Fisher scoring parameter updates
#'
#' The updates here are based on a transformation to use a simpler
#' updating rule in the case of the Probit regression. The used
#' form is (XWX)^{-1} * lambda where lambda contains the Probit
#' regression specific information used for the update.
#'
#' @param beta (`numeric()`) Current parameter vector.
#' @param X (`matrix()`) Data matrix.
#' @param lambda (`numeric()`) Vector of values containing the Probit information.
#' @param w (`numeric()`) Weights.
#' @return (`numeric()`) New parameter value.
updateParam = function(beta, X, lambda, w = NULL) {
  checkmate::assertNumeric(x = beta, len = ncol(X))
  checkmate::assertMatrix(x = X, mode = "numeric")
  checkmate::assertNumeric(x = lambda, len = nrow(X))
  checkmate::assertNumeric(x = w, len = nrow(X), null.ok = TRUE)

  W = diag(as.vector(lambda * (X %*% beta + lambda)))
  if (! is.null(w[1])) {
    X = diag(sqrt(w)) %*% X
    lambda = diag(sqrt(w)) %*% lambda
  }
  return(beta + solve(t(X) %*% W %*% X) %*% t(X) %*% lambda)
}

#' Conduct Probit regression
#'
#' Note that the values for this custom Probit regression are the
#' same as when calling the `glm` function.
#'
#' @param y (`numeric()`) Response vector with 0-1 entries.
#' @param X (`matrix()`) Data matrix.
#' @param w (`numeric()`) Weights.
#' @param beta_start (`numeric()`) Initial parameter value (default = 0).
#' @param stop_tol (`numeric(1)`) Value used to stop the Fisher scoring (default = 1e-8).
#' @param iter_max (`integer(1)`) Maximal number of Fisher scoring iterations.
#' @param trace (`logical(1)`) Flag to indicate whether to print fitting information of not.
#' @return (`numeric()`) Parameter estimates.
probitRegr = function(y, X, w = NULL, beta_start = 0, stop_tol = 1e-8, iter_max = 25L, trace = FALSE) {
  checkmate::assertIntegerish(x = y, lower = 0, upper = 1, any.missing = FALSE, len = nrow(X))
  checkmate::assertMatrix(x = X, mode = "numeric")
  checkmate::assertNumeric(x = w, len = length(y), null.ok = TRUE)
  
  if (length(beta_start) == 1) beta_start = rep(beta_start, ncol(X))
  
  checkmate::assertNumeric(x = beta_start, len = ncol(X))
  checkmate::assertNumeric(x = stop_tol, len = 1L, lower = 0)
  checkmate::assertCount(x = iter_max)
  checkmate::assertLogical(x = trace, len = 1L)

  if (length(beta_start) == 1) beta_start = rep(beta_start, ncol(X))
  if (is.vector(beta_start)) beta_start = cbind(beta_start)

  beta = beta_start
  iter = 0L

  deviance = numeric(iter_max)
  deviance[1] = probitDeviance(y, X, beta)

  if (trace) cat("\n")

  while (iter <= iter_max) {

    beta_start = beta
    beta = updateParam(beta_start, X = X, lambda = calculateLambda(y = y, X = X, beta = beta_start), w = w)

    iter = iter + 1L
    deviance[iter + 1L] = probitDeviance(y, X, beta)

    if (trace) cat("Deviance of iter", iter, "=", round(deviance[iter + 1L], digits = 4L), "\n")

    if (probitDevianceStop(y = y, X = X, beta = beta, beta_old = beta_start) < stop_tol) {
      if (trace) { cat("\n"); break; }
    }
  }
  out = list(iter = iter, parameter = beta, deviance = deviance[seq_len(iter + 1L)])
  return(out)
}

#' Predict the estimated Probit regression on a regular grid between 0 and 1.
#'
#' @param mod (Object returned from `probitRegr`) Estimated Probit regression parameter.
#' @return (`numeric()`) Predicted values on the regular grid using the binormal form.
predictProbit = function(mod) {
  x = seq(0, 1, 0.01)
  y = pnorm(mod$parameter[1] + mod$parameter[2] * qnorm(x))

  return(list(x = x, y = y))
}

#' Calculate the stopping criteria of the Probit regression based on the deviance
#'
#' @param y (`numeric()`) Response vector with 0-1 entries.
#' @param X (`matrix()`) Data matrix.
#' @param beta (`numeric()`) New parameter value.
#' @param beta_old (`numeric()`) Old parameter value.
#' @return (`numeric(1)`) Deviance stop criteria.
probitDevianceStop = function(y, X, beta, beta_old) {
  checkmate::assertIntegerish(x = y, lower = 0, upper = 1, any.missing = FALSE, len = nrow(X))
  checkmate::assertMatrix(x = X, mode = "numeric")
  checkmate::assertNumeric(x = beta, len = ncol(X))
  checkmate::assertNumeric(x = beta_old, len = ncol(X))

  dev = -2 * log(probitLikelihood(y, X, beta))
  dev_old = -2 * log(probitLikelihood(y, X, beta_old))

  # from ?glm.control:
  out = abs(dev - dev_old) / (abs(dev) + 0.1)
  return(out)
}

#' Calculate likelihood for the Probit regression
#'
#' @param y (`numeric()`) Response vector with 0-1 entries.
#' @param X (`matrix()`) Data matrix.
#' @param beta (`numeric()`) Parameter value.
#' @param w (`numeric()`) Weights.
#' @return (`numeric(1)`) Log-likelihood.
probitLikelihood = function(y, X, beta, w = NULL) {
  checkmate::assertIntegerish(x = y, lower = 0, upper = 1, any.missing = FALSE, len = nrow(X))
  checkmate::assertMatrix(x = X, mode = "numeric")
  checkmate::assertNumeric(x = beta, len = ncol(X))
  checkmate::assertNumeric(x = w, len = length(y), null.ok = TRUE)

  eta = X %*% beta

  if (is.null(w)) {
    w = rep(1, times = nrow(X))
  }
  lh = pnorm(eta)^y * (1 - pnorm(eta))^(1 - y)
  return(prod(lh^w))
}

#' Calculate the deviance for the Probit regression
#'
#' @param y (`numeric()`) Response vector with 0-1 entries.
#' @param X (`matrix()`) Data matrix.
#' @param beta (`numeric()`) Parameter value.
#' @return (`numeric(1)`) Deviance.
probitDeviance = function(y, X, beta) {
  checkmate::assertIntegerish(x = y, lower = 0, upper = 1, any.missing = FALSE, len = nrow(X))
  checkmate::assertMatrix(x = X, mode = "numeric")
  checkmate::assertNumeric(x = beta, len = ncol(X))

  return(-2 * log(probitLikelihood(y, X, beta)))
}

#' Calculate the lambda for the Probit regression
#'
#' @param y (`numeric()`) Response vector with 0-1 entries.
#' @param X (`matrix()`) Data matrix.
#' @param beta (`numeric()`) Parameter value.
#' @return (`numeric()`) Vector of lambda values.
calculateLambda = function(y, X, beta) {
  checkmate::assertIntegerish(x = y, lower = 0, upper = 1, any.missing = FALSE, len = nrow(X))
  checkmate::assertMatrix(x = X, mode = "numeric")
  checkmate::assertNumeric(x = beta, len = ncol(X))

  eta = X %*% beta
  q = 2 * y - 1
  qeta = q * eta

  return((dnorm(qeta) * q) / (pnorm(qeta)))
}

# ============================================================================ #
#                                  ROC GLM
# ============================================================================ #

#' Calculate the U matrix for the ROC-GLM
#'
#' The U matrix is the response used for the Probit regression to
#' fit the ROG-GLM.
#'
#' @param tset (`numeric()`) Threshold values.
#' @param placement_values (`numeric()`) Placement values.
#' @return (`matrix()`) 0-1-matrix containing the response for the Probit regression.
calcU = function(tset, placement_values) {
  tset_sorted = sort(tset)
  out = vapply(X = tset, FUN.VALUE = integer(length(placement_values)),
    FUN = function(t) { ifelse(placement_values <= t, 1L, 0L) })
    return(out)
}

#' Create the data matrix for the Probit regression
#'
#' The ROC-GLM is basically a Probit regression on specific data. This
#' function creates this data frame.
#'
#' @param U (`matrix()`) 0-1-matrix of response values.
#' @param tset (`numeric()`) Threshold values.
#' @return (`data.frame()`) Data frame with columns `y` (response) `x` (covariate
#'   based on the binormal form), and `w` (weights).
dataRocGLM = function(U, tset) {
  roc_glm_data = data.frame(
    y = rep(c(0, 1), times = length(tset)),
    x = rep(qnorm(tset), each = 2L),
    w = as.vector(apply(U, 2, function(x) c(sum(x == 0), sum(x == 1)))))
  return(roc_glm_data)
}

#' Integrate over the binormal form to get AUC estimate
#'
#' @param mod (Object returned from `probitRegr`) Estimated Probit regression parameter.
#' @return (`numeric(1)`) Estimated AUC based on the ROC-GLM.
integrateBinormal = function(params) {
  temp = function(x) pnorm(params[1] + params[2] * qnorm(x))
  int = integrate(f = temp, lower = 0, upper = 1)
  return(int$value)
}

#' Calculate AUC based on the ROC-GLM
#'
#' @param data (`data.frame()`) Data frame containing columns `score` and `truth`.
#' @param ind (`integer()`) Indices for subsetting the data.
#' @param unlogit (`logical(1)`) If `TRUE` the AUC is given as it is on the [0,1]
#'   scale. If `FALSE`, the AUC is transformed with log(AUC / (1 - AUC)).
#' @return (`numeric(1)`) The value of AUC (based on the ROC-GLM) either with or
#'   without logit transformation.
rocGLMlogitAUC = function(data, ind = NULL, unlogit = FALSE) {
  if (is.null(ind[1])) ind = seq_len(nrow(data))

  scores = data$score[ind]
  truth = data$truth[ind]

  Fn_global = ecdf(scores[truth == 0])
  Sn_global = function(x) 1 - Fn_global(x)

  thresh_set = seq(0, 1, length.out = 30L)
  pv_global = Sn_global(scores[truth == 1])
  U_global = calcU(thresh_set, pv_global)

  roc_data_global = dataRocGLM(U = U_global, tset = thresh_set)
  roc_data_global = roc_data_global[is.finite(roc_data_global$x),]

  y_global = roc_data_global$y
  X_global = model.matrix(y ~ x, data = roc_data_global)
  w_global = roc_data_global$w

  my_roc_glm_global = tryCatch(
    expr = {probitRegr(y = y_global, X = X_global, w = w_global)},
    error = function(e) return("fail")
  )
  if (! is.character(my_roc_glm_global)) {
    auc = integrateBinormal(my_roc_glm_global$parameter)
    attr(auc, "params") = my_roc_glm_global$parameter
    if (unlogit) return(auc)
    return(log(auc / (1 - auc)))
  } else {
    return(NA)
  }
}


# ============================================================================ #
#                                  Gaussian Noise
# ============================================================================ #

#' Calculates the error function used for the analytic Gaussian mechanism
#'
#' @param x (`numeric()`) A (vector of) real number(s)
#' @return (`numeric()`) The evaluated error function
#' @author Raphael Rehms
erf = function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
}

#' Calculates the analytic Gaussian mechanism for given privacy parameters (see https://arxiv.org/abs/1805.06530)
#'
#' @param epsilon (`numeric()`) Epsilon > 0
#' @param delta (`numeric()`) Delta that is between (0,1)
#' @param sens (`numeric()`) Sensitivity of the algorithm
#' @param tol (`numeric()`) Tolerance for binary search
#' @return (`numeric()`) Sigma that can be used to generate calibrated Gaussian noise
#' @author Raphael Rehms
analyticGaussianMechanism = function(epsilon, delta, sens, tol = 1e-12) {
  phi = function(t) {
    0.5 * (1.0 + erf(t / sqrt(2.0)))
  }
  
  caseA = function(epsilon, s) {
    phi(sqrt(epsilon * s)) - exp(epsilon) * phi(-sqrt(epsilon * (s + 2.0)))
  }
  
  caseB = function(epsilon, s) {
    phi(-sqrt(epsilon * s)) - exp(epsilon) * phi(-sqrt(epsilon * (s + 2.0)))
  }
  
  doubling_trick = function(predicate_stop, s_inf, s_sup) {
    while (!predicate_stop(s_sup)) {
      s_inf = s_sup
      s_sup = 2.0 * s_inf
    }
    return(c(s_inf, s_sup))
  }
  
  binary_search = function(predicate_stop, predicate_left, s_inf, s_sup) {
    s_mid = s_inf + (s_sup - s_inf) / 2.0
    while (!predicate_stop(s_mid)) {
      if (predicate_left(s_mid)) {
        s_sup = s_mid
      } else {
        s_inf = s_mid
      }
      s_mid = s_inf + (s_sup - s_inf) / 2.0
    }
    return(s_mid)
  }
  
  delta_thr = caseA(epsilon, 0.0)
  
  if (delta == delta_thr) {
    alpha = 1.0
  } else {
    if (delta > delta_thr) {
      predicate_stop_DT = function(s) caseA(epsilon, s) >= delta
      function_s_to_delta = function(s) caseA(epsilon, s)
      predicate_left_BS = function(s) function_s_to_delta(s) > delta
      function_s_to_alpha = function(s) sqrt(1.0 + s / 2.0) - sqrt(s / 2.0)
    } else {
      predicate_stop_DT = function(s) caseB(epsilon, s) <= delta
      function_s_to_delta = function(s) caseB(epsilon, s)
      predicate_left_BS = function(s) function_s_to_delta(s) < delta
      function_s_to_alpha = function(s) sqrt(1.0 + s / 2.0) + sqrt(s / 2.0)
    }
    
    predicate_stop_BS = function(s) abs(function_s_to_delta(s) - delta) <= tol
    
    s_ = doubling_trick(predicate_stop_DT, 0.0, 1.0)
    s_inf = s_[1]; s_sup = s_[2]
    
    s_final = binary_search(predicate_stop_BS, predicate_left_BS, s_inf, s_sup)
    alpha = function_s_to_alpha(s_final)
  }
  
  sigma = alpha * sens / sqrt(2.0 * epsilon)
  
  return(sigma)
}
