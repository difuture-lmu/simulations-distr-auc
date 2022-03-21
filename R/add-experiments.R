fAlgo = function(job, data, instance, base_seed, l2sens, epsilon, delta, reps) {

  if ((delta == 0) || (epsilon == 0)) {
    noise = 0
  } else {
    cs = 2 * log(1.25 / delta)
    noise = sqrt(cs) * l2sens / epsilon
  }

  ll = list()
  for (i in seq_len(reps)) {
    alpha = 0.05
    auc_emp = 1

    df_roc  = generateROCData(i, base_seed)
    auc_emp = logitAUC(df_roc, unlogit = TRUE)
    nsim    = nrow(df_roc)
    npos    = sum(df_roc$truth)
    #method  = c("Bootstrap Empirical", "Bootstrap ROC-GLM", "DeLong Empirical", "DeLong ROC-GLM")
    method  = c("DeLong Empirical", "DeLong ROC-GLM")

    e = try({
      df_roc_glm = df_roc

      x = sort(df_roc_glm$score)

      # Generate noisy scores:
      df_roc_glm$score = rnorm(nsim, df_roc_glm$score, noise)

      # Truncate if values are bigger than 1 or smaller than 0:
      df_roc_glm$score = ifelse(df_roc_glm$score < 0, 0, df_roc_glm$score)
      df_roc_glm$score = ifelse(df_roc_glm$score > 1, 1, df_roc_glm$score)

      # Calculate ROC-GLM on noisy scores:
      auc_roc = rocGLMlogitAUC(df_roc_glm, ind = seq_len(nsim), unlogit = TRUE)
      auc_roc_params = attr(auc_roc, "params")

      # Transform AUC to log scale:
      lauc_emp = log(auc_emp / (1 - auc_emp))
      lauc_roc = log(auc_roc / (1 - auc_roc))

      # Calculate CI using the variance based on pooled data without noisy score values:
      var_auc = deLongVar(scores = df_roc$score, truth = df_roc$truth)
      ci_emp_log = pepeCI(lauc_emp, alpha, var_auc)
      ci_emp_auc = 1 / (1 + exp(-ci_emp_log))

      # Calculate CI using the variance based on the noisy score values:
      var_auc_app = deLongVar(scores = df_roc_glm$score, truth = df_roc_glm$truth)
      ci_app_log = pepeCI(lauc_emp, alpha, var_auc_app)
      ci_app_auc = 1 / (1 + exp(-ci_app_log))

      # Calculate the error of the approximated CI:
      delta_ci = sum(abs(ci_app_auc - ci_emp_auc))# / sum(diff(ci_emp_auc))

      # Gather results:
      df_aucs = data.frame(auc_emp = auc_emp, auc_roc = auc_roc, auc_roc_param1 = auc_roc_params[1],
        auc_roc_param2 = auc_roc_params[2], n = nsim, npos = npos, threshold = 0.5, noise = noise,
        epsilon = epsilon, delta = delta, l2sens = l2sens, delta_ci = delta_ci)

      df_cis = data.frame(
        method    = method,
        log_lower = c(ci_emp_log[1], ci_app_log[1]),
        log_upper = c(ci_emp_log[2], ci_app_log[2]),
        log_auc   = c(lauc_emp, lauc_roc),
        lower     = c(ci_emp_auc[1], ci_app_auc[1]),
        upper     = c(ci_emp_auc[2], ci_app_auc[2]),
        auc       = c(auc_emp, auc_roc))

      list(aucs = df_aucs, cis = df_cis)
    }, silent = TRUE)

    if (class(e) == "try-error") {
      df_aucs = data.frame(auc_emp = auc_emp, auc_roc = NA, auc_roc_param1 = NA, auc_roc_param2 = NA,
        n = nsim, npos = npos, threshold = 0.5, noise = noise, epsilon = epsilon, delta = delta,
        l2sens = l2sens, delta_ci = NA)

      df_cis = data.frame(
        method   = method,
        log_lower = rep(NA, 2),
        log_upper = rep(NA, 2),
        log_auc   = rep(NA, 2),
        lower     = rep(NA, 2),
        upper     = rep(NA, 2),
        auc       = rep(NA, 2)
      )

      out = list(aucs = df_aucs, cis = df_cis)
    } else {
      out = e
    }
    ll[[i]] = out
  }
  return(ll)
}

addProblem("dummy")
addAlgorithm(name = "auc-values", fun = fAlgo)
addExperiments(algo.design = list('auc-values' = rbind(

  # This is the base configuration without any noise. We use
  # this to get an idea of the accuracy of the ROC-GLM.
  data.frame(
    base_seed = BASE_SEED,
    l2sens    = 0,
    epsilon   = 0,
    delta     = 0,
    reps      = REPETITIONS
  ),
  # This is the grid of all tried out values:
  expand.grid(
    base_seed  = BASE_SEED,
    l2sens     = L2SENS,
    epsilon    = EPSILON,
    delta      = DELTA,
    reps       = REPETITIONS
  )
)))
