#' Estimates the conditional effect for binary outcome MRT using
#' truncated weight
#'
#' This estimator returns the estimates for the semiparametric locally efficient estimator
#' using truncated weight and provides the estimated variance and standard error for
#' the estimators, with small sample correction for an MRT with binary outcome using the
#' "Hat" matrix in the variance estimate and t-distribution or F-distribution critical
#' value with corrected degrees of freedom.
#'
#' This estimator should theoretically output approximately the same
#' results as the efficient_ee estimator and the efficient_ee_twostep
#' estimator. If this estimator does not return appropriate estimates,
#' please try out the other two estimators.
#'
#' @param dta a data set in long format
#' @param id_varname the variable name that specifies the subject IDs
#' @param decision_time_varname the variable name that specifies the decision points
#' @param treatment_varname the variable name that specifies the assigned treatments for subjects
#' @param outcome_varname the variable name that specifies the outcomes for subjects
#' @param control_varname a vector of variable names used to reduce noise, this could be NULL
#' @param moderator_varname a vector of variable names of the effect modifiers, this could be NULL
#' @param rand_prob_varname the variable name that specifies the treatment randomizing probability in the data set
#' @param avail_varname the variable name that specifies the availability of the subjects,
#'                      default to be always available at any decision points using NULL
#' @param estimator_initial_value a numeric vector of the initial value for the estimator,
#'                                its length should be the sum of the length of control and moderator variables plus 2
#'                                default to be all 0's using NULL
#' @param weight_threshold a numerical value that specifies the truncated weight used in the estimating equation,
#'                         default to be 0.95 using NULL
#'
#' @return the estimated beta with its intercept and alpha with its intercept,
#'         the standard error of the estimated beta with its intercept and alpha with its intercept,
#'         the adjusted standard error of the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the value of the estimating function at the estimated beta and alpha
#' @export
#'
#' @examples efficient_ee_modified_weight(dta = dgm_sample,
#'                          id_varname = "userid",
#'                          decision_time_varname = "day",
#'                          treatment_varname = "A",
#'                          outcome_varname = "Y",
#'                          control_varname = c("time_var1", "time_var2"),
#'                          moderator_varname = "time_var1",
#'                          rand_prob_varname = "prob_A",
#'                          weight_threshold = 0.90)
#'
efficient_ee_modified_weight <- function(
  dta,
  id_varname,
  decision_time_varname,
  treatment_varname,
  outcome_varname,
  control_varname,
  moderator_varname,
  rand_prob_varname,
  avail_varname = NULL,
  estimator_initial_value = NULL,
  weight_threshold = 0.95
)
{
  message("\nTheoretically, efficient_ee(), efficient_ee_twostep(), and efficient_ee_modified_weight() should return approximately the same estimates.")
  message("\nIf this estimator does not return appropriate outputs, please try out the other two estimators.")
  # weight in EE is modified to make it more stable

  lambda1 <- weight_threshold
  lambda2 <- weight_threshold

  ### 1. preparation ###

  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)

  A <- dta[, treatment_varname]
  p_t <- dta[, rand_prob_varname]
  cA <- A - p_t # centered A
  Y <- dta[, outcome_varname]
  Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added

  if (is.null(avail_varname)) {
    avail <- rep(1, total_person_decisionpoint)
  } else {
    avail <- dta[, avail_varname]
  }

  p <- length(moderator_varname) + 1 # dimension of beta
  q <- length(control_varname) + 1 # dimension of alpha

  Xnames <- c("Intercept", moderator_varname)
  Znames <- c("Intercept", control_varname)

  ### 2. estimation ###

  estimating_equation <- function(theta) {
    alpha <- as.matrix(theta[1:q])
    beta <- as.matrix(theta[(q+1):(q+p)])

    exp_Zdm_alpha <- exp(Zdm %*% alpha)
    exp_Xdm_beta <- exp(Xdm %*% beta)
    exp_AXdm_beta <- exp(A * (Xdm %*% beta))
    exp_negAXdm_beta <- exp_AXdm_beta^(-1)
    exp_negXdm_beta <- exp_Xdm_beta^(-1)

    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta

    # modifying weights

    xi1_vector <- ifelse(exp_Zdm_alpha <= lambda1, 1 - exp_Zdm_alpha, 1 - lambda1)
    xi2_vector <- ifelse(exp_Zdm_alpha * exp_Xdm_beta <= lambda2, 1 - exp_Zdm_alpha * exp_Xdm_beta, 1 - lambda2)

    weight <- exp_negAXdm_beta / ( xi1_vector * p_t + exp_negXdm_beta * xi2_vector * (1 - p_t) )

    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:q) {
      ef[i] <- sum( weight * residual * avail * Zdm[, i])
    }
    for (i in 1:p) {
      ef[q + i] <- sum( weight * residual * avail * cA * Xdm[, i])
    }

    ef <- ef / sample_size
    return(ef)
  }

  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p + q)
  }

  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee():")
      message("\nThe program cannot find a numerical solution to the estimating eqaution.")
      message(cond)
      return(list(root = rep(NaN, p + q), msg = cond,
                  f.root = rep(NaN, p + q)))
    })

  estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
  alpha_hat <- as.vector(estimator$alpha)
  beta_hat <- as.vector(estimator$beta)

  ### 3. asymptotic variance ###

  ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###

  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
  # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
  # See note 2018.08.06 about small sample correction

  r_term_collected <- rep(NA, total_person_decisionpoint)
  D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
  partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)

  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }
    if (q == 1) {
      Zalpha <- Zdm[it, ] * alpha_hat
    } else {
      Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
    }

    eta1 <- as.numeric(exp(Zalpha) <= lambda1)
    xi1 <- ifelse(eta1, 1 - exp(Zalpha), 1 - lambda1)
    eta2 <- as.numeric(exp(Zalpha + Xbeta) <= lambda2)
    xi2 <- ifelse(eta2, 1 - exp(Zalpha + Xbeta), 1 - lambda2)

    denom <- xi1 * p_t[it] + exp(-Xbeta) * xi2 * (1 - p_t[it])

    W1 <- exp(- A[it] * Xbeta) / (denom^2)
    W2 <- p_t[it] * eta1 + (1 - p_t[it]) * eta2
    W3 <- (1 - p_t[it]) * (exp(-Xbeta) * xi2 + exp(Zalpha) * eta2) - A[it] * denom

    # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
    partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
    partialD_partialtheta[1:q, 1:q] <- W1 * W2 * (Zdm[it, ] %o% Zdm[it, ])
    partialD_partialtheta[1:q, (q+1):(q+p)] <- W1 * W3 * (Zdm[it, ] %o% Xdm[it, ])
    partialD_partialtheta[(q+1):(q+p), 1:q] <- W1 * W2 * cA[it] * (Xdm[it, ] %o% Zdm[it, ])
    partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- W1 * W3 * cA[it] * (Xdm[it, ] %o% Xdm[it, ])

    # r_term = r^(t) (scalar)
    r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
    r_term_collected[it] <- r_term

    # D_term = D^{(t),T} (vector of length (p+q))
    D_term <- exp(- A[it] * Xbeta) / denom * c(Zdm[it, ], cA[it] * Xdm[it, ])
    D_term_collected[, it] <- D_term

    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T} (vector of length (p+q))
    partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
    partialr_partialtheta_collected[it, ] <- partialr_partialtheta

    Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)

  ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###

  Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
  # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
  # See note 2018.08.06 about small sample correction

  person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)

  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
    r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)

    Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
  }
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size

  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
  alpha_se <- sqrt(diag(varcov)[1:q])
  beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])


  ### 4. small sample correction ###

  Sigman_tilde <- 0
  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
    r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
    partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
    H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
    Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)

    Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
  }
  Sigman_tilde <- Sigman_tilde / sample_size

  varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
  alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
  beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])

  ### 6. return the result with variable names ###

  names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
  names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames

  return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
              beta_se = beta_se, alpha_se = alpha_se,
              beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
              # test_result_t = test_result_t,
              # test_result_f = test_result_f,
              varcov = varcov,
              varcov_adjusted = varcov_adjusted,
              dims = list(p = p, q = q),
              f.root = solution$f.root))
}
