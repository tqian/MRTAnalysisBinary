#' Estimates the conditional effect for binary outcome MRT using
#' a two-step estimator
#'
#' This estimator returns the estimates for the semiparametric locally efficient estimator
#' using a two-step estimator and provides the estimated variance and standard error for
#' the estimators, with small sample correction for an MRT with binary outcome using the
#' "Hat" matrix in the variance estimate and t-distribution or F-distribution critical
#' value with corrected degrees of freedom.
#'
#' This estimator should theoretically output approximately the same
#' results as the efficient_ee estimator and the efficient_ee_modified_weight
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
#'
#' @return the estimated beta with its intercept and alpha with its intercept,
#'         the standard error of the estimated beta with its intercept and alpha with its intercept,
#'         the adjusted standard error of the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the 95 percent confidence interval for beta_hat, and the adjusted 95 percent confidence interval for beta_hat,
#'         the dimension of the moderated variables, and the dimension of the control variables,
#'         the value of the estimating function at the estimated beta and alpha
#' @export
#' @import rootSolve
#' @examples efficient_ee_twostep(dta = dgm_sample,
#'                          id_varname = "userid",
#'                          decision_time_varname = "day",
#'                          treatment_varname = "A",
#'                          outcome_varname = "Y",
#'                          control_varname = c("time_var1", "time_var2"),
#'                          moderator_varname = "time_var1",
#'                          rand_prob_varname = "prob_A")
#'
efficient_ee_twostep <- function(
  dta,
  id_varname,
  decision_time_varname,
  treatment_varname,
  outcome_varname,
  control_varname,
  moderator_varname,
  rand_prob_varname,
  avail_varname = NULL
)
{
  message("\nTheoretically, efficient_ee(), efficient_ee_twostep(), and efficient_ee_modified_weight() should return approximately the same estimates.")
  message("\nIf this estimator does not return appropriate outputs, please try out the other two estimators.")
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

  ### 1.5. get initial values of beta and alpha to be used in the efficient estimating equation ###

  # This is by assuming E(Y_{t+1} | H_t, A_t) = exp(Zdm %*% alpha + A_t * Xdm %*% beta)
  # So the parameters can be found using glm(family = binomial(link = "log"))

  # browser()

  # glm_formula <- as.formula(paste0(outcome_varname, "~", paste(control_varname, collapse = "+"), "+", treatment_varname, "*(",
  #       paste(moderator_varname, collapse = "+"), ")"))
  # glm_fit <- glm(glm_formula, family = binomial(link = "log"), data = dta)

  # glm_fit <- glm.fit(x = cbind(Zdm, A * Xdm), y = Y, family = binomial(link = "log"))
  # geese.fit(x = cbind(Zdm, A * Xdm), y = Y, id = dta[, id_varname], family = binomial(link = "log"))

  # alphabeta_init <- glm_fit$coefficients # initial values of alpha and beta
  # alpha_init <- alphabeta_init[1:q]
  # beta_init <- alphabeta_init[(q+1):(q+p)]


  ## First, use an estimating equation without weights (from the derivative) to obtain the initial parameters

  ee_init <- function(theta) {
    alpha <- as.matrix(theta[1:q])
    beta <- as.matrix(theta[(q+1):(q+p)])

    exp_Zdm_alpha <- exp(Zdm %*% alpha)
    exp_AXdm_beta <- exp(A * Xdm %*% beta)

    # first time starting value through the EE without weight (from the derivative)
    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta

    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:q) {
      ef[i] <- sum( residual * avail * Zdm[, i])
    }
    for (i in 1:p) {
      ef[q + i] <- sum( residual * avail * A * Xdm[, i])
    }

    ef <- ef / sample_size
    return(ef)
  }

  solution_init <- tryCatch(
    {
      multiroot(ee_init, rep(0, p + q)) # initial value is all zero's
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee_twostep():")
      message("\nThe program cannot find a numerical solution to the first step estimating equation.")
      message(cond)
      return(list(root = rep(NaN, p+q), msg = cond,
                  f.root = rep(NaN, p+q)))
    })

  estimator_init <- get_alpha_beta_from_multiroot_result(solution_init, p, q)
  alpha_init <- as.vector(estimator_init$alpha)
  beta_init <- as.vector(estimator_init$beta)


  ## Second, use the score equation from MLE

  ee_mle <- function(theta) {
    alpha <- as.matrix(theta[1:q])
    beta <- as.matrix(theta[(q+1):(q+p)])

    exp_Zdm_alpha <- exp(Zdm %*% alpha)
    exp_AXdm_beta <- exp(A * Xdm %*% beta)

    residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
    weight <- 1 / (1 - exp_Zdm_alpha * exp_AXdm_beta)

    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:q) {
      ef[i] <- sum( weight * residual * avail * Zdm[, i])
    }
    for (i in 1:p) {
      ef[q + i] <- sum( weight * residual * avail * A * Xdm[, i])
    }

    ef <- ef / sample_size
    return(ef)
  }

  solution_mle <- tryCatch(
    {
      multiroot(ee_mle, c(alpha_init, beta_init)) # initial value is from ee_init
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee_twostep():")
      message("\nThe program cannot find a numerical solution to the MLE score equation.")
      message(cond)
      return(list(root = rep(NaN, p+q), msg = cond,
                  f.root = rep(NaN, p+q),
                  iter = NaN, estim.precis = NaN))
    })

  # summary(1 / (1 - exp(Zdm %*% alpha_init) * exp(A * Xdm %*% beta_init)))

  if (!is.nan(solution_mle$estim.precis) & solution_mle$estim.precis < 1e-2) {
    estimator_mle <- get_alpha_beta_from_multiroot_result(solution_mle, p, q)
    alpha_mle <- as.vector(estimator_mle$alpha)
    beta_mle <- as.vector(estimator_mle$beta)
  } else {
    alpha_mle <- alpha_init
    beta_mle <- beta_init
  }

  ### 2. estimation ###

  exp_Zdm_alpha <- exp(Zdm %*% alpha_mle)
  exp_Xdm_beta <- exp(Xdm %*% beta_mle)
  exp_negXdm_beta <- exp_Xdm_beta^(-1)
  weight <- ( (1 - exp_Zdm_alpha) * p_t + (exp_negXdm_beta - exp_Zdm_alpha) * (1 - p_t) )^(-1)

  estimating_equation <- function(beta) {

    # only the blipping-down part uses the new beta; all other parts uses old beta
    residual <- exp( - A * (Xdm %*% beta) ) * Y - exp_Zdm_alpha

    ef <- rep(NA, length(beta)) # value of estimating function
    for (i in 1:p) {
      ef[i] <- sum( weight * residual * avail * cA * Xdm[, i])
    }

    ef <- ef / sample_size
    return(ef)
  }

  solution <- tryCatch(
    {
      multiroot(estimating_equation, beta_mle) # initial value here is from ee_mle
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee_twostep():")
      message("\nThe program cannot find a numerical solution to the second step estimating equation.")
      message(cond)
      return(list(root = rep(NaN, p), msg = cond,
                  f.root = rep(NaN, p),
                  iter = NaN, estim.precis = NaN))
    })

  if (!is.nan(solution$estim.precis) & solution$estim.precis < 1e-2) {
    alpha_hat <- alpha_mle
    beta_hat <- solution$root
  } else {
    alpha_hat <- alpha_mle
    beta_hat <- beta_mle
  }


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

    #  old version, seems incorrect
    # W1 <- ( (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it]) ) ^ (-2)
    # W2 <- exp(Zalpha - A[it] * Xbeta)
    # W3 <- - ( exp(A[it] * Xbeta) * p_t[it] * A[it] - exp(Zalpha + A[it] * Xbeta) * A[it]
    #           + exp((A[it] - 1) * Xbeta) * (1 - p_t[it]) * (A[it] - 1) )

    denom <- (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it])
    W1 <- exp(- A[it] * Xbeta) / (denom^2)
    W2 <- exp(Zalpha)
    W3 <- exp(-Xbeta) * (1 - p_t[it]) - A[it] * denom


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
    D_term <- exp(- A[it] * Xbeta) / ( (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it]) ) *
      c(Zdm[it, ], cA[it] * Xdm[it, ])
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

  ### 7. calculate confidence interval

  conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
  c <- qt(1 - 0.05/2, df = sample_size - p - q)
  conf_int_adjusted <- cbind(beta_hat - c * beta_se_adjusted,
                             beta_hat + c * beta_se_adjusted)
  colnames(conf_int) <- colnames(conf_int_adjusted) <- c("2.5 %", "97.5 %")

  return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
              beta_se = beta_se, alpha_se = alpha_se,
              beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
              # test_result_t = test_result_t,
              # test_result_f = test_result_f,
              varcov = varcov,
              varcov_adjusted = varcov_adjusted,
              conf_int = conf_int, conf_int_adjusted = conf_int_adjusted,
              dims = list(p = p, q = q),
              f.root = solution$f.root))
}
