#' Estimates the marginal excursion effect for binary outcome MRT
#' where the treatment indicator is always centered
#'
#' This estimator assumes that the randomization probability does not
#' depend on H_t, and that the treatment indicator A is always centered.
#' This estimator returns the estimates for the marginal excursion effect estimator
#' with the above two assumptions and provides the estimated variance and
#' standard error for the estimators, with small sample correction for
#' an MRT with binary outcome using the "Hat" matrix in the variance estimate
#' and t-distribution or F-distribution critical value with corrected degrees of freedom.
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
#'
#' @return Returns the estimated beta with its intercept and alpha with its intercept,
#'         the standard error of the estimated beta with its intercept and alpha with its intercept,
#'         the adjusted standard error of the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the 95 percent confidence interval for beta_hat, and the adjusted 95 percent confidence interval for beta_hat,
#'         the dimension of the moderated variables, and the dimension of the control variables,
#'         the value of the estimating function at the estimated beta and alpha
#' @import rootSolve
#' @export
#'
#' @examples estimator_EMEE_alwaysCenterA(dta = dgm_sample,
#'                          id_varname = "userid",
#'                          decision_time_varname = "day",
#'                          treatment_varname = "A",
#'                          outcome_varname = "Y",
#'                          control_varname = c("time_var1", "time_var2"),
#'                          moderator_varname = "time_var1",
#'                          rand_prob_varname = "prob_A",
#'                          avail_varname = "avail")
#'
estimator_EMEE_alwaysCenterA <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname, # intercept will be automatically added
    moderator_varname, # intercept will be automatically added
    rand_prob_varname,
    avail_varname = NULL,
    estimator_initial_value = NULL
)
{

    ### 1. preparation ###

    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)

    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }

    A <- dta[, treatment_varname]
    # checking for NA in treatment indicator
    if (any(is.na(A[avail == 1]))) {
        stop("Treatment indicator is NA where availability = 1.")
    }
    A[avail == 0] <- 0

    p_t <- dta[, rand_prob_varname]
    cA <- A - p_t # centered A
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added

    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha

    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)

    ### 2. estimation ###

    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])

        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_negcAXdm_beta <- exp(- cA * (Xdm %*% beta))

        residual <- exp_negcAXdm_beta * Y - exp_Zdm_alpha

        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( residual * avail * cA * Xdm[, i])
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
            message("\nCatched error in multiroot inside estimator_EMEE_alwaysCenterA():")
            message("\nThe program cannot find a numerical solution to the estimating eqaution.")
            message(cond)
            return(list(root = rep(NaN, p + q), msg = cond,
                        f.root = rep(NaN, p + q)))
        })

    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)

    ### 3. asymptotic variance and small sample correction ###

    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###

    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction

    for (it in 1:total_person_decisionpoint) {

        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.
        if (p == 1) {
            Xbeta_it <- Xdm[it, ] * beta_hat # a scalar
        } else {
            Xbeta_it <- as.numeric(Xdm[it, ] %*% beta_hat) # a scalar
        }
        if (q == 1) {
            Zalpha_it <- Zdm[it, ] * alpha_hat # a scalar
        } else {
            Zalpha_it <- as.numeric(Zdm[it, ] %*% alpha_hat) # a scalar
        }

        Zdm_it <- as.numeric(Zdm[it, ])
        Xdm_it <- as.numeric(Xdm[it, ])

        stopifnot(class(Zalpha_it) == "numeric")
        stopifnot(class(Zdm_it) == "numeric")
        stopifnot(length(Zdm_it) == q)

        stopifnot(class(Xbeta_it) == "numeric")
        stopifnot(class(Xdm_it) == "numeric")
        stopifnot(length(Xdm_it) == p)

        Mn_summand[it, 1:q, 1:q] <- - avail[it] * exp(Zalpha_it) * (Zdm_it %o% Zdm_it)
        Mn_summand[it, 1:q, (q+1):(q+p)] <- - avail[it] * exp(- cA[it] * Xbeta_it) * Y[it] * cA[it] * (Zdm_it %o% Xdm_it)
        Mn_summand[it, (q+1):(q+p), 1:q] <- - avail[it] * exp(Zalpha_it) * cA[it] * (Xdm_it %o% Zdm_it)
        Mn_summand[it, (q+1):(q+p), (q+1):(q+p)] <- - avail[it] * exp(- cA[it] * Xbeta_it) * Y[it] * cA[it]^2 * (Xdm_it %o% Xdm_it)
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)

    ### 3.2 Compute \Sigma_n matrix and \tilde{\Sigma}_n matrix ###

    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    Sigman_tilde_summand <- array(NA, dim = c(sample_size, p+q, p+q))

    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)

    for (i in 1:sample_size) {

        T_i <- person_first_index[i+1] - person_first_index[i] # number of time points for individual i
        index_i <- person_first_index[i] : (person_first_index[i+1] - 1) # row number of individual i's data in dta

        Zdm_i <- Zdm[index_i, ] # each row in Zdm_i corresponds to a time point
        Xdm_i <- Xdm[index_i, ] # each row in Xdm_i corresponds to a time point
        cA_i <- cA[index_i]     # each entry in cA_i corresponds to a time point
        Y_i <- Y[index_i]

        if (p == 1) {
            Xdm_i <- matrix(Xdm_i, ncol = 1)
        }
        if (q == 1) {
            Zdm_i <- matrix(Zdm_i, ncol = 1)
        }

        stopifnot(class(Zdm_i)[1] == "matrix")
        stopifnot(class(Xdm_i)[1] == "matrix")
        stopifnot(nrow(Zdm_i) == T_i & ncol(Zdm_i) == q)
        stopifnot(nrow(Xdm_i) == T_i & ncol(Xdm_i) == p)

        D_i <- cbind(Zdm_i, cA_i * Xdm_i)
        r_i <- matrix(exp(- cA_i * Xdm_i %*% beta_hat) * Y_i - exp(Zdm_i %*% alpha_hat), nrow = T_i, ncol = 1)
        I_i <- diag(avail[index_i])

        stopifnot(nrow(D_i) == T_i & ncol(D_i) == (q+p))
        stopifnot(nrow(r_i) == T_i & ncol(r_i) == 1)
        stopifnot(nrow(I_i) == T_i & ncol(I_i) == T_i)

        Sigman_summand[i, , ] <- t(D_i) %*% I_i %*% r_i %*% t(r_i) %*% t(I_i) %*% D_i

        # deriv_r_i is \partial r(\theta) / \partial \theta^T for the i-th individual
        deriv_r_i <- cbind( - as.numeric(exp(Zdm_i %*% alpha_hat)) * Zdm_i,
                            - as.numeric(exp(- cA_i * Xdm_i %*% beta_hat)) * Y_i * cA_i * Xdm_i)
        stopifnot(nrow(deriv_r_i) == T_i & ncol(deriv_r_i) == (q+p))

        H_ii <- deriv_r_i %*% Mn_inv %*% t(D_i) / sample_size
        stopifnot(nrow(H_ii) == T_i & ncol(H_ii) == T_i)

        I_minus_H_i <- diag(1, nrow = T_i, ncol = T_i)
        I_minus_H_i_inv <- solve(I_minus_H_i - H_ii)

        Sigman_tilde_summand[i, , ] <- t(D_i) %*% I_i %*% I_minus_H_i_inv %*% r_i %*%
            t(r_i) %*% t(I_minus_H_i_inv) %*% t(I_i) %*% D_i
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size

    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])

    Sigman_tilde <- apply(Sigman_tilde_summand, c(2,3), sum) / sample_size

    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])

    ### 4. return the result with variable names ###

    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames

    ### 5. calculate confidence interval

    conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
    c <- qt(1 - 0.05/2, df = sample_size - p - q)
    conf_int_adjusted <- cbind(beta_hat - c * beta_se_adjusted,
                               beta_hat + c * beta_se_adjusted)
    colnames(conf_int) <- colnames(conf_int_adjusted) <- c("2.5 %", "97.5 %")

    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                conf_int = conf_int, conf_int_adjusted = conf_int_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}
