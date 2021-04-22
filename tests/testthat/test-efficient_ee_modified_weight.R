# generating data

set.seed(123)

dgm_demo <- function(sample_size, total_T) {

  alpha_0 <- - 1.5
  alpha_1 <- 0.5
  alpha_2 <- 0.3

  beta_0 <- 0.1
  beta_1 <- 0.3

  # With the above parameter specification, the range of success probability of Y would be
  # [exp(-1.5), exp(-1.5 + 0.5 + 0.3 + 0.1 + 0.3)] = [0.223, 0.741]

  df_names <- c("userid", "time", "time_var1", "time_var2", "Y", "A", "avail", "prob_Y", "prob_Y_A0", "prob_A")
  # time_var1 is time / total_T
  # time_var2 is 1(time > total_T/2)

  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names

  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  dta$time_var1 <- dta$time / total_T
  dta$time_var2 <- as.numeric(dta$time > (total_T/2))

  for (t in 1:total_T) {
    # row index for the rows corresponding to time t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)

    dta$avail[row_index] <- rbinom(sample_size, 1, 0.8) # 0.8 probability to be available

    dta$prob_A[row_index] <- ifelse(t %% 3 == 1, 0.3, ifelse(t %% 3 == 2, 0.5, 0.7))
    dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index]) * dta$avail[row_index] # A can only be 1 if avail = 1
    # We keep prob_A as-is for those observations with avail = 0. It's OK because those prob_A won't be used in the estimation.

    dta$prob_Y_A0[row_index] <- exp(alpha_0 + alpha_1 * dta$time_var1[row_index] + alpha_2 * dta$time_var2[row_index])
    dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * (beta_0 + beta_1 * dta$time_var1[row_index]))
    dta$Y[row_index] <- rbinom(sample_size, 1, dta$prob_Y[row_index])
  }

  return(dta)
}

dgm_sample <- dgm_demo(sample_size = 100, total_T = 30)

# tests

test_that(
  "check beta_hat",
  {
    expect_equal(
      as.numeric(efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$beta_hat),
      as.vector(c(0.07473133, 0.43884629)),
      tolerance = 1e-7
    )})

test_that(
  "check alpha_hat",
  {
    expect_equal(
      as.numeric(efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$alpha_hat),
      as.vector(c(-1.4841301, 0.4835868, 0.2641755)),
      tolerance = 1e-6
    )})

test_that(
  "check beta_se",
  {
    expect_equal(
      as.numeric(efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$beta_se),
      as.vector(c(0.1156131, 0.1596739)),
      tolerance = 1e-6
    )})

test_that(
  "check alpha_se",
  {
    expect_equal(
      as.numeric(efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$alpha_se),
      as.vector(c(0.08276097, 0.20129953, 0.09133819)),
      tolerance = 1e-7
    )})

test_that(
  "check beta_se_adjusted",
  {
    expect_equal(
      as.numeric(efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$beta_se_adjusted),
      as.vector(c(0.1168487, 0.1613801)),
      tolerance = 1e-6
    )})

test_that(
  "check alpha_se_adjusted",
  {
    expect_equal(
      as.numeric(efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$alpha_se_adjusted),
      as.vector(c(0.08364861, 0.20347557, 0.09226770)),
      tolerance = 1e-7
    )})

test_that(
  "check varcov",
  {
    expect_equal(
      efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$varcov,
      matrix(c(0.006849378, -0.01421519, 0.003620799, -0.006387413, 0.009333666,
               -0.014215190, 0.04052150, -0.015126125, 0.011912498, -0.020217459,
               0.003620799, -0.01512613, 0.008342665, -0.002321266, 0.004017219,
               -0.006387413, 0.01191250, -0.002321266, 0.013366383, -0.017095715,
               0.009333666, -0.02021746, 0.004017219, -0.017095715, 0.025495758),
             nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

test_that(
  "check varcov_adjusted",
  {
    expect_equal(
      efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90
        )$varcov_adjusted,
      matrix(c(0.006997089, -0.01452552, 0.003700014, -0.006522308, 0.009532709,
               -0.014525523, 0.04140231, -0.015446323, 0.012164067, -0.020653046,
               0.003700014, -0.01544632, 0.008513328, -0.002370222, 0.004104497,
               -0.006522308, 0.01216407, -0.002370222, 0.013653617, -0.017461601,
               0.009532709, -0.02065305, 0.004104497, -0.017461601, 0.026043543),
             nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

# extract the confidence interval from the output and drop its column and row names

conf_int <- efficient_ee_modified_weight(
  dta = dgm_sample,
  id_varname = "userid",
  decision_time_varname = "day",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("time_var1", "time_var2"),
  moderator_varname = "time_var1",
  rand_prob_varname = "prob_A",
  weight_threshold = 0.90)$conf_int
dimnames(conf_int) <- c()

test_that(
  "check conf_int",
  {
    expect_equal(conf_int,
                 matrix(c(-0.1518703, 0.3013330, 0.1258854, 0.7518072),
                        nrow = 2, ncol = 2, byrow = TRUE),
                 tolerance = 1e-6
    )})

# extract the adjusted confidence interval from the output and drop its column and row names

conf_int_adjusted <- efficient_ee_modified_weight(
  dta = dgm_sample,
  id_varname = "userid",
  decision_time_varname = "day",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("time_var1", "time_var2"),
  moderator_varname = "time_var1",
  rand_prob_varname = "prob_A",
  weight_threshold = 0.90)$conf_int_adjusted
dimnames(conf_int_adjusted) <- c()

test_that(
  "check conf_int_adjusted",
  {
    expect_equal(conf_int_adjusted,
                 matrix(c(-0.1572427, 0.3067053, 0.1184662, 0.7592263),
                        nrow = 2, ncol = 2, byrow = TRUE),
                 tolerance = 1e-6
    )})

dgm_sample$A_test <- dgm_sample$A
dgm_sample$A_test[dgm_sample$avail == 1] <- rep(NA, 2420)

test_that(
  "check error when treatment indicator is NA where availability = 1",
  {
    expect_error(
      efficient_ee_modified_weight(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A_test",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        weight_threshold = 0.90),
      "Treatment indicator is NA where availability = 1."
    )
  }
)
