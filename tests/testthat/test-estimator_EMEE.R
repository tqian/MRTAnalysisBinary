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
    as.numeric(estimator_EMEE(
      dta = dgm_sample,
      id_varname = "userid",
      decision_time_varname = "day",
      treatment_varname = "A",
      outcome_varname = "Y",
      control_varname = c("time_var1", "time_var2"),
      moderator_varname = "time_var1",
      rand_prob_varname = "prob_A",
      avail_varname = "avail"
    )$beta_hat),
    as.vector(c(0.08114495, 0.42931332))
)})

test_that(
  "check alpha_hat",
  {
    expect_equal(
      as.numeric(estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"
      )$alpha_hat),
      as.vector(c(-1.454236899, 0.456590517, 0.247983626))
    )})

test_that(
  "check beta_se",
  {
    expect_equal(
      as.numeric(estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"
      )$beta_se),
      as.vector(c(0.1301265007, 0.1882974587))
    )})

test_that(
  "check alpha_se",
  {
    expect_equal(
      as.numeric(estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"
      )$alpha_se),
      as.vector(c(0.09699208, 0.24074140, 0.11167475))
    )})

test_that(
  "check beta_se_adjusted",
  {
    expect_equal(
      as.numeric(estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"
      )$beta_se_adjusted),
      as.vector(c(0.13164362, 0.19054544))
    )})

test_that(
  "check alpha_se_adjusted",
  {
    expect_equal(
      as.numeric(estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"
      )$alpha_se_adjusted),
      as.vector(c(0.09813404, 0.24373700, 0.11301692))
    )})

test_that(
  "check varcov",
  {
    expect_equal(
      estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"
      )$varcov,
      matrix(c(0.009407463, -0.01848050, 0.0036293844, -0.0088951364, 0.013687095,
                  -0.018480501, 0.05795642, -0.0213308200, 0.0142698457, -0.028903632,
                0.003629384, -0.02133082, 0.0124712494, -0.0009844053, 0.004400843,
                -0.008895136, 0.01426985, -0.0009844053, 0.0169329062, -0.022639280,
                0.013687095, -0.02890363, 0.0044008428, -0.0226392795, 0.035455933),
                nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

test_that(
  "check varcov_adjusted",
  {
    expect_equal(
      estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"
      )$varcov_adjusted,
      matrix(c(0.009630289, -0.01893458, 0.0037223411, -0.0090997258, 0.01401250,
               -0.018934578, 0.05940773, -0.0218615441, 0.0145917501, -0.02960477,
               0.003722341, -0.02186154, 0.0127728233, -0.0009971158, 0.00450630,
               -0.009099726, 0.01459175, -0.0009971158, 0.0173300425, -0.02317355,
               0.014012503, -0.02960477, 0.0045063005, -0.0231735513, 0.03630756),
             nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

rand_prob_tilde_test <- rep(0.5, length = 20)

test_that(
  "check error when p_tilde is of incorrect length",
  {
    expect_error(
      estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail",
        rand_prob_tilde = rand_prob_tilde_test),
      "rand_prob_tilde is of incorrect length."
    )
  }
)

dgm_sample$A_test <- dgm_sample$A
dgm_sample$A_test[dgm_sample$avail == 1] <- rep(NA, 2420)
test_that(
  "check error when treatmeng indicator is NA where availability = 1",
  {
    expect_error(
      estimator_EMEE(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A_test",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        avail_varname = "avail"),
      "Treatment indicator is NA where availability = 1."
    )
  }
)

v_test <- {}
test_that(
  "check error for find_change_location",
  {
    expect_error(
      find_change_location(v_test),
      "The vector need to have length > 1."
    )
  }
)
