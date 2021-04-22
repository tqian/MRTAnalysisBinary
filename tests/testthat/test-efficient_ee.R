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
      as.numeric(efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$beta_hat),
      as.vector(c(0.07473133, 0.43884629)),
      tolerance = 1e-7
    )})

test_that(
  "check alpha_hat",
  {
    expect_equal(
      as.numeric(efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$alpha_hat),
      as.vector(c(-1.4841301, 0.4835868, 0.2641755)),
      tolerance = 1e-6
    )})

test_that(
  "check beta_se",
  {
    expect_equal(
      as.numeric(efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$beta_se),
      as.vector(c(0.1161303, 0.1613436)),
      tolerance = 1e-6
    )})

test_that(
  "check alpha_se",
  {
    expect_equal(
      as.numeric(efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$alpha_se),
      as.vector(c(0.07988023, 0.19605660, 0.09279685)),
      tolerance = 1e-7
    )})

test_that(
  "check beta_se_adjusted",
  {
    expect_equal(
      as.numeric(efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$beta_se_adjusted),
      as.vector(c(0.1173742, 0.1630781)),
      tolerance = 1e-6
    )})

test_that(
  "check alpha_se_adjusted",
  {
    expect_equal(
      as.numeric(efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$alpha_se_adjusted),
      as.vector(c(0.08070538, 0.19811555, 0.09375403)),
      tolerance = 1e-7
    )})

test_that(
  "check varcov",
  {
    expect_equal(
      efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$varcov,
      matrix(c(0.006380851, -0.01308811, 0.003305393, -0.006204126, 0.009073781,
               -0.013088110, 0.03843819, -0.014921059, 0.011267004, -0.019398075,
               0.003305393, -0.01492106, 0.008611255, -0.001976675, 0.003608133,
               -0.006204126, 0.01126700, -0.001976675, 0.013486249, -0.017351647,
               0.009073781, -0.01939807, 0.003608133, -0.017351647, 0.026031756),
             nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

test_that(
  "check varcov_adjusted",
  {
    expect_equal(
      efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE
        )$varcov_adjusted,
      matrix(c(0.006513359, -0.01336179, 0.003374504, -0.006332914, 0.009263806,
               -0.013361786, 0.03924977, -0.015233734, 0.011497775, -0.019806389,
               0.003374504, -0.01523373, 0.008789818, -0.002014988, 0.003682553,
               -0.006332914, 0.01149778, -0.002014988, 0.013776692, -0.017724414,
               0.009263806, -0.01980639, 0.003682553, -0.017724414, 0.026594464),
             nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

# extract the confidence interval from the output and drop its column and row names

conf_int <- efficient_ee(
  dta = dgm_sample,
  id_varname = "userid",
  decision_time_varname = "day",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("time_var1", "time_var2"),
  moderator_varname = "time_var1",
  rand_prob_varname = "prob_A",
  EMEE_initial_value = TRUE)$conf_int
dimnames(conf_int) <- c()

test_that(
  "check conf_int",
  {
    expect_equal(conf_int,
                 matrix(c(-0.1528841, 0.3023467, 0.1226128, 0.7550797),
                        nrow = 2, ncol = 2, byrow = TRUE),
                 tolerance = 1e-6
    )})

# extract the adjusted confidence interval from the output and drop its column and row names

conf_int_adjusted <- efficient_ee(
  dta = dgm_sample,
  id_varname = "userid",
  decision_time_varname = "day",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("time_var1", "time_var2"),
  moderator_varname = "time_var1",
  rand_prob_varname = "prob_A",
  EMEE_initial_value = TRUE)$conf_int_adjusted
dimnames(conf_int_adjusted) <- c()

test_that(
  "check conf_int_adjusted",
  {
    expect_equal(conf_int_adjusted,
                 matrix(c(-0.1582858, 0.3077485, 0.1150953, 0.7625972),
                        nrow = 2, ncol = 2, byrow = TRUE),
                 tolerance = 1e-6
    )})

dgm_sample$A_test <- dgm_sample$A
dgm_sample$A_test[dgm_sample$avail == 1] <- rep(NA, 2420)

test_that(
  "check error when treatment indicator is NA where availability = 1",
  {
    expect_error(
      efficient_ee(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A_test",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        EMEE_initial_value = TRUE),
      "Treatment indicator is NA where availability = 1."
    )
  }
)
