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
      as.numeric(efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"
        )$beta_hat),
      as.vector(c(0.07588889, 0.43685503)),
      tolerance = 1e-7
    )})

test_that(
  "check alpha_hat",
  {
    expect_equal(
      as.numeric(efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"
        )$alpha_hat),
      as.vector(c(-1.4813222, 0.4747662, 0.2676603)),
      tolerance = 1e-6
    )})

test_that(
  "check beta_se",
  {
    expect_equal(
      as.numeric(efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"
        )$beta_se),
      as.vector(c(0.1161012, 0.1612776)),
      tolerance = 1e-6
    )})

test_that(
  "check alpha_se",
  {
    expect_equal(
      as.numeric(efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"
        )$alpha_se),
      as.vector(c(0.07987289, 0.19625116, 0.09282339)),
      tolerance = 1e-7
    )})

test_that(
  "check beta_se_adjusted",
  {
    expect_equal(
      as.numeric(efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"
        )$beta_se_adjusted),
      as.vector(c(0.1173458, 0.1630112)),
      tolerance = 1e-6
    )})

test_that(
  "check alpha_se_adjusted",
  {
    expect_equal(
      as.numeric(efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"
        )$alpha_se_adjusted),
      as.vector(c(0.08069838, 0.19831366, 0.09378108)),
      tolerance = 1e-7
    )})

test_that(
  "check varcov",
  {
    expect_equal(
      efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A",
        )$varcov,
      matrix(c(0.006379679, -0.01310027, 0.003309269, -0.006205946, 0.009075697,
               -0.013100266, 0.03851452, -0.014943879, 0.011278757, -0.019418665,
               0.003309269, -0.01494388, 0.008616181, -0.001983497, 0.003620061,
               -0.006205946, 0.01127876, -0.001983497, 0.013479499, -0.017339188,
               0.009075697, -0.01941866, 0.003620061, -0.017339188, 0.026010469),
             nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

test_that(
  "check varcov_adjusted",
  {
    expect_equal(
      efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"
        )$varcov_adjusted,
      matrix(c(0.006512228, -0.01337437, 0.003378510, -0.006334869, 0.009265842,
               -0.013374374, 0.03932831, -0.015257219, 0.011509878, -0.019827472,
               0.003378510, -0.01525722, 0.008794891, -0.002021982, 0.003694789,
               -0.006334869, 0.01150988, -0.002021982, 0.013770028, -0.017711846,
               0.009265842, -0.01982747, 0.003694789, -0.017711846, 0.026572636),
             nrow = 5, ncol = 5, byrow = TRUE),
      tolerance = 1e-7
    )})

# extract the confidence interval from the output and drop its column and row names

conf_int <- efficient_ee_twostep(
  dta = dgm_sample,
  id_varname = "userid",
  decision_time_varname = "day",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("time_var1", "time_var2"),
  moderator_varname = "time_var1",
  rand_prob_varname = "prob_A")$conf_int
dimnames(conf_int) <- c()

test_that(
  "check conf_int",
  {
    expect_equal(conf_int,
                 matrix(c(-0.1516696, 0.3034473, 0.1207509, 0.7529592),
                        nrow = 2, ncol = 2, byrow = TRUE),
                 tolerance = 1e-6
    )})

# extract the adjusted confidence interval from the output and drop its column and row names

conf_int_adjusted <- efficient_ee_twostep(
  dta = dgm_sample,
  id_varname = "userid",
  decision_time_varname = "day",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("time_var1", "time_var2"),
  moderator_varname = "time_var1",
  rand_prob_varname = "prob_A")$conf_int_adjusted
dimnames(conf_int_adjusted) <- c()

test_that(
  "check conf_int_adjusted",
  {
    expect_equal(conf_int_adjusted,
                 matrix(c(-0.1570719, 0.3088497, 0.1132370, 0.7604731),
                        nrow = 2, ncol = 2, byrow = TRUE),
                 tolerance = 1e-6
    )})

dgm_sample$A_test <- dgm_sample$A
dgm_sample$A_test[dgm_sample$avail == 1] <- rep(NA, 2420)

test_that(
  "check error when treatment indicator is NA where availability = 1",
  {
    expect_error(
      efficient_ee_twostep(
        dta = dgm_sample,
        id_varname = "userid",
        decision_time_varname = "day",
        treatment_varname = "A_test",
        outcome_varname = "Y",
        control_varname = c("time_var1", "time_var2"),
        moderator_varname = "time_var1",
        rand_prob_varname = "prob_A"),
      "Treatment indicator is NA where availability = 1."
    )
  }
)
