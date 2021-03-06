% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator_EMEE_alwaysCenterA.R
\name{estimator_EMEE_alwaysCenterA}
\alias{estimator_EMEE_alwaysCenterA}
\title{Estimates the marginal excursion effect for binary outcome MRT
where the treatment indicator is always centered}
\usage{
estimator_EMEE_alwaysCenterA(
  dta,
  id_varname,
  decision_time_varname,
  treatment_varname,
  outcome_varname,
  control_varname,
  moderator_varname,
  rand_prob_varname,
  avail_varname = NULL,
  estimator_initial_value = NULL
)
}
\arguments{
\item{dta}{a data set in long format}

\item{id_varname}{the variable name that specifies the subject IDs}

\item{decision_time_varname}{the variable name that specifies the decision points}

\item{treatment_varname}{the variable name that specifies the assigned treatments for subjects}

\item{outcome_varname}{the variable name that specifies the outcomes for subjects}

\item{control_varname}{a vector of variable names used to reduce noise, this could be NULL}

\item{moderator_varname}{a vector of variable names of the effect modifiers, this could be NULL}

\item{rand_prob_varname}{the variable name that specifies the treatment randomizing probability in the data set}

\item{avail_varname}{the variable name that specifies the availability of the subjects,
default to be always available at any decision points using NULL}

\item{estimator_initial_value}{a numeric vector of the initial value for the estimator,
its length should be the sum of the length of control and moderator variables plus 2
default to be all 0's using NULL}
}
\value{
Returns the estimated beta with its intercept and alpha with its intercept,
        the standard error of the estimated beta with its intercept and alpha with its intercept,
        the adjusted standard error of the estimated beta with its intercept and alpha with its intercept for small sample,
        the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept,
        the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept for small sample,
        the 95 percent confidence interval for beta_hat, and the adjusted 95 percent confidence interval for beta_hat,
        the dimension of the moderated variables, and the dimension of the control variables,
        the value of the estimating function at the estimated beta and alpha
}
\description{
This estimator assumes that the randomization probability does not
depend on H_t, and that the treatment indicator A is always centered.
This estimator returns the estimates for the marginal excursion effect estimator
with the above two assumptions and provides the estimated variance and
standard error for the estimators, with small sample correction for
an MRT with binary outcome using the "Hat" matrix in the variance estimate
and t-distribution or F-distribution critical value with corrected degrees of freedom.
}
\examples{
estimator_EMEE_alwaysCenterA(dta = dgm_sample,
                         id_varname = "userid",
                         decision_time_varname = "day",
                         treatment_varname = "A",
                         outcome_varname = "Y",
                         control_varname = c("time_var1", "time_var2"),
                         moderator_varname = "time_var1",
                         rand_prob_varname = "prob_A",
                         avail_varname = "avail")

}
