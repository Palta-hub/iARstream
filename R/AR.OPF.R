# DEFINICION CLASE --------------------------------------------------------

AR_OPF_class <- S7::new_class(
  "AR.OPF",
  parent = Bayesian,
  package = "iARstream",
  properties = list(

    y     = S7::new_property(S7::class_numeric),
    p     = S7::new_property(S7::class_numeric),
    error_distribution  = S7::new_property(S7::class_character, default = "gaussian"),

    phi     = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted  = S7::new_property(S7::class_any, default = NULL),
    sigma_est     = S7::new_property(S7::class_any, default = NULL)),

  validator = function(self) {
    if (length(self@p) != 1 || self@p <= 0) {
      return("Property @p must be a positive integer indicating the model order.")
    }
    if (length(self@y) <= self@p) {
      return("The length of time series @y must be greater than order @p.")
    }
    valid_distribution <- c("gaussian", "laplace", "student")
    if (!(self@error_distribution %in% valid_distribution)) {
      return("Property @error_distribution must be 'gaussian', 'laplace' or 'student'.")
    }
    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, AR_OPF_class) <- function(x, ...) {
  header    <- sprintf("Online Particle Filter applied to AR(%s)", x@p)
  line_distribution <- paste("Error Distribution =", tools::toTitleCase(x@error_distribution))
  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", paste(round(x@phi, 3), collapse = ", "))
  }

  lines_vector <- c(header, line_distribution, line_phi)
  max_width    <- max(nchar(lines_vector))
  sep_double <- strrep("=", max_width)
  sep_single <- strrep("-", max_width)
  cat(sep_double, "\n")
  cat(header, "\n")
  cat(sep_double, "\n")
  cat(line_distribution, "\n")
  cat(line_phi, "\n")
  cat(sep_single, "\n")
  invisible(x)
}

method(est, AR_OPF_class) <- function(x) {
  res <- AR.OPF(
    y     = x@y,
    p     = x@p,
    error_distribution  = x@error_distribution
  )
  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted
  x@sigma_est <- res@sigma_est
  return(x)
}

# FUNCION AR.OPF ----------------------------------------------------------

#' Online Particle Filter for AR(p) Models
#'
#' Estimates the parameters and dynamic volatility of an Autoregressive AR(p) model using
#' an Online Particle Filter (Sequential Monte Carlo) algorithm. The method performs
#' sequential observation-by-observation updates and incorporates Liu & West shrinkage
#' to prevent particle degeneracy and maintain parameter stability. It supports
#' Gaussian, Laplace, and Student-t error distributions for robust estimation under
#' heavy-tailed noise.
#'
#' @details
#' The algorithm utilizes \code{Np = 1000 * p} particles to approximate the posterior
#' distribution of both the AR parameters and the log-volatility. Parameter evolution is
#' controlled via the Liu & West shrinkage method with a discount factor of \code{a = 0.99}.
#'
#' To initialize the dynamic variance, the algorithm computes the empirical variance of a
#' minimum warm-up window of \code{max(p + 1, 15)} observations. For robust estimation,
#' if the Student-t distribution is selected, the degrees of freedom are fixed at \code{nu = 5}.
#' Similar to classical optimization methods, the estimated AR parameters are box-constrained
#' within \code{[-0.999, 0.999]} at each resampling step to ensure numerical stability.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param p An integer specifying the order of the AR(p) model.
#' @param error_distribution A character string specifying the likelihood error distribution.
#'   Options are \code{"gaussian"} (default), \code{"laplace"}, or \code{"student"}.
#'
#' @return An object of class \code{AR_OPF_class} containing the following slots:
#' \item{y}{The original input time series.}
#' \item{p}{The specified AR order.}
#' \item{error_distribution}{The error distribution utilized.}
#' \item{phi}{A numeric vector with the final estimated AR parameters.}
#' \item{history}{A numeric matrix of dimension \code{n x p} tracking the sequential AR parameter estimates.}
#' \item{fitted}{A numeric vector containing the 1-step ahead in-sample predicted values.}
#' \item{sigma_est}{A numeric vector containing the dynamically estimated sequential volatility (\code{sigma}).}
#'
#' @importFrom utils tail
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(model = list(ar = c(0.6, -0.3)), n = 100)
#' fit <- AR.OPF(y = y, p = 2, error_distribution = "laplace")
#' fit
AR.OPF <- function(y, p, error_distribution = c("gaussian", "laplace", "student")) {
  if (p < 1) {
    stop("Property @p must be a positive integer indicating the model order.")
  }
  if (length(y) <= p) {
    stop("The length of time series @y must be greater than order @p.")
  }
  error_distribution <- match.arg(error_distribution)
  n <- length(y)
  nu <- 5
  Np <- 1000 * p
  phi0 <- matrix(runif(Np * p, min = -0.99, max = 0.99), Np, p)
  min_window <- max(p + 1, 15)
  init_var <- var(y[1:min_window], na.rm = TRUE)
  if (is.na(init_var) || init_var <= 0) init_var <- 1
  log_sigma2_init <- rnorm(Np, mean = log(init_var), sd = 1)
  theta <- cbind(phi0, log_sigma2_init)
  phi <- matrix(NA, n, p)
  sigma_est <- numeric(n)
  fitted_vals <- numeric(n)
  a <- 0.99
  h_sq <- 1 - a^2

  for (t in (min_window + 1):n) {
    lags <- rev(y[(t - p):(t - 1)])
    mbar <- colMeans(theta)
    theta_contr <- a * theta + (1 - a) * matrix(mbar, nrow = Np, ncol = p + 1, byrow = TRUE)
    phi_contr <- theta_contr[, 1:p, drop = FALSE]
    log_sigma2_contr <- theta_contr[, p + 1]
    sigma_y <- exp(0.5 * log_sigma2_contr)
    mu <- as.vector(phi_contr %*% lags)
    fitted_vals[t] <- mean(mu)
    err <- y[t] - mu

    if (error_distribution == "gaussian") {
      log_w <- -0.5 * (err / sigma_y)^2 - log(sigma_y) - 0.5 * log(2 * pi)
    } else if (error_distribution == "laplace") {
      log_w <- -abs(err) / sigma_y - log(2 * sigma_y)
    } else if (error_distribution == "student") {
      log_w <- dt(err / sigma_y, df = nu, log = TRUE) - log(sigma_y)
    }
    w <- exp(log_w - max(log_w))
    sum_w <- sum(w)
    if (is.na(sum_w) || sum_w < 1e-12) {
      w <- rep(1 / Np, Np)
    } else {
      w <- w / sum_w
    }

    theta_est_t <- colSums(theta_contr * w)
    phi[t, ] <- theta_est_t[1:p]
    sigma_est[t] <- sum(sigma_y * w)
    theta_centered <- t(t(theta) - theta_est_t)
    cov_theta <- t(theta_centered) %*% (theta_centered * w)
    cov_theta <- cov_theta + diag(1e-6, p + 1)
    u1 <- runif(1, min = 0, max = 1/Np)
    u <- u1 + seq(0, Np - 1) / Np
    cw <- cumsum(w)
    idx <- findInterval(u, cw) + 1
    idx[idx > Np] <- Np
    L <- chol(h_sq * cov_theta)
    noise <- matrix(rnorm(Np * (p + 1)), Np, p + 1) %*% L
    theta <- theta_contr[idx, , drop = FALSE] + noise
    theta[, 1:p] <- pmax(pmin(theta[, 1:p], 0.999), -0.999)
  }
  phi_final <- as.numeric(tail(phi, 1))

  return(AR_OPF_class(
    y             = y,
    p             = p,
    error_distribution = error_distribution,
    phi           = phi_final,
    history       = phi,
    fitted        = fitted_vals,
    sigma_est = sigma_est
  ))
}
