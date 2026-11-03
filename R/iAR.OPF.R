# DEFINICION CLASE --------------------------------------------------------

iAR_OPF_class <- S7::new_class(
  "iAR.OPF",
  parent = Bayesian,
  package = "iARstream",
  properties = list(

    y      = S7::new_property(S7::class_numeric),
    times  = S7::new_property(S7::class_numeric),
    error_distribution   = S7::new_property(S7::class_character, default = "gaussian"),

    phi = S7::new_property(S7::class_any, default = NULL),
    history       = S7::new_property(S7::class_any, default = NULL),
    fitted        = S7::new_property(S7::class_any, default = NULL)),

  validator = function(self) {
    if (length(self@y) < 2) {
      return("The length of time series @y must be at least 2.")
    }
    if (length(self@times) != length(self@y)) {
      return("The length of @times must be the same as the length of @y")
    }
    valid_distribution <- c("gaussian", "student") # <-- Modificado
    if (!(self@error_distribution %in% valid_distribution)) {
      return("Property @error_distribution must be 'gaussian' or 'student'.") # <-- Modificado
    }
    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, iAR_OPF_class) <- function(x, ...) {
  header    <- "Online Particle Filter applied to iAR"
  line_distribution <- paste("Error Distribution =", tools::toTitleCase(x@error_distribution))
  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", round(x@phi, 3))
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

method(est, iAR_OPF_class) <- function(x) {
  res <- iAR.OPF(
    y      = x@y,
    times  = x@times,
    error_distribution   = x@error_distribution
  )
  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted
  return(x)
}

# FUNCION OPF -------------------------------------------------------------

#' Online Particle Filter for Irregularly Observed Autoregressive (iAR) Models
#'
#' Performs online parameter estimation for irregularly observed Autoregressive (iAR) models
#' using an Online Particle Filter (Sequential Monte Carlo) algorithm.
#'
#' The algorithm accounts for irregular observation intervals using the continuous-time
#' prediction rule \eqn{\hat{y}_j = \phi^{\Delta t_j} y_{j-1}}, where \eqn{\Delta t_j = t_j - t_{j-1}},
#' with a dynamically scaled conditional variance \eqn{\text{Var}(y)(1 - \phi^{2\Delta t_j})}.
#' To accommodate diverse noise structures and outliers, it supports Gaussian, Laplace, and
#' Student-t likelihood formulations.
#'
#' @details
#' The algorithm approximates the posterior distribution of the \eqn{\phi} parameter using
#' \code{Np = 500} particles, initialized from a truncated normal distribution on \eqn{[0, 1]}.
#' At each time step, the sequential estimate is extracted as the weighted median of the
#' particle population.
#'
#' To prevent particle degeneracy, systematic resampling is dynamically triggered whenever
#' the Effective Sample Size (ESS) falls below the threshold \code{Np / 2}. Following resampling,
#' particles undergo an artificial evolution step via Liu & West shrinkage (discount factor
#' \code{a = 0.98}) combined with a reflection boundary mechanism to strictly enforce the
#' stationarity and invertibility constraint \eqn{\phi \in (10^{-4}, 1 - 10^{-4})}.
#' For robust estimation under the Student-t distribution, the degrees of freedom are fixed
#' at \code{nu = 5}.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param times A numeric vector of observation times corresponding to `y`. Must have the same length as `y`.
#' @param error_distribution A character string specifying the assumed noise distribution family.
#'   Options are `"gaussian"` (default) or `"student"`.
#'
#' @return An object of class `iAR_OPF_class` containing the following slots:
#' \item{y}{The original input time series.}
#' \item{times}{The vector of observation times.}
#' \item{error_distribution}{The error distribution utilized.}
#' \item{phi}{A numeric scalar with the final estimated iAR parameter.}
#' \item{history}{A numeric vector tracking the sequential iAR parameter estimates at each time step.}
#' \item{fitted}{A numeric vector containing the in-sample fitted values.}
#'
#' @importFrom utils tail
#'
#' @export
#'
#' @examples
#' # ------------------------------------------------------------------
#' # Example 1: iAR Normal model
#' # ------------------------------------------------------------------
#' set.seed(123)
#' o <- iAR::utilities()
#' o <- iAR::gentime(o, n = 200, distribution = "uniform")
#' times <- o@times
#' y_norm <- rnorm(200)
#'
#' fit_norm <- iAR.OPF(y = y_norm, times = times, error_distribution = "gaussian")
#' fit_norm
#'
#' # ------------------------------------------------------------------
#' # Example 2: iAR Student-t model
#' # ------------------------------------------------------------------
#' set.seed(123)
#' o <- iAR::utilities()
#' o <- iAR::gentime(o, n = 200, distribution = "uniform")
#' times <- o@times
#' y_t <- rt(200, df = 5)
#'
#' fit_t <- iAR.OPF(y = y_t, times = times, error_distribution = "student")
#' fit_t
iAR.OPF <- function(y, times,
                    error_distribution = c("gaussian", "student")) {
  if (length(y) < 2) {
    stop("The length of time series @y must be at least 2.")
  }
  if (length(y) != length(times)) {
    stop("The length of @times must be the same as the length of @y")
  }
  error_distribution <- match.arg(error_distribution)
  n <- length(y)
  var_y <- var(y)
  nu <- 5
  Np <- 500
  m0 <- 0.5
  P0 <- 0.1
  a <- 0.98
  mu_process <- 0

  phi <- truncnorm::rtruncnorm(Np, a = 0, b = 1, mean = m0, sd = sqrt(P0))
  log_w <- rep(0, Np)

  phi_est <- numeric(n)
  phi_est[1] <- m0

  # Inicializamos el vector para los valores ajustados (fitted values)
  fitted_vals <- numeric(n)
  fitted_vals[1] <- mean(y) # Una inicialización razonable para el primer t

  h_sq <- 1 - a^2
  ESS_threshold <- Np / 2

  for (t in 2:n) {
    delta_t <- times[t] - times[t-1]
    y_prev <- y[t-1]
    phi_delta <- phi^delta_t

    mu_cond <- phi_delta * y_prev + mu_process * (1 - phi_delta)
    # Guardamos la predicción media esperada (antes de la actualización de pesos)
    fitted_vals[t] <- mean(mu_cond)

    var_cond <- var_y * pmax(0.0001, 1 - phi^(2 * delta_t))
    sigma_t <- sqrt(var_cond)
    err <- y[t] - mu_cond

    if (error_distribution == "gaussian") {
      log_lik <- -0.5 * (err / sigma_t)^2 - log(sigma_t)
    } else if (error_distribution == "student") {
      s_t <- sigma_t * sqrt((nu - 2) / nu)
      log_lik <- dt(err / s_t, df = nu, log = TRUE) - log(s_t)
    }

    log_w <- log_w + log_lik
    max_log_w <- max(log_w, na.rm = TRUE)
    log_w <- log_w - max_log_w
    w <- exp(log_w)

    sum_w <- sum(w, na.rm = TRUE)
    if (is.na(sum_w) || sum_w < 1e-12) {
      w_norm <- rep(1 / Np, Np)
    } else {
      w_norm <- w / sum_w
    }

    ord <- order(phi)
    cum_w <- cumsum(w_norm[ord])
    med_idx <- which(cum_w >= 0.5)[1]
    if (is.na(med_idx)) med_idx <- 1
    phi_est[t] <- phi[ord[med_idx]]

    ESS <- 1 / sum(w_norm^2)
    if (ESS < ESS_threshold) {
      idx <- sample.int(Np, size = Np, replace = TRUE, prob = w_norm)
      phi_res <- phi[idx]
      mpar <- mean(phi_res)
      vpar <- var(phi_res)
      if (is.na(vpar) || vpar < 1e-6) vpar <- 1e-6
      phi_ms <- a * phi_res + (1 - a) * mpar
      noise_sd <- sqrt(max(0, h_sq * vpar))
      phi <- phi_ms + rnorm(Np, mean = 0, sd = noise_sd)
      phi <- abs(phi)
      phi <- 1 - abs(1 - phi)
      phi <- pmax(1e-4, pmin(1 - 1e-4, phi))
      log_w <- rep(0, Np)
    }
  }

  # Instanciamos la clase S7
  phi_final <- as.numeric(tail(phi_est, 1))

  return(iAR_OPF_class(
    y                  = y,
    times              = times,
    error_distribution = error_distribution,
    phi                = phi_final,
    history            = phi_est,
    fitted             = fitted_vals
  ))
}
