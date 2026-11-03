# DEFINICION CLASE --------------------------------------------------------

iAR_OBR_class <- S7::new_class(
  "iAR.OBR",
  parent = Bayesian,
  package = "iARstream",
  properties = list(

    y     = S7::new_property(S7::class_numeric),
    times = S7::new_property(S7::class_numeric),

    phi = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted  = S7::new_property(S7::class_any, default = NULL)),

  validator = function(self) {
    if (!is.null(dim(self@times)) || !is.null(dim(self@y))) {
      return("Propertys @times and @y must be numeric vectors")
    }
    if (length(self@times) != length(self@y)) {
      return("The length of @times must be the same as the length of @y")
    }
    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, iAR_OBR_class) <- function(x, ...) {
  header    <- sprintf("Online Bayesian Regression applied to iAR")
  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", paste(round(x@phi, 3), collapse = ", "))
  }

  lines_vector <- c(header, line_phi)
  max_width    <- max(nchar(lines_vector))
  sep_double <- strrep("=", max_width)
  sep_single <- strrep("-", max_width)
  cat(sep_double, "\n")
  cat(header, "\n")
  cat(sep_double, "\n")
  cat(line_phi, "\n")
  cat(sep_single, "\n")
  invisible(x)
}

method(est, iAR_OBR_class) <- function(x) {
  res <- iAR.OBR(
    y        = x@y,
    times    = x@times
    )
    x@phi     <- res@phi
    x@history <- res@history
    x@fitted  <- res@fitted
  return(x)
}

# FUNCION iAR.OBR -------------------------------------------------------------
#' Online Bayesian Regression for Irregularly Observed Autoregressive (iAR) Models
#'
#' @description
#' Performs online parameter estimation for irregularly observed Autoregressive (iAR) models
#' using a discrete-grid Online Bayesian Regression (OBR) approach.
#'
#' The algorithm accounts for irregular observation intervals using the exact continuous-time
#' prediction rule \eqn{\hat{y}_j = \phi^{\Delta t_j} y_{j-1}}, where \eqn{\Delta t_j = t_j - t_{j-1}},
#' combined with the theoretical conditional variance \eqn{1 - \phi^{2\Delta t_j}}.
#'
#' @details
#' The function approximates the continuous posterior distribution by discretizing the
#' parameter space of \eqn{\phi} into a fine grid of 1000 equally spaced candidate values
#' ranging from \eqn{0.001} to \eqn{0.999}.
#'
#' Starting from a uniform prior on the logarithmic scale, the algorithm sequentially
#' updates the posterior distribution utilizing the exact Gaussian log-likelihood of the
#' iAR process. To guarantee numerical stability over long sequences and prevent floating-point
#' underflow, it employs the Log-Sum-Exp trick before converting the log-posteriors back
#' to normalized probabilities. At each time step, the sequential point estimate is
#' extracted as the posterior mean of the distribution.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param times A numeric vector of observation times corresponding to \code{y}.
#'   Must have the same length as \code{y}.
#'
#' @return An object of class \code{iAR_OBR_class} containing the following slots:
#' \item{y}{The original input time series.}
#' \item{times}{The vector of observation times.}
#' \item{phi}{A numeric scalar with the final estimated iAR parameter (the posterior mean at the final step).}
#' \item{history}{A numeric vector tracking the sequence of posterior means at each time step.}
#' \item{fitted}{A numeric vector containing the in-sample fitted values (\eqn{\hat{y}_t}).}
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
#'
#' y_norm <- rnorm(200)
#'
#' fit_norm <- iAR.OBR(y = y_norm, times = times)
#' fit_norm
#'
#' # ------------------------------------------------------------------
#' # Example 2: iAR Student-t model
#' # ------------------------------------------------------------------
#' set.seed(123)
#' o <- iAR::utilities()
#' o <- iAR::gentime(o, n = 200, distribution = "uniform")
#' times <- o@times
#'
#' y_t <- rt(200, df = 5)
#'
#' fit_t <- iAR.OBR(y = y_t, times = times)
#' fit_t
iAR.OBR <- function(y, times) {
  y <- as.vector(as.numeric(y))
  times <- as.vector(as.numeric(times))
  if (length(y) < 2) {
    stop("La serie de tiempo 'y' debe tener al menos 2 observaciones.")
  }
  if (length(y) != length(times)) {
    stop("Los vectores 'y' y 'times' deben tener la misma longitud.")
  }
  n <- length(y)
  phi_est <- rep(NA, n)

  # 1. Definimos el espacio de búsqueda (Prior Discreto)
  # Creamos 1000 posibles candidatos para phi entre 0.001 y 0.999
  phis <- seq(0.001, 0.999, length.out = 1000)

  # Prior uniforme inicial (en escala logarítmica para estabilidad)
  log_prior <- rep(0, length(phis))

  for (t in 2:n) {
    delta_t <- times[t] - times[t-1]
    if (delta_t < 1e-6) delta_t <- 1e-6

    # 2. Likelihood exacto del proceso iAR
    # Para cada posible phi, calculamos su media y varianza teórica exacta
    mu <- (phis^delta_t) * y[t-1]
    var_cond <- 1 - phis^(2 * delta_t)

    # Prevenir varianzas cero por redondeo numérico
    var_cond[var_cond < 1e-8] <- 1e-8

    # Calculamos la Log-Verosimilitud Gaussiana vectorizada
    log_lik <- -0.5 * log(var_cond) - 0.5 * (y[t] - mu)^2 / var_cond

    # 3. Actualización Bayesiana (Posterior = Prior + Likelihood)
    log_post <- log_prior + log_lik

    # Estabilización numérica: restamos el máximo (Log-Sum-Exp trick)
    log_post <- log_post - max(log_post)

    # Pasamos a probabilidades reales y normalizamos (que sumen 1)
    post_probs <- exp(log_post)
    post_probs <- post_probs / sum(post_probs)

    # 4. Estimación: La Media Posterior de la distribución
    phi_est[t] <- sum(phis * post_probs)

    # El posterior de hoy es el prior de mañana
    log_prior <- log_post
  }

  phi_final <- phi_est[n] # Tomamos la estimación del último paso
  phi_hist_matrix <- phi_est

  # Calculamos los valores ajustados teóricos
  y_hat <- rep(NA, n)
  for (t in 2:n) {
    delta_t <- times[t] - times[t-1]
    y_hat[t] <- (phi_final^delta_t) * y[t-1]
  }

  iAR_OBR_class(
    y       = as.numeric(y),
    times   = as.numeric(times),
    phi     = as.numeric(phi_final),
    history = as.numeric(phi_hist_matrix),
    fitted  = as.numeric(y_hat)
  )
}


