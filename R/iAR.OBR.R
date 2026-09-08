# DEFINICION CLASE --------------------------------------------------------

iAR_OBR_class <- S7::new_class(
  "iAR.OBR",
  parent = Bayesian,
  package = "iARest",
  properties = list(
    # Entradas
    y     = S7::new_property(S7::class_numeric),
    times = S7::new_property(S7::class_numeric),

    #Salidas
    phi = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted  = S7::new_property(S7::class_any, default = NULL)
  ),
  validator = function(self) {
    # y, times deben ser vectores unidimensionales
    if (!is.null(dim(self@times)) || !is.null(dim(self@y))) {
      return("Propertys @times and @y must be numeric vectors")
    }
    # deben tener la misma longitud
    if (length(self@times) != length(self@y)) {
      return("The length of @times must be the same as the length of @y")
    }

    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, iAR_OBR_class) <- function(x, ...) {
  # Armar cadenas de texto de cada linea
  header    <- sprintf("Online Bayesian Regression applied to iAR")

  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", paste(round(x@phi, 3), collapse = ", "))
  }

  # Calcular la longitud maxima del contenido
  lines_vector <- c(header, line_loss, line_phi)
  max_width    <- max(nchar(lines_vector))

  # Generar las separaciones con la funcion strrep()
  sep_double <- strrep("=", max_width)
  sep_single <- strrep("-", max_width)

  # Imprimir en consola
  cat(sep_double, "\n")
  cat(header, "\n")
  cat(sep_double, "\n")
  cat(line_loss, "\n")
  cat(line_phi, "\n")
  cat(sep_single, "\n")

  invisible(x)
}

method(est, iAR_OBR_class) <- function(x) {
  res <- iAR.ORB(
    y        = x@y,
    times    = x@times
    )

  #devoler mateiz
    x@phi     <- res@phi
    x@history <- res@history
    x@fitted  <- res@fitted
  return(x)
}

# FUNCION iAR.OBR -------------------------------------------------------------

#' Online Bayesian Regression for irregular Autoregressive (iAR) Models
#'
#' @description
#' Fits an irregularly observed autoregressive time series model (iAR) using an
#' Online Bayesian Regression (OBR) approach. It utilizes a discrete grid filter
#' to approximate the exact iAR likelihood and updates the posterior distribution
#' of the autoregressive parameter sequentially. Prior parameters (uniform grid)
#' are fixed internally.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param times A numeric vector of observation times corresponding to \code{y}.
#'   Must have the same length as \code{y}.
#'
#' @return An S7 object of class \code{iAR.OBR} containing:
#' \itemize{
#'   \item \code{y}: The input numeric time series.
#'   \item \code{times}: Observation times.
#'   \item \code{phi}: A numeric value representing the final estimated iAR coefficient.
#'   \item \code{history}: A numeric vector containing the sequence of posterior means at each time step.
#'   \item \code{fitted}: A numeric vector of fitted values (\eqn{\hat{y}_t}).
#' }
#'
#' @export
#'
#' @examples
#' # Simulate irregular times and an iAR process
#' set.seed(123)
#' times <- cumsum(rexp(100, rate = 1))
#' y <- rnorm(100) # Replace with actual iAR simulation if available
#'
#' # Fit the iAR model via OBR
#' fit <- iAR.OBR(y = y, times = times)
#' fit@phi


iAR.OBR <- function(y, times) {
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


