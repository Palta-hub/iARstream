# DEFINICION CLASE --------------------------------------------------------

iAR_ONS_class <- S7::new_class(
  "iAR.ONS",
  parent = Classical,
  package = "iARstream",
  properties = list(

    y      = S7::new_property(S7::class_numeric),
    times  = S7::new_property(S7::class_numeric),
    phi_0  = S7::new_property(S7::class_any, default = NULL),
    loss_function   = S7::new_property(S7::class_character, default = "quadratic"),

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
    valid_losses <- c("quadratic", "absolute", "logcosh", "huber")
    if (!(self@loss_function %in% valid_losses)) {
      return("Property @loss_function must be 'quadratic', 'absolute', 'logcosh', or 'huber'.")
    }
    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, iAR_ONS_class) <- function(x, ...) {
  header    <- "Online Newton Step applied to iAR"
  line_loss <- paste("Loss Function =", tools::toTitleCase(x@loss_function))
  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", round(x@phi, 3))
  }

  lines_vector <- c(header, line_loss, line_phi)
  max_width    <- max(nchar(lines_vector))
  sep_double <- strrep("=", max_width)
  sep_single <- strrep("-", max_width)
  cat(sep_double, "\n")
  cat(header, "\n")
  cat(sep_double, "\n")
  cat(line_loss, "\n")
  cat(line_phi, "\n")
  cat(sep_single, "\n")
  invisible(x)
}

method(est, iAR_ONS_class) <- function(x) {
  res <- iAR.ONS(
    y      = x@y,
    times  = x@times,
    phi_0  = x@phi_0,
    loss_function   = x@loss_function
  )
  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted
  return(x)
}

# FUNCION ONS -------------------------------------------------------------

#' Online Newton Step for Irregularly Observed Autoregressive (iAR) Models
#'
#' Performs online parameter estimation for irregularly observed Autoregressive (iAR) models
#' via the Online Newton Step (ONS) optimization algorithm.
#'
#' The algorithm accounts for irregular observation intervals using the continuous-time
#' prediction rule \eqn{\hat{y}_j = \phi^{\Delta t_j} y_{j-1}}, where \eqn{\Delta t_j = t_j - t_{j-1}}.
#' Optimization is conducted in an unconstrained parameter space via the logistic link
#' \eqn{\gamma = \text{logit}(\phi)} to strictly enforce \eqn{\phi \in (0, 1)}. It supports
#' four customizable robust loss functions (quadratic, absolute, Huber, and log-cosh) to
#' accommodate diverse noise structures and outliers.
#'
#' @details
#' Parameter updates are computed using a Damped Online Newton Step mechanic, which maintains
#' an exponentially decayed pseudo-Hessian running sum \eqn{A_t = \beta A_{t-1} + g_t^2}
#' initialized with `epsilon = 1e-3` and a decay rate of `beta = 0.99`. The optimization utilizes
#' a Cosine Annealing learning rate schedule across 20 epochs, decaying smoothly from `eta = 0.2`
#' down to `eta_min = 0.01`.
#'
#' If no initial parameter `phi_0` is provided, it is smartly initialized using the empirical
#' autocorrelation scaled by the mean time gap. The final parameter estimate `phi` is computed
#' using Polyak-Ruppert averaging over the last 5 epochs, removing a 20% burn-in period from the
#' time series trajectory to ensure robust convergence.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param times A numeric vector of observation times corresponding to `y`. Must have the same length as `y`.
#' @param phi_0 An optional numeric value specifying the initial estimate for the iAR parameter in `(0, 1)`.
#'   If `NULL` (default), it is initialized using the empirical autocorrelation.
#' @param loss_function A character string specifying the loss function to minimize.
#'   Options are `"quadratic"` (default), `"absolute"`, `"logcosh"`, or `"huber"`.
#'
#' @return An object of class `iAR_ONS_class` containing the following slots:
#' \item{y}{The original input time series.}
#' \item{times}{The vector of observation times.}
#' \item{phi_0}{The initial parameter value used (either user-provided or auto-initialized).}
#' \item{loss_function}{The loss function utilized.}
#' \item{phi}{A numeric scalar with the final robust estimated iAR parameter.}
#' \item{history}{A numeric matrix tracking parameter estimates during the last 5 epochs (Polyak-Ruppert averaging window).}
#' \item{fitted}{A numeric vector containing the in-sample fitted values.}
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
#' fit_norm <- iAR.ONS(y = y_norm, times = times, loss_function = "quadratic")
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
#' fit_t <- iAR.ONS(y = y_t, times = times, loss_function = "logcosh")
#' fit_t
iAR.ONS <- function(y, times,
                    phi_0 = NULL,
                    loss_function = c("quadratic", "absolute", "logcosh", "huber")) {
  if (length(y) != length(times)) {
    stop("The length of 'times' must be equal to the length of 'y'.")
  }

  epochs <- 20
  eta <- 0.2
  epsilon <- 1e-3
  beta <- 0.99
  delta <- 1.345
  loss_function <- match.arg(loss_function)
  n <- length(y)

  # Inicialización de phi_0 si no se provee
  if (is.null(phi_0)) {
    r1 <- abs(cor(y[2:n], y[1:(n-1)], use = "complete.obs"))
    r1 <- pmax(r1, 1e-4)
    mean_delta <- mean(diff(times))
    phi_0_smart <- r1^(1 / mean_delta)
    phi_0 <- pmax(pmin(phi_0_smart, 0.999), 0.001)
  }

  phi_0_safe <- pmax(pmin(phi_0, 0.999), 0.001)
  gamma_prev <- log(phi_0_safe / (1 - phi_0_safe))
  y_hat <- rep(NA, n)

  # Configuración para promediar las últimas épocas (Polyak-Ruppert Averaging)
  avg_epochs <- 5
  start_avg <- epochs - avg_epochs + 1
  gamma_hist_matrix <- matrix(NA, nrow = n, ncol = avg_epochs)

  # Inicialización del pseudo-Hessiano
  A_t <- epsilon

  # <-- NUEVO: Valor mínimo al que llegará el learning rate
  eta_min <- 0.01

  for (e in 1:epochs) {

    # === NUEVO: Decaimiento del Learning Rate (Cosine Annealing) ===
    if (epochs > 1) {
      fraccion_avance <- (e - 1) / (epochs - 1)
      eta_e <- eta_min + 0.5 * (eta - eta_min) * (1 + cos(fraccion_avance * pi))
    } else {
      eta_e <- eta
    }

    for (j in 2:n) {
      delta_t <- times[j] - times[j-1]

      # Mapeo de gamma (inversamente Logit) a phi
      phi_actual <- 1 / (1 + exp(-gamma_prev))
      phi_actual <- pmax(pmin(phi_actual, 0.999), 1e-6)

      # Predicción y error
      y_hat_j <- (phi_actual^delta_t) * y[j-1]
      error_j <- y[j] - y_hat_j
      if (e == epochs) y_hat[j] <- y_hat_j

      # Derivadas (Regla de la cadena)
      d_yhat_d_phi <- delta_t * (phi_actual^(delta_t - 1)) * y[j-1]
      d_phi_d_gamma <- max(phi_actual * (1 - phi_actual), 0.01)

      # Selección de la función de pérdida
      if (loss_function == "quadratic") {
        grad_phi <- -error_j * d_yhat_d_phi
      } else if (loss_function == "absolute") {
        grad_phi <- -sign(error_j) * d_yhat_d_phi
      } else if (loss_function == "logcosh") {
        grad_phi <- -tanh(error_j) * d_yhat_d_phi
      } else if (loss_function == "huber") {
        grad_phi <- if (abs(error_j) <= delta) -error_j * d_yhat_d_phi else -delta * sign(error_j) * d_yhat_d_phi
      }

      # Gradiente final respecto a gamma
      grad_gamma <- grad_phi * d_phi_d_gamma

      # === ACTUALIZACIÓN ONS CON DESCUENTO (DAMPED ONS) ===
      # El acumulador olvida gradientes antiguos de forma exponencial
      A_t <- beta * A_t + (grad_gamma^2)

      # Paso de Newton puro escalado por la inversa de la curvatura reciente
      gamma_prev <- gamma_prev - eta_e * (grad_gamma / A_t)

      # Guardar historial para el promedio final
      if (e >= start_avg) {
        col_idx <- e - start_avg + 1
        gamma_hist_matrix[j, col_idx] <- gamma_prev
      }
    }
  }

  # Descartar el "burn-in" de la serie de tiempo (20% inicial)
  burn_in <- max(2, floor(0.2 * n))
  gamma_final <- mean(gamma_hist_matrix[burn_in:n, ], na.rm = TRUE)

  # Transformar de vuelta a phi
  phi_final <- 1 / (1 + exp(-gamma_final))
  phi_hist_matrix <- 1 / (1 + exp(-gamma_hist_matrix)) # Para que el usuario vea valores entre 0 y 1

  # Devolver directamente la clase estructurada
  iAR_ONS_class(
    y       = as.numeric(y),
    times   = as.numeric(times),
    phi_0   = phi_0,
    loss_function    = loss_function,
    phi     = as.numeric(phi_final),
    history = phi_hist_matrix,
    fitted  = as.numeric(y_hat)
  )
}
