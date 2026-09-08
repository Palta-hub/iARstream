# DEFINICION CLASE --------------------------------------------------------

iAR_OGD_class <- S7::new_class(
  "iAR.OGD",
  parent = Classical, # Asumiendo que Classical ya se definió antes
  package = "iARest",
  properties = list(
    ## Entradas del modelo ##
    y             = S7::class_numeric,
    times         = S7::class_numeric,
    family        = S7::new_property(S7::class_character, default = "norm"),
    loss_function = S7::new_property(S7::class_character, default = "quadratic"),
    phi0          = S7::new_property(S7::class_any, default = NULL),

    ## Salidas del modelo ##
    phi = S7::new_property(S7::class_any, default = NULL),
    history       = S7::new_property(S7::class_any, default = NULL),
    fitted        = S7::new_property(S7::class_any, default = NULL)
  ),
  validator = function(self) {
    # Validar tamaño de y
    if (length(self@y) < 2) {
      return("The length of time series @y must be at least 2.")
    }

    # Validar que y y times tengan la misma longitud
    if (length(self@times) != length(self@y)) {
      return("The length of @times must be the same as the length of @y")
    }

    # Validar opciones de 'family'
    valid_families <- c("norm", "t")
    if (!(self@family %in% valid_families)) {
      return("Property @family must be 'norm' or 't'.")
    }

    # Validar opciones de 'loss_function'
    valid_losses <- c("quadratic", "absolute", "logcosh", "huber")
    if (!(self@loss_function %in% valid_losses)) {
      return("Property @loss_function must be 'quadratic', 'absolute', 'logcosh', or 'huber'.")
    }

    return(NULL) # Todo correcto
  }
)

# METODO PARA LA CLASE ----------------------------------------------------
method(print, iAR_OGD_class) <- function(x, ...) {
  # Armar cadenas de texto de cada linea
  header    <- "Online Gradient Descent applied to iAR "
  line_fam  <- paste("Family        =", x@family)
  line_loss <- paste("Loss Function =", x@loss_function)

  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", round(x@phi, 3))
  }

  # Calcular la longitud maxima del contenido
  lines_vector <- c(header, line_fam, line_loss, line_phi)
  max_width    <- max(nchar(lines_vector))

  # Generar las separaciones con la funcion strrep()
  sep_double <- strrep("=", max_width)
  sep_single <- strrep("-", max_width)

  # Imprimir en consola
  cat(sep_double, "\n")
  cat(header, "\n")
  cat(sep_double, "\n")
  cat(line_fam, "\n")
  cat(line_loss, "\n")
  cat(line_phi, "\n")
  cat(sep_single, "\n")

  invisible(x)
}


method(est, iAR_OGD_class) <- function(x) {
  res <- iAR.OGD(
    y             = x@y,
    times         = x@times,
    family        = x@family,
    loss_function = x@loss_function,
    phi0          = x@phi0
  )

  # Tu función matemática devuelve un elemento llamado "phi_estimate"
  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted

  return(x)
}
# FUNCION OGD -------------------------------------------------------------

#' Online Gradient Descent for iAR Models
#'
#' Estimates the parameters of an irregularly observed Autoregressive (iAR) model
#' using Online Gradient Descent (OGD) coupled with the Adam optimizer. Supports
#' adaptive learning rates, different distribution families, and robust loss functions.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param times A numeric vector of observation times corresponding to \code{y}.
#'   Must have the same length as \code{y}.
#' @param family A character string specifying the assumed distribution family.
#'   Options are \code{"norm"} (default) or \code{"t"}.
#' @param loss_function A character string specifying the loss function to minimize.
#'   Options are \code{"quadratic"} (default), \code{"absolute"}, \code{"logcosh"}, or \code{"huber"}.
#' @param phi0 An optional numeric value containing the initial parameter estimate.
#'   If \code{NULL}, defaults to 0.5.
#'
#' @return A list containing:
#' \item{phi_estimate}{A numeric value with the final robust estimated iAR parameter.}
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' times <- cumsum(rexp(100, rate = 1))
#' y <- rnorm(100) # Replace with actual iAR simulation
#' fit <- iAR.OGD(y = y, times = times)
#' fit@phi


iAR.OGD <- function(y, times,
                    family = c("norm", "t"),
                    loss_function = c("quadratic", "absolute", "logcosh", "huber"),
                    phi0 = NULL) {

  if (length(y) != length(times)) {
    stop("The length of 'times' must be equal to the length of 'y'.")
  }

  family <- match.arg(family)
  loss_function <- match.arg(loss_function)
  delta <- 1.345

  n <- length(y)
  if (is.null(phi0)) phi0 <- 0.5

  ## CONFIGURACION SEGUN FAMILIA ##
  if (family == "norm") {
    eta_max <- 0.1
    epochs      <- 20
    use_mle <- loss_function %in% c("quadratic", "huber")
  } else {
    eta_max <- 0.1
    epochs      <- 20
    use_mle <- FALSE
  }

  eta_min <- 0.01

  sigma2 <- if (use_mle) var(y, na.rm = TRUE) else NULL

  theta_hist <- rep(NA, n)
  y_hat      <- rep(NA, n)
  theta      <- log(phi0 / (1 - phi0))

  ## PARAMETROS DE ADAM ##
  beta1   <- 0.9
  beta2   <- 0.999
  epsilon <- 1e-8
  m       <- 0
  v       <- 0
  global_t <- 0

  ## BUCLE DE EPOCHS ##
  for (e in 1:epochs) {
    fraccion_avance <- (e - 1) / (epochs - 1)
    eta_e <- eta_min + 0.5 * (eta_max - eta_min) * (1 + cos(fraccion_avance * pi))

    for (j in 2:n) {
      global_t <- global_t + 1
      delta_t  <- times[j] - times[j-1]

      phi <- 1 / (1 + exp(-theta))

      # ACOTACION Y PREDICCION DE PHI ##
      if (family == "norm") {
        phi_safe <- min(max(phi, 1e-4), 1 - 1e-4)
        y_hat_j  <- (phi_safe^delta_t) * y[j-1]
      } else {
        phi_safe <- max(phi, 1e-7)
        y_hat_j  <- (phi_safe^delta_t) * y[j-1]
      }

      if (e == epochs) y_hat[j] <- y_hat_j

      error_j <- y[j] - y_hat_j

      ## FUNCIONES DE PERDIDA ##
      if (loss_function == "quadratic") {
        if (use_mle) val_loss <- (error_j^2) / 2
        d_loss <- error_j
      } else if (loss_function == "logcosh") {
        if (use_mle) val_loss <- log(cosh(error_j))
        d_loss <- tanh(error_j)
      } else if (loss_function == "absolute") {
        if (use_mle) val_loss <- abs(error_j)
        d_loss <- sign(error_j)
      } else if (loss_function == "huber") {
        if (abs(error_j) <= delta) {
          if (use_mle) val_loss <- (error_j^2) / 2
          d_loss <- error_j
        } else {
          if (use_mle) val_loss <- delta * abs(error_j) - (delta^2) / 2
          d_loss <- delta * sign(error_j)
        }
      }

      ## GRADIENTE ##
      if (!use_mle) {
        d_yhat_d_phi <- delta_t * (phi_safe^(delta_t - 1)) * y[j-1]
        grad_phi     <- -d_yhat_d_phi * d_loss
      } else {
        var_t         <- 1 - phi_safe^(2 * delta_t)
        d_error_d_phi <- -delta_t * (phi_safe^(delta_t - 1)) * y[j-1]
        d_var_d_phi   <- -2 * delta_t * (phi_safe^(2 * delta_t - 1))

        term1     <- d_var_d_phi / (2 * var_t)
        term2_num <- (d_loss * d_error_d_phi * var_t) - (val_loss * d_var_d_phi)
        term2     <- term2_num / (sigma2 * (var_t^2))

        grad_phi  <- term1 + term2
      }

      ## THETA ##
      if (family == "norm") {
        grad_theta <- grad_phi * (phi_safe * (1 - phi_safe))
      } else {
        grad_theta <- grad_phi * (phi_safe * (1 - phi_safe))
      }

      ## ESTABILIZACION ##
      if (is.na(grad_theta) || is.infinite(grad_theta)) grad_theta <- 0
      grad_theta <- max(min(grad_theta, 5), -5)

      ## ACTUALIZACION ADAM ##
      m <- beta1 * m + (1 - beta1) * grad_theta
      v <- beta2 * v + (1 - beta2) * (grad_theta^2)

      m_hat <- m / (1 - beta1^global_t)
      v_hat <- v / (1 - beta2^global_t)

      theta <- theta - (eta_e / (sqrt(v_hat) + epsilon)) * m_hat

      if (e == epochs) {
        theta_hist[j] <- theta
      }
    }
  }

  ## PROMEDIO FINAL ##
  burn_in       <- max(2, floor(n * 0.2))
  theta_robusto <- mean(theta_hist[burn_in:n], na.rm = TRUE)
  phi_robusto   <- 1 / (1 + exp(-theta_robusto))

  phi_hist <- 1 / (1 + exp(-theta_hist))

  iAR_OGD_class(
    y             = as.numeric(y),
    times         = as.numeric(times),
    family        = family,
    loss_function = loss_function,
    phi0          = phi0,
    phi           = as.numeric(phi_robusto),
    history       = as.vector(phi_hist),
    fitted        = as.vector(y_hat)
  )
}
