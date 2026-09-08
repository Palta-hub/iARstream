# DEFINICION CLASE --------------------------------------------------------

AR_OGD_class <- S7::new_class(
  "AR.OGD",
  parent = Classical,
  package = "iARest",
  properties = list(
    ## Entradas del modelo ##
    y = S7::class_numeric,
    p = S7::class_numeric,
    loss_function = S7::new_property(S7::class_character, default = "quadratic"),
    phi0 = S7::new_property(S7::class_any, default = NULL),

    ## Salidas del modelo ##
    phi = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted = S7::new_property(S7::class_any, default = NULL)
  ),
  validator = function(self) {
    # 'p' debe ser un unico numero positivo
    if (length(self@p) != 1 || self@p <= 0) {
      return("Property @p must be a positive integer indicating the model order.")
    }

    # 'y' debe tener mas datos que el orden 'p' para poder estimar
    if (length(self@y) <= self@p) {
      return("The length of time series @y must be greater than order @p.")
    }

    # Validar que la funcion de perdida ingresada sea correcta
    valid_losses <- c("quadratic", "absolute", "logcosh", "huber")
    if (!(self@loss_function %in% valid_losses)) {
      return("Property @loss_function must be 'quadratic', 'absolute', 'logcosh', or 'huber'.")
    }

    return(NULL) # NULL = pasó todas las pruebas
  }
)
# METODO PARA LA CLASE ----------------------------------------------------

method(print, AR_OGD_class) <- function(x, ...) {
  # Armar cadenas de texto de cada linea
  header    <- sprintf(" Online Gradient Descent applied to AR(%s) ", x@p)
  line_loss <- paste("Loss Function =", x@loss_function)

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

method(est, AR_OGD_class) <- function(x) {

  # Ejecutamos la función base usando SOLO los parámetros que esta acepta
  res <- AR.OGD(
    y             = x@y,
    p             = x@p,
    loss_function = x@loss_function,
    phi0          = x@phi0
  )

  # Guardar resultados (res) en las propiedades vacías del objeto (x)
  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted

  # Devolver objeto actualizado
  return(x)
}


# FUNCION AR.OGD ----------------------------------------------------------

#' Online Gradient Descent for AR(p) Models
#'
#' Estimates the parameters of an Autoregressive model of order p, AR(p),
#' using Online Gradient Descent (OGD) with adaptive learning rate and customizable loss functions.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param p An integer specifying the order of the AR(p) model.
#' @param loss_function A character string specifying the loss function to minimize.
#'   Options are \code{"quadratic"} (default), \code{"absolute"}, \code{"logcosh"}, or \code{"huber"}.
#' @param phi0 An optional numeric vector of length \code{p} containing initial parameter values.
#'   If \code{NULL}, defaults to a vector of 0.1 values.
#'
#' @return A list containing:
#' \item{phi}{A numeric vector with the final estimated AR parameters.}
#' \item{history}{A numeric matrix of dimension \code{n x p} tracking parameter estimates during the last epoch.}
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(model = list(ar = c(0.5, -0.2)), n = 100)
#' fit <- AR.OGD(y = y, p = 2, loss_function = "quadratic")
#' fit@phi

AR.OGD <- function(y, p,
                   loss_function = c("quadratic", "absolute", "logcosh", "huber"),
                   phi0 = NULL) {

  if (length(y) <= p) {
    stop("The length of time series @y must be greater than order @p.")
  }
  n <- length(y)
  loss_function <- match.arg(loss_function)
  delta <- 1.345
  epochs <- 30
  eta <- 0.1
  eta_min <- 0.001
  phi <- if (is.null(phi0)) rep(0.1, p) else phi0
  phi_hist <- matrix(NA, nrow = n, ncol = p)

  for (e in 1:epochs) {
    fraccion_avance <- (e - 1) / max(1, (epochs - 1))
    eta_e <- eta_min + 0.5 * (eta - eta_min) * (1 + cos(fraccion_avance * pi))

    for (t in (p + 1):n) {
      x_t <- matrix(y[(t - 1):(t - p)], ncol = 1)
      y_hat_t <- sum(phi * x_t)
      error_t <- y[t] - y_hat_t

      if (loss_function == "quadratic") {
        d_loss <- error_t
      } else if (loss_function == "logcosh") {
        d_loss <- tanh(error_t)
      } else if (loss_function == "absolute") {
        d_loss <- sign(error_t)
      } else if (loss_function == "huber") {
        d_loss <- if (abs(error_t) <= delta) error_t else delta * sign(error_t)
      }
      grad_t <- -x_t * d_loss

      phi <- phi - eta_e * as.vector(grad_t)
      phi <- pmax(pmin(phi, 0.999), -0.999)
      if (e == epochs) {
        phi_hist[t, ] <- phi
      }
    }
  }
  phi_final <- colMeans(phi_hist[(p + 1):n, , drop = FALSE], na.rm = TRUE)

  fitted_vals <- rep(NA, n)
  for (t in (p + 1):n) {
    x_t <- y[(t - 1):(t - p)]
    fitted_vals[t] <- sum(phi_final * x_t)
  }

  AR_OGD_class(
    y             = as.numeric(y),
    p             = as.integer(p),
    loss_function = loss_function,
    phi0          = phi0,
    phi           = as.vector(phi_final),
    history       = phi_hist,
    fitted        = fitted_vals
  )
}
