# DEFINICION CLASE --------------------------------------------------------

AR_ONS_class <- S7::new_class(
  "AR.ONS",
  parent = Classical,
  package = "iARest",
  properties = list(
    ## Entradas del modelo ##
    y     = S7::new_property(S7::class_numeric),
    p     = S7::new_property(S7::class_numeric),
    phi_0 = S7::new_property(S7::class_any, default = NULL),
    loss_function  = S7::new_property(S7::class_character, default = "quadratic"),

    ## Salidas del modelo ##
    phi     = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted  = S7::new_property(S7::class_any, default = NULL)
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

    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, AR_ONS_class) <- function(x, ...) {
  # Armar cadenas de texto de cada linea
  header    <- sprintf("Online Newton Step applied to AR(%s) ", x@p)
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

method(est, AR_ONS_class) <- function(x) {

  # Correr función matemática
  res <- AR.ONS(
    y     = x@y,
    p     = x@p,
    phi_0 = x@phi_0,
    loss_function  = x@loss_function
  )

  # Guardar TODOS los resultados que entrega la función en el objeto
  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted

  return(x)
}


# FUNCION AR.ONS ----------------------------------------------------------

#' Online Newton Step for Autoregressive (AR) Models
#'
#' @description
#' Fits an autoregressive time series model of order \eqn{p} using the Online Newton
#' Step (ONS) optimization algorithm. It supports multiple loss functions and uses
#' a learning rate schedule with cosine annealing.
#'
#' @param y A numeric vector or univariate time series (\code{ts}).
#' @param p A positive integer indicating the autoregressive order \eqn{p}.
#' @param phi_0 An optional numeric vector of length \code{p} with the initial
#'   values for the AR coefficients. If \code{NULL} (default), coefficients are
#'   initialized to \code{0.1}.
#' @param loss_function A character string specifying the loss function to minimize.
#'   Must be one of \code{"quadratic"} (default), \code{"absolute"},
#'   \code{"logcosh"}, or \code{"huber"}.
#'
#' @return An S7 object of class \code{AR.ONS} containing:
#' \itemize{
#'   \item \code{y}: The input numeric time series.
#'   \item \code{p}: The autoregressive order.
#'   \item \code{phi_0}: Initial coefficients vector.
#'   \item \code{loss_function}: The loss function used.
#'   \item \code{phi}: A numeric vector of length \code{p} with the estimated autoregressive parameters (averaged over the series).
#'   \item \code{history}: An \eqn{n \times p} numeric matrix containing the coefficient trajectory at the final epoch.
#'   \item \code{fitted}: A numeric vector of fitted values (\eqn{\hat{y}_t}).
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(model = list(ar = c(0.6, -0.3)), n = 500)
#' fit <- AR.ONS(y = as.numeric(y), p = 2, loss_function = "quadratic")
#' fit@phi

AR.ONS <- function(y, p,
                   phi_0 = NULL,
                   loss_function = c("quadratic", "absolute", "logcosh", "huber")) {
  if (length(y) <= p) {
    stop("The length of time series @y must be greater than order @p.")
  }
  epochs  <- 30
  eps     <- 1e-4
  delta   <- 1.345
  beta    <- 0.99
  eta     <- 0.1
  eta_min <- 0.001
  loss_function <- match.arg(loss_function)
  n <- length(y)

  phi_prev <- if (is.null(phi_0)) matrix(0.1, nrow = p, ncol = 1)
  else matrix(phi_0, nrow = p, ncol = 1)
  phi_hist <- matrix(NA, nrow = n, ncol = p)
  y_hat <- rep(NA, n)
  A <- diag(eps, p)

  for (e in 1:epochs) {

    fraccion_avance <- (e - 1) / max(1, (epochs - 1))
    eta_e <- eta_min + 0.5 * (eta - eta_min) * (1 + cos(fraccion_avance * pi))

    for (t in (p + 1):n) {
      x_t <- matrix(y[(t - 1):(t - p)], ncol = 1)
      y_hat_t <- as.numeric(t(phi_prev) %*% x_t)
      err <- y[t] - y_hat_t

      if (e == epochs) {
        y_hat[t] <- y_hat_t
      }

      if (loss_function == "quadratic") {
        grad_t <- -x_t * err
      } else if (loss_function == "logcosh") {
        grad_t <- -x_t * tanh(err)
      } else if (loss_function == "absolute") {
        grad_t <- -x_t * sign(err)
      } else if (loss_function == "huber") {
        grad_t <- if (abs(err) <= delta) -x_t * err else -delta * x_t * sign(err)
      }

      A <- beta * A + (1 - beta) * (grad_t %*% t(grad_t))

      tryCatch({
        phi_new <- phi_prev - eta_e * solve(A + diag(1e-8, p), grad_t)
        phi_prev <- phi_new
      }, error = function(e) { })

      phi_prev <- pmin(pmax(phi_prev, -0.999), 0.999)

      if (e == epochs) {
        phi_hist[t, ] <- as.vector(phi_prev)
      }
    }
  }

  phi_final <- colMeans(phi_hist[(p + 1):n, , drop = FALSE], na.rm = TRUE)

  fitted_vals <- rep(NA, n)
  for (t in (p + 1):n) {
    x_t <- y[(t - 1):(t - p)]
    fitted_vals[t] <- sum(phi_final * x_t)
  }

  AR_ONS_class(
    y             = as.numeric(y),
    p             = as.integer(p),
    loss_function = loss_function,
    phi_0         = phi_0,
    phi           = as.vector(phi_final),
    history       = phi_hist,
    fitted        = fitted_vals
  )
}
