# DEFINICION CLASE --------------------------------------------------------

AR_OBR_class <- S7::new_class(
  "AR.OBR",
  parent = Bayesian,
  package = "iARest",
  properties = list(
    # Entradas
    y = S7::new_property(S7::class_numeric),
    p = S7::new_property(S7::class_numeric),

    # Salida
    phi = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted  = S7::new_property(S7::class_any, default = NULL)
  ),
  validator = function(self) {
    # p debe ser un unico número (el rezago), no un vector del tamaño de y
    if (length(self@p) != 1 || self@p < 1) {
      return("Property @p must be a positive integer indicating the model order.")
    }
    # y debe ser lo suficientemente largo para el rezago p
    if (length(self@y) <= self@p) {
      return("The length of time series @y must be greater than order @p.")
    }

    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, AR_OBR_class) <- function(x, ...) {
  # Armar cadenas de texto de cada linea
  header    <- sprintf("Online Bayesian Regression applied to AR(%s) ", x@p)

  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", paste(round(x@phi, 3), collapse = ", "))
  }

  # Calcular la longitud maxima del contenido
  lines_vector <- c(header, line_phi)
  max_width    <- max(nchar(lines_vector))

  # Generar las separaciones con la funcion strrep()
  sep_double <- strrep("=", max_width)
  sep_single <- strrep("-", max_width)

  # Imprimir en consola
  cat(sep_double, "\n")
  cat(header, "\n")
  cat(sep_double, "\n")
  cat(line_phi, "\n")
  cat(sep_single, "\n")

  invisible(x)
}

method(est, AR_OBR_class) <- function(x) {
  res <- AR.OBR(
    y = x@y,
    p = x@p
  )

  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted

  return(x)
}

# FUNCION OBR --------------------------------------------------------------

#' Online Bayesian Regression for Autoregressive (AR) Models
#'
#' @description
#' Fits an autoregressive time series model of order \eqn{p} using an Online
#' Bayesian Regression (OBR) approach. It updates the posterior distribution
#' of the AR coefficients sequentially as new observations arrive.
#'
#' @param y A numeric vector or univariate time series.
#' @param p A positive integer indicating the autoregressive order \eqn{p}.
#'
#' @return An S7 object of class \code{AR.OBR} containing:
#' \itemize{
#'   \item \code{y}: The input numeric time series.
#'   \item \code{p}: The autoregressive order.
#'   \item \code{phi}: A numeric vector of length \code{p} with the final estimated AR parameters.
#'   \item \code{history}: An \eqn{n \times p} numeric matrix containing posterior mean trajectories.
#'   \item \code{fitted}: A numeric vector of fitted values (\eqn{\hat{y}_t}).
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(model = list(ar = c(0.4, -0.6)), n = 500)
#' fit <- AR.OBR(y = as.numeric(y), p = 2)
#' fit@phi

AR.OBR <- function(y, p){

  if (length(y) <= p) {
    stop("The length of time series @y must be greater than order @p.")
  }

  n <- length(y)
  phi_est <- matrix(NA, n, p)

  ## INICIALIZACIÓN (Creación del Prior) ##
  # Asumimos que los pesos parten de 0 y su matriz de precisión inicial es alpha * I
  alpha <- 2.0
  beta <- 1.0
  m <- matrix(0, nrow = p, ncol = 1)
  S_inv <- alpha * diag(p)

  ## AQUI PASA LA MAGIA (Actualización Bayesiana Secuencial) ##
  for (t in (p + 1):n) {

    # Vector de características: los rezagos del modelo AR (vector columna p x 1)
    phi_t <- as.matrix(rev(y[(t - p):(t - 1)]))

    ## ACTUALIZACIÓN DE LA PRECISIÓN (Covarianza inversa posterior) ##
    S_inv_new <- S_inv + beta * (phi_t %*% t(phi_t))

    # Invertimos para obtener la covarianza real actual
    S <- solve(S_inv_new)

    ## ACTUALIZACIÓN DE LA MEDIA (Estimación de los parámetros) ##
    m_new <- S %*% (S_inv %*% m + beta * phi_t * y[t])

    ## ESTIMACIÓN DEL ESTADO ##
    phi_est[t, ] <- as.vector(m_new)

    ## EL POSTERIOR DE HOY ES EL PRIOR DE MAÑANA (Sequential Learning) ##
    m <- m_new
    S_inv <- S_inv_new
  }

  phi_final <- as.numeric(phi_est[n, ])

  # Calcular valores ajustados
  fitted_vals <- rep(NA, n)
  for (t in (p + 1):n) {
    x_t <- y[(t - 1):(t - p)]
    fitted_vals[t] <- sum(phi_final * x_t)
  }

  # Devolver el objeto S7
  AR_OBR_class(
    y       = as.numeric(y),
    p       = as.integer(p),
    phi     = as.vector(phi_final),
    history = phi_est,
    fitted  = as.numeric(fitted_vals)
  )
}
