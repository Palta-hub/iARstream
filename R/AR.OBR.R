# DEFINICION CLASE --------------------------------------------------------

AR_OBR_class <- S7::new_class(
  "AR.OBR",
  parent = Bayesian,
  package = "iARstream",
  properties = list(

    y = S7::new_property(S7::class_numeric),
    p = S7::new_property(S7::class_numeric),

    phi = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted  = S7::new_property(S7::class_any, default = NULL)),

  validator = function(self) {
    if (length(self@p) != 1 || self@p < 1) {
      return("Property @p must be a positive integer indicating the model order.")
    }
    if (length(self@y) <= self@p) {
      return("The length of time series @y must be greater than order @p.")
    }
    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, AR_OBR_class) <- function(x, ...) {
  header    <- sprintf("Online Bayesian Regression applied to AR(%s)", x@p)
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

# FUNCION AR.OBR ----------------------------------------------------------

#' Online Bayesian Regression for AR(p) Models
#'
#' Estimates the parameters of an Autoregressive AR(p) model using an Online Bayesian
#' Regression (OBR) approach. The method performs sequential observation-by-observation
#' updates of the parameter posterior distribution, leveraging mathematical conjugacy
#' under a Gaussian likelihood framework to provide exact, closed-form sequential inference.
#'
#' @details
#' The algorithm initializes the AR parameters with a zero-mean isotropic Gaussian prior.
#' The prior precision matrix is initialized as \eqn{S_0^{-1} = \alpha I_p}, where
#' the prior precision hyperparameter is fixed at \code{alpha = 2.0}. The observation
#' noise precision is assumed fixed at \code{beta = 1.0}.
#'
#' At each time step \eqn{t}, the algorithm performs sequential learning (where the
#' posterior of time \eqn{t-1} becomes the prior for time \eqn{t}) by analytically
#' updating the precision matrix and the mean vector:
#' \itemize{
#'   \item Precision update: \eqn{S_{t}^{-1} = S_{t-1}^{-1} + \beta (x_t x_t^T)}
#'   \item Mean update: \eqn{m_t = S_t (S_{t-1}^{-1} m_{t-1} + \beta x_t y_t)}
#' }
#' This recursive formulation allows the model to continuously adapt to streaming data
#' without needing to store past observations in memory.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param p An integer specifying the order of the AR(p) model.
#'
#' @return An object of class \code{AR_OBR_class} containing the following slots:
#' \item{y}{The original input time series.}
#' \item{p}{The specified AR order.}
#' \item{phi}{A numeric vector with the final estimated AR parameters.}
#' \item{history}{A numeric matrix of dimension \code{n x p} tracking the sequential posterior mean trajectories.}
#' \item{fitted}{A numeric vector containing the in-sample fitted values.}
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(model = list(ar = c(0.4, -0.6)), n = 100)
#' fit <- AR.OBR(y = as.numeric(y), p = 2)
#' fit
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
