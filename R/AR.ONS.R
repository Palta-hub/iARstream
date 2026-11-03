# DEFINICION CLASE --------------------------------------------------------

AR_ONS_class <- S7::new_class(
  "AR.ONS",
  parent = Classical,
  package = "iARstream",
  properties = list(

    ## Entradas del modelo ##

    y     = S7::new_property(S7::class_numeric),
    p     = S7::new_property(S7::class_numeric),
    phi_0 = S7::new_property(S7::class_any,
                             default = NULL),
    loss_function  = S7::new_property(S7::class_character,
                                      default = "quadratic"),

    ## Salidas del modelo ##

    phi     = S7::new_property(S7::class_any, default = NULL),
    history = S7::new_property(S7::class_any, default = NULL),
    fitted  = S7::new_property(S7::class_any, default = NULL)),

  validator = function(self) {
    if (length(self@p) != 1 || self@p <= 0) {
      return("Property @p must be a positive integer indicating the model order.")
    }
    if (length(self@y) <= self@p) {
      return("The length of time series @y must be greater than order @p.")
    }
    valid_losses <- c("quadratic", "absolute", "logcosh", "huber")
    if (!(self@loss_function %in% valid_losses)) {
      return("Property @loss_function must be 'quadratic', 'absolute', 'logcosh', or 'huber'.")
    }
    return(NULL)
  }
)

# METODO PARA LA CLASE ----------------------------------------------------

method(print, AR_ONS_class) <- function(x, ...) {
  header    <- sprintf("Online Newton Step applied to AR(%s)", x@p)
  line_loss <- paste("Loss Function =", tools::toTitleCase(x@loss_function))
  if (is.null(x@phi)) {
    line_phi  <- "Phi estimated = Model not estimated yet"
  } else {
    line_phi  <- paste("Phi estimated =", paste(round(x@phi, 3), collapse = ", "))
  }

  # Calcular la longitud maxima del contenido y Generar las separaciones

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

method(est, AR_ONS_class) <- function(x) {
  res <- AR.ONS(
    y     = x@y,
    p     = x@p,
    phi_0 = x@phi_0,
    loss_function  = x@loss_function
  )
  x@phi     <- res@phi
  x@history <- res@history
  x@fitted  <- res@fitted
  return(x)
}


# FUNCION AR.ONS ----------------------------------------------------------
#' Online Newton Step for AR(p) Models
#'
#' Estimates the parameters of an Autoregressive AR(p) model using Online Newton Step (ONS)
#' optimization with sequential observation-by-observation updates across multiple epochs.
#' The method incorporates second-order optimization mechanics by dynamically tracking a
#' regularized second-moment matrix alongside quadratic, absolute, Huber, and log-cosh
#' loss functions for robust estimation under non-Gaussian noise. Parameter updates are
#' controlled via a Cosine Annealing learning rate schedule, with box-constraints to maintain
#' stability and final estimates computed via epoch-averaged trajectories.
#'
#' @details
#' The optimization process uses a Cosine Annealing learning rate schedule across 30 epochs,
#' smoothly decaying the learning rate from `eta = 0.1` down to `eta_min = 0.001`.
#' The Online Newton Step algorithm maintains a running second-moment matrix
#' \eqn{A_t = \beta A_{t-1} + (1 - \beta) g_t g_t^T} initialized with \eqn{A_0 = \epsilon I_p}
#' (`eps = 1e-4`, `beta = 0.99`). To ensure numerical stability, the estimated AR parameters
#' are bounded within `[-0.999, 0.999]` at each iteration. For the Huber loss function,
#' the clipping threshold is set to `delta = 1.345`.
#'
#' @param y A numeric vector representing the univariate time series.
#' @param p An integer specifying the order of the AR(p) model.
#' @param loss_function A character string specifying the loss function to minimize.
#'   Options are `"quadratic"` (default), `"absolute"`, `"logcosh"`, or `"huber"`.
#' @param phi_0 An optional numeric vector of length `p` containing initial parameter values.
#'   If `NULL`, defaults to a vector with values 0.1.
#'
#' @return An object of class `AR_ONS_class` containing the following slots:
#' \item{y}{The original input time series.}
#' \item{p}{The specified AR order.}
#' \item{loss_function}{The loss function utilized.}
#' \item{phi0}{The initial parameter values used.}
#' \item{phi}{A numeric vector with the final estimated AR parameters (averaged over the last epoch).}
#' \item{history}{A numeric matrix of dimension `n x p` tracking parameter estimates during the last epoch.}
#' \item{fitted}{A numeric vector containing the in-sample fitted values.}
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(model = list(ar = c(0.6, -0.3)), n = 100)
#' fit <- AR.ONS(y = y, p = 2, loss_function = "quadratic")
#' fit

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
