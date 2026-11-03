# Definición clases

Bayesian <- S7::new_class("Bayesian", package = "iARstream")
Classical <- S7::new_class("Classical",package = "iARstream")

# Definición principal
est <- S7::new_generic("est", "x")

#' @importFrom stats cor var dt rnorm runif
#' @import S7
NULL
