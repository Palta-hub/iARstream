# Definición clases

Bayesian <- S7::new_class("Bayesian", package = "iARest")
Classical <- S7::new_class("Classical",package = "iARest")

# Definición principal
est <- S7::new_generic("est", "x")

#' @importFrom stats cor var
#' @import S7
NULL
