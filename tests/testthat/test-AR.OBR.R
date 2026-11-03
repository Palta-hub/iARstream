# TEST DE VALIDACION DE ENTRADAS Y ERRORES --------------------------------

test_that("AR.OBR y su clase validan entradas incorrectas", {
  set.seed(123)
  sim_y <- stats::rnorm(50)
  expect_error(AR.OBR(y = sim_y, p = 0))
  expect_error(AR.OBR(y = sim_y[1:2], p = 2))
  expect_error(AR_OBR_class(y = sim_y[1:2], p = 2))
  expect_error(AR_OBR_class(y = sim_y, p = 0))
})

# TEST DE ESTRUCTURA Y DIMENSIONES DE SALIDA ------------------------------

test_that("AR.OBR retorna las dimensiones, tipos y slots correctos", {
  set.seed(123)
  n_obs <- 100
  sim_y <- stats::rnorm(n_obs)
  p_order <- 2
  fit <- AR.OBR(y = sim_y, p = p_order)

  expect_true(any(grepl("AR.OBR", class(fit))))
  expect_length(fit@phi, p_order)
  expect_type(fit@phi, "double")

  expect_true(is.matrix(fit@history))
  expect_equal(nrow(fit@history), n_obs)
  expect_equal(ncol(fit@history), p_order)
  expect_length(fit@fitted, n_obs)
  expect_type(fit@fitted, "double")

  expect_true(all(is.na(fit@history[1:p_order, ])))
  expect_true(all(is.na(fit@fitted[1:p_order])))
  expect_false(any(is.na(fit@history[(p_order + 1):n_obs, ])))
  expect_false(any(is.na(fit@fitted[(p_order + 1):n_obs])))
})

# TEST DE METODOS DE LA CLASE ---------------------------------------------

test_that("Los metodos print y est de AR_OBR_class funcionan correctamente", {
  set.seed(123)
  sim_y <- stats::rnorm(100)
  obj <- AR_OBR_class(y = sim_y, p = 2)
  expect_null(obj@phi)

  out_unest <- capture.output(print(obj))
  expect_true(any(grepl("Model not estimated yet", out_unest)))

  expect_no_error(
    obj_est <- est(obj)
  )
  expect_false(is.null(obj_est@phi))
  out_est <- capture.output(print(obj_est))
  expect_false(any(grepl("Model not estimated yet", out_est)))
  expect_true(any(grepl("Phi estimated =", out_est)))
})

# TEST DE CONVERGENCIA Y ESTIMACION ANALITICA ----------------------------

test_that("AR.OBR se aproxima con precision a los parametros reales", {
  set.seed(123)
  n <- 800
  phi_true <- c(0.4, -0.6)
  y <- numeric(n)
  y[1:2] <- stats::rnorm(2)

  for (t in 3:n) {
    y[t] <- phi_true[1] * y[t-1] + phi_true[2] * y[t-2] + stats::rnorm(1, sd = 0.5)
  }
  fit <- AR.OBR(y = y, p = 2)
  expect_equal(fit@phi, phi_true, tolerance = 0.15)
})
