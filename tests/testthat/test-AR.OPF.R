# TEST DE VALIDACION DE ENTRADAS Y ERRORES --------------------------------

test_that("AR.OPF y su clase validan entradas incorrectas", {
  set.seed(123)
  sim_y <- stats::rnorm(50)
  expect_error(AR.OPF(y = sim_y, p = 0))
  expect_error(AR.OPF(y = sim_y, p = 2, error_distribution = "inexistente"))
  expect_error(AR_OPF_class(y = 1, p = 2))
  expect_error(AR_OPF_class(y = sim_y, p = 2, error_distribution = "invalida"))
})

# TEST DE ESTRUCTURA Y DIMENSIONES DE SALIDA ------------------------------

test_that("AR.OPF retorna las dimensiones, tipos y slots correctos", {
  set.seed(123)
  n_obs <- 100
  sim_y <- stats::rnorm(n_obs)
  p_order <- 2
  fit <- AR.OPF(y = sim_y, p = p_order, error_distribution = "gaussian")

  expect_true(any(grepl("AR.OPF", class(fit))))
  expect_length(fit@phi, p_order)
  expect_type(fit@phi, "double")

  expect_true(is.matrix(fit@history))
  expect_equal(nrow(fit@history), n_obs)
  expect_equal(ncol(fit@history), p_order)
  expect_length(fit@fitted, n_obs)
  expect_length(fit@sigma_est, n_obs)

  min_win <- max(p_order + 1, 15)
  expect_true(all(is.na(fit@history[1:min_win, ])))
  expect_false(any(is.na(fit@history[(min_win + 1):n_obs, ])))
})

# TEST DE DISTRIBUCIONES DE ERROR ----------------------------------------

test_that("AR.OPF funciona con todas las distribuciones de error soportadas", {
  set.seed(123)
  n_obs <- 100
  sim_y <- stats::rt(n_obs, df = 3)

  distributions <- c("gaussian", "laplace", "student")
  for (d in distributions) {
    expect_no_error(
      fit <- AR.OPF(y = sim_y, p = 1, error_distribution = d)
    )
  }
})

# TEST DE METODOS DE LA CLASE ---------------------------------------------

test_that("Los metodos print y est de AR_OPF_class funcionan correctamente", {
  set.seed(123)
  sim_y <- stats::rnorm(100)
  obj <- AR_OPF_class(y = sim_y, p = 2, error_distribution = "laplace")

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

# TEST DE CONVERGENCIA (LOGICA DE ESTIMACION) -----------------------------

test_that("AR.OPF se aproxima a los parametros reales", {
  set.seed(123)
  n <- 1000
  phi_true <- c(0.5, -0.3)
  y <- numeric(n)
  y[1:2] <- stats::rnorm(2)

  for (t in 3:n) {
    y[t] <- phi_true[1] * y[t-1] + phi_true[2] * y[t-2] + stats::rnorm(1, sd = 0.5)
  }
  fit <- AR.OPF(y = y, p = 2, error_distribution = "gaussian")
  expect_equal(fit@phi, phi_true, tolerance = 0.35)
})
