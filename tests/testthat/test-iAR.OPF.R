# TEST DE VALIDACION DE ENTRADAS Y ERRORES --------------------------------

test_that("iAR.OPF y su clase validan entradas incorrectas", {
  set.seed(123)
  sim_times <- c(1, 3, 5, 8, 10)
  sim_y <- stats::rnorm(5)
  expect_error(iAR.OPF(y = sim_y, times = sim_times[-1]))
  expect_error(iAR.OPF(y = sim_y[1], times = sim_times[1]))
  expect_error(iAR.OPF(y = sim_y, times = sim_times, error_distribution = "inexistente"))
  expect_error(iAR_OPF_class(y = sim_y, times = sim_times[-1]))
  expect_error(iAR_OPF_class(y = sim_y, times = sim_times, error_distribution = "invalida"))
})

# TEST DE ESTRUCTURA Y DIMENSIONES DE SALIDA ------------------------------

test_that("iAR.OPF retorna las dimensiones, tipos y slots correctos", {
  set.seed(123)
  n_obs <- 100
  sim_times <- cumsum(sample(1:3, n_obs, replace = TRUE))
  sim_y <- stats::rnorm(n_obs)
  fit <- iAR.OPF(y = sim_y, times = sim_times, error_distribution = "gaussian")

  expect_true(any(grepl("iAR.OPF", class(fit))))
  expect_length(fit@phi, 1)
  expect_type(fit@phi, "double")

  expect_type(fit@history, "double")
  expect_length(fit@history, n_obs)
  expect_type(fit@fitted, "double")
  expect_length(fit@fitted, n_obs)

  expect_false(any(is.na(fit@history[2:n_obs])))
  expect_false(any(is.na(fit@fitted[2:n_obs])))
})

# TEST DE DISTRIBUCIONES DE ERROR ----------------------------------------

test_that("iAR.OPF funciona con todas las distribuciones de error soportadas", {
  set.seed(123)
  n_obs <- 50
  sim_times <- cumsum(sample(1:2, n_obs, replace = TRUE))
  sim_y <- stats::rt(n_obs, df = 4)

  distributions <- c("gaussian", "student")
  for (d in distributions) {
    expect_no_error(
      fit <- iAR.OPF(y = sim_y, times = sim_times, error_distribution = d)
    )
  }
})

# TEST DE METODOS DE LA CLASE ---------------------------------------------

test_that("Los metodos print y est de iAR_OPF_class funcionan correctamente", {
  set.seed(123)
  n_obs <- 50
  sim_times <- cumsum(sample(1:2, n_obs, replace = TRUE))
  sim_y <- stats::rnorm(n_obs)
  obj <- iAR_OPF_class(y = sim_y, times = sim_times, error_distribution = "gaussian")

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

# TEST DE ESTIMACION LOGICA Y CONVERGENCIA --------------------------------

test_that("iAR.OPF se aproxima al parametro real en datos simulados", {
  set.seed(123)
  n <- 600
  times <- cumsum(sample(1:3, n, replace = TRUE))
  y <- numeric(n)
  y[1] <- stats::rnorm(1)
  phi_true <- 0.7

  for (t in 2:n) {
    dt <- times[t] - times[t-1]
    y[t] <- (phi_true^dt) * y[t-1] + stats::rnorm(1, sd = 0.4)
  }
  fit <- iAR.OPF(y = y, times = times, error_distribution = "gaussian")
  expect_true(fit@phi > 0 && fit@phi < 1)
  expect_equal(fit@phi, phi_true, tolerance = 0.25)
})

