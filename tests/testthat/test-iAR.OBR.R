# TEST DE VALIDACION DE ENTRADAS Y ERRORES --------------------------------

test_that("iAR.OBR y su clase validan entradas incorrectas", {
  set.seed(123)
  sim_times <- c(1, 3, 5, 8, 10)
  sim_y <- stats::rnorm(5)
  expect_error(iAR.OBR(y = sim_y, times = sim_times[-1]))
  expect_error(iAR.OBR(y = sim_y[1], times = sim_times[1]))
  expect_error(iAR_OBR_class(y = sim_y, times = sim_times[-1]))
  expect_error(iAR_OBR_class(y = matrix(sim_y, ncol = 1), times = sim_times))
})

# TEST DE ESTRUCTURA Y DIMENSIONES DE SALIDA ------------------------------

test_that("iAR.OBR retorna las dimensiones, tipos y slots correctos", {
  set.seed(123)
  n_obs <- 100
  sim_times <- cumsum(sample(1:3, n_obs, replace = TRUE))
  sim_y <- stats::rnorm(n_obs)
  fit <- iAR.OBR(y = sim_y, times = sim_times)

  expect_true(any(grepl("iAR.OBR", class(fit))))
  expect_length(fit@phi, 1)
  expect_type(fit@phi, "double")
  expect_type(fit@history, "double")
  expect_length(fit@history, n_obs)

  expect_false(is.matrix(fit@history))
  expect_length(fit@fitted, n_obs)
  expect_type(fit@fitted, "double")
  expect_true(is.na(fit@history[1]))
  expect_true(is.na(fit@fitted[1]))

  expect_false(any(is.na(fit@history[2:n_obs])))
  expect_false(any(is.na(fit@fitted[2:n_obs])))
})

# TEST DE METODOS DE LA CLASE ---------------------------------------------

test_that("Los metodos print y est de iAR_OBR_class funcionan correctamente", {
  set.seed(123)
  n_obs <- 50
  sim_times <- cumsum(sample(1:2, n_obs, replace = TRUE))
  sim_y <- stats::rnorm(n_obs)
  obj <- iAR_OBR_class(y = sim_y, times = sim_times)

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

# TEST DE CONVERGENCIA Y ESTIMACION BAYESIANA ----------------------------

test_that("iAR.OBR se aproxima con precision al parametro real en series irregulares", {
  set.seed(123)
  n <- 600
  times <- cumsum(sample(1:3, n, replace = TRUE))
  y <- numeric(n)
  y[1] <- stats::rnorm(1)
  phi_true <- 0.65

  for (t in 2:n) {
    dt <- times[t] - times[t-1]
    var_cond <- 1 - phi_true^(2 * dt)
    y[t] <- (phi_true^dt) * y[t-1] + stats::rnorm(1, sd = sqrt(var_cond))
  }
  fit <- iAR.OBR(y = y, times = times)
  expect_true(fit@phi > 0 && fit@phi < 1)
  expect_equal(fit@phi, phi_true, tolerance = 0.1)
})
