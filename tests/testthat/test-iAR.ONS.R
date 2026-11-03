# TEST DE VALIDACION DE ENTRADAS Y ERRORES --------------------------------

test_that("iAR.ONS y su clase validan entradas incorrectas", {
  set.seed(123)
  sim_times <- c(1, 3, 5, 8, 10)
  sim_y <- stats::rnorm(5)
  expect_error(iAR.ONS(y = sim_y, times = sim_times[-1]))
  expect_error(iAR.ONS(y = sim_y, times = sim_times, loss_function = "inexistente"))
})

# TEST DE ESTRUCTURA Y DIMENSIONES DE SALIDA ------------------------------

test_that("iAR.ONS retorna las dimensiones, tipos y clases correctas", {
  set.seed(123)
  n_obs <- 200
  sim_times <- cumsum(sample(1:3, n_obs, replace = TRUE))
  sim_y <- stats::rnorm(n_obs)
  fit <- iAR.ONS(y = sim_y, times = sim_times, loss_function = "quadratic")
  fit <- est(fit)

  expect_true(any(grepl("iAR.ONS", class(fit))))
  expect_length(fit@phi, 1)
  expect_type(fit@phi, "double")

  expect_false(is.null(fit@history))
  expect_true(length(fit@history) > 0)
  expect_false(is.null(fit@fitted))
  expect_true(length(fit@fitted) > 0)

  expect_true(any(!is.na(fit@history)))
  expect_true(any(!is.na(fit@fitted)))
})

# TEST DE FUNCIONES DE PERDIDA Y PHI_0 -------------------------------------

test_that("iAR.ONS funciona con todas las perdidas y admite phi_0 auto u opcional", {
  set.seed(123)
  n_obs <- 50
  sim_times <- cumsum(sample(1:2, n_obs, replace = TRUE))
  sim_y <- stats::rnorm(n_obs)

  losses <- c("quadratic", "absolute", "logcosh", "huber")
  for (l in losses) {
    expect_no_error(
      fit <- iAR.ONS(y = sim_y, times = sim_times, loss_function = l)
    )
  }

  fit_auto <- iAR.ONS(y = sim_y, times = sim_times, phi_0 = NULL)
  expect_true(is.null(fit_auto@phi_0) || is.numeric(fit_auto@phi_0))

  fit_custom <- iAR.ONS(y = sim_y, times = sim_times, phi_0 = 0.8)
  expect_equal(fit_custom@phi_0, 0.8)
})

# TEST DE METODOS DE LA CLASE ---------------------------------------------

test_that("Los metodos print y est de iAR_ONS_class funcionan correctamente", {
  set.seed(123)
  n_obs <- 50
  sim_times <- cumsum(sample(1:2, n_obs, replace = TRUE))
  sim_y <- stats::rnorm(n_obs)
  obj <- iAR_ONS_class(y = sim_y, times = sim_times, loss_function = "absolute")

  expect_null(obj@phi)
  out_unest <- capture.output(print(obj))
  expect_true(any(grepl("Model not estimated yet", out_unest)))

  expect_no_error(
    obj_est <- est(obj)
  )
  expect_false(is.null(obj_est@phi))
  out_est <- capture.output(print(obj_est))
  expect_false(any(grepl("@ phi\\s*:\\s*NULL", out_est)))
})

# TEST DE ESTIMACION LOGICA -----------------------------------------------

test_that("iAR.ONS devuelve valores razonables en el intervalo esperado", {
  set.seed(123)
  n <- 500
  times <- cumsum(sample(1:3, n, replace = TRUE))
  y <- numeric(n)
  y[1] <- stats::rnorm(1)
  phi_true <- 0.7

  for (t in 2:n) {
    dt <- times[t] - times[t-1]
    y[t] <- (phi_true^dt) * y[t-1] + stats::rnorm(1, sd = 0.5)
  }
  fit <- iAR.ONS(y = y, times = times, loss_function = "quadratic")
  fit <- est(fit)
  expect_true(fit@phi > -1 && fit@phi < 1)
})
