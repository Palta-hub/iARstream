# TEST DE VALIDACION DE ENTRADAS Y ERRORES --------------------------------

test_that("AR.ONS y su clase validan entradas incorrectas", {
  set.seed(123)
  sim_y <- stats::rnorm(50)
  expect_error(AR.ONS(y = sim_y, p = 0))
  expect_error(AR.ONS(y = sim_y, p = 2, loss_function = "inexistente"))
  expect_error(AR_ONS_class(y = 1, p = 2))
})

# TEST DE ESTRUCTURA Y DIMENSIONES DE SALIDA ------------------------------

test_that("AR.ONS retorna las dimensiones y tipos correctos", {
  set.seed(123)
  n_obs <- 100
  sim_y <- stats::rnorm(n_obs)
  p_order <- 2
  fit <- AR.ONS(y = sim_y, p = p_order, loss_function = "quadratic")
  fit <- est(fit)

  expect_true(any(grepl("AR.ONS", class(fit))))
  expect_length(fit@phi, p_order)
  expect_type(fit@phi, "double")

  expect_true(is.matrix(fit@history))
  expect_equal(nrow(fit@history), n_obs)
  expect_equal(ncol(fit@history), p_order)
  expect_length(fit@fitted, n_obs)

  expect_true(all(is.na(fit@history[1:p_order, ])))
  expect_true(all(is.na(fit@fitted[1:p_order])))
  expect_false(any(is.na(fit@history[(p_order + 1):n_obs, ])))
  expect_false(any(is.na(fit@fitted[(p_order + 1):n_obs])))
})

# TEST DE TODAS LAS FUNCIONES DE PERDIDA Y PHI_0 --------------------------

test_that("AR.ONS funciona con todas las perdidas y admite phi_0 automatico u opcional", {
  set.seed(123)
  n_obs <- 100
  sim_y <- stats::rt(n_obs, df = 3)

  losses <- c("quadratic", "absolute", "logcosh", "huber")
  for (l in losses) {
    expect_no_error(
      fit <- AR.ONS(y = sim_y, p = 1, loss_function = l)
    )
  }

  fit_auto <- AR.ONS(y = sim_y, p = 2, phi_0 = NULL)
  expect_true(is.null(fit_auto@phi_0) || is.numeric(fit_auto@phi_0))

  custom_phi0 <- c(0.5, 0.2)
  fit_custom <- AR.ONS(y = sim_y, p = 2, phi_0 = custom_phi0)

  expect_equal(fit_custom@phi_0, custom_phi0)
})

# TEST DE METODOS DE LA CLASE ---------------------------------------------

test_that("Los metodos print y est funcionan correctamente", {
  set.seed(123)
  sim_y <- stats::rnorm(100)
  obj <- AR_ONS_class(y = sim_y, p = 2, loss_function = "absolute")

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

# TEST DE CONVERGENCIA (LOGICA DE ESTIMACION) -----------------------------

test_that("AR.ONS se aproxima a los parametros reales", {
  set.seed(123)
  n <- 1000
  phi_true <- c(0.5, -0.3)
  y <- numeric(n)
  y[1:2] <- stats::rnorm(2)

  for (t in 3:n) {
    y[t] <- phi_true[1] * y[t-1] + phi_true[2] * y[t-2] + stats::rnorm(1, sd = 0.5)
  }
  fit <- AR.ONS(y = y, p = 2, loss_function = "quadratic")
  expect_equal(fit@phi, phi_true, tolerance = 0.3)
})
