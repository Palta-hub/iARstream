# TEST DE VALIDACION DE ENTRADAS Y ERRORES --------------------------------

test_that("AR.OGD y su clase validan entradas incorrectas", {
  set.seed(123)
  sim_y <- stats::rnorm(50)
  expect_error(AR.OGD(y = sim_y, p = 0))
  expect_error(AR.OGD(y = sim_y, p = 2, loss_function = "inventada"))
  expect_error(AR_OGD_class(y = 1, p = 2))
  expect_error(AR_OGD_class(y = sim_y, p = 2, loss_function = "inventada"))
})

# TEST DE ESTRUCTURA Y DIMENSIONES DE SALIDA ------------------------------

test_that("AR.OGD retorna las dimensiones y tipos correctos", {
  set.seed(123)
  n_obs <- 100
  sim_y <- stats::rnorm(n_obs)
  p_order <- 2
  fit <- AR.OGD(y = sim_y, p = p_order, loss_function = "quadratic")
  fit <- est(fit)

  expect_true(any(grepl("AR.OGD", class(fit))))
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

test_that("AR.OGD funciona con todas las perdidas y admite phi0", {
  set.seed(123)
  n_obs <- 100
  sim_y <- stats::rt(n_obs, df = 3)

  losses <- c("quadratic", "absolute", "logcosh", "huber")
  for (l in losses) {
    expect_no_error(
      fit <- AR.OGD(y = sim_y, p = 1, loss_function = l)
    )
  }

  custom_phi0 <- c(0.5, 0.2)
  fit_custom <- AR.OGD(y = sim_y, p = 2, phi0 = custom_phi0)

  expect_equal(fit_custom@phi0, custom_phi0)
})

# TEST DE METODOS DE LA CLASE ---------------------------------------------

test_that("Los métodos est y print funcionan correctamente", {
  set.seed(123)
  sim_y <- stats::rnorm(100)
  obj <- AR_OGD_class(y = sim_y, p = 2, loss_function = "absolute")

  expect_null(obj@phi)
  out_unestimated <- capture.output(print(obj))
  expect_true(any(grepl("Model not estimated yet", out_unestimated)))

  expect_no_error(
    obj_est <- est(obj)
  )
  expect_false(is.null(obj_est@phi))
  out_estimated <- capture.output(print(obj_est))
  expect_false(any(grepl("@ phi\\s*:\\s*NULL", out_estimated)))
})

# TEST DE CONVERGENCIA (LOGICA DE ESTIMACION) -----------------------------

test_that("AR.OGD se aproxima a los parametros reales", {
  set.seed(42)
  n <- 1000
  phi_true <- c(0.5, -0.3)
  y <- numeric(n)
  y[1:2] <- stats::rnorm(2)

  for (t in 3:n) {
    y[t] <- phi_true[1] * y[t-1] + phi_true[2] * y[t-2] + stats::rnorm(1, sd = 0.5)
  }
  fit <- AR.OGD(y = y, p = 2, loss_function = "quadratic")
  expect_equal(fit@phi, phi_true, tolerance = 0.25)
})
