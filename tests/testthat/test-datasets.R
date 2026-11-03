test_that("Base de datos IGBC es correcta", {
  expect_true("IGBC" %in% data(package = "iARstream")$results[, "Item"])

  expect_equal(nrow(IGBC), 2361)
  expect_equal(ncol(IGBC), 2)

  expect_named(IGBC, c("DATE", "IGBC"))

  expect_true(is.numeric(IGBC$IGBC))
  expect_true(inherits(IGBC$DATE, "Date") || is.character(IGBC$DATE))
  expect_false(any(is.na(IGBC$IGBC)))
})

test_that("Base de datos IC20 es correcta", {
  expect_true("IC20" %in% data(package = "iARstream")$results[, "Item"])

  expect_equal(nrow(IC20), 770)
  expect_equal(ncol(IC20), 2)

  expect_named(IC20, c("DATE", "IC20"))

  expect_true(is.numeric(IC20$IC20))
  expect_true(inherits(IC20$DATE, "Date") || is.character(IC20$DATE))
  expect_false(any(is.na(IC20$IC20)))
})

test_that("Base de datos ICAP es correcta", {
  expect_true("ICAP" %in% data(package = "iARstream")$results[, "Item"])

  expect_equal(nrow(ICAP), 770)
  expect_equal(ncol(ICAP), 2)

  expect_named(ICAP, c("DATE", "ICAP"))

  expect_true(is.numeric(ICAP$ICAP))
  expect_true(inherits(ICAP$DATE, "Date") || is.character(ICAP$DATE))
  expect_false(any(is.na(ICAP$ICAP)))
})
