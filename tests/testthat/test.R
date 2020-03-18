library("ridge")

tol <- 0.0001

context("Basic tests")

test_that("test linearRidge agrees with lm when lambda = 0", {
  model1 <- lm(mpg ~ wt + cyl, data = mtcars)
  model2 <- linearRidge(mpg ~ wt + cyl, data = mtcars, lambda = 0)
  
  expect_equal(coef(model1), coef(model2), tolerance = tol, label = "coefficients")
  expect_equal(predict(model1), predict(model2), tolerance = tol, label = "predictions")
})

test_that("test linearRidge near lm when lambda = 0.01", {
  model1 <- lm(mpg ~ wt + cyl, data = mtcars)
  model2 <- linearRidge(mpg ~ wt + cyl, data = mtcars, lambda = 0.01)
  
  # at least one coefficient is greater than 0.01 different:
  expect(any(abs(coef(model1) - coef(model2)) > 0.01),
         "coefficients agree too strongly")
  # all coefficients are less than 0.11 different:
  expect(all(abs(coef(model1) - coef(model2)) < 0.11),
         "coefficients too different")
})

test_that("test linearRidge vcov method for lambda = 0, two terms", {
  model1 <- lm(mpg ~ wt + cyl, data = mtcars)
  model2 <- linearRidge(mpg ~ wt + cyl, data = mtcars, lambda = 0.0)
  
  vc1 <- vcov(model1)
  vc2 <- vcov(model2)
  # vcov agrees between model1 and model2 when lambda is 0
  expect_equal(as.vector(vc1), as.vector(vc2), tolerance = tol, label = "vcov lambda 0")
})

test_that("test linearRidge vcov method for lambda = 0, three terms", {
  model1 <- lm(mpg ~ wt + cyl + disp, data = mtcars)
  model2 <- linearRidge(mpg ~ wt + cyl + disp, data = mtcars, lambda = 0.0)
  
  vc1 <- vcov(model1)
  vc2 <- vcov(model2)
  # vcov agrees between model1 and model2 when lambda is 0
  expect_equal(as.vector(vc1), as.vector(vc2), tolerance = tol, label = "vcov lambda 0")
})

test_that("test linearRidge vcov method for lambda = 0.01", {
  model1 <- lm(mpg ~ wt + cyl, data = mtcars)
  model2 <- linearRidge(mpg ~ wt + cyl, data = mtcars, lambda = 0.01)
  
  vc1 <- vcov(model1)
  vc2 <- vcov(model2)
  # at least one coefficient is greater than 0.01 different:
  expect(any(abs(vc1 - vc2) > 0.01), "coefficients agree too strongly")
  # all coefficients are less than 0.013 different:
  expect(all(abs(vc1 - vc2) < 0.013), "coefficients too different")
})

test_that("test linearRidge vcov method for lambda = 0, with a factor", {
  foo <- mtcars
  foo$cyl_factor <- as.factor(paste0("cyl", foo$cyl))
  model1 <- lm(mpg ~ wt + cyl_factor, data = foo)
  model2 <- linearRidge(mpg ~ wt + cyl_factor, data = foo, lambda = 0.0)

  # coef agrees between model1 and model2 when lambda is 0, with factor
  expect_equal(coef(model1), coef(model2), tolerance = tol, label = "coefficients")
  # predict agrees between model1 and model2 when lambda is 0, with factor
  expect_equal(predict(model1), predict(model2), tolerance = tol, label = "predictions")

  # vcov agrees between model1 and model2 when lambda is 0, with factor
  expect_equal(as.vector(vcov(model1)), as.vector(vcov(model2)),
               tolerance = tol, label = "vcov lambda 0, factors")
})
