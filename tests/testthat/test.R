library("ridge")
library("datasets")

tol <- 0.0001

context("Basic tests")

test_that("test linearRidge agrees with lm when lambda = 0", {
  model1 <- lm(mpg ~ wt + cyl, data = mtcars)
  model2 <- linearRidge(mpg ~ wt + cyl, data = mtcars, lambda = 0)
  
  expect_equal(coef(model1), coef(model2), tolerance = tol, label = "coefficients")
  expect_equal(predict(model1), predict(model2), tolerance = tol, label = "predictions")
})

test_that("test linearRidge with formula variable that goes out of scope", {
  the_formula <- formula('mpg ~ wt + cyl')
  model1 <- lm(the_formula, data = mtcars)
  model2 <- linearRidge(the_formula, data = mtcars, lambda = 0)
  # suppose the_formula goes away:
  rm(the_formula)
  
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

test_that("test linearRidge predict method for lambda = 0, with a factor and newdata", {
  foo <- mtcars
  foo$cyl_factor <- as.factor(paste0("cyl", foo$cyl))
  model1 <- lm(mpg ~ wt + cyl_factor, data = foo)
  model2 <- linearRidge(mpg ~ wt + cyl_factor, data = foo, lambda = 0.0)

  # predict agrees between model1 and model2 when lambda is 0, with factor
  newdata <- data.frame(wt=c(1.0), cyl_factor=c("cyl4"))
  preds1 <- predict(model1, newdata)
  preds2 <- predict(model2, newdata)
  expect_equal(preds1, preds2, tolerance = tol, label = "predictions")
})


context("Simple run tests for different datasets")

test_that("Wrong argument -  Hair is not numeric or logical", {
  data(HairEyeColor)
  expect_warning(linearRidge(Hair ~ ., data = as.data.frame(HairEyeColor)))
})

test_that("Runs model + predict with HairEyeColor dataset", {
  
  data(HairEyeColor)
  model <- linearRidge(Freq ~ ., data = as.data.frame(HairEyeColor))
  pred <- predict(model, as.data.frame(HairEyeColor))
  expect_more_than(pred[1],18)
})

test_that("Runs model + predict with HairEyeColor dataset - version with different formula", {
  
  data(HairEyeColor)
  model <- linearRidge(Freq ~ Eye, data = as.data.frame(HairEyeColor))
  pred <- predict(model, as.data.frame(HairEyeColor))
  expect_more_than(pred[1],22)
})


test_that("Runs model + predict with GenBin dataset", {
  
  data(GenBin)
  model <- logisticRidge(Phenotypes ~ ., data = as.data.frame(GenBin))
  pred <- predict(model, as.data.frame(GenBin))
  expect_less_than(pred[1],0)
})


test_that("Runs model + predict with Hald dataset", {
  
  data(Hald)
  model <- linearRidge(y ~ ., data = as.data.frame(Hald))
  pred <- predict(model, as.data.frame(Hald))
  expect_more_than(pred[1],70)
})

test_that("Runs model + predict with Hald dataset", {
  
  data(Hald)
  model <- linearRidge(X1 ~ ., data = as.data.frame(Hald), scaling="none")
  pred <- predict(model, as.data.frame(Hald))
  expect_more_than(pred[1],7)
})

test_that("Runs model + predict with Hald dataset", {
  
  data(Hald)
  model <- linearRidge(y ~ X1 + X2 + X3, data = as.data.frame(Hald), scaling="none")
  pred <- predict(model, as.data.frame(Hald))
  expect_more_than(pred[1],78)
})

test_that("Runs model + predict with Hald dataset", {
  
  data(Hald)
  model <- linearRidge(y ~ X1 + X2 + X3, data = as.data.frame(Hald), lambda = 0.01, scaling="none")
  pred <- predict(model, as.data.frame(Hald))
  expect_more_than(pred[1],78)
})


test_that("Runs model + predict with Gorman dataset", {
  
  data(Gorman)
  model <- linearRidge(logY ~ ., data = as.data.frame(Gorman))
  pred <- predict(model, as.data.frame(Gorman))
  expect_more_than(pred[1],2)
})

test_that("Runs model + predict with iris dataset", {
  
  data(iris)
  model <- linearRidge(Sepal.Length ~ ., data = as.data.frame(iris))
  pred <- predict(model, as.data.frame(iris))
  expect_more_than(pred[1],4)
})


test_that("Runs model + predict with ToothGrowth dataset", {
  
  data(ToothGrowth)
  model <- linearRidge(len ~ ., data = as.data.frame(ToothGrowth))
  pred <- predict(model, as.data.frame(ToothGrowth))
  expect_more_than(pred[1],10)
})
                          




