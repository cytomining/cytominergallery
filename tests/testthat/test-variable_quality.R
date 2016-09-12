context("variable_quality")

test_that("`variable_quality` measure replicate reproducibility of variables", {

  doParallel::registerDoParallel()

  x1 <- rnorm(50)
  x2 <- x1 + rnorm(50)/10
  y1 <- rnorm(50)
  y2 <- y1 + rnorm(50)

  sample <- data.frame(x=c(x1, x2), y=c(y1, y2), class=rep(1:50, times = 2))

  variable_quality_df <-
  variable_quality(sample = sample,
                   variables = c("x", "y"),
                   strata = c("class"),
                   replicates = 2)

  expect_equal(
    variable_quality_df %>%
      dplyr::filter(variable == "x") %>%
      magrittr::extract2("median"),
    cor(x1, x2)
  )

  expect_equal(
    variable_quality_df %>%
      dplyr::filter(variable == "y") %>%
      magrittr::extract2("median"),
    cor(y1, y2)
  )

})
