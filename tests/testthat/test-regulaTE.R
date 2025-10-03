test_that("regulaTE runs without error on synthetic example data", {
  rawdata <- data.frame(
    X = c(1, 2, 3, 4),
    probTreat = c(0.1, 0.5, 0.5, 0.9),
    probX = c(0.25, 0.25, 0.25, 0.25),
    Y0mean = c(0, 0, 0, 0),
    tauX = c(0, 5, 5, 0)
  )

  sigma <- 5

  set.seed(1)
  sim_size <- 500
  simdata <- rawdata[sample(seq_len(nrow(rawdata)), size = sim_size, replace = TRUE, prob = rawdata$probX), ]
  simdata$x <- as.integer(runif(sim_size) < simdata$probTreat)

  Xvalues <- unique(simdata$X)
  for (i in seq_along(Xvalues)) {
    dummy <- as.integer(simdata$X == Xvalues[i])
    if (sum(dummy) > 0) {
      simdata[[paste0("Z", i)]] <- dummy
    }
  }
  Z_uc <- as.matrix(simdata[, grep("^Z", names(simdata))])
  Z_uc <- Z_uc[, -ncol(Z_uc)]  # drop last column to avoid collinearity

  simdata$Y1mean <- simdata$Y0mean + simdata$tauX
  C0 <- sd(simdata$tauX)
  simdata$Y1 <- rnorm(sim_size, mean = simdata$Y1mean, sd = sigma)
  simdata$Y0 <- rnorm(sim_size, mean = simdata$Y0mean, sd = sigma)
  simdata$Y <- simdata$x * simdata$Y1 + (1 - simdata$x) * simdata$Y0

  expect_error(
      regulaTE(
        outcome = simdata$Y,
        treatment = simdata$x,
        covariates = Z_uc,
        confounders = NULL,
        parameter = "ATE",
        C = C0,
        se = "hom",
        cluster = NULL,
        sig_level = 0.05,
        df_corr = TRUE,
        digits = 4,
        trimmed_data = NULL
      ),
      regexp = NA  
    )
})
