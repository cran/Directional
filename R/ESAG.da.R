ESAG.da <- function(y, ina, fraction = 0.2, R = 100, seed = FALSE) {
  runtime <- proc.time()
  ina <- as.numeric(ina)
  n <- dim(y)[1]  
  frac <- round(fraction * n)
  g <- max(ina)
  mesi <- matrix(nrow = g, ncol = 3)
  if (seed)  set.seed(1234567)
  mat <- matrix(0, frac, g)
  per <- numeric(R)
  for (i in 1:R) {
    nu <- sample(1:n, frac)
    id <- ina[-nu]
    ytrain <- y[-nu, ]
    ytest <- y[nu, ]
    for (j in 1:g) {
      mod <- Directional::ESAGmle( ytrain[id == j, ] )
      mat[, j] <- Directional::ESAGdensity(ytest, c(mod$mu, mod$gam), logden = TRUE )
    }
    est <- Rfast::rowMaxs(mat)
    per[i] <- sum(est == ina[nu])/frac
  }
  percent <- mean(per)
  s1 <- sd(per)
  s2 <- sqrt(percent * (1 - percent)/R)
  conf1 <- c(percent - 1.96 * s1, percent + 1.96 * s1)
  conf2 <- c(percent - 1.96 * s2, percent + 1.96 * s2)
  if (conf1[2] > 1)  conf1[2] <- 1
  if (conf1[1] < 0)  conf1[1] <- 0
  if (conf2[2] > 1)  conf2[2] <- 1
  if (conf2[1] < 0)  conf2[1] <- 0
  conf3 <- quantile(per, probs = c(0.025, 0.975))
  ci <- rbind(conf1, conf2, conf3)
  runtime <- proc.time() - runtime
  colnames(ci) <- c("2.5%", "97.5%")
  rownames(ci) <- c("standard", "binomial", "empirical")
  percent <- c(percent, s1, s2)
  names(percent) <- c("percent", "sd1", "sd2")
  list(percent = percent, ci = ci, runtime = runtime)
}


