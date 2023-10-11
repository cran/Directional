pc.test <- function(x, B = 1, tol = 1e-6) {
  mod <- Directional::sespc.mle(x, tol = tol)
  stat <- 2 * mod$loglik - 2 * mod$sipc.loglik

  if (B == 1) {
    p.value <- pchisq(stat, 2, lower.tail = FALSE)
    parameter <- 2     ;   names(parameter) <- "df"
    statistic <- stat  ;   names(statistic) <- "Test statistic"
    alternative <- "SESPC is prefered to SIPC"
    method <- "Asymptotic rotational symmetry (SIPC versus SESPC distribution)"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"

  } else {
    tb <- numeric(B)
    n <- dim(x)[1]
    for (i in 1:B) {
      nu <- Rfast2::Sample.int(n, n, replace = TRUE)
      mod <- Directional::sespc.mle(x[nu, ], tol = tol)
      tb[i] <- 2 * mod$loglik - 2 * mod$sipc.loglik
    }
    p.value <- ( sum(tb > stat) + 1) / (B + 1)
    parameter <- "NA"     ;   names(parameter) <- "df"
    statistic <- stat  ;   names(statistic) <- "Test statistic"
    alternative <- "SESPC is prefered to SIPC"
    method <- "Bootstrap rotational symmetry (SIPC versus SESPC distribution)"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"
  }

  result
}
