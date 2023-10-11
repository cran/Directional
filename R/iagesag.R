### IAG versus ESAG
iagesag <- function(x, B = 1, tol = 1e-07) {
  mod <- Directional::esag.mle(x, tol = tol)
  stat <- 2 * mod$loglik - 2 * mod$iag.loglik

  if (B == 1) {
    p.value <- pchisq(stat, 2, lower.tail = FALSE)
    parameter <- 2     ;   names(parameter) <- "df"
    statistic <- stat  ;   names(statistic) <- "Test statistic"
    alternative <- "ESAG is prefered to IAG"
    method <- "Asymptotic rotational symmetry (IAG versus ESAG distribution)"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"

  } else {
    tb <- numeric(B)
    n <- dim(x)[1]
    for (i in 1:B) {
      nu <- Rfast2::Sample.int(n, n, replace = TRUE)
	    mod <- Directional::esag.mle(x[nu, ], tol = tol)
      tb[i] <- 2 * mod$loglik - 2 * mod$iag.loglik
    }
    p.value <- ( sum(tb > stat) + 1) / (B + 1)
    parameter <- "NA"     ;   names(parameter) <- "df"
    statistic <- stat  ;   names(statistic) <- "Test statistic"
    alternative <- "ESAG is prefered to IAG"
    method <- "Bootstrap rotational symmetry (IAG versus ESAG distribution)"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"
  }

  result
}
