vmf2test <- function(y1, y2, B = 1) {
  mod0 <- vmf2(y1, y2)
  lik0 <- mod0$loglik
  mod1 <- Directional::vmf.mle(y1, fast = TRUE)
  mod2 <- Directional::vmf.mle(y2, fast = TRUE)
  lik1 <- mod1$loglik + mod2$loglik
  stat <- 2 * lik1 - 2 * lik0

  if ( B == 1 ) {
    pvalue <- pchisq(stat, dim(y1)[2] - 1, lower.tail = FALSE)

  } else {
    n1 <- dim(y1)[1]  ;  n2 <- dim(y2)[1]
    rot1 <- t( Directional::rotation(mod1$mu, mod0$mu) )
    rot2 <- t( Directional::rotation(mod2$mu, mod0$mu) )
    x1 <- y1 %*% rot1
    x2 <- y2 %*% rot2
    bstat <- numeric(B)
    for (i in 1:B) {
      y1b <- x1[sample(n1, n1, replace = TRUE), ]
      y2b <- x2[sample(n2, n2, replace = TRUE), ]
      lik0 <- vmf2(y1b, y2b)$loglik
      lik1 <- Directional::vmf.mle(y1b, fast = TRUE)$loglik + Directional::vmf.mle(y2b, fast = TRUE)$loglik
      bstat[i] <- lik1 - lik0
    }
    pvalue <- ( sum(bstat >= 0.5 * stat) + 1 ) / (B + 1)
  }

  statistic <- stat  ;   names(statistic) <- "Test statistic"
  parameter <- dim(y1)[2] - 1    ;   names(parameter) <- "df"
  alternative <- "The 2 location parameters differ"
  method <- "Location parameter testing using the von Mises-Fisher distribution"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = pvalue,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
