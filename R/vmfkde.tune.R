vmfkde.tune <- function (x, low = 0.1, up = 1) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    d <- tcrossprod(x)
    diag(d) <- NA
    con <- (2 * pi)^(p/2)
    funa <- function(h) {
        A <- d/h^2
        cpk <- (1/h^2)^(p/2 - 1)/besselI(1/h^2, p/2 - 1)
        f <- colSums(exp(A), na.rm = TRUE)
        n * log(cpk) + sum( log(f) ) 
    }
    a <- optimize(funa, c(low, up), maximum = TRUE)
    res <- c(a$maximum, a$objective/n - log(con) - log(n - 1) )
    names(res) <- c("Optimal h", "cv")
    res
}

