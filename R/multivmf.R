#[export]
multivmf <- function(x, ina, tol = 1e-07, ell = FALSE) {
   Rfast::multivmf.mle(x = x, ina = ina, tol = tol, ell = ell)
}



















