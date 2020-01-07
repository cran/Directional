spml.fbed <- function(y, x, alpha = 0.05, K = 0, 
    backward = FALSE, parallel = FALSE, tol = 1e-07, maxiters = 100) {
	Rfast2::fbed.reg(y, x, alpha = alpha, type = "spml", K = K, 
    backward = backward, parallel = parallel, tol = tol, maxiters = maxiters)
}   
