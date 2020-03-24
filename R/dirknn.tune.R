################################
#### Discriminant analysis for directional data
#### using the k-NN alorithm, tuning the k neighbours
#### Tsagris Michail 01/2016
#### mtsagris@yahoo.gr
################################
dirknn.tune <- function(x, k = 2:10, ina, type = "S", mesos = TRUE, nfolds = 10, folds = NULL,
                        parallel = FALSE, stratified = TRUE, seed = FALSE, rann = FALSE) {
  ## x is the matrix containing the data
  ## nfolds is the number of folds, set to 10 by default
  ## A is the maximum number of neighbours to use
  ## ina indicates the groups, numerical variable
  ## type is either 'S' or 'NS'. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distange of the k nearest
  ## points will be used.
  ## If not, then the harmonic mean will be used. Both of these apply for
  ## the non-standard algorithm, that is when type='NS'
  runtime <- proc.time()
  n <- dim(x)[1]  ## sample size
  ina <- as.numeric(ina) ## makes sure ina is numeric
  if ( is.null(folds) )  folds <- Directional::makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  nfolds <- length(folds)

  per <- matrix( nrow = nfolds, ncol = length(k) )
  for (vim in 1:nfolds) {
    id <- ina[ folds[[ vim ]] ]  ## groups of test sample
    ina2 <- ina[ -folds[[ vim ]] ]   ## groups of training sample
    aba <- as.vector( folds[[ vim ]] )
    aba <- aba[aba > 0]
    g <- Directional::dirknn(x = x[-aba, ], xnew = x[aba, ,drop = FALSE], k = k, ina = ina2,
                             type = type, mesos = mesos, parallel = parallel, rann = rann)
    be <- g - id
    per[vim, ] <- Rfast::colmeans(be == 0)
  }

  ela <- Rfast::colmeans(per)
  runtime <- proc.time() - runtime
  names(ela) <- paste("k=", k, sep = "")
  plot(k, ela, type = "b", xlab = "k nearest neighbours", pch = 9,
       ylab = "Estimated percentage of correct classification")
  percent <- max(ela)
  names(percent) <- c("Estimated percentage")
  list( per = ela, percent = percent, runtime = runtime )
}
