
####################         corpcor vs corshrink         ################################

cor.shrink = function(x, lambda, w, verbose=TRUE)
{
  return ( powcor.shrink(x=x, alpha=1, lambda=lambda, w=w, verbose=verbose) )
}

powcor.shrink = function(x, alpha, lambda, w, verbose=TRUE)
{
  if (missing(alpha)) stop("Please specify the exponent alpha!")

  x = as.matrix(x)

  # matrix power of shrinkage correlation
  powr = pvt.powscor(x=x, alpha=alpha, lambda=lambda, w=w, verbose=verbose)

  return(powr)
}

pvt.powscor = function(x, alpha, lambda, w, verbose)
{
  #### determine correlation shrinkage intensity
  if (missing(lambda))
  {
    lambda = estimate.lambda(x, w, verbose)
    lambda.estimated=TRUE
  }
  else
  {
    if (lambda < 0) lambda = 0
    if (lambda > 1) lambda = 1
    if (verbose)
    {
      cat(paste("Specified shrinkage intensity lambda (correlation matrix):", round(lambda, 4), "\n"))
    }
    lambda.estimated=FALSE
  }
  #####


  n = nrow(x)
  w = pvt.check.w(w, n)
  xs = wt.scale(x, w, center=TRUE, scale=TRUE) # standardize data matrix

  # bias correction factor
  h1 = 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)

  p = ncol(xs)

  if (lambda == 1 | alpha == 0) # result in both cases is the identity matrix
  {
    powr = diag(p)    # return identity matrix
    rownames(powr) = colnames(xs)
    colnames(powr) = colnames(xs)
  }
  else if (alpha == 1) # don't do SVD in this simple case
  {
    # unbiased empirical estimator
    # for w=1/n  the following  would simplify to:  r = 1/(n-1)*crossprod(xs)
    #r0 = h1 * t(xs) %*% diag(w) %*% xs
    #r0 = h1 * t(xs) %*% sweep(xs, 1, w, "*") # sweep requires less memory
    r0 = h1 * crossprod( sweep(xs, 1, sqrt(w), "*") ) # even faster

    # shrink off-diagonal elements
    powr = (1-lambda)*r0
    diag(powr) = 1
  }
  else
  {
    # number of zero-variance variables
    zeros = (attr(xs, "scaled:scale")==0.0)

    svdxs = fast.svd(xs)
    m = length(svdxs$d)  # rank of xs

    UTWU = t(svdxs$u) %*% sweep(svdxs$u, 1, w, "*") #  t(U) %*% diag(w) %*% U
    C = sweep(sweep(UTWU, 1, svdxs$d, "*"), 2, svdxs$d, "*") # D %*% UTWU %*% D
    C = (1-lambda) * h1 * C

    C = (C + t(C))/2  # symmetrize for numerical reasons (mpower() checks symmetry)

    # note: C is of size m x m, and diagonal if w=1/n

    if (lambda==0.0) # use eigenvalue decomposition computing the matrix power
    {
      if (m < p-sum(zeros))
        warning(paste("Estimated correlation matrix doesn't have full rank",
                      "- pseudoinverse used for inversion."), call. = FALSE)

      powr =  svdxs$v %*% tcrossprod( mpower(C, alpha), svdxs$v)
    }
    else # use a special identity for computing the matrix power
    {
      F = diag(m) - mpower(C/lambda + diag(m), alpha)
      powr = (diag(p) - svdxs$v %*% tcrossprod(F, svdxs$v))*(lambda)^alpha
    }

    # set all diagonal entries corresponding to zero-variance variables to 1
    diag(powr)[zeros] = 1

    rownames(powr) = colnames(xs)
    colnames(powr) = colnames(xs)
  }
  rm(xs)

  attr(powr, "lambda") = lambda
  attr(powr, "lambda.estimated") = lambda.estimated
  attr(powr, "class") = "shrinkage"

  if (verbose) cat("\n")

  return( powr )
}

