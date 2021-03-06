#' Orthogonal to Groups analysis
#' @description Perform matrix decomposition under group constraints.
#' @param X n x p data matrix to preprocess.
#' @param z n x 1 vector with group information. 
#' @param K approximation rank. The defaul is \code{max(2, round(NCOL(X)/100))}
#' @param rescaled Should the matrix \code{X} be rescaled? Default is \code{TRUE}.
#' @return  \itemize{
#' \item \verb{U} n x K matrix of scores.
#' \item \verb{S} K x p matrix of sparse loadings
#' }
#' @details The function performs a matrix decomposition of the input matrix \code{X} with constraints on the left singular vectors and the group variable \code{z}.
#' The output is a matrix \code{U} of basis and a matrix \code{S} of scores such that \eqn{\tilde X = S * U^T}, where \eqn{\tilde X} is the rank-K approximation of \code{X}
#' @references Aliverti, Lum, Johndrow and Dunson (2018). Removing the influence of a group variable in high-dimensional predictive modelling (https://arxiv.org/abs/1810.08255).


OG = function(X, z, K = max(2, round(NCOL(X)/10)), rescale = T) {

	if(!is.matrix(X)) stop("X must be a matrix")
	if(rescale) X = scale(X)
	SVD = svd(X, nu = K, nv = K)
	temp = lm(SVD$u %*% diag(SVD$d[1:K]) ~ z)
	S = SVD$u %*% diag(SVD$d[1:K]) - cbind(1,z)%*%temp$coef
	return(list(S = S, U = t(SVD$v)))
}
