#' Control parameters for the sOG algorithm
#' 
#' @param maxiter An integer indicating the maximum number of iterations (within each loop). Default: 2000.
#' @param tol  A (small) positive number controlling the convergence of the algorith. The method converges when the abs difference between two consecutive iterations is smaller than this value.
#' @param lim A list with three elements. \code{lim_std} set the default upper limit for lambda_star. If \code{nlminb} cannot find a solution in \code{(0, lim_std)}, the upper bound is increased by \code{lim_step}, until a maximum of \code{lim_max}
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' @export
#' 

control_sOG = function(maxiter = 2e3, tol = .Machine$double.eps^.25, lim_std = 60, lim_max = 500, lim_step = 5) {
  if(lim_max < lim_std) stop("lim max must be greater or equal than lim_std")
       list(maxiter = maxiter, tol = tol, lim = list(lim_std = lim_std, lim_max = lim_max, lim_step = lim_step))
}


#' sparse Orthogonal to Groups analysis
#' @description Perform sparse matrix decomposition under group constraints.
#' @param X n x p data matrix to preprocess.
#' @param z n x 1 vector with group information. 
#' @param t_val penalty parameter. Must be a positive number. Control the strength of the l-1 penalty.
#' @param K approximation rank. The defaul is \code{max(2, round(NCOL(X)/100))}
#' @param rescaled Should the matrix \code{X} be rescaled? Default is \code{TRUE}.
#' @param control Control parameters for the fitting process. A list as returned by control_sOG.
#' @return  \itemize{
#' \item \verb{U} n x K matrix of scores.
#' \item \verb{S} K x p matrix of sparse loadings
#' }
#' @details The function performs a sparse matrix decomposition of the input matrix \code{X} with constraints on the left singular vectors and the group variable \code{z}.
#' The output is a matrix \code{U} of basis and a matrix \code{S} of scores such that \eqn{\tilde X = S * U^T}, where \eqn{\tilde X} is the rank-K approximation of \code{X}
#' @examples 
#' k.rid = 10
#' n = 5000
#' p = 200
#' W = matrix(rnorm(p*k.rid), k.rid)
#' S = matrix(rnorm(n*k.rid), n)
#' z = sample(rep(0:1, each=n/2))
#' lambda = rnorm( k.rid, mean = 0, sd = 1)
#' A = jitter( (S - lambda * z ) %*% W)
#' res = sOG(X = A, z = z, t_val = 10, K = 30)
#' @references Aliverti, Lum, Johndrow and Dunson (2018). Removing the influence of a group variable in high-dimensional predictive modelling (https://arxiv.org/abs/1810.08255).


sOG = function(X = warning('X missing'), 
	       z = stop('z missing'), 
	       t_val = stop('t_val missing'),
	       K = max(2, round(NCOL(X)/100)),
	       rescale = T,
	       control = control_sOG()
			  ){

  
  
	#++++++++++++
	# Controls
	#++++++++++++

	if(!is.matrix(X)) stop("X MUST BE A MATRIX")
	if(!is.numeric(t_val) | t_val <=0) stop("t_val must be a positive number")

	Nit_max = control$maxit
	eps = control$tol
	lim = control$lim

	# require(Matrix)
	if(rescale) X = scale(X)

	n = nrow(X)
	p = ncol(X)

	# initialise u and s at random
	u = rnorm(p)
	#u = numeric(p)+1
	u0 = u = u/sqrt(sum(u^2))

	s = rnorm(n)
	#s = numeric(n)+1
	s0 = s = s/sqrt(sum(s^2))

	#Fixed quantities
	nit=0
	err=10
	U = matrix(data=0, nrow = K, ncol = p)
	S = matrix(nrow = n, ncol = K)


	cat("############################################## \n")
	cat('## Cycle over other coord_j, j = 1, ...',K,  '##\n')
	cat("############################################## \n")

	P = matrix(0,n,n)
	for(j in (1:K))
	{
		# Re initialize quantities
		nit=0
		err=10
		u = u0
		s = s0

		# Slow part, projection matrix. Pretty unefficient.
		if( j !=1){
			PI = diag(n) - P
			px = PI %*% X
		}
		else{ px = X}

		cat(j,'-th starting!','\n')
		while (nit < control$maxiter & err > control$tol) {
			sold = s
			uold = u
			s = px %*% u
			s = s/sqrt(sum(s^2))
			temp = lm(s~z)$coef
			s = s - temp[1] - temp[2]*z
			s = s/sqrt(sum(s^2))

			Xts = t(X) %*% s
			lambda_grid = seq(0,20,l=100)

			if(sum(abs(Xts)) < t_val) 
			{
				u = Xts/sqrt(sum(Xts^2))
				# cat("--",' small u\n')
			}
			else
			{
				lambda_star = try(suppressWarnings(uniroot(function(x) sum(abs(softV(Xts,x))) - t_val, interval = c(0, lim$lim_std))$root))
				# Increase the bound

					new_lim = lim$lim_step + lim$lim_std
				while(!is.numeric(lambda_star) & new_lim < lim$lim_max){
					lambda_star = try(uniroot(function(x) 
								  sum(abs(softV(Xts,x))) - t_val, interval = c(0,new_lim))$root)
					new_lim = new_lim + lim$lim_step
				}
				u = softV(Xts,lambda_star)
			}

			nit = nit + 1
			err = sum(abs(s-sold)) + sum(abs(u-uold))
			if(nit %% 100 == 0) cat("-----", j,"-th coordinate", 'nit:',nit, " with err=", err,'\n')
		}

		U[j,] = u
		S[,j] = s
		P = P + s %*% t(s)
	}
	d = numeric(K)
	for(i in 1:K) d[i] = S[,i] %*% X %*% U[i,]
	S = S %*% diag(d)
	return(list(U=U, S=S))
}

