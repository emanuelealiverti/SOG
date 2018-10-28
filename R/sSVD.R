###################
# Computes spars PCA (adapting sFPML removing fairnsess constraint)
###################
# rm(list = ls())

sSVD = function(X = warning('X missing'), 
                  t_val = stop('t_val missing'), 
                  K = round(NCOL(X)/100), 
                  ctr = list(Nit_max = 2e3, 
                             eps = .Machine$double.eps^.25)){
  if(!is.matrix(X)) stop("X MUST BE A MATRIX")
  Nit_max = ctr$Nit_max
  eps = ctr$eps
  require(Matrix)
  ###################
  # Soft thresholding
  ###################
  require(Rcpp)
  sourceCpp("MAIN_FUNCTIONS/softV.cpp")

  n = nrow(X)
  p = ncol(X)
  
  # initialise u and s at random
  # u = rnorm(p)
  u = numeric(p)+1
  u0 = u = u/sqrt(sum(u^2))
  
  # s = rnorm(n)
  s = numeric(n)+1
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
    while (nit < Nit_max & err > eps) {
      sold = s
      uold = u
      s = px %*% u
      s = s/sqrt(sum(s^2))
      # temp = lm(s~z)$coef
      # s = s - temp[1] - temp[2]*z
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
        lambda_star = try(suppressWarnings(uniroot(function(x) sum(abs(softV(Xts,x))) - t_val, interval = c(0,60))$root))
        # Increase the bound
        
          new_lim = 70
        while(!is.numeric(lambda_star) & new_lim < 200){
          lambda_star = try(uniroot(function(x) sum(abs(softV(Xts,x))) - t_val, interval = c(0,new_lim))$root)
          new_lim = new_lim + 5
        }
        # lambda_star = uniroot.all(function(x) sum(abs(softV(Xts,x))) - t_val, interval = c(0,20))
        # lambda_star = nlminb(function(x) {abs(sum(abs(softV(Xts,x))) - t_val)}, start = 0, lower=0)$par
        
        # val_grid = softV(Xts, lambda_grid)
        # val_grid_abs = colSums(abs(val_grid))
        # best = which.min(abs(val_grid_abs - t_val))
        # cat(min(abs(val_grid_abs - t_val)),'\n')
        # u = val_grid[,best]
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
  for(i in 1:K) d[i] = S[,i]%*% X %*% U[i,]
  S = S %*% diag(d)
  return(list(U=Matrix(U), S=S))
}
