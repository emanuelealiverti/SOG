#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Reproduce simulations for the paper.
# The seed to create the data is always the same, while changing the seed before train/test splitting allows to evaluate results over different splits, as suggested by an anonymous referee
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = ls())

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This funcion computes the peformacnce of different approaches
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
err_matrix = function(predicted, truth)
{
	require(Metrics)
	mm = c("rmse","mae","mdae")
	sapply(mm, function(m) unlist(lapply(predicted, function(p) get(m)(p,truth))))
}

#+++++++++++++++++++++++++
# Takes roughly 30 minutes
#+++++++++++++++++++++++++
nsim = 50 
out = array(0, dim = c(nsim,12,5))
seeds = 1:nsim
for(i in 1:nsim){
	#this seed will be used only for the splitting
	seed = seeds[i]
	set.seed(2711)

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# Scenario 1. Low rank structures, Y depend on X and Z.
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++
	k.rid = 10
	#n = 500
	#p = 1000
	n = 1000
	p = 200

	W = matrix(rnorm(p*k.rid), k.rid)
	S = matrix(rnorm(n*k.rid), n)
	z=rep(0:1, each=n/2)

	lambda = rnorm(k.rid)
	A = (S - lambda * z ) %*% W
	X = scale(A)
	beta = runif(k.rid,-5,5)

	y = rnorm(n, mean = (S - lambda * z ) %*% beta )

	set.seed(seed)
	id.train = sample(1:n,.75*n)
	id.test = setdiff(1:n, id.train) 

	require(sOG)
	#Estimation
	K=10
	res_sfpl = sOG(X[id.train,], z=z[id.train], K=10,  t_val = 7,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl = OG(X=X[id.train,], z=z[id.train], K=10)
	res_comp = svd(t(sva::ComBat(dat = t(X[id.train,]),batch = as.vector(z[id.train]))),nu=K,nv=K)
	?sva::sva
	res_sva = svd(t(sva::psva(dat = t(X[id.train,]),batch = as.factor(z[id.train]),n.sv=K)),nu=K,nv=K)
	SVD = svd(X[id.train,], nu = K, nv = K)

	#Linear models
	y_train = y[id.train]
	y_test = y[id.test]

	m1 = lm(y_train ~ res_fpl$S)
	m2 = lm(y_train ~ res_sfpl$S)
	m_comp = lm(y_train ~ res_comp$u)
	m_sva = lm(y_train ~ res_sva$u)
	m_SVD = lm(y_train ~ SVD$u)

	#Test matrices
	res_sfpl_test =  sOG(X[id.test,], z=z[id.test], K=10,  t_val = 5,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl_test  =  OG(X=X[id.test,], z=z[id.test], K=10)
	res_comp_test =  svd(t(sva::ComBat(dat = t(X[id.test,]),batch = as.vector(z[id.test]))),nu=K,nv=K)
	res_sva_test = svd(t(sva::psva(dat = t(X[id.test,]),batch = as.factor(z[id.test]),n.sv=K)),nu=K,nv=K)
	SVD_test      =  svd(X[id.test,], nu = K, nv = K)


	preds = list("ComBat" = cbind(1,res_comp_test$u) %*% coef(m_comp),
		     "sva" = cbind(1, res_sva_test$u) %*% coef(m_sva),
		     "OG" = cbind(1, res_fpl_test$S) %*% coef(m1),
		     "sOG" = cbind(1, res_sfpl_test$S) %*% coef(m2),
		     "unADJ" = cbind(1, SVD_test$u) %*% coef(m_SVD)
	)


	res1 = err_matrix(predicted = preds,truth = y_test)

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# Scenario 2. Low rank structures, Y depend on X and not on Z.
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	set.seed(2711)

	k.rid = 10
	n = 1000
	p = 200

	W = matrix(rnorm(p*k.rid), k.rid)
	S = matrix(rnorm(n*k.rid), n)
	z=rep(0:1, each=n/2)

	lambda = rnorm(k.rid)
	A = (S - lambda * z ) %*% W
	#A = matrix(rnorm(n*p,sd=2),n) 
	#A = (S + lambda*z)%*%W
	X = scale(A)
	#lattice::levelplot(cor(X))
	beta = runif(k.rid,-5,5)

	y = rnorm(n, mean = S %*% beta )

	set.seed(seed)
	id.train = sample(1:n,.75*n)
	id.test = setdiff(1:n, id.train) 

	require(sOG)
	#Estimation
	K=10
	res_sfpl = sOG(X[id.train,], z=z[id.train], K=10,  t_val = 7,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl = OG(X=X[id.train,], z=z[id.train], K=10)
	res_comp = svd(t(sva::ComBat(dat = t(X[id.train,]),batch = as.vector(z[id.train]))),nu=K,nv=K)
	res_sva = svd(t(sva::psva(dat = t(X[id.train,]),batch = as.factor(z[id.train]),n.sv=K)),nu=K,nv=K)
	SVD = svd(X[id.train,], nu = K, nv = K)

	#Linear models
	y_train = y[id.train]
	y_test = y[id.test]

	m1 = lm(y_train ~ res_fpl$S)
	m2 = lm(y_train ~ res_sfpl$S)
	m_comp = lm(y_train ~ res_comp$u)
	m_sva = lm(y_train ~ res_sva$u)
	m_SVD = lm(y_train ~ SVD$u)

	#Test matrices
	res_sfpl_test =  sOG(X[id.test,], z=z[id.test], K=10,  t_val = 5,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl_test  =  OG(X=X[id.test,], z=z[id.test], K=10)
	res_comp_test =  svd(t(sva::ComBat(dat = t(X[id.test,]),batch = as.vector(z[id.test]))),nu=K,nv=K)
	res_sva_test = svd(t(sva::psva(dat = t(X[id.test,]),batch = as.factor(z[id.test]),n.sv=K)),nu=K,nv=K)
	SVD_test      =  svd(X[id.test,], nu = K, nv = K)


	preds = list("ComBat" = cbind(1,res_comp_test$u) %*% coef(m_comp),
		     "sva" = cbind(1, res_sva_test$u) %*% coef(m_sva),
		     "OG" = cbind(1, res_fpl_test$S) %*% coef(m1),
		     "sOG" = cbind(1, res_sfpl_test$S) %*% coef(m2),
		     "unADJ" = cbind(1, SVD_test$u) %*% coef(m_SVD)
	)


	(res2 = err_matrix(predicted = preds,truth = y_test))


	#+++++++++++++++++++
	# Scenario 3: p > n 
	#+++++++++++++++++++

	#rm(list=ls())
	set.seed(2711)

	k.rid = 10
	p = 1000
	n = 200

	W = matrix(rnorm(p*k.rid), k.rid)
	S = matrix(rnorm(n*k.rid), n)
	z=rep(0:1, each=n/2)

	lambda = rnorm(k.rid)
	A = (S - lambda * z ) %*% W
	#A = matrix(rnorm(n*p,sd=2),n) 
	#A = (S + lambda*z)%*%W
	X = scale(A)
	#lattice::levelplot(cor(X))
	beta = rnorm(k.rid, sd = .2)
	beta = runif(k.rid,-5,5)

	y = rnorm(n, mean = (S - lambda * z ) %*% beta )

	set.seed(seed)
	id.train = sample(1:n,.75*n)
	id.test = setdiff(1:n, id.train) 

	require(sOG)
	#Estimation
	K=10
	res_sfpl = sOG(X[id.train,], z=z[id.train], K=10,  t_val = 7,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl = OG(X=X[id.train,], z=z[id.train], K=10)
	res_comp = svd(t(sva::ComBat(dat = t(X[id.train,]),batch = as.vector(z[id.train]))),nu=K,nv=K)
	res_sva = svd(t(sva::psva(dat = t(X[id.train,]),batch = as.factor(z[id.train]),n.sv=K)),nu=K,nv=K)
	SVD = svd(X[id.train,], nu = K, nv = K)

	#Linear models
	y_train = y[id.train]
	y_test = y[id.test]

	m1 = lm(y_train ~ res_fpl$S)
	m2 = lm(y_train ~ res_sfpl$S)
	m_comp = lm(y_train ~ res_comp$u)
	m_sva = lm(y_train ~ res_sva$u)
	m_SVD = lm(y_train ~ SVD$u)

	#Test matrices
	res_sfpl_test =  sOG(X[id.test,], z=z[id.test], K=10,  t_val = 5,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl_test  =  OG(X=X[id.test,], z=z[id.test], K=10)
	res_comp_test =  svd(t(sva::ComBat(dat = t(X[id.test,]),batch = as.vector(z[id.test]))),nu=K,nv=K)
	res_sva_test = svd(t(sva::psva(dat = t(X[id.test,]),batch = as.factor(z[id.test]),n.sv=K)),nu=K,nv=K)
	SVD_test      =  svd(X[id.test,], nu = K, nv = K)


	preds = list("ComBat" = cbind(1,res_comp_test$u) %*% coef(m_comp),
		     "sva" = cbind(1, res_sva_test$u) %*% coef(m_sva),
		     "OG" = cbind(1, res_fpl_test$S) %*% coef(m1),
		     "sOG" = cbind(1, res_sfpl_test$S) %*% coef(m2),
		     "unADJ" = cbind(1, SVD_test$u) %*% coef(m_SVD)
	)

	(res3 = err_matrix(predicted = preds,truth = y_test))


	#+++++++++++++++++++++++++++++++++++++++
	# Scenario 4: no low rank structure
	#+++++++++++++++++++++++++++++++++++++++

	#rm(list=ls())
	set.seed(2711)

	k.rid = 200
	#n = 500
	#p = 1000
	n = 1000
	p = 200

	W = matrix(rnorm(p*k.rid), k.rid)
	S = matrix(rnorm(n*k.rid), n)
	z=rep(0:1, each=n/2)

	lambda = rnorm(k.rid)
	A = (S - lambda * z ) %*% W
	X = scale(A)
	beta = runif(k.rid,-5,5)

	y = rnorm(n, mean = (S - lambda * z ) %*% beta )

	set.seed(seed)
	id.train = sample(1:n,.75*n)
	id.test = setdiff(1:n, id.train) 

	require(sOG)
	#Estimation
	K=10
	res_sfpl = sOG(X[id.train,], z=z[id.train], K=10,  t_val = 7,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl = OG(X=X[id.train,], z=z[id.train], K=10)
	res_comp = svd(t(sva::ComBat(dat = t(X[id.train,]),batch = as.vector(z[id.train]))),nu=K,nv=K)
	res_sva = svd(t(sva::psva(dat = t(X[id.train,]),batch = as.factor(z[id.train]),n.sv=K)),nu=K,nv=K)
	SVD = svd(X[id.train,], nu = K, nv = K)

	#Linear models
	y_train = y[id.train]
	y_test = y[id.test]

	m1 = lm(y_train ~ res_fpl$S)
	m2 = lm(y_train ~ res_sfpl$S)
	m_comp = lm(y_train ~ res_comp$u)
	m_sva = lm(y_train ~ res_sva$u)
	m_SVD = lm(y_train ~ SVD$u)

	#Test matrices
	res_sfpl_test =  sOG(X[id.test,], z=z[id.test], K=10,  t_val = 5,control = control_sOG(lim_step = 1,lim_max=5000))
	res_fpl_test  =  OG(X=X[id.test,], z=z[id.test], K=10)
	res_comp_test =  svd(t(sva::ComBat(dat = t(X[id.test,]),batch = as.vector(z[id.test]))),nu=K,nv=K)
	res_sva_test = svd(t(sva::psva(dat = t(X[id.test,]),batch = as.factor(z[id.test]),n.sv=K)),nu=K,nv=K)
	SVD_test      =  svd(X[id.test,], nu = K, nv = K)


	preds = list("ComBat" = cbind(1,res_comp_test$u) %*% coef(m_comp),
		     "sva" = cbind(1, res_sva_test$u) %*% coef(m_sva),
		     "OG" = cbind(1, res_fpl_test$S) %*% coef(m1),
		     "sOG" = cbind(1, res_sfpl_test$S) %*% coef(m2),
		     "unADJ" = cbind(1, SVD_test$u) %*% coef(m_SVD)
	)


	res4 = err_matrix(predicted = preds,truth = y_test)
	resFULL =rbind(t(res1),t(res2),t(res3),t(res4))
	out[i,,] = resFULL
	cat(i,'\n')
}

#++++++++++++++
# print results
#++++++++++++++

tot = apply(out, c(2,3),mean)
sds = apply(out, c(2,3), sd)
sds
dimnames(tot) = dimnames(resFULL)

make_tab_v = function(cc, v) {
	out_ic  = cc*0
	for(i in 1:NROW(out_ic)) {
		for(j in 1:NCOL(out_ic)) {
			out_ic[i,j] = paste0(sprintf('%.2f',cc[i,j]), " (",sprintf('%.2f',v[i,j]), ")")
		}
	}
	return(out_ic)
}


xtable::xtable(make_tab_v(tot,sds))
