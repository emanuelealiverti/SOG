#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Reproduce simulations for the paper.
# The seed to create the data is always the same, while changing the seed before train/test splitting allows to evaluate results over different splits, as suggested by an anonymous referee
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = ls())

#+++++++++++++++++++++++++
# Takes roughly 40 minutes
#+++++++++++++++++++++++++
nsim = 50
r_max = 11
out = array(0, dim = c(nsim,12,5))
seeds = 1:nsim
#this seed will be used only for the splitting
out = array(NA, dim = c(nsim, r_max-1,5))
dimnames(out) = list(NULL, paste0("rank ", 2:r_max),c("COMBAT", "PSVA",  "OG", "SOG", "SVD"))
for(i in 1:nsim) {
	seed = seeds[i]
	#++++++++++++++++++++++++++++++++++++++
	# Focus only on scenario 1 of the paper
	#++++++++++++++++++++++++++++++++++++++
	set.seed(2711)
	k.rid = 10
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

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# We run everything for the maximum rank, and then consider all the lowest version.
	# The trick does not woek with sva since it always provides X as a result
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	set.seed(seed)
	id.train = sample(1:n,.75*n)
	id.test = setdiff(1:n, id.train) 
	require(sOG)
	res_sfpl = sOG(X[id.train,], z=z[id.train], K=r_max,  
		       t_val = 7,control = control_sOG(lim_step = 1,maxiter = 300))
	res_fpl = OG(X=X[id.train,], z=z[id.train], K=r_max)
	res_comp = svd(t(sva::ComBat(dat = t(X[id.train,]),batch = as.vector(z[id.train]))), nu = r_max, nv = r_max)
	SVD = svd(X[id.train,], nu = r_max, nv = r_max)

	#res_sva = t(sva::psva(dat = t(X[id.train,]),batch = as.factor(z[id.train]),n.sv=r_max))
	for(k_curr in 2:r_max) {
		K=k_curr

		# Need to rerun
		res_sva = t(sva::psva(dat = t(X[id.train,]),batch = as.factor(z[id.train]),n.sv=K))

		out[seed,k_curr-1,1] = norm(X[id.train,] - residuals(lm(res_comp$u[,1:K]~z[id.train])) %*% diag(res_comp$d[1:K])%*%t(res_comp$v[,1:K]),"F")
		out[seed,k_curr-1,2] =norm(X[id.train,] - res_sva,"F")
		out[seed,k_curr-1,3] =norm(X[id.train,] - res_fpl$S[,1:K] %*% res_fpl$U[1:K,] ,"F")
		out[seed,k_curr-1,4] =norm(X[id.train,] - res_sfpl$S[,1:K] %*% res_sfpl$U[1:K,] ,"F")
		out[seed,k_curr-1,5] =norm(X[id.train,] - SVD$u[,1:K] %*% diag(SVD$d[1:K])%*%t(SVD$v[,1:K]),"F")
		cat("\n", "#################", '\n')
		cat("split ", seed, " rank ", K, '\n', "COMPLETED",'\n')
		cat("#################", '\n')
	}
}
#save(out, file = "Ranks50.RData")
out_rel = out/norm(X,"F")
#save(out_rel, file = "Ranks_rel50.RData")
# Original data


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The plot is very confusing so I have preferred to show a table
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

df_complete = reshape2::melt(apply(out, c(2,3), mean))
df_up = reshape2::melt(apply(out, c(2,3), quantile, .9))
df_down = reshape2::melt(apply(out, c(2,3), quantile, .1))
colnames(df_var) = colnames(df_down) = colnames(df_up) = colnames(df_complete) = c("K", "Method","value")
df_error = cbind(df_up, df_down$value)
names(df_error)[3:4] = c("up","low")
head(df_error)
require(ggplot2)
rank_pl = ggplot(df_complete) + 
	geom_point( aes(x = K, y =(value), group = Method, colour = Method, shape=Method),size=3) + 
	geom_line( aes(x = K, y =(value), group = Method, colour = Method, linetype=Method)) +

	geom_errorbar(data = df_error, aes(x = K, ymin =(low), ymax =(up),colour = Method)) +
	#scale_color_manual(values = c('#000000', '#808080', '#BEBEBE'))+
	#facet_wrap(~w) + 
	theme_bw() +
	ylab('Frobenius norm of reconstruction error')+
	xlab('Approximation rank')+
	ylim(c(0, 400))+
	#xlim(c(1, 100))+
	theme(
	      # axis.title.x=element_blank(),
	      axis.title = element_text(size=18),
	      # axis.text.x=element_blank(),
	      axis.ticks.x=element_blank(),
	      # axis.text.y=element_text('Reconstruction error X - X'),
	      # axis.text.y=element_blank(),
	      axis.ticks.y=element_blank(),
	      legend.title = element_blank(),
	      legend.key = element_rect(colour = 'black'),
	      legend.text = element_text(size=18),
	      legend.key.size = unit(2, 'lines'),
	      legend.position = 'top',
	      panel.spacing=unit(.5, "lines"),
	      panel.border = element_rect(color = "black", fill = NA, size = .5), 
	      strip.background = element_rect(color = "black", size = .5),
	      strip.text = element_text(size=18, face ='bold')) + 
guides(points.shape = guide_legend(override.aes = list(size = 2)))
rank_pl

# at this point I think a table is better
# This function makes sth like value (ic1, ic2)
make_tab = function(cc, l, u, dig = 2) {
	out_ic  = cc*0
	for(i in 1:NROW(out_ic)) {
		for(j in 1:NCOL(out_ic)) {
			out_ic[i,j] = paste0(sprintf('%.2f',cc[i,j]), " (",sprintf('%.2f',l[i,j]), ",", sprintf('%.2f',u[i,j]), ")")
		}
	}
	return(out_ic)
}

make_tab_v = function(cc, v) {
	out_ic  = cc*0
	for(i in 1:NROW(out_ic)) {
		for(j in 1:NCOL(out_ic)) {
			out_ic[i,j] = paste0(sprintf('%.2f',cc[i,j]), " (",sprintf('%.2f',v[i,j]), ")")
		}
	}
	return(out_ic)
}

cc = apply(out, c(2,3), mean)
l = apply(out, c(2,3), min)
u= apply(out, c(2,3), max)
(v= apply(out, c(2,3), sd))

make_tab(cc,l,u,2)
make_tab_v(cc,v)
xtable::xtable(make_tab(cc,l,u,2))
xtable::xtable(make_tab_v(cc,v))
