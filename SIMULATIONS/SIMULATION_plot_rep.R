#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Reproduce simulations for the paper.
# The seed to create the data is always the same, while changing the seed before train/test splitting allows to evaluate results over different splits, as suggested by an anonymous referee
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = ls())

#+++++++++++++++++++++++++
# Takes roughly 10 minutes
#+++++++++++++++++++++++++
nsim = 50 
out = array(0, dim = c(nsim,12,5))
seeds = 1:nsim
preds = list()
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


	preds[[i]] = data.frame("COMBAT" = cbind(1,res_comp_test$u) %*% coef(m_comp),
		     "PSVA" = cbind(1, res_sva_test$u) %*% coef(m_sva),
		     "OG" = cbind(1, res_fpl_test$S) %*% coef(m1),
		     "SOG" = cbind(1, res_sfpl_test$S) %*% coef(m2),
		     "SVD" = cbind(1, SVD_test$u) %*% coef(m_SVD)
		     #"z" = z[id.test]
	)

	cat(i,'\n')
}

#If z is added later
predsN = list()
for(i in 1:length(preds)){
	seed = seeds[i]
	set.seed(seed)
	id.train = sample(1:n,.75*n)
	id.test = setdiff(1:n, id.train) 
	tmp = preds[[i]]
	tmp$z = z[id.test]
	predsN[[i]] = tmp
}

cor_res = matrix(unlist(lapply(predsN, function(x) cor(x)[6,-6])),nrow = nsim,byrow = T)
colnames(cor_res) = colnames(predsN[[1]])[-6]
boxplot(cor_res)

df = reshape2::melt(cor_res)
df$Var2 = factor(df$Var2, labels = c("COMBAT", "PSVA", "OG", "SOG", "SVD"))
require(latex2exp)
require(ggplot2)
ggplot(df, aes(x = Var2, y = value)) + geom_boxplot(size = 1.5) + theme_bw(base_size = 18) +
	#geom_jitter(aes(x = Var2, y = value), size = 2, alpha = .6, position=position_jitter(0.2)) +
		#theme(panel.grid.major = element_blank(), 
		      #panel.grid.minor = element_blank(), 
		      #panel.background = element_blank()) +
		scale_y_continuous(expand = c(0, 0.05)) +

		  #stat_summary(fun.y = mean, geom = "errorbar", 
               #aes(ymax = ..y.., ymin = ..y.., group = factor(Var2)),
               #width = 0.75, linetype = "dashed", position = position_dodge()) +

		xlab('')+
		ylab(TeX('Correlation between \\hat{Y} and Z'))+
		theme(strip.text = element_text(size=16),
		      #axis.title = element_text(size=18),
		      axis.text =  element_text(size=18),
		      axis.text.x =  element_text(angle = 0),
		      legend.position = 'none',
		      panel.spacing=unit(.5, "lines"),
		      panel.border = element_rect(color = "black", fill = NA, size = .5), 
		      strip.background = element_rect(color = "black", size = .5),
		      ) + guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(file = "boxplot.pdf", width = 20, height = 4.5)


